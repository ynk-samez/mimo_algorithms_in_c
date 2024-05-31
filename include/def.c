
#include "def.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <time.h>
#include <unistd.h>
#include <math.h>
#include "RFunction.h"
#include "const.h"
#include "complexArith.h"

void printInfo()
{
    printf("\e[32m");
    printf("Eb/N0 : %d --> %d\n", ENMIN, ENMAX);
    printf("NT = %d\n", NT);
    printf("NR = %d\n", NR);
    printf("Send = %e\n", LOOP);
    printf("\e[m\n");
}

void MatTrans(complex **A, complex **ans, int A_row, int A_column)
{
    int i, j;

    for (i = 0; i < A_row; i++)
    {
        for (j = 0; j < A_column; j++)
        {
            A[i][j].im = -A[i][j].im;
            ans[j][i] = A[i][j];
            A[i][j].im = -A[i][j].im;
        }
    }
}

void VectTrans(complex *A, complex **ans, int value)
{
    int i;

    for (i = 0; i < value; i++)
    {
        A[i].im = -A[i].im;
        ans[0][i] = A[i];
        A[i].im = -A[i].im;
    }
}

double ComAbsSqu(complex A)
{
    double abs;
    abs = pow(A.re, 2) + pow(A.im, 2);
    return abs;
}

double ComAbsSqu2(complex A, int *com_mul_count, int *re_add_count)
{
    double abs;
    abs = pow(A.re, 2) + pow(A.im, 2);
    com_mul_count[0] += 2;
    re_add_count[0] += 1;
    return abs;
}

double STV(double en) // 標準偏差を算出
{
    double std;
    std = 2 * pow(10, en / 10);
    std = sqrt(1 / std);
    return std;
}

double STV_2(double en)
{
    double q;
    q = pow(10, en / 10);
    q = sqrt(1 / q);
    return q;
}

double STV_3(double en) // 雑音の標準偏差を算出
{
    double std2;
    std2 = pow(10, en / 10);
    std2 = sqrt(1 / std2);
    return std2;
}

complex wfading(complex h, long int t) // Hの作成
{
    double *alpha;
    double omega, beta = 0.0;
    double xc, xs = 0.0;
    double xc_sum, xs_sum;
    long int phase_shift = F_WAVES * 4 + 2;
    long int j;

    alpha = (double *)calloc(F_WAVES, sizeof(double));
    // ランダム位相生成

    for (j = 0; j < F_WAVES; ++j)
    {
        alpha[j] = 2 * M_PI * (double)rand() / (RAND_MAX + 1.0);
    }
    // Fading係数生成
    xc_sum = 0.0;
    xs_sum = 0.0;

    for (j = 0; j < F_WAVES; ++j)
    {
        omega = 2 * M_PI * (double)DOPPLER * cos((2 * M_PI * (j + 1)) / phase_shift);
        beta = (M_PI * (j + 1)) / F_WAVES;
        xc = cos(beta) * cos(omega * t + alpha[j]);
        xs = sin(beta) * sin(omega * t + alpha[j]);
        xc_sum += xc;
        xs_sum += xs;
    }

    h.re = (2 * xc_sum + cos(2 * M_PI * DOPPLER * t)) / sqrt(2 * F_WAVES + 1);
    h.im = (2 * xs_sum + sin(2 * M_PI * DOPPLER * t)) / sqrt(2 * F_WAVES + 1);

    free(alpha);

    return h;
}

complex gaussianRand(double vd_stdev, double vd_meanValue) //  標準偏差σ  平均
{
    complex n;

    double randN_1; //  乱数1
    double randN_2; //  乱数2
    double wgn_1;   //  ホワイトガウスノイズ1
    double wgn_2;   //  ホワイトガウスノイズ2

    randN_1 = ((double)rand()) / ((double)RAND_MAX); //  (0,9]の乱数取得
    randN_2 = ((double)rand()) / ((double)RAND_MAX); //  (0,9]の乱数取得

    wgn_1 = vd_stdev * sqrt(-2.0 * log(randN_1)) * cos(2.0 * M_PI * randN_2) + vd_meanValue;
    wgn_2 = vd_stdev * sqrt(-2.0 * log(randN_1)) * sin(2.0 * M_PI * randN_2) + vd_meanValue;

    n.re = wgn_1;
    n.im = wgn_2;

    return n;
}

complex rand_complex(complex x) // ±1をランダムで出力する
{
    x.re = (rand() % 2) * 2 - 1;
    x.im = (rand() % 2) * 2 - 1;
    return (x);
}

double count(complex x, complex y, double k)
{
    if (x.re != y.re)
    {
        k++;
    }
    if (x.im != y.im)
    {
        k++;
    }
    return (k);
}

double Int_Count(int x, int y, double k)
{ // エラーカウント
    if (x != y)
    {
        k++;
    }
    return (k);
}

double scalar(complex A)
{
    double sca;
    sca = sqrt(pow(A.re, 2) + pow(A.im, 2));
    return sca;
}

void IdentitystvComMat(complex **imat, double stv, int n) // 単位複素数行列と雑音の標準偏差の乗算
{
    int i, j;
    complex one, zero;
    one.re = stv;
    one.im = 0;
    zero.re = 0;
    zero.im = 0;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j) // 行と列が同じとき標準偏差を代入
                imat[i][j] = one;
            else
                imat[i][j] = zero;
        }
    }
}

void Identitystv2ComMat(complex **imat, double stv2, int n) // 単位複素数行列と雑音の標準偏差の2乗の乗算
{
    int i, j;
    complex one, zero;
    one.re = stv2;
    one.im = 0;
    zero.re = 0;
    zero.im = 0;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j) // 行と列が同じとき雑音の2乗を代入
                imat[i][j] = one;
            else
                imat[i][j] = zero;
        }
    }
}

void MatMatAddCom(complex **A, complex **B, complex **ans, int row, int col)
{
    int i, j;

    InitComMat(ans, row, col); // 複素数行列の初期化
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            ans[i][j].re = A[i][j].re + B[i][j].re;
            ans[i][j].im = A[i][j].im + B[i][j].im;
        }
    }
}

void MatMatSubCom(complex **A, complex **B, complex **ans, int row, int col)
{
    int i, j;

    InitComMat(ans, row, col);

    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            ans[i][j].re = A[i][j].re - B[i][j].re;
            ans[i][j].im = A[i][j].im - B[i][j].im;
        }
    }
}

void MatMatMulCom(complex **A, complex **B, complex **ans, int Arow, int Acol, int Bcol)
{
    int i, j, k;

    InitComMat(ans, Arow, Bcol);

    for (i = 0; i < Arow; i++)
    {
        for (j = 0; j < Bcol; j++)
        {
            for (k = 0; k < Acol; k++)
            {
                ans[i][j] = complex_add(ans[i][j], complex_mul(A[i][k], B[k][j]));
            }
        }
    }
}
void Mat4MatMulCom(complex ****A, complex **B, complex **ans, int Arow, int Acol, int A3, int Bcol)
{ // 4次元の配列と2次元配列の掛け算(A3はA[][][A3][]の箇所)
    int i, j, k, l;
    InitComMat(ans, Arow, A3);

    for (i = 0; i < Arow; i++)
    {
        for (j = 0; j < A3; j++)
        {
            for (k = 0; k < Acol; k++)
            {
                for (l = 0; l < Bcol; l++)
                {
                    ans[i][j] = complex_add(ans[i][j], complex_mul(A[i][k][j][l], B[k][l]));
                }
            }
        }
    }
}

void InitComMat(complex **init, int row, int col) // 複素数行列の初期化
{
    int i, j;
    complex zero = {0.0, 0.0};

    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            init[i][j] = zero;
        }
    }
}

void InitComVect(complex *init, int value)
{ // 複素数ベクトル行列の初期化
    int i;
    complex zero = {0.0, 0.0};

    for (i = 0; i < value; i++)
        init[i] = zero;
}

void InitMat(double **init, int col, int row)
{ // 行列の初期化
    int i, j;

    for (i = 0; i < col; i++)
    {
        for (j = 0; j < row; j++)
        {
            init[i][j] = 0.0;
        }
    }
}

void InitDMat(double **init, int row, int col) // double行列の初期化
{
    int i, j;

    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            init[i][j] = 0;
        }
    }
}

void InitDVect(double *init, int value)
{ // doubleベクトル行列の初期化
    int i;

    for (i = 0; i < value; i++)
        init[i] = 0;
}

complex Hard(complex x) // 硬判定
{
    complex z;
    if (x.re > 0)
        z.re = 1;
    if (x.re < 0)
        z.re = -1;
    if (x.im > 0)
        z.im = 1;
    if (x.im < 0)
        z.im = -1;
    return z;
}

void MatVectMulCom(complex **A, complex *B, complex *ans, int row, int col) // ベクトル行列と行列の乗算
{
    int i, j;
    InitComVect(ans, row);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            ans[i] = complex_add(complex_mul(A[i][j], B[j]), ans[i]);
        }
    }
}

void VectAddCom(complex *A, complex *B, complex *ans, int value) // ベクトル行列の加算
{
    int i;
    InitComVect(ans, value);
    for (i = 0; i < value; i++)
    {
        ans[i] = complex_add(A[i], B[i]);
    }
}

void VectSubCom(complex *A, complex *B, complex *ans, int value) // ベクトル行列の減算
{
    int i;
    InitComVect(ans, value);
    for (i = 0; i < value; i++)
    {
        ans[i] = complex_sub(A[i], B[i]);
    }
}

complex ****calloc_com4d(int value1, int value2, int value3, int value4) // 2列以上の行列のメモリ領域の確保
{
    complex ****c;
    int i, j, k;

    c = (complex ****)calloc(value1, sizeof(complex ***));
    if (c == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    for (i = 0; i < value1; i++)
    {
        c[i] = (complex ***)calloc(value2, sizeof(complex **));
        if (c[i] == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
        for (j = 0; j < value2; j++)
        {
            c[i][j] = (complex **)calloc(value3, sizeof(complex *));
            if (c[i][j] == NULL)
            {
                perror("MemError.\n");
                exit(1);
            }
            for (k = 0; k < value3; k++)
            {
                c[i][j][k] = (complex *)calloc(value4, sizeof(complex));
                if (c[i][j][k] == NULL)
                {
                    perror("MemError.\n");
                    exit(1);
                }
            }
        }
    }
    return c;
}

complex ***calloc_com3d(int value1, int value2, int value3) // 2列以上の行列のメモリ領域の確保
{
    complex ***c;
    int i, j;

    c = (complex ***)calloc(value1, sizeof(complex **));
    if (c == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    for (i = 0; i < value1; i++)
    {
        c[i] = (complex **)calloc(value2, sizeof(complex *));
        if (c[i] == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
        for (j = 0; j < value2; j++)
        {
            c[i][j] = (complex *)calloc(value3, sizeof(complex));
            if (c[i][j] == NULL)
            {
                perror("MemError.\n");
                exit(1);
            }
        }
    }
    return c;
}

complex **recalloc_com2d(complex **arr2d, int newCol, int newRow) // 2列以上の行列のメモリ領域の確保
{
    int x, y;
    arr2d = realloc(arr2d, sizeof(complex *) * newRow);
    for (y = 0; y < newRow; y++)
    {
        arr2d[y] = realloc(arr2d[y], sizeof(complex) * newCol);
        if (arr2d[y] == NULL)
        {
            // TODO: handle memory allocation fail error
            exit(1);
        }
    }
    return arr2d;
}

complex **calloc_com2d(int w, int h) // 2列以上の行列のメモリ領域の確保
{
    int i, j;
    complex **c = (complex **)calloc(w, sizeof(complex *));
    for (i = 0; i < w; i++)
    {
        c[i] = (complex *)calloc(h, sizeof(complex));
    }

    return c;
}

complex *calloc_com(int value) // 1列の行列のメモリ領域の確保
{
    complex *ptr;
    ptr = (complex *)calloc(value, sizeof(complex));
    if (ptr == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    return ptr;
}

double *calloc_d(int value)
{
    double *d;

    d = (double *)calloc(value, sizeof(double));
    if (d == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    return d;
}

double **calloc_d2d(int value1, int value2)
{
    double **d;
    int i;

    d = (double **)calloc(value1, sizeof(double *));
    if (d == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    for (i = 0; i < value1; i++)
    {
        d[i] = (double *)calloc(value2, sizeof(double));
        if (d[i] == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
    }
    return d;
}

double ***calloc_d3d(int value1, int value2, int value3)
{
    double ***d;
    int j, k;

    d = (double ***)calloc(value1, sizeof(double **));
    if (d == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    for (j = 0; j < value1; j++)
    {
        d[j] = (double **)calloc(value2, sizeof(double *));
        if (d[j] == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
        for (k = 0; k < value2; k++)
        {
            d[j][k] = (double *)calloc(value3, sizeof(double));
            if (d[j][k] == NULL)
            {
                perror("MemError.\n");
                exit(1);
            }
        }
    }
    return d;
}

double ****calloc_d4d(int value1, int value2, int value3, int value4)
{
    double ****d;
    int i, j, k;

    d = (double ****)calloc(value1, sizeof(double ***));
    if (d == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }

    for (i = 0; i < value1; i++)
    {
        d[i] = (double ***)calloc(value2, sizeof(double **));
        if (d == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
        for (j = 0; j < value2; j++)
        {
            d[i][j] = (double **)calloc(value3, sizeof(double *));
            if (d[i][j] == NULL)
            {
                perror("MemError.\n");
                exit(1);
            }
            for (k = 0; k < value3; k++)
            {
                d[i][j][k] = (double *)calloc(value4, sizeof(double));
                if (d[i][j][k] == NULL)
                {
                    perror("MemError.\n");
                    exit(1);
                }
            }
        }
    }
    return d;
}

int **calloc_i2d(int value1, int value2)
{
    int **i;
    int j;

    i = (int **)calloc(value1, sizeof(int *));
    if (i == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    for (j = 0; j < value1; j++)
    {
        i[j] = (int *)calloc(value2, sizeof(int));
        if (i[j] == NULL)
        {
            perror("MemError.\n");
            exit(1);
        }
    }
    return i;
}

int *calloc_i(int value)
{
    int *i;

    i = (int *)calloc(value, sizeof(int));
    if (i == NULL)
    {
        perror("MemError.\n");
        exit(1);
    }
    return i;
}

void ShinFree2d(complex ***c, int w) // メモリの開放
{
    int i;

    for (i = 0; i < w; i++)
    {
        free((*c)[i]);
    }
    free(*c);
    *c = NULL;
}

void Free_com2d(complex **c, int w) // メモリの開放
{
    int i;

    for (i = 0; i < w; i++)
    {
        free(c[i]);
    }
    free(c);
    c = NULL;
}

void Free_com3d(complex ***c, int value1, int value2)
{
    int i;

    for (i = 0; i < value1; i++)
    {
        Free_com2d(c[i], value2);
    }
    free(c);
    c = NULL;
}

void Free_com4d(complex ****c, int value1, int value2, int value3)
{
    int i;

    for (i = 0; i < value1; i++)
    {
        Free_com3d(c[i], value2, value3);
    }
    free(c);
    c = NULL;
}

void Free_d2d(double **d, int value)
{
    int i;

    for (i = 0; i < value; i++)
    {
        free(d[i]);
    }
    free(d);
}

void Free_d3d(double ***d, int value1, int value2)
{
    int i;

    for (i = 0; i < value1; i++)
    {
        Free_d2d(d[i], value2);
    }
    free(d);
    d = NULL;
}

void Free_d4d(double ****d, int value1, int value2, int value3)
{
    int i;

    for (i = 0; i < value1; i++)
    {
        Free_d3d(d[i], value2, value3);
    }
    free(d);
    d = NULL;
}

void Free_i2d(int **i, int value)
{
    int j;

    for (j = 0; j < value; j++)
    {
        free(i[j]);
    }
    free(i);
}

void IdentityComMat(complex **imat, int n) // 単位複素数行列の生成
{
    int i, j;
    complex one = {1.0, 0.0};
    complex zero = {0.0, 0.0};

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
                imat[i][j] = one;
            else
                imat[i][j] = zero;
        }
    }
}

void InvComMat(complex **A, complex **inv_A, int n) // 複素数逆行列の生成
{
    complex p;
    int i, j, k;

    IdentityComMat(inv_A, n); // 単位行列生成

    for (i = 0; i < n; i++)
    { // 掃き出し法
        p = A[i][i];
        for (j = 0; j < n; j++)
        {
            A[i][j] = complex_div(A[i][j], p);         // 対角成分を1にする
            inv_A[i][j] = complex_div(inv_A[i][j], p); // 単位行列もpで割る
        }
        for (j = 0; j < n; j++)
        {
            if (i != j)
            { // 対角成分以外の計算
                p = A[j][i];
                for (k = 0; k < n; k++)
                {
                    A[j][k] = complex_sub(A[j][k], complex_mul(A[i][k], p));
                    inv_A[j][k] = complex_sub(inv_A[j][k], complex_mul(inv_A[i][k], p));
                }
            }
        }
    }
}

void MMSE(complex **h, complex *y, complex *ans, double stv2, int nt, int nr)
{
    complex **h_h; // hの逆行列 hの共役転置行列
    complex **h3, **Ih3, **Ih3_IM, **w, *wy;
    complex **I; // σ^2と単位行列の積
    int i;

    h_h = calloc_com2d(NT, NR);    // hの共役転置行列
    h3 = calloc_com2d(NR, NR);     // hとhの共役転置行列の積
    Ih3 = calloc_com2d(NR, NR);    // h3とσ^2Iの和
    Ih3_IM = calloc_com2d(NR, NR); // 逆行列
    w = calloc_com2d(NT, NR);
    wy = calloc_com(NT);
    I = calloc_com2d(NR, NR);

    // 共役転置行列H^Hの作成
    MatTrans(h, h_h, NR, NT);

    // hとh_hの積
    MatMatMulCom(h, h_h, h3, NR, NT, NR);

    // σ^2Iの生成
    Identitystv2ComMat(I, stv2, NR);

    // Ih3 = H*H^H + σ^2*I
    MatMatAddCom(h3, I, Ih3, NR, NR);

    // Ih3の逆行列
    InvComMat(Ih3, Ih3_IM, NR);

    // wの生成
    MatMatMulCom(h_h, Ih3_IM, w, NT, NR, NR);

    // wyの作成
    MatVectMulCom(w, y, ans, NT, NR);

    Free_com2d(h_h, NT);
    Free_com2d(h3, NR);
    Free_com2d(Ih3, NR);
    Free_com2d(Ih3_IM, NR);
    Free_com2d(w, NT);
    free(wy);
    Free_com2d(I, NR);
}

void PrintComVect(complex *A, int value)
{
    int i, j;
    for (i = 0; i < value; i++)
    {
        if (A[i].im >= 0)
        {
            printf("z[%d] = %lf+%lfi  ", i, A[i].re, A[i].im);
        }
        else
        {
            printf("z[%d] = %lf %lfi  ", i, A[i].re, A[i].im);
        }
        printf("\n");
    }
}

void PrintComMat(complex **A, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            if (A[i][j].im >= 0)
            {
                printf("%lf+%lfi  ", A[i][j].re, A[i][j].im);
            }

            else
            {
                printf("%lf %lfi  ", A[i][j].re, A[i][j].im);
            }
        }
        printf("\n\n");
    }
    printf("\n");
}

void PrintComReMat(complex **A, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%lf  ", A[i][j].re);
        }
        printf("\n\n");
    }
    printf("\n");
}

void print_d_Vect(double *A, int value)
{
    int i;
    for (i = 0; i < value; i++)
    {
        printf("%f ", A[i]);
    }
    printf("\n");
}

void Print_d_Mat(double **A, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%lf  ", A[i][j]);
        }
        printf("\n\n");
    }
    printf("\n");
}

void Print_i_Mat(int **A, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%d  ", A[i][j]);
        }
        printf("\n\n");
    }
    printf("\n");
}

void printVect_i(int *A, int value)
{
    int i;
    for (i = 0; i < value; i++)
    {
        printf("%d ", A[i]);
    }
    printf("\n");
}

void excel_com_Vect(complex *z, int n, int a)
{ // aでファイル名を数字で変える
    int i, j;

    FILE *fp;
    char filename[50];
    sprintf(filename, "%d.csv", a);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    for (i = 0; i < n; i++)
    {

        fprintf(fp, "%f,%f\n", z[i].re, z[i].im);
    }
    fclose(fp);
}

void excel_com_Mat(complex **z, int row, int col, int a)
{ // aでファイル名を数字で変える
    int i, j;

    FILE *fp;
    char filename[50];
    sprintf(filename, "%d.csv", a);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    for (i = 0; i < col; i++)
    {
        for (j = 0; j < row; j++)
        {

            fprintf(fp, "%f,%f\n", z[j][i].re, z[j][i].im);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void excel_d_Vect3(double *x, int x_n, double *y, int y_n, double *z, int z_n, int a)
{ // aでファイル名を数字で変える
    int i, j;

    FILE *fp;
    char filename[50];
    sprintf(filename, "%d.csv", a);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    for (i = 0; i < x_n; i++)
    {

        fprintf(fp, "%f\n", x[i]);
    }

    for (i = 0; i < y_n; i++)
    {

        fprintf(fp, ",%f\n", y[i]);
    }

    for (i = 0; i < z_n; i++)
    {

        fprintf(fp, "%f\n", z[i]);
    }
    fclose(fp);
}

void excel_d_Vect(double *z, int n, int a)
{ // aでファイル名を数字で変える
    int i, j;

    FILE *fp;
    char filename[50];
    sprintf(filename, "%d.csv", a);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    for (i = 0; i < n; i++)
    {

        fprintf(fp, "%f\n", z[i]);
    }
    fclose(fp);
}

void ChangeCol_i(int *A, int x, int y)
{
    int tmp;

    tmp = A[x];
    A[x] = A[y];
    A[y] = tmp;
}

void ChangeCol_d(double *A, int x, int y)
{
    double tmp;

    tmp = A[x];
    A[x] = A[y];
    A[y] = tmp;
}

void ChangeCol_com(complex *A, int x, int y)
{
    complex tmp;

    tmp = A[x];
    A[x] = A[y];
    A[y] = tmp;
}

void ChangeCol_MatCom(complex **A, int A_row, int x, int y)
{
    int i, j;
    complex tmp;

    for (i = 0; i < A_row; i++)
    {
        tmp = A[i][x];
        A[i][x] = A[i][y];
        A[i][y] = tmp;
    }
}

void Cp_Mat(double **A, double **cp_A, int col, int row)
{
    int i, j;
    for (i = 0; i < col; i++)
    {
        for (j = 0; j < row; j++)
        {
            cp_A[i][j] = A[i][j];
        }
    }
}

void Cp_Mat_i(int **A, int **cp_A, int col, int row)
{
    int i, j;
    for (i = 0; i < col; i++)
    {
        for (j = 0; j < row; j++)
        {
            cp_A[i][j] = A[i][j];
        }
    }
}

void Cp_Vect(double *A, double *cp_A, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        cp_A[i] = A[i];
    }
}
