#include "complexArith.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "def.h"
#include "RFunction.h"
#include "const.h"

void printComplex(complex a)
{
    printf("%f + j %f\n", a.re, a.im);
}

double norm(complex *vec, int h)
{
    double nm = 0;
    for (int i = 0; i < h; i++)
        nm += pow(vec[h].re, 2) + pow(vec[h].im, 2);

    return sqrt(nm);
}

complex **exchangeColmn(complex **A, int w, int h, int col_a, int col_b)
{
    complex *tmp;
    tmp = calloc_com(h);

    for (int j = 0; j < h; j++)
    {
        tmp[j] = A[col_a][j];
    }

    for (int j = 0; j < h; j++)
    {
        A[col_a][j] = A[col_b][j];
        A[col_b][j] = tmp[j];
    }
    free(tmp);
    return A;
}

int argmin(double *norms, int h)
{
    int i = 0;
    double _norms[NT];

    for (i = 0; i < NT; i++)
    {
        _norms[i] = norms[i];
    }

    qsort(_norms, NT, sizeof(double), cmp);

    for (i = 0; i < NT; i++)
    {
        if (_norms[0] == norms[i])
        {
            return i;
        }
    }
}

complex **makeZeroMatrix(int w, int h)
{
    complex **Z = calloc_com2d(w, h);
    return Z;
}

void printPotinter(char *error_pointer, complex **A)
{
    printf("%p\n", A);
}

void copyMatrix(complex **source, complex **destination, int w, int h)
{
    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < h; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

complex *complexVecSub(complex *A, complex *B, int Size)
{
    complex *Ans = calloc_com(Size);
    for (int cnt = 0; cnt < Size; cnt++)
        Ans[cnt] = complex_sub(A[cnt], B[cnt]);

    return Ans;
}

complex **calcWeightMatrix(complex **H, int w, int h, double stv)
{

    complex **Ans = calloc_com2d(h, w);
    complex **HH = calloc_com2d(h, w);
    complex **HxHH = calloc_com2d(h, h);
    complex **Inv = calloc_com2d(h, h);
    complex **s2E = calloc_com2d(h, h);

    HH = hermitianMatrix(H, w, h);

    HxHH = mulComplexMat(H, HH, w, h, h, w);
    s2E = scalarMat(genE(h), h, h, pow(stv, 2.0));

    Inv = invMatrix(addComplexMat(HxHH, s2E, h, h), h);

    Ans = mulComplexMat(HH, Inv, h, w, h, h);

    ShinFree2d(&HH, h);
    ShinFree2d(&HxHH, h);
    ShinFree2d(&Inv, h);
    ShinFree2d(&s2E, h);
    return Ans; // (colSize x NR ) x (NR x NR)
}

complex **hermitianMatrix(complex **Matrix, int w, int h)
{

    DEBUG(MODE)
    printf("hermitianMatrix : recieved [ %d x %d ]\n", w, h);
    complex **Conj;
    complex **Ans;

    Conj = calloc_com2d(w, h);
    Ans = calloc_com2d(h, w);

    Conj = conjugateMatrix(Matrix, w, h);
    // printf("hermitianMatrix: w=%d, h=%d の転置を依頼\n", w, h);
    Ans = transposeMatrix(Conj, w, h);

    ShinFree2d(&Conj, w);
    return Ans;
}

complex **deleteCol(complex **H, int w, int h, int index)
{
    int x, y;
    for (x = 0; x < w - 1; x++)
    {
        for (y = 0; y < h; y++)
        {
            if (x < index)
                H[x][y] = H[x][y];
            if (x >= index)
                H[x][y] = H[x + 1][y];
        }
    }
    H = recalloc_com2d(H, w - 1, h);

    return H;
}

complex *mulEach(complex *Vec, complex a)
{
    complex *Ans;

    Ans = calloc_com(NR);
    int y;
    for (y = 0; y < NR; y++)
    {
        Ans[y] = complex_mul(Vec[y], a);
    }

    return Ans;
}

complex *makeVecFromMat(complex **Mat, int h, int index)
{

    complex *Vec = calloc_com(h);
    for (int cnt = 0; cnt < h; cnt++)
    {
        //   printf("%d\n", cnt);
        Vec[cnt] = Mat[index][cnt];
    }
    return Vec;
}

complex **ShinMulComplexMat(complex **A, complex **B, int Size)
{
    complex **Ans;
    complex tmp, pre;
    Ans = calloc_com2d(NT, NR);

    int x, y, m;
    for (x = 0; x < NT; x++)
    {
        for (y = 0; y < NR; y++)
        {
            pre.re = 0;
            pre.im = 0;
            for (m = 0; m < Size; m++)
            {
                pre = complex_add(complex_mul(A[m][y], B[x][m]), pre);
            }
            // printf("( %f+j%f )\t", pre.re, pre.im);
            Ans[x][y] = pre;
        }
        // printf("\n");
    }
    return Ans;
}

complex *
ShinMulMatVec(complex **Mat, complex *Vec, int row, int col)
{
    complex *Ans = calloc_com(NR);
    MatVectMulCom(Mat, Vec, Ans, row, col);
    return Ans;
}

complex **updateChannelMatrix(complex **H, int minIndex)
{

    int i, j;
    for (i = minIndex; i < NT - 1; i++)
    {
        for (j = 0; j < NR; j++)
        {
            H[j][i] = H[j][i + 1];
        }
    }
    return H;
}

int cmp(const void *a, const void *b)
{
    double diff = *(double *)a - *(double *)b;
    if (diff < 0)
        return -1;
    else if (diff > 0)
        return 1;
    else
        return 0;
}

int minDiag(complex **A, int preColSize)
{
    double diagArray[preColSize];
    double originArray[preColSize];
    complex tmp = {0, 0};
    complex min = {DBL_MAX, DBL_MAX};

    for (int i = 0; i < preColSize; i++)
    {
        originArray[i] = A[i][i].re;
        diagArray[i] = A[i][i].re;
    }
    qsort(diagArray, preColSize, sizeof(double), cmp);
    for (int i = 0; i < preColSize; i++)
    {
        if (originArray[i] == diagArray[0])
            return i;
    }
    return -1;
}

complex **mulComplexVecVec(complex *A, complex *B, int sizeA, int sizeB)
{
    complex **Ans;
    Ans = calloc_com2d(sizeA, sizeB);
    int x, y;
    for (x = 0; x < sizeB; x++)
    {
        for (y = 0; y < sizeA; y++)
        {
            Ans[x][y] = complex_mul(A[y], B[x]);
            // printf("( %f+j%f )\t", pre.re, pre.im)
        }
        // printf("\n");
    }

    return Ans;
}

complex **complexMatrixAdd(complex **A, complex **B, int w, int h)
{
    complex **Ans;
    Ans = calloc_com2d(w, h);
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y] = complex_add(A[x][y], B[x][y]);
        }
    }
    return Ans;
}

complex **complexMatrixSub(complex **A, complex **B, int w, int h)
{
    complex **Ans;
    Ans = calloc_com2d(w, h);
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y] = complex_sub(A[x][y], B[x][y]);
        }
    }
    return Ans;
}

complex **complexMatrixAve(complex **A, int w, int h)
{

    complex **Ans;
    Ans = calloc_com2d(w, h);
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y].re = A[x][y].re / Pilots;
            Ans[x][y].im = A[x][y].im / Pilots;
        }
    }
    return Ans;
}

complex **genE(int Size)
{
    int x, y;
    complex **Ans;
    Ans = calloc_com2d(Size, Size);
    for (y = 0; y < Size; y++)
    {
        for (x = 0; x < Size; x++)
        {
            Ans[x][y].im = 0;
            Ans[x][y].re = (x == y) ? 1 : 0;
        }
    }
    return Ans;
}

complex **scalarMat(complex **A, int w, int h, double lamda)
{
    complex **Ans;
    Ans = calloc_com2d(w, h);
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y].re = A[x][y].re * lamda;
            Ans[x][y].im = A[x][y].im * lamda;
        }
    }

    return Ans;
}

complex **addComplexMat(complex **A, complex **B, int w, int h)
{
    complex **Ans;
    Ans = calloc_com2d(h, w);
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y] = complex_add(A[x][y], B[x][y]);
        }
    }

    return Ans;
}

complex *conjugateVec(complex *A, int Size)
{
    int y;
    complex *Ans;
    Ans = calloc_com(Size);
    for (y = 0; y < Size; y++)
    {
        Ans[y].re = A[y].re;
        Ans[y].im = -1 * A[y].im;
    }
    return Ans;
}

complex **conjugateMatrix(complex **A, int w, int h)
{
    /*共役　*/
    int x, y;
    complex **Ans;
    Ans = calloc_com2d(w, h);
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            Ans[x][y].re = A[x][y].re;
            Ans[x][y].im = -1 * A[x][y].im;
        }
    }
    return Ans;
}

complex **transposeMatrix(complex **A, int w, int h)
{
    complex **T;
    //  puts("転置関数: 定義成功");
    T = calloc_com2d(h, w);

    // puts("転置関数: メモリ確保成功");
    int x, y;
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            T[y][x] = A[x][y];
        }
    }

    // printAnyComplexMatrix(T, h, w);
    return T;
}

complex **invMatrix(complex **A, int Size)
{
    complex **Ans;
    Ans = calloc_com2d(Size, Size);
    InvComMat(A, Ans, Size);

    return Ans;
}

complex *genRandomInput()
{
    complex *Ans;
    Ans = calloc_com(NT);
    for (int i = 0; i < NT; i++)
    {
        Ans[i].re = (rand() % 2) * 2 - 1;
        Ans[i].im = (rand() % 2) * 2 - 1;
    }
    return Ans;
}

complex **genChannelComplex()
{
    complex **Ans;
    complex h = {0, 0};
    Ans = calloc_com2d(NT, NR);

    for (int y = 0; y < NR; y++)
    {
        for (int x = 0; x < NT; x++)
        {
            Ans[x][y] = wfading(h, 0);
        }
    }

    return Ans;
}

complex *genNoiseComplex(double stv)
{
    complex *Ans;
    Ans = calloc_com(NR);
    for (int i = 0; i < NR; i++)
    {
        Ans[i] = gaussianRand(stv, 0);
    }
    return Ans;
}

complex *addComplexVec(complex *A, complex *B)
{
    complex *Ans;
    Ans = calloc_com(NR);
    for (int i = 0; i < NR; i++)
    {
        Ans[i] = complex_add(A[i], B[i]);
    }
    return Ans;
}

complex *mulComplexMatVec(complex **A, complex *B, int mat_w, int mat_h)
{
    complex *Ans;
    complex pre;

    Ans = calloc_com(mat_h);
    int y, m;
    for (y = 0; y < mat_h; y++)
    {
        pre.re = 0;
        pre.im = 0;
        for (m = 0; m < mat_w; m++)
        {
            Ans[y] = complex_mul(A[m][y], B[m]);
            pre = complex_add(Ans[y], pre);
        }
        Ans[y] = pre;
    }
    return Ans;
}

/*Check*/
complex **calcJmatrix(complex **W, complex **H, int w, int h)
{
    DEBUG(MODE)
    printf("calcJmatrix : recieved [ %d x %d ]\n", w, h);
    complex **Ans;
    complex **HH, **WW, **E, **HHWW;
    Ans = calloc_com2d(w, w);
    WW = calloc_com2d(w, h);
    HH = calloc_com2d(h, w);
    E = calloc_com2d(w, w);

    HHWW = calloc_com2d(w, w);
    WW = hermitianMatrix(W, h, w);
    HH = hermitianMatrix(H, w, h);

    E = genE(w);

    HHWW = mulComplexMat(HH, WW, h, w, w, h);
    Ans = complexMatrixSub(E, HHWW, w, w);

    ShinFree2d(&WW, w);
    ShinFree2d(&HH, h);
    ShinFree2d(&E, w);
    ShinFree2d(&HHWW, w);

    return Ans;
}

complex **mulComplexMat(complex **A, complex **B, int input1_w, int input1_h, int input2_w, int input2_h)
{

    complex **Ans;
    complex tmp, pre;

    Ans = calloc_com2d(input2_w, input1_h);

    int x, y, m;
    for (x = 0; x < input2_w; x++)
    {
        for (y = 0; y < input1_h; y++)
        {
            pre.re = 0;
            pre.im = 0;
            for (m = 0; m < input1_w; m++)
            {
                Ans[x][y] = complex_mul(A[m][y], B[x][m]);
                pre = complex_add(Ans[x][y], pre);
            }
            Ans[x][y] = pre;
        }
    }

    return Ans;
}

void printAnyComplexVec(complex *A, int W)
{
    int y;
    printf("----------------- Vector : matrix (%d) -----------------------\n", W);
    for (y = 0; y < W; y++)
    {
        if (A[y].im < 0)
            printf("%+1.1f -j %.1f  ", A[y].re, fabs(A[y].im));
        else
            printf("%+1.1f +j %.1f  ", A[y].re, fabs(A[y].im));
        printf("\n");
    }
    puts("");
}

void printAnyComplexMatrix(complex **A, int W, int H)
{
    int x, y;
    printf("----------------- Matrix : complex (%d x %d) -----------------------\n", W, H);
    for (y = 0; y < H; y++)
    {
        for (x = 0; x < W; x++)
        {
            if (A[x][y].im < 0)
                printf("%+1.1f -j %.1f  ", A[x][y].re, fabs(A[x][y].im));
            else
                printf("%+1.1f +j %.1f  ", A[x][y].re, fabs(A[x][y].im));
        }
        printf("\n");
    }
    puts("");
}

void printAnyComplexMatrixColored(complex **A, int W, int H)
{
    int x, y;
    printf("----------------- Matrix : complex (%d x %d) -----------------------\n", W, H);
    for (y = 0; y < H; y++)
    {
        for (x = 0; x < W; x++)
        {
            if ((floor(A[x][y].re * 10.0) / 10.0 != 0.0L || floor(A[x][y].im * 10.0) / 10.0 != 0.0L) || (round(A[x][y].re) == 1.0))
                printf("\e[46m");

            if (A[x][y].im < 0)
                printf("%+1.1f -j %.1f  ", A[x][y].re, fabs(A[x][y].im));
            else
                printf("%+1.1f +j %.1f  ", A[x][y].re, fabs(A[x][y].im));

            printf("\e[m");
        }
        printf("\n");
    }
    puts("");
}

complex complex_add(complex x, complex y)
{
    complex z;

    z.re = x.re + y.re;
    z.im = x.im + y.im;

    return z;
}

complex complex_sub(complex x, complex y)
{
    complex z;
    z.re = x.re - y.re;
    z.im = x.im - y.im;

    return z;
}

complex complex_mul(complex x, complex y)
{
    complex z;
    z.re = x.re * y.re - x.im * y.im;
    z.im = x.re * y.im + x.im * y.re;
    return z;
}

complex complex_div(complex x, complex y)
{
    complex z = {0.0, 0.0};
    double d;

    d = y.re * y.re + y.im * y.im;

    z.re = (double)(x.re * y.re + x.im * y.im) / (double)d;
    z.im = (double)(-x.re * y.im + x.im * y.re) / (double)d;

    return z;
}

complex complex_add2(complex x, complex y, int *com_add_count)
{
    complex z;

    z.re = x.re + y.re;
    z.im = x.im + y.im;

    com_add_count[0] += 2;

    return z;
}

complex complex_sub2(complex x, complex y, int *com_add_count)
{
    complex z;
    z.re = x.re - y.re;
    z.im = x.im - y.im;
    com_add_count[0] += 2;
    return z;
}

complex complex_mul2(complex x, complex y, int *com_add_count, int *com_mul_count)
{
    complex z;
    z.re = x.re * y.re - x.im * y.im;
    z.im = x.re * y.im + x.im * y.re;

    com_add_count[0] += 2;
    com_mul_count[0] += 4;

    return z;
}

complex complex_div2(complex x, complex y, int *com_add_count, int *com_mul_count)
{
    complex z = {0.0, 0.0};
    double d;

    d = y.re * y.re + y.im * y.im;

    z.re = (x.re * y.re + x.im * y.im) / d;
    z.im = (-x.re * y.im + x.im * y.re) / d;

    com_add_count[0] += 3;
    com_mul_count[0] += 8;

    return z;
}