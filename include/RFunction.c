#include "RFunction.h"
#include "const.h"
#include "def.h"
#include <math.h>
#include <stdlib.h>
void EnseAve(complex **A, complex **ans, int row, int col, int n)
{
    int i, j;

    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            ans[i][j].re = A[i][j].re / n;
            ans[i][j].im = A[i][j].im / n;
        }
    }
}

void EstH(complex **y, complex **x_p, complex **ans, int symbol, int packl)
{
    int symbit;
    complex **x_h;      // X^H
    complex **xx_h;     // X*X^H
    complex **e_xx_h;   // X*X^H
    complex **inv_xx_h; // X*X^Hの逆行列
    complex **yx_h;     // Y*X^H
    complex **e_yx_h;   // Y*X^H

    x_h = calloc_com2d(symbol, NT);
    xx_h = calloc_com2d(NT, NT);
    e_xx_h = calloc_com2d(NT, NT);
    inv_xx_h = calloc_com2d(NT, NT);
    yx_h = calloc_com2d(NR, NT);
    e_yx_h = calloc_com2d(NR, NT);

    // X^Hの作成
    MatTrans(x_p, x_h, NT, symbol);

    // Y*X^H
    MatMatMulCom(y, x_h, yx_h, NR, symbol, NT);

    // Y*X^Hのアンサンブルアベレージ
    EnseAve(yx_h, e_yx_h, NR, NT, symbol);

    // X*X^H
    MatMatMulCom(x_p, x_h, xx_h, NT, symbol, NT);

    // X*X^Hのアンサンブルアベレージ計算
    EnseAve(xx_h, e_xx_h, NT, NT, symbol);

    // X*X^Hのアンサンブルアベレージの逆行列
    InvComMat(e_xx_h, inv_xx_h, NT);

    // Hの推定
    MatMatMulCom(e_yx_h, inv_xx_h, ans, NR, NT, NT);

    Free_com2d(x_h, symbol);
    Free_com2d(xx_h, NT);
    Free_com2d(e_xx_h, NT);
    Free_com2d(inv_xx_h, NT);
    Free_com2d(yx_h, NR);
    Free_com2d(e_yx_h, NR);
}

double VectEuclideanDis(complex *v1, complex *v2, int dim)
{

    complex *ans;
    ans = calloc_com(dim);
    double c[dim], phi;
    phi = 0;
    int i;

    VectSubCom(v1, v2, ans, dim);
    for (i = 0; i < dim; i++)
    {
        c[i] = sqrt(pow((ans[i].re), 2.0) + pow((ans[i].im), 2.0));
        phi += c[i];
    }
    free(ans);
    return phi;
}

// 送信信号レプリカの作成
void MakeRep(complex **X_rep, int N_T)
{

    int i, j, k, tmp;
    complex QPSK_pt[4];

    QPSK_pt[0].re = 1.0;
    QPSK_pt[0].im = 1.0;

    QPSK_pt[1].re = 1.0;
    QPSK_pt[1].im = -1.0;

    QPSK_pt[2].re = -1.0;
    QPSK_pt[2].im = -1.0;

    QPSK_pt[3].re = -1.0;
    QPSK_pt[3].im = 1.0;

    for (i = 0; i < pow(4, N_T); i++)
    {
        for (j = 0; j < N_T; j++)
        {
            tmp = i / pow(4, j);
            X_rep[i][j] = QPSK_pt[tmp % 4];
        }
    }
}

void MLD(complex **x_rep, complex **h, complex *y, complex *ans, int N_T, int N_R)
{
    int i, j;
    double phi;
    double d;
    d = 10000;
    complex **y_rep = calloc_com2d(pow(4, N_T), N_R); // yのレプリカ

    for (i = 0; i < pow(4, N_T); i++)
        MatVectMulCom(h, x_rep[i], y_rep[i], N_R, N_T); // yのレプリカ作成

    for (i = 0; i < pow(4, N_T); i++)
    {
        phi = VectEuclideanDis(y, y_rep[i], N_R);
        if (d > phi)
        {
            d = phi;
            j = i;
        }
    }

    for (i = 0; i < N_T; i++)
    {
        ans[i] = x_rep[j][i];
    }
    Free_com2d(y_rep, pow(4, N_T));
}

void MMSE_est(complex **h, complex **y, complex **ans, double stv2, int packl)
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
    MatMatMulCom(w, y, ans, NT, NR, packl);

    Free_com2d(h_h, NT);
    Free_com2d(h3, NR);
    Free_com2d(Ih3, NR);
    Free_com2d(Ih3_IM, NR);
    Free_com2d(w, NT);
    free(wy);
    Free_com2d(I, NR);
}

void ZF(complex **h, complex *y, complex *ans)
{
    complex **inv_h;
    inv_h = calloc_com2d(NR, NT);

    // 逆行列作成
    InvComMat(h, inv_h, NR);
    // 復調
    MatVectMulCom(inv_h, y, ans, NR, NT);
    Free_com2d(inv_h, NR);
}

void ZF_est(complex **h, complex **y, complex **ans, int packl)
{
    complex **inv_h;
    inv_h = calloc_com2d(NR, NT);

    // 逆行列作成
    InvComMat(h, inv_h, NR);
    // 復調
    MatMatMulCom(inv_h, y, ans, NR, NT, packl);
    Free_com2d(inv_h, NR);
}

void MMSE_w(complex **h, complex **w, double stv2, int nt, int nr)
{
    complex **h_h; // hの逆行列 hの共役転置行列
    complex **h3, **Ih3, **Ih3_IM;
    complex **I; // σ^2と単位行列の積
    int i;

    h_h = calloc_com2d(NT, NR);    // hの共役転置行列
    h3 = calloc_com2d(NR, NR);     // hとhの共役転置行列の積
    Ih3 = calloc_com2d(NR, NR);    // h3とσ^2Iの和
    Ih3_IM = calloc_com2d(NR, NR); // 逆行列
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

    Free_com2d(h_h, NT);
    Free_com2d(h3, NR);
    Free_com2d(Ih3, NR);
    Free_com2d(Ih3_IM, NR);
    Free_com2d(I, NR);
}
void MMSE_SIC(complex **h, complex *y, complex *ans, double stv2)
{
    complex **h_h; // hの逆行列 hの共役転置行列
    complex **h3, **Ih3, **Ih3_IM, **w, *wy, **trans_w, **J, **HHWW, k1, *t, *hz, *Z, *Zbar;
    complex **I;        // σ^2と単位行列の積と単位行列
    complex **h1, **h2; // hの複製
    complex *Y;         // yの複製
    int i, j, nt, minIndex;
    double *t_real;
    double *m1;
    double tmp;

    Z = calloc_com(NT);
    Zbar = calloc_com(NT);      // 復調信号
    h_h = calloc_com2d(NT, NR); // hの共役転置行列
    w = calloc_com2d(NT, NR);
    trans_w = calloc_com2d(NR, NT);
    wy = calloc_com(NT);
    I = calloc_com2d(NT, NT);
    J = calloc_com2d(NT, NT);
    HHWW = calloc_com2d(NT, NT);
    t = calloc_com(NT);
    hz = calloc_com(NR);
    h1 = calloc_com2d(NR, NT);
    Y = calloc_com(NR);
    t_real = calloc_d(NT);
    m1 = calloc_d(NT);

    // 引数hの複製

    for (i = 0; i < NR; i++)
    {
        for (j = 0; j < NT; j++)
        {
            h1[i][j] = h[i][j];
        }
    }

    // 引数yの複製
    for (i = 0; i < NR; i++)
    {
        Y[i] = y[i];
    }

    for (nt = NT; nt > 0; nt = nt - 1)
    {
        // wの作成
        MMSE_w(h1, w, stv2, nt, NR);

        MatVectMulCom(w, Y, wy, nt, NR); // wyの作成

        for (i = 0; i < nt; i++)
        {
            Z[i] = Hard(wy[i]); // 硬判定
        }

        // Jの作成
        IdentityComMat(I, nt);                        // 単位行列の作成
        MatTrans(w, trans_w, nt, NR);                 // wの共役転置行列
        MatTrans(h1, h_h, NR, nt);                    // Hの共役転置行列
        MatMatMulCom(h_h, trans_w, HHWW, nt, NR, nt); // HHWW=H_H*W_H
        MatMatSubCom(I, HHWW, J, nt, nt);             // J=I-H_H*W_H

        // Jの体格成分の最小値を求める
        for (i = 0; i < nt; i++)
        {
            t_real[i] = J[i][i].re;
        }
        // 最小値格納。。。関数にしたかったけどメモリの関係でなぜかできなかった。。。

        for (i = 0; i < nt; i++)
            m1[i] = t_real[i];

        for (i = 0; i < nt; i++)
        {
            for (j = i + 1; j < nt; j++)
            {
                if (m1[i] < m1[j])
                {
                    tmp = m1[i];
                    m1[i] = m1[j];
                    m1[j] = tmp;
                }
            }
        }
        for (i = 0; i < nt; i++)
        {
            if (t_real[i] == m1[nt - 1])
            {
                minIndex = i;
                break;
            }
        }
        // minIndex = Array_SmallNumber(t_real,nt);

        // Yの更新
        for (i = 0; i < NR; i++)
        {
            hz[i] = complex_mul(h1[i][minIndex], Z[minIndex]);
            Y[i] = complex_sub(Y[i], hz[i]);
        }

        // Z_にZの保存
        for (i = 0; i < NT; i++)
        {
            if (h1[0][minIndex].re == h[0][i].re)
            {
                Zbar[i] = Z[minIndex];
            }
        }

        // hの更新
        for (i = minIndex; i < NT - 1; i++)
        {
            for (j = 0; j < NR; j++)
            {
                h1[j][i] = h1[j][i + 1];
            }
        }
    }

    for (i = 0; i < NT; i++)
    {
        ans[i] = Zbar[i];
    }

    free(Z);
    free(Zbar);
    Free_com2d(h_h, NT);
    Free_com2d(w, NT);
    Free_com2d(trans_w, NR);
    free(wy);
    Free_com2d(I, NT);
    Free_com2d(J, NT);
    Free_com2d(HHWW, NT);
    free(t);
    free(hz);
    Free_com2d(h1, NR);
    free(t_real);
    free(Y);
    free(m1);
}

void MMSE_SQRD(complex **stvI, complex **h, complex *y, complex *z)
{
    int i, j, k, tmp, *p;
    complex **q_hq, **qr, *yr, **hr, **r, **q, **qt, *s, tmp2, tmp3;

    yr = calloc_com(NR + NT);
    hr = calloc_com2d(NR + NT, NT);
    p = calloc_i(NT);
    r = calloc_com2d(NT, NT);
    q_hq = calloc_com2d(NT, NT);
    q = calloc_com2d(NR + NT, NT);
    qr = calloc_com2d(NR + NT, NT);
    qt = calloc_com2d(NT, NR + NT);
    s = calloc_com(NT);

    for (i = 0; i < NR; i++)
    {
        yr[i] = y[i];
    }
    for (i = NR; i < NR + NT; i++)
    {
        yr[i].re = 0;
        yr[i].im = 0;
    }

    for (i = 0; i < NR; i++)
    {
        for (j = 0; j < NT; j++)
        {
            hr[i][j] = h[i][j];
        }
    }
    /*
      for(i = NR;i < NR+NT;i++){
        for(j = 0;j < NT;j++){
          if(i-NR == j){
            hr[i][j] = stvI[j][j];
          }
        }
      }
    */

    for (i = NR; i < NR + NT; i++)
    {
        for (j = 0; j < NT; j++)
        {
            for (k = 0; k < NT; k++)
            {
                if (i == k + NR)
                {
                    hr[i][j] = stvI[k][j];
                }
            }
        }
    }

    for (i = 0; i < NR + NT; i++)
    {
        for (j = 0; j < NT; j++)
        {
            q[i][j] = hr[i][j];
        }
    }

    for (i = 0; i < NT; i++)
        p[i] = i;
    /*
      InitComMat(r, NT,NT);
    */

    for (i = 0; i < NT; i++)
    {
        for (j = 0; j < NT; j++)
        {
            r[i][j].re = 0;
            r[i][j].im = 0;
        }
    }
    SQRD(p, q, r);
    MatTrans(q, qt, NR + NT, NT);

    /*
      printf("q\n");
      PrintComMat(q, NR+NT, NT);

      printf("q_h\n");
      PrintComMat(qt, NT, NR+NT);
    */

    MatMatMulCom(qt, q, q_hq, NT, NR + NT, NT);
    /*
     printf("Q_H*Q\n");
     PrintComMat(q_hq, NT, NT);

     printf("p\n");
     for(i = 0;i<NT;i++){
       printf("%d",p[i]);
     }
     printf("\n\nr\n");
     PrintComMat(r, NT, NT);
   */

    MatVectMulCom(qt, yr, s, NT, NR + NT);

    for (i = NT - 1; i >= 0; i--)
    {
        z[i] = complex_div(s[i], r[i][i]);
        z[i] = Hard(z[i]);
        for (j = 0; j < NT; j++)
        {
            s[j] = complex_sub(s[j], complex_mul(r[j][i], z[i]));
        }
    }

    // 復調信号zの整列
    for (i = 0; i < NT; i++)
    {
        for (j = i + 1; j < NT; j++)
        {

            if (p[i] > p[j])
            {

                tmp = p[i];
                p[i] = p[j];
                p[j] = tmp;

                tmp2 = z[i];
                z[i] = z[j];
                z[j] = tmp2;

                for (k = 0; k < NT; k++)
                {
                    tmp3 = r[k][i];
                    r[k][i] = r[k][j];
                    r[k][j] = tmp3;
                }
            }
        }
    }

    /*
      printf("r(After Sorting)\n");
      PrintComMat(r, NT, NT);
      MatMatMulCom(q, r, qr, NR+NT, NT, NT);
      printf("h\n");
      PrintComMat(hr, NR+NT, NT);
      printf("qr\n");
      PrintComMat(qr, NR+NT, NT);
    */

    free(yr);
    Free_com2d(hr, NR + NT);
    free(p);
    Free_com2d(r, NT);
    Free_com2d(q_hq, NT);
    Free_com2d(q, NR + NT);
    Free_com2d(qr, NR + NT);
    Free_com2d(qt, NT);
    free(s);
}

void SQRD(int *p, complex **q, complex **r)
{
    int i, j, k;
    int minnum;
    double *norm, min;
    complex *qi, **qt, *qt_q, *qk;

    norm = calloc_d(NT);
    qi = calloc_com(NR + NT);
    qt = calloc_com2d(1, NR + NT);
    qt_q = calloc_com(1);
    qk = calloc_com(NR + NT);

    for (i = 0; i < NT; i++)
    {
        norm[i] = 0;
        for (j = 0; j < NR + NT; j++)
            norm[i] += pow(scalar(q[j][i]), 2.0);
    }

    for (i = 0; i < NT; i++)
    {
        min = norm[i];
        minnum = i;

        for (j = i + 1; j < NT; j++)
        {
            if (min > norm[j])
            {
                min = norm[j];
                minnum = j;
            }
        }

        ChangeCol_i(p, i, minnum);
        ChangeCol_d(norm, i, minnum);
        ChangeCol_MatCom(q, NR + NT, i, minnum);
        ChangeCol_MatCom(r, NT, i, minnum);
        r[i][i].re = sqrt(norm[i]);

        for (j = 0; j < NR + NT; j++)
            q[j][i] = complex_div(q[j][i], r[i][i]);

        for (j = i + 1; j < NT; j++)
        {
            for (k = 0; k < NR + NT; k++)
                qi[k] = q[k][i];

            VectTrans(qi, qt, NR + NT);

            for (k = 0; k < NR + NT; k++)
                qk[k] = q[k][j];

            MatVectMulCom(qt, qk, qt_q, 1, NR + NT);

            r[i][j] = qt_q[0];

            for (k = 0; k < NR + NT; k++)
            {
                qk[k] = complex_sub(qk[k], complex_mul(r[i][j], q[k][i]));
                q[k][j] = qk[k];
            }
            norm[j] = norm[j] - pow(scalar(r[i][j]), 2.0);
        }
    }

    Free_com2d(qt, 1);
    free(qi);
    free(qt_q);
    free(qk);
    free(norm);
}

/*↓ビタビ復号用の関数*/
/*拘束長可変，符号化率(出力は2bit)固定？*/
/*トレリス線図を並び替えて，input=0のときは上半分の状態，じゃないときは下半分の状態に遷移するようにしている*/

void QPSK(int *x, complex *ans, int Trabit)
{
    int i;
    for (i = 0; i < Trabit; i++)
    {
        if (x[i * 2] == 0)
            ans[i].re = -1.0;
        else if (x[i * 2] == 1)
            ans[i].re = 1.0;

        if (x[i * 2 + 1] == 0)
            ans[i].im = -1.0;
        else
            ans[i].im = 1.0;
    }
}

void DeQPSK(int *x, complex *ans, int Trabit)
{
    int i;
    for (i = 0; i < Trabit; i++)
    {
        if (ans[i].re == -1.0)
            x[i * 2] = 0;
        else if (ans[i].re == 1.0)
            x[i * 2] = 1;
        if (ans[i].im == -1.0)
            x[i * 2 + 1] = 0;
        else
            x[i * 2 + 1] = 1;
    }
}

/*畳み込み符号回路生成*/
/*結合関数がなんたら, main関数内で，8進数のst1,st2をグローバル変数とかで定義する*/
/*g[outputnum][cl], g[何番目の出力ビットか][何番目の内部レジスタか] = 0 or 1(線が繋がっていれば1)*/
void convon_init_s(int BITRA, int KLEN, int **g, int **t, int STATE, int St1, int St2)
{
    int i, j, x;

    for (i = 0; i < BITRA; i++)
    {
        for (j = 0; j < KLEN; j++)
        {
            g[i][j] = 0;
        }
    }

    /*状態生成*/
    for (i = 0; i < STATE; i++)
    {
        x = i;
        for (j = KLEN - 2; j > -1; j--)
        {
            t[i][j] = x % 2;
            x = x / 2;
        }
    }

    /*畳み込み符号回路生成*/
    for (i = 0; i < BITRA; i++)
    {
        for (j = 0; j < KLEN; j++)
        {
            if (i == 0)
            {
                g[i][j] = St1 % 2;
                St1 = St1 / 2;
            }

            else
            {
                g[i][j] = St2 % 2;
                St2 = St2 / 2;
            }
        }
    }
}

void convon_code(int n, int cl, int **g, int *bin, int *result, int Trabit)
{
    int i, j, k;
    int **x_i;
    int *m;

    x_i = calloc_i2d(Trabit, n);
    m = calloc_i(cl);

    for (i = 0; i < Trabit; i++)
    {
        for (j = 0; j < n; j++)
        {
            x_i[i][j] = 0;
        }
    }

    for (i = 0; i < cl; i++)
    {
        m[i] = 0;
    }

    for (i = 0; i < Trabit * 2; i++)
    {
        result[i] = 0;
    }

    /* 符号化 */
    for (i = 0; i < Trabit; i++)
    {
        /* 入力ビット */
        m[0] = bin[i];
        /* 排他的論理和 */
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < cl; k++)
            {
                if (g[j][k] == 1)
                    x_i[i][j] ^= m[k];
            }
        }
        /* シフト */
        for (k = cl - 1; k >= 1; k--)
            m[k] = m[k - 1];
    }

    for (i = 0; i < Trabit; i++)
    {
        result[2 * i] = x_i[i][0];
        result[2 * i + 1] = x_i[i][1];
    }
    /*
    // 結果出力
      printf("\n\n| The Code Vectors Are: |\n");
      for (i = 0; i < Trabit; i++) {
        printf("| \n| ");

        for (j = 0; j < n; j++)
          printf("%d ",x_i[i][j]);
      }

      printf("\n\n");
      for(i = 0;i < Trabit*2;i++){
        printf("%d ",result[i]);
       }
         printf("\n");
    */
    Free_i2d(x_i, Trabit);
    free(m);
}

/*ビタビ復号*/
void Viterbi(complex *y, int **t, double **eucdis, double *dis_buf, int *vi_result, int **route, int LEN, int KLEN, int STATE)
{
    int i, j;
    double **eucdis_buf = calloc_d2d(STATE, 2);
    int **out = calloc_i2d(2 * STATE, 2);

    dis_buf[0] = 0;
    for (i = 1; i < STATE; i++)
    {
        dis_buf[i] = 100;
    }

    /*ブランチメトリック計算*/
    for (i = 0; i < LEN; i++)
    {
        for (j = 0; j < STATE; j++)
        {
            ViterbiDecoding(y[i], t[j], out[2 * j], out[2 * j + 1], eucdis[j], KLEN);
        }

        /*パスメトリック計算*/
        for (j = 0; j < STATE; j++)
        {
            eucdis_buf[j][0] = eucdis[j][0] + dis_buf[j];
            eucdis_buf[j][1] = eucdis[j][1] + dis_buf[j];
        }
        /*パスメトリックの比較*/
        PathComparison(eucdis_buf, dis_buf, route[i], STATE);
    }

    /*経路で復号*/
    FollowRoute(route, vi_result, LEN, STATE);
    free(eucdis_buf);
    Free_i2d(out, STATE);
}

/**/
void ViterbiDecoding(complex y, int *t, int *out0, int *out1, double *eucdis, int KLEN)
{
    int i;
    complex com0, com1, dif;

    Convolution(t, out0, out1, KLEN);

    for (i = 0; i < 2; i++)
    {
        if (out0[i] == 0)
            out0[i] = -1;
    }
    com0.re = out0[0];
    com0.im = out0[1];
    for (i = 0; i < 2; i++)
    {
        if (out1[i] == 0)
            out1[i] = -1;
    }
    com1.re = out1[0];
    com1.im = out1[1];

    dif = complex_sub(y, com0);
    eucdis[0] = pow(dif.re, 2.0) + pow(dif.im, 2.0);
    dif = complex_sub(y, com1);
    eucdis[1] = pow(dif.re, 2.0) + pow(dif.im, 2.0);
}

/*畳み込み計算*/
void Convolution(int *t, int *out0, int *out1, int KLEN)
{
    int i;
    switch (KLEN)
    {
    case 3:
        out0[0] = 0 ^ t[1];
        out0[1] = 0 ^ t[0] ^ t[1];
        out1[0] = 1 ^ t[1];
        out1[1] = 1 ^ t[0] ^ t[1];
        break;
    case 4:
        out0[0] = 0 ^ t[0] ^ t[2];
        out0[1] = 0 ^ t[0] ^ t[1] ^ t[2];
        out1[0] = 1 ^ t[0] ^ t[2];
        out1[1] = 1 ^ t[0] ^ t[1] ^ t[2];
        break;
    case 5:
        out0[0] = 0 ^ t[2] ^ t[3];
        out0[1] = 0 ^ t[0] ^ t[1] ^ t[3];
        out1[0] = 1 ^ t[2] ^ t[3];
        out1[1] = 1 ^ t[0] ^ t[1] ^ t[3];
        break;
    case 6:
        out0[0] = 0 ^ t[1] ^ t[3] ^ t[4];
        out0[1] = 0 ^ t[0] ^ t[1] ^ t[2] ^ t[4];
        out1[0] = 1 ^ t[1] ^ t[3] ^ t[4];
        out1[1] = 1 ^ t[0] ^ t[1] ^ t[2] ^ t[4];
        break;
    case 7:
        out0[0] = 0 ^ t[0] ^ t[1] ^ t[2] ^ t[5];
        out0[1] = 0 ^ t[1] ^ t[2] ^ t[4] ^ t[5];
        out1[0] = 1 ^ t[0] ^ t[1] ^ t[2] ^ t[5];
        out1[1] = 1 ^ t[1] ^ t[2] ^ t[4] ^ t[5];
        break;
    }
}

/*パスメトリック比較*/
void PathComparison(double **eucdis_buf, double *sum, int *route, int STATE)
{
    int i, j;
    for (i = 0; i < STATE / 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            /*小さいほうを*/
            if (eucdis_buf[2 * i][j] < eucdis_buf[2 * i + 1][j])
            {
                sum[i + STATE / 2 * j] = eucdis_buf[2 * i][j];
                route[i + STATE / 2 * j] = 2 * i;
            }
            else if (eucdis_buf[2 * i + 1][j] <= eucdis_buf[2 * i][j])
            {
                sum[i + STATE / 2 * j] = eucdis_buf[2 * i + 1][j];
                route[i + STATE / 2 * j] = 2 * i + 1;
            }
        }
    }
}

/**/
void FollowRoute(int **route, int *result, int LEN, int STATE)
{
    int i;
    result[LEN - 1] = 0;
    for (i = 1; i < LEN; i++)
    {
        result[LEN - i - 1] = route[LEN - i][result[LEN - i]];
    }

    /*トレリス線図を並び替えて，input=0のときは上半分の状態，じゃないときは下半分の状態に遷移するようにしている*/
    for (i = 0; i < LEN; i++)
    {
        if (result[i] < STATE / 2)
            result[i] = 0;
        else
            result[i] = 1;
    }
}