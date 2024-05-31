#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#include "../include/def.h"
#include "../include/RFunction.h"
#include "../include/const.h"
#include "../include/complexArith.h"

int main(void)
{
    int x, y, *p;
    complex **H, *X, *Y, *N;
    complex *Yb, **Hb, **Q, **R, **E, *S, *Z, **QH;
    complex *qi, *qk, *qiTrans, *qtq;
    double *norms;
    complex rik;

    FILE *fp;
    char filename[50], fullpath[PATH_SIZE], *dir;

    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf("%s\n", dir);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    } // ok

    printInfo();

    // fprintf(fp, "Eb/N0[db]\tBER\n");
    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
    int i, j, k, loop;

    // srand((unsigned)time(NULL));
    srand(1);
    for (en = ENMIN; en <= ENMAX; en++)
    {
        stv = STV(en); // 標準偏差
        error_sum = 0.0;
        H = calloc_com2d(NT, NR);
        X = calloc_com(NT);
        Y = calloc_com(NR);
        N = calloc_com(NR);

        norms = (double *)calloc(NT, sizeof(double));
        p = (int *)calloc(NT, sizeof(int));

        Hb = calloc_com2d(NT, NR + NT);
        R = calloc_com2d(NT, NT);
        E = calloc_com2d(NT, NT);
        Yb = calloc_com(NR + NT);
        Q = calloc_com2d(NT, NR + NT);
        QH = calloc_com2d(NR + NT, NT);
        S = calloc_com(NT);
        Z = calloc_com(NT);

        qi = calloc_com(NT + NR);
        qk = calloc_com(NT + NR);
        qiTrans = calloc_com(NT + NR);
        qtq = calloc_com(1);
        // qi = calloc_com(NR + NT);
        // qk = calloc_com(NR + NT);
        E = genE(NT);

        for (y = 0; y < NT; y++)
            p[y] = y;

        int L = LOOP;

        for (loop = 0; loop < L; loop++)
        {
            X = genRandomInput();
            H = genChannelComplex();
            N = genNoiseComplex(stv);
            Y = mulComplexMatVec(H, X, NT, NR);
            Y = addComplexVec(Y, N);

            for (x = 0; x < NT; x++)
            {
                for (y = 0; y < NR + NT; y++)
                {
                    if (y < NR)
                    {
                        Hb[x][y] = H[x][y];
                    }
                    else
                    {
                        Hb[x][y].re = stv * E[x][y - NR].re;
                        Hb[x][y].im = stv * E[x][y - NR].im;
                    }
                }
            }

            Yb = mulComplexMatVec(Hb, X, NT, NR + NT);
            R = makeZeroMatrix(NT, NT);
            copyMatrix(Hb, Q, NT, NR + NT);

            // printAnyComplexMatrix(Q, NT, NR + NT);
            //  SQRD:Sorted QR Decomposition
            for (i = 0; i < NT; i++)
            {
                norms[i] = 0;
                for (j = 0; j < NT + NR; j++)
                {
                    norms[i] += pow(Q[i][j].re, 2.0) + pow(Q[i][j].im, 2.0);
                }
            }

            for (i = 0; i < NT; i++)
            {
                int ki = i;
                double min = norms[i];

                for (j = i + 1; j < NT; j++)
                {
                    if (min > norms[j])
                    {
                        min = norms[j];
                        ki = j;
                    }
                }

                //-------------  EXCHANGE  ----------------

                R = exchangeColmn(R, NT, NT, i, ki);

                int itmp = p[i];
                p[i] = p[ki];
                p[ki] = itmp;

                double dtmp = norms[i];
                norms[i] = norms[ki];
                norms[ki] = dtmp;

                for (j = 0; j < NR + i; j++)
                {
                    complex ctmp = Q[i][j];
                    Q[i][j] = Q[ki][j];
                    Q[ki][j] = ctmp;
                }

                //-----------------------------------
                R[i][i].re = sqrt(norms[i]);
                R[i][i].im = 0.0;

                qi = makeVecFromMat(Q, (NT + NR), i);
                for (j = 0; j < NT + NR; j++)
                {
                    qi[j] = complex_div(qi[j], R[i][i]);
                    Q[i][j] = qi[j];
                }

                for (k = i + 1; k < NT; k++)
                {
                    qk = makeVecFromMat(Q, (NT + NR), k);
                    for (j = 0; j < NT + NR; j++)
                    {
                        qiTrans[j].im = -1.0 * qi[j].im;
                        qiTrans[j].re = qi[j].re;
                    }

                    rik.re = 0.0;
                    rik.im = 0.0;
                    for (j = 0; j < NT + NR; j++)
                    {
                        rik = complex_add(rik, complex_mul(qiTrans[j], qk[j]));
                    }

                    R[k][i] = rik; // do NOT change this

                    for (j = 0; j < NT + NR; j++)
                    {
                        qk[j] = complex_sub(qk[j], complex_mul(rik, qi[j]));
                        Q[k][j] = qk[j];
                    }

                    norms[k] = norms[k] - pow(scalar(rik), 2.0);
                }
            }
            QH = hermitianMatrix(Q, NT, NT + NR);
            // 上三角行列になってるか
            printAnyComplexMatrixColored(R, NT, NT);
            //   単位行列になってるか
            // printAnyComplexMatrixColored(mulComplexMat(QH, Q, NT + NR, NT, NT, NT + NR), NT, NT);

            /*---------- 復調 -----------*/

            S = mulComplexMatVec(QH, Yb, NR + NT, NT);
            for (i = NT - 1; i >= 0; i--)
            {
                Z[i] = complex_div(S[i], R[i][i]);
                for (j = 0; j < NT; j++)
                {
                    S[j] = complex_sub(S[j], complex_mul(R[i][j], Z[i]));
                }
            }

            //-----------------並び替え:わからなかったので借りた-------------
            for (i = 0; i < NT; i++)
            {
                for (j = i + 1; j < NT; j++)
                {
                    if (p[i] > p[j])
                    {
                        int tmp = p[i];
                        p[i] = p[j];
                        p[j] = tmp;

                        complex ctmp = Z[i];
                        Z[i] = Z[j];
                        Z[j] = ctmp;

                        for (k = 0; k < NT; k++)
                        {
                            complex ctmp2 = R[i][k];
                            R[i][k] = R[j][k];
                            R[j][k] = ctmp2;
                        }
                    }
                }
            }

            //--------------硬判定------------
            for (int i = 0; i < NT; i++)
                Z[i] = Hard(Z[i]);

            for (int i = 0; i < NT; i++)
                (X[i].re != Z[p[i]].re) ? error_sum++ : error_sum;
            for (int i = 0; i < NT; i++)
                (X[i].im != Z[p[i]].im) ? error_sum++ : error_sum;
        }

        free(X);
        Free_com2d(Hb, NT);
        Free_com2d(H, NT);
        Free_com2d(E, NT);
        Free_com2d(Q, NT);
        Free_com2d(QH, NR + NT);
        Free_com2d(R, NT);

        free(S);
        free(N);
        free(Y);
        free(Yb);
        free(p);
        free(norms);
        free(Z);

        free(qi);
        free(qiTrans);
        free(qk);

        free(qtq);

        ber = error_sum / ((L) * 2 * NT);
        printf("EN:%.2f  ber:%.18f\n", en, ber);
        fprintf(fp, "%.2f\t%.18f\n", en, ber);
    }
    fclose(fp);
    return 0;
}