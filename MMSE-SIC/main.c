#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <unistd.h>
#include <string.h>
#include "../include/def.h"
#include "../include/RFunction.h"
#include "../include/const.h"
#include "../include/complexArith.h"
/*

__attribute__((destructor)) static void destructor()
{
    system("leaks -q main.app");
}
*/

int main(void)
{
    complex **H, *X, *Y, *N;
    complex **E, **W, **HH, **HHxWW, **nextH, **orgnH;
    complex *mmseX;
    complex *Z;
    complex zi;
    // complex **tmp;
    complex **J, **WW;
    FILE *fp;

    char fullpath[PATH_SIZE];
    char *dir, filename[50];
    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf(" %s (%d -> %d)\n", dir, NT, NR);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;

    srand((unsigned)time(NULL));
    // srand(1);
    int i, j, minIndex, loop, x, y;

    int colSize, cnt;
    int L = LOOP;

    for (en = ENMIN; en <= ENMAX; en++)
    {
        mmseX = calloc_com(NT);
        H = calloc_com2d(NT, NR);
        orgnH = calloc_com2d(NT, NR);
        X = calloc_com(NT);
        Y = calloc_com(NR);
        N = calloc_com(NR);
        Z = calloc_com(NT);

        Z = calloc_com(NT);
        W = calloc_com2d(NR, NT);
        J = calloc_com2d(NT, NT);
        nextH = calloc_com2d(NT - 1, NR);
        stv = STV(en); // 標準偏差
        error_sum = 0.0;

        for (loop = 0; loop < L; loop++)
        {

            X = genRandomInput();
            H = genChannelComplex();
            N = genNoiseComplex(stv);

            Y = mulComplexMatVec(H, X, NT, NR);
            Y = addComplexVec(Y, N);

            for (i = 0; i < NT; i++)
            {
                for (j = 0; j < NR; j++)
                {
                    orgnH[i][j] = H[i][j];
                }
            }

            // printAnyComplexMatrix(H, NT, NR);

            /*-------------------- MMSE-SIC : body------------------*/
            colSize = NT;
            for (cnt = NT; cnt > 0; cnt--)
            {
                Z = (complex *)realloc(Z, colSize);
                W = recalloc_com2d(W, NR, colSize);
                J = recalloc_com2d(J, colSize, colSize);
                W = calcWeightMatrix(H, colSize, NR, stv);
                Z = mulComplexMatVec(W, Y, NR, colSize);

                for (i = 0; i < colSize; i++)
                    Z[i] = Hard(Z[i]);

                //----------------- SIC ------------------
                J = calcJmatrix(W, H, colSize, NR);
                minIndex = minDiag(J, colSize);

                for (i = 0; i < NR; i++)
                    Y[i] = complex_sub(Y[i], complex_mul(H[minIndex][i], Z[minIndex]));

                for (i = 0; i < NT; i++)
                {
                    if (H[minIndex][0].re == orgnH[i][0].re)
                    {
                        mmseX[i] = Z[minIndex];
                    }
                }

                H = deleteCol(H, colSize, NR, minIndex);

                colSize--;
            }
            /* calc error */
            for (i = 0; i < NT; i++)
                (X[i].re != mmseX[i].re) ? error_sum++ : error_sum;
            for (i = 0; i < NT; i++)
                (X[i].im != mmseX[i].im) ? error_sum++ : error_sum;
        }
        //  printf("----> L=%d\n", L);

        ber = (double)error_sum / (double)(L * 2.0 * NT); // 実数と虚数があるので分母2倍

        printf("%.2f\t%.20f\n", en, ber);
        fprintf(fp, "%.2f\t%.20f\n", en, ber);

        free(Z);
        Free_com2d(nextH, NT);
        Free_com2d(W, NR);
        Free_com2d(J, NT);
        Free_com2d(H, NT);
        Free_com2d(orgnH, NT);
        free(mmseX);
        free(X);
        free(N);
        free(Y);
    }
    fclose(fp);
    PASS;
    return 0;
}
