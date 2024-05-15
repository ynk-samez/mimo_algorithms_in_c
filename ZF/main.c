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

    complex **H, *X, *Y, *N;
    complex **invH, *Xd;
    // complex **estH, *estX, *estY, **tmp1, **tmp2;
    FILE *fp;

    char filename[50];
    char fullpath[PATH_SIZE];
    char *dir;
    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf("%s\n", dir);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);
    printf("%s\n", filename);
    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }
    // fprintf(fp, "Eb/N0[db]\tBER\n");
    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
    int i, j, k, loop;

    if ((int)NT != (int)NR)
    {
        puts("error");
        pause();
    }
    int NN = NT;
    srand((unsigned)time(NULL));
    for (en = ENMIN; en <= ENMAX; en++)
    {

        H = calloc_com2d(NN, NN);
        X = calloc_com(NN);
        Y = calloc_com(NN);
        N = calloc_com(NN);
        Xd = calloc_com(NN);
        invH = calloc_com2d(NN, NN);

        stv = STV(en); // 標準偏差
        error_sum = 0.0;
        for (loop = 0; loop < LOOP; loop++)
        {
            /*
                        estH = calloc_com2d(NT, NR);
                        estX = calloc_com(NT);
                        estY = calloc_com(NR);

                        tmp1 = calloc_com2d(NT, NR);
                        tmp2 = calloc_com2d(NT, NR);
            */

            /* ------------- estimation start ------------- */

            /*
                        H = genChannelComplex(NN, NN); // generate H
                        for (int i = 0; i < S; i++)
                        {
                            estX = genRandomInput();
                            estY = addComplexVec(mulComplexMatVec(H, estX), genNoiseComplex(stv)); // Y=HX+N

                            // printAnyComplexVec(estY, NN);
                            tmp1 = complexMatrixAdd(tmp1, mulComplexVecVec(estY, conjugateVec(estX)));
                            tmp2 = complexMatrixAdd(tmp2, mulComplexVecVec(estX, conjugateVec(estX)));
                        }
                        estH = mulComplexMat(complexMatrixAve(tmp1), invMatrix(complexMatrixAve(tmp2)));
                        */

            /*-------------  estimation end ------------- */
            N = genNoiseComplex(stv);
            H = genChannelComplex();
            X = genRandomInput();
            Y = mulComplexMatVec(H, X, NN, NN);
            Y = addComplexVec(Y, N);

            invH = invMatrix(H, NN);
            Xd = mulComplexMatVec(invH, Y, NN, NN);

            for (int f = 0; f < NN; f++)
            {
                Xd[f] = Hard(Xd[f]);
            }

            for (int i = 0; i < NN; i++)
                (X[i].re != Xd[i].re) ? error_sum++ : error_sum;
            for (int i = 0; i < NN; i++)
                (X[i].im != Xd[i].im) ? error_sum++ : error_sum;
        }

        free(X);

        free(N);

        free(Y);
        Free_com2d(invH, NN);
        Free_com2d(H, NN);
        free(Xd);
        ber = error_sum / ((LOOP) * 2 * NN); // 実数と虚数があるので分母2倍
        printf("EN:%.2f  ber:%.8f\n", en, ber);
        // PASS;
        fprintf(fp, "%.2f\t%.8f\n", en, ber);
    }

    fclose(fp);
    return 0;
}
