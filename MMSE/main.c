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

#define PASS puts("\e[32m PASS \e[m\n")

int main(void)
{
    complex **H, *X, *Y, *N;
    complex **invH, *Nd, **E, *Xd, **W, **HH;
    complex **estH, *estX, *estY, **tmp1, **tmp2;
    FILE *fp;
    char filename[50];
    char fullpath[PATH_SIZE];
    char *dir;

    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf("%s\n", dir);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);
    printInfo();
    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    // fprintf(fp, "Eb/N0[db]\tBER\n");
    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
    int i, j, k, loop;

    H = calloc_com2d(NT, NR);
    X = calloc_com(NT);
    Y = calloc_com(NR);
    N = calloc_com(NR);
    Xd = calloc_com(NT);

    W = calloc_com2d(NR, NT);
    srand((unsigned)time(NULL));
    // srand(1);
    for (en = ENMIN; en <= ENMAX; en++)
    {
        stv = STV(en); // 標準偏差
        error_sum = 0.0;
        for (loop = 0; loop < LOOP; loop++)
        {
            X = genRandomInput();
            N = genNoiseComplex(stv);
            H = genChannelComplex();

            Y = mulComplexMatVec(H, X, NT, NR);
            Y = addComplexVec(Y, N);

            W = calcWeightMatrix(H, NT, NR, stv);
            Xd = mulComplexMatVec(W, Y, NR, NT);

            // printAnyComplexVec(Y, NR);
            //  printAnyComplexMatrix(H, NR, NT);

            for (int f = 0; f < NT; f++)
            {
                Xd[f] = Hard(Xd[f]);
            }

            for (int i = 0; i < NT; i++)
                (X[i].re != Xd[i].re) ? error_sum++ : error_sum;
            for (int i = 0; i < NT; i++)
                (X[i].im != Xd[i].im) ? error_sum++ : error_sum;
        }

        ber = error_sum / ((LOOP) * 2 * NT); // 実数と虚数があるので分母2倍
        printf("EN:%.2f  ber:%.8f\n", en, ber);
        fprintf(fp, "%.2f\t%.8f\n", en, ber);
    }
    free(X);
    free(H);
    free(N);
    free(Y);
    free(Xd);
    Free_com2d(W, NR);

    fclose(fp);
    return 0;
}
