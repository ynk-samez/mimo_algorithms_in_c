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
    complex **H, *X, *Y, *N, **XRep, **YRep;
    FILE *fp;
    char filename[50];
    char fullpath[PATH_SIZE];
    char *dir;

    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf("%s\n", dir);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);

    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

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

        /* NR x NT x NT x pow(4,NT)*/
        XRep = calloc_com2d(pow(4, NT), NT);
        YRep = calloc_com2d(pow(4, NT), NR);
        for (loop = 0; loop < LOOP; loop++)
        {
            X = genRandomInput();
            H = genChannelComplex();
            N = genNoiseComplex(stv);
            Y = mulComplexMatVec(H, X, NT, NR);
            Y = addComplexVec(Y, N);
        }

        for (int i = 0; i < NT; i++)
            (X[i].re != XRep[c][i].re) ? error_sum++ : error_sum;
        for (int i = 0; i < NT; i++)
            (X[i].im != XRep[c][i].im) ? error_sum++ : error_sum;
        free(X);
        Free_com2d(H, NT);
        free(N);
        free(Y);
        ber = error_sum / ((LOOP) * 2 * NT); // 実数と虚数があるので分母2倍
        printf("EN:%.2f  ber:%.8f\n", en, ber);
        fprintf(fp, "%.2f\t%.8f\n", en, ber);
    }
    fclose(fp);
    return 0;
}