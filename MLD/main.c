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
    printf("%s\n", filename);
    if ((fp = fopen(filename, "w+")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

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

            MakeRep(XRep, NT);
            YRep = mulComplexMat(H, XRep, NT, NR, pow(4, NT), NT);
            // printAnyComplexMatrix(YRep, pow(4, NT), NR);

            /*----------- MLD ------------*/
            double pre = 9999,
                   phi;
            int c = 0;

            for (int cnt = 0; cnt < pow(4, NT); cnt++)
            {
                phi = 0.0;
                for (int i = 0; i < NR; i++)
                {
                    complex z = complex_sub(YRep[cnt][i], Y[i]);
                    phi += sqrt(pow(z.re, 2) + pow(z.im, 2));
                }
                if (pre > phi)
                {
                    c = cnt;
                    pre = phi;
                }
            }
            // printf("min = %d\n", c);

            for (int i = 0; i < NT; i++)
                (X[i].re != XRep[c][i].re) ? error_sum++ : error_sum;
            for (int i = 0; i < NT; i++)
                (X[i].im != XRep[c][i].im) ? error_sum++ : error_sum;

            // printf("[%1.4f+j %1.4f ]\t[%1.4f+j %1.4f ]\n ", X[0].re, X[0].im, X[1].re, X[1].im);
        }
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