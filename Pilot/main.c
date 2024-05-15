#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "def.h"
#include "RFunction.h"
#include "const.h"
#include "complexArith.h"

clock_t cpu_time_start;
clock_t cpu_time_end;
double sec_start;
double sec_end;
double result_time;

int main(int argc, char **argv)
{

    int L = atoi(*++argv);

    /* スレッド数を取得 */
    n = omp_get_max_threads();
    complex **H, *X, *Y, *N, **XRep, **YRep;
    FILE *fp;
    char filename[50];
    sprintf(filename, "../Results/MLD-mimo%d-%d.d", NT, NR);
    if ((fp = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    // fprintf(fp, "Eb/N0[db]\tBER\n");
    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
    int i, j, k, loop;
    double idl;
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
            X = genRandomInput(NT);
            H = genChannelComplex(NT, NR);
            N = genNoiseComplex(stv);
            Y = mulComplexMatVec(H, X);
            Y = addComplexVec(Y, N);

            MakeRep(XRep, NT);
            YRep = mulComplexMat(H, XRep);
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

            for (int i = 0; i < NT; i++)
                (X[i].re != XRep[c][i].re) ? error_sum++ : error_sum;
            for (int i = 0; i < NT; i++)
                (X[i].im != XRep[c][i].im) ? error_sum++ : error_sum;
        }
        free(X);
        free(H);
        free(N);
        free(Y);
        ber = error_sum / ((LOOP) * 2 * NT); // 実数と虚数があるので分母2倍
        printf("EN:%.2f  ber:%.8f \n", en, ber);
        fprintf(fp, "%.2f\t%.8f\n", en, ber);
    }
    fclose(fp);
    return 0;
}
