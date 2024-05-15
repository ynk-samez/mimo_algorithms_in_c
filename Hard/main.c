#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "def.h"
#include "RFunction.h"
#include <unistd.h>
#include "const.h"

int main(void)
{
  FILE *fp;
  char filename[50];
  sprintf(filename, "rinko.csv");

  if ((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "can not open write file.\n");
    exit(1);
  }

  fprintf(fp, "Eb/N0[db]\tBER\n");
  double error_sum = 0.0;
  double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
  int loop;
  complex X, N, Y, Xd;

  srand((unsigned)time(NULL)); // 乱数を時間をもとに生成
  for (en = ENMIN; en <= ENMAX; en++)
  {
    stv = STV(en); // 標準偏差
    error_sum = 0.0;
    for (loop = 0; loop < LOOP; loop++)
    {
      X = rand_complex(X);
      N = gaussianRand(stv, 0);

      Y = complex_add(X, N);
      Xd = Hard(Y);
      (X.re != Xd.re) ? error_sum++ : error_sum;
      (X.im != Xd.im) ? error_sum++ : error_sum;
    }
    ber = error_sum / ((LOOP) * 2 * NT); // 実数と虚数があるので分母2倍
    printf("EN:%.2f  ber:%.8f\n", en, ber);
    fprintf(fp, "%.2f,%.8f\n", en, ber);
  }
  fclose(fp);
  return 0;
}