#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "../include/def.h"
#include "../include/RFunction.h"
#include <unistd.h>
#include "../include/const.h"
#include "../include/error_detection.h" // ← new
#include "../include/econst.h"

int main(void)
{
    FILE *fp;
    char filename[50];
    char fullpath[PATH_SIZE];
    char *dir;

    getcwd(fullpath, PATH_SIZE);
    dir = strrchr(fullpath, '/');
    dir = dir + 1;
    printf("%s\n", dir);
    sprintf(filename, "../results/%s_%dx%d.d", dir, NT, NR);

    if ((fp = fopen(filename, "w+")) == NULL)
    {
        fprintf(stderr, "can not open write file.\n");
        exit(1);
    }

    printInfo();
    double error_sum = 0.0;
    double error_counter = 0.0, ber, en, stv, vd_meanValue = 0.0;
    int loop;
    complex X, N, Xd;
    complex Y[BEFORE_DECORDE_LEN + 1];
    //   Node *node = (Node *)malloc((NODE_SIZE) * sizeof(Node));
    int input[BEFORE_DECORDE_LEN];
    int *convoluted_input;
    int *ans;
    complex *qpsk_input;
    Node *node;
    int cnt;
    int *bits;
    int *ans_nodes;
    srand((unsigned)time(NULL));

    for (en = ENMIN; en <= ENMAX; en++)
    {
        stv = STV(en); // 標準偏差
        error_sum = 0.0;

        for (int loop = 0; loop < LOOP; loop++)
        {
            node = (Node *)malloc((NODE_SIZE) * sizeof(Node));
            initNode(node, NODE_SIZE);

            // 信号系列の作成 : BEFORE_DECORDE_LEN
            for (int j = 0; j < BEFORE_DECORDE_LEN; j++)
                input[j] = rand() % 2;
            input[0] = 0;
            // 畳み込み : BEFORE_DECORDE_LEN --> BEFORE_DECORDE_LEN *　2 + TAIL_BITS
            convoluted_input = encorder(input, BEFORE_DECORDE_LEN, 3);
            // QPSK :  BEFORE_DECORDE_LEN *　2 + TAIL_BITS -->  (BEFORE_DECORDE_LEN *2 + TAIL_BITS) /2
            qpsk_input = convo2QPSK(convoluted_input, AFTER_DECORDE_LEN);
            // 送信　：　(BEFORE_DECORDE_LEN *2 + TAIL_BITS) /2　bits
            for (int j = 0; j < BEFORE_DECORDE_LEN; j++)
            {
                N = gaussianRand(stv, 0);
                Y[j] = complex_add(qpsk_input[j], N);
                Y[j] = Hard(Y[j]);
            }
            // viterbi
            //  for (int j = 0; j < BEFORE_DECORDE_LEN; j++)
            //      printf("(%.1f+j%.1f)\t", Y[j].re, Y[j].im);
            //  printf("\n");

            bits = QPSK2bit(Y, BEFORE_DECORDE_LEN);
            //  printSignals(bits, AFTER_DECORDE_LEN);
            node = trellis(node, bits);
            ans_nodes = viterbi(node);

            // printResults(node);
            //  printSignals(ans_nodes, BEFORE_DECORDE_LEN);
            ans = getSurviversPath(ans_nodes);
            // printSurvivorPath(ans, convoluted_input, "answer  ");

            // printSignals(ans, AFTER_DECORDE_LEN);
            cnt = 0;
            for (int j = 0; j < AFTER_DECORDE_LEN; j++)
            {
                (convoluted_input[j] = ans[j]) ? error_sum++ : error_sum;
                cnt++;
            }
            free(qpsk_input);
            free(convoluted_input);
            free(node);
            free(ans);
            free(bits);
            free(ans_nodes);
            // PASS;
            // printf("%d\n", loop);
        }

        ber = error_sum / (LOOP * AFTER_DECORDE_LEN * 2); // 実数と虚数があるので分母2倍
        printf("EN:%.2f  ber:%.16f\n", en, ber);
        fprintf(fp, "%.2f\t%.16f\n", en, ber);
        fflush(fp); // バッファフラッシュ
    }

    fclose(fp);
    return 0;
}