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
#include "../include/error_detection.h"

int main(void)
{
    int *ans;
    Node *node = (Node *)malloc((NODE_SIZE) * sizeof(Node));
    initNode(node, NODE_SIZE);

    signals input[] = {0, 1, 0, 0, 1, 0, 0};
    int input_len = sizeof(input) / sizeof(signals);
    signals *encorded_signals = encorder(input, input_len, 3);
    printSignals(encorded_signals, input_len * 2);
    node = trellis(node, encorded_signals);

    printResults(node);
    ans = viterbi(node);
    printSurvivorPath(ans, encorded_signals, "encorded");

    //  encorded_signals[5] ^= encorded_signals[5];

    // node = trellis(node, encorded_signals);
    // ans = viterbi(node);
    // printSurvivorPath(ans, encorded_signals, "error   ");

    // string msg = "hello world";
    // asciiToBinSeries(msg, 12);
    return 0;
}