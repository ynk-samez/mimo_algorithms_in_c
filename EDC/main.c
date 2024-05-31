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
    Node *node = (Node *)malloc((BIT_LENGTH * 4 * 2 - 3) * sizeof(Node));
    initNode(node, BIT_LENGTH * 2 * 4 - 3);

    signals input[] = {0, 1, 0, 0, 1, 0, 0};
    signals *encorded_signals = encorder(input, BIT_LENGTH);

    printSignals(encorded_signals, BIT_LENGTH * 2);
    node = trellis(node, input, encorded_signals);

    printResults(node);
    viterbi(node);
    return 0;
}