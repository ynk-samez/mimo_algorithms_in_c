#ifndef ERROR_DETECTION_H
#define ERROR_DETECTION_H
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include "def.h"
#include "const.h"
#include "complexArith.h"
#include "econst.h"

typedef int signals;
typedef char *string;

typedef struct
{
    int state;
    int output;
} DFF;

struct node
{
    int pre_id;
    int id;
    int pm; // path metric
    struct node *pre;
};
typedef struct node Node;

DFF initRegistor(DFF D);
signals *encorder(signals *input, int len, int k);

signals viterbi_encorder(signals *input, int len);
void printSignals(signals *input, int len);

int hamming_distance(signals *x, signals *y, int len);
int transitionsCheck(string src, string dst);
string getLabels(string src, string dst);

int calcBranchMetric(string input, string labels);
void initNode(Node *nd, int size);

int search_min_id(Node *node, int dst_start_block, int dst_end_block);
void printResults(Node *node);
Node *trellis(Node *node, signals *encorded_signals);
int *viterbi(Node *node);
void printSurvivorPath(int *ans_nodes, signals *orign, string title);
char *ascii2bin(char ascii_code);
int *asciiToBinSeries(string msg, int len);
char binary2Ascii(const char *binary);

complex *convo2QPSK(signals *recived_signals, int len);
int *QPSK2bit(complex *qpsked, int len);

int *getSurviversPath(int *ans_nodes);
#endif // ERROR_DETECTION_H
