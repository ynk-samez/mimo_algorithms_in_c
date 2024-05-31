#ifndef ERROR_DETECTION_H
#define ERROR_DETECTION_H
#include <stdlib.h>
#include <limits.h>
#include "const.h"
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
signals *encorder(signals *input, int len);

signals viterbi_encorder(signals *input, int len);
void printSignals(signals *input, int len);

int hamming_distance(signals *x, signals *y, int len);
int transitionsCheck(string src, string dst);
string getLabels(string src, string dst);

int calcBranchMetric(string input, string labels);
void initNode(Node *nd, int size);
Node decorder(Node *nd, int size);

int search_min_id(Node *node, int dst_start_block, int dst_end_block);
void printResults(Node *node);
Node *trellis(Node *node, signals *input, signals *encorded_signals);
string viterbi(Node *node);
#endif // ERROR_DETECTION_H
