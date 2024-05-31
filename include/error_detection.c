#include "error_detection.h"
#include <stdio.h>
#include <string.h>
#include <string.h>
#include "../include/const.h"
DFF initRegistor(DFF D)
{
    D.state = 0;
    D.output = 0;
    return D;
}

void initNode(Node *nd, int size)
{
    for (int i = 0; i < size; i++)
    {
        nd[i].id = i;
        nd[i].pre = VOIDS;
        nd[i].pm = VOIDS;
        nd[i].pre_id = VOIDS;
    }
}

int calcBranchMetric(string input, string labels)
{
    int bm = 0;
    for (int i = 0; i < 2; i++)
    {
        if ((input[i] - '0') != (labels[i] - '0'))
            bm++;
    }
    return bm;
}

string getLabels(string src, string dst)
{
    string ans = "NG";
    if (0 == strcmp(src, "00") && 0 == strcmp(dst, "00"))
        ans = "00";
    if (0 == strcmp(src, "00") && 0 == strcmp(dst, "01"))
        ans = "11";

    if (0 == strcmp(src, "01") && 0 == strcmp(dst, "10"))
        ans = "01";
    if (0 == strcmp(src, "01") && 0 == strcmp(dst, "11"))
        ans = "10";

    if (0 == strcmp(src, "11") && 0 == strcmp(dst, "11"))
        ans = "01";
    if (0 == strcmp(src, "11") && 0 == strcmp(dst, "10"))
        ans = "10";

    if (0 == strcmp(src, "10") && 0 == strcmp(dst, "00"))
        ans = "11";
    if (0 == strcmp(src, "10") && 0 == strcmp(dst, "01"))
        ans = "00";

    return ans;
}

int transitionsCheck(string src, string dst)
{
    int ans = 0;
    // 00 --> 00 || 00 --> 01
    if (!strcmp(src, "00"))
    {
        if (!strcmp(dst, "00") || !strcmp(dst, "01"))
            ans = 1;
    }
    // 11 --> 11 or 11 --> 10
    if (!strcmp(src, "11") == 1)
    {
        if (!strcmp(dst, "11") || !strcmp(dst, "10"))
            ans = 1;
    }
    // 01 --> 10 or 01 --> 11
    if (!strcmp(src, "01"))
    {
        if (!strcmp(dst, "10") || !strcmp(dst, "11"))
            ans = 1;
    }
    // 10 --> 00 or 10 --> 00
    if (!strcmp(src, "10"))
    {
        if (!strcmp(dst, "00") || !strcmp(dst, "01"))
            ans = 1;
    }
    return ans;
}

signals *encorder(signals *input, int len)
{
    signals g0, g1;
    signals *ans = (signals *)calloc(len * 2, sizeof(signals));
    DFF D1, D2;

    D1 = initRegistor(D1);
    D2 = initRegistor(D2);

    int cnt = 0;
    for (int i = 0; i < len; i++)
    {
        D1.state = input[i];
        D2.state = D1.output;

        g0 = input[i] ^ D2.output;
        g1 = input[i] ^ D2.output ^ D1.output;

        ans[cnt] = g0;
        cnt++;
        ans[cnt] = g1;
        cnt++;

        D1.output = D1.state;
        D2.output = D2.state;
    }
    return ans;
}

signals viterbi_encorder(signals *input, int len)
{
    return 0;
}

void printSignals(signals *input, int len)
{
    printf("\e[32m ──────────── signals ──────────── \e[m\n");
    for (int i = 0; i < len; i++)
        printf("%d ", input[i]);
    printf("\n\n");
}

int hamming_distance(signals *x, signals *y, int len)
{
    int dist = 0;

    for (int i = 0; i < len; i++)
        if (x[i] != y[i])
            dist++;
    return dist;
}

Node decorder(Node *nd, int size)
{
    Node snd;
    return snd;
}
int search_min_id(Node *node, int dst_start_block, int dst_end_block)
{
    int i;
    int minval, minId;
    for (i = dst_start_block; VOIDS == node[i].pm; i++)
        ;
    minval = node[i].pm;
    minId = i;

    for (int i = dst_start_block; i < dst_end_block + 1; i++)
    {
        if (node[i].pm != VOIDS && node[i].pm < minval && i != minId)
        {
            minval = node[i].pm;
            minId = i;
        }
    }
    return minId;
}

void printResults(Node *node)
{
    printf("\e[32m \uf0a9  node[%02d]\e[m  pm = %d\n", 0, node[0].pm);
    for (int i = 1; i < BIT_LENGTH * 2 * 4 - 3; i++)
    {
        if (i % 4 == 1)
            printf("\n");
        if (node[i].pm == VOIDS)
            printf("\x1b[1m\e[32m \uf0a9  node[%02d] \e[m\x1b[1m  (null)  \n", i);
        else
            printf("\x1b[1m\e[32m \uf0a9  node[%02d] \e[m \x1b[1m\x1b[4m pm = %d\x1b[m  pre_id = %d\n", i, node[i].pm, node[i].pre_id);
    }
}

Node *trellis(Node *node, signals *input, signals *encorded_signals)
{
    string states[4] = {"00", "01", "10", "11"};
    int src, dst, src_id = 0, dst_id;
    char input_set[32];
    int cnt = 0, bm = 0, dst_start_block = 1, dst_end_block = dst_start_block + STATES - 1;
    src = 0;

    node[src].pm = 0;
    sprintf(input_set, "%d%d", encorded_signals[cnt], encorded_signals[cnt + 1]);
    cnt += 2;
    for (dst_id = 1; dst_id < STATES + 1; dst_id++)
    {
        dst = (dst_id - 1) % 4;
        if (transitionsCheck(input_set, states[dst]))
        {
            bm = calcBranchMetric(input_set, getLabels(states[src], states[dst]));
            if (node[dst_id].pm == VOIDS || node[dst_id].pm >= node[src_id].pm + bm)
            {
                node[dst_id].pm = node[src_id].pm + bm;
                node[dst_id].pre_id = src_id;
            }
        }
    }

    while (dst_end_block < BIT_LENGTH * 4 * 2 - 3)
    {
        sprintf(input_set, "%d%d", encorded_signals[cnt], encorded_signals[cnt + 1]);
        cnt += 2;
        // 過去のdstの中でpath metricが最小のnodeを探す.
        dst_start_block += STATES;
        dst_end_block = dst_start_block + STATES - 1;
        for (src_id = dst_start_block - STATES; src_id < dst_end_block - STATES + 1; src_id++)
        {
            if (src_id != VOIDS)
            {
                //---------- branch metricの計算 -------------
                for (dst_id = dst_start_block; dst_id < dst_end_block + 1; dst_id++)
                {
                    src = (src_id - 1) % STATES;
                    dst = (dst_id - 1) % STATES;

                    if (transitionsCheck(states[src], states[dst]) && node[src_id].pm != VOIDS) // その状態遷移がありえるか
                    {
                        bm = calcBranchMetric(input_set, getLabels(states[src], states[dst]));

                        if (node[dst_id].pm == VOIDS || node[dst_id].pm >= node[src_id].pm + bm)
                        {
                            node[dst_id].pm = node[src_id].pm + bm;
                            node[dst_id].pre_id = src_id;
                        }
                    }
                }
            }
        }
    }
    return node;
}

string viterbi(Node *node)
{
    int ans[7];
    string states[4] = {"00", "01", "10", "11"};
    int index, minId, minVal, i, j, cnt = 0;
    for (index = BIT_LENGTH * 4 * 2 - 3 - 1; index >= STATES; index -= STATES)
    {
        for (j = index; node[j].pm == VOIDS; j--)
            ;
        minVal = node[j].pm;
        minId = index;
        for (i = index; i > index - STATES; i--)
        {
            if (minVal > node[i].pm && node[i].pm != VOIDS)
            {
                minVal = node[i].pm;
                minId = i;
            }
        }
        ans[cnt++] = node[minId].pre_id;
        printf("%d\n", node[minId].pre_id);
    }

    //  for (cnt = 7 - 1; cnt >= 1; cnt--)
    //    printf("%s ・・・・　state(%s) -> state(%s)\n", getLabels(states[(ans[cnt] - 1) % 4], states[(ans[cnt - 1] - 1) % 4]), states[(ans[cnt] - 1) % 4], states[(ans[cnt - 1] - 1) % 4]);
}

// 00 00 00 NG 00 NG