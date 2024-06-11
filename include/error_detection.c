#include "error_detection.h"

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
        nd[i].pre = NULL;
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
    // 00 -> 00
    if (!strcmp(src, "00") && !strcmp(dst, "00"))
        ans = "00";
    // 00 -> 11
    if (!strcmp(src, "00") && !strcmp(dst, "01"))
        ans = "11";

    // 01 -> 10
    if (!strcmp(src, "01") && !strcmp(dst, "10"))
        ans = "01";
    // 01 -> 11
    if (!strcmp(src, "01") && !strcmp(dst, "11"))
        ans = "10";

    // 11 -> 11
    if (!strcmp(src, "11") && !strcmp(dst, "11"))
        ans = "01";
    // 11 -> 10
    if (!strcmp(src, "11") && !strcmp(dst, "10"))
        ans = "10";

    // 10 -> 00
    if (!strcmp(src, "10") && !strcmp(dst, "00"))
        ans = "11";
    // 10 -> 01
    if (!strcmp(src, "10") && !strcmp(dst, "01"))
        ans = "00";

    return ans;
}

int transitionsCheck(string src, string dst)
{
    int ans = 0; // NG
    // 00 --> 00 || 00 --> 01
    if (!strcmp(src, "00"))
    {
        if (!strcmp(dst, "00") || !strcmp(dst, "01"))
            ans = 1;
    }
    // 11 --> 11 or 11 --> 10
    if (!strcmp(src, "11"))
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

signals *encorder(signals *input, int len, int k)
{
    signals g0, g1;
    signals *ans = (signals *)calloc(len * 2 + TAIL_BITS, sizeof(signals));
    DFF D1, D2;
    D1 = initRegistor(D1);
    D2 = initRegistor(D2);

    int cnt = -1;
    for (int i = 0; i < len; i++)
    {
        D1.state = input[i];
        D2.state = D1.output;

        g0 = input[i] ^ D2.output;
        g1 = input[i] ^ D1.output ^ D2.output;

        ans[++cnt] = g0;
        ans[++cnt] = g1;

        D1.output = D1.state;
        D2.output = D2.state;
    }
    ans[len * 2] = ans[len * 2 + 1] = 0;

    return ans;
}

void printSignals(signals *input, int len)
{
    printf("\e[32m ──────────── signals ──────────── \e[m\n");
    for (int i = 0; i < len; i++)
        printf("%d ", input[i]);

    printf("\e[33m + %d ", input[len]);
    printf("%d \e[m", input[len + 1]);
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
    printf("\x1b[1m\e[32m \uf0a9  node[%02d]\e[m  \x1b[1m\x1b[4m pm = %d\e[m\n", 0, node[0].pm);
    for (int i = 1; i < NODE_SIZE - 1; i++)
    {
        if (i % 4 == 1)
            printf("\n");
        if (node[i].pm == VOIDS)
            printf("\x1b[1m\e[32m \uf0a9  node[%02d] \e[m\x1b[1m  (null)  \n", i);
        else
            printf("\x1b[1m\e[32m \uf0a9  node[%02d] \e[m \x1b[1m\x1b[4m pm = %d\x1b[m  pre_id = %d\n", i, node[i].pm, node[i].pre_id);
    }
    printf("\n\n\x1b[1m\e[32m \uf0a9  node[%02d]\e[m  \x1b[1m\x1b[4m pm = %d\e[m pre_id = %d\n\n", NODE_SIZE - 1, node[NODE_SIZE - 1].pm, node[NODE_SIZE - 1].pre_id);
}

Node *trellis(Node *node, signals *encorded_signals)
{
    string states[4] = {"00", "01", "10", "11"};
    int src, dst, src_id = 0, dst_id;
    char input_set[32];
    int cnt = 0, bm = 0;
    int src_start_block, src_end_block;
    int dst_start_block = 1, dst_end_block = dst_start_block + STATES;
    int min;
    src = 0;

    node[src].pm = 0;
    sprintf(input_set, "%d%d", encorded_signals[cnt], encorded_signals[cnt + 1]);
    cnt += 2;
    // printf("inputs set: %s\n", input_set);

    for (dst_id = 1; dst_id < STATES + 1; dst_id++)
    {
        dst = (dst_id - 1) % 4;
        if (transitionsCheck(input_set, states[dst]))
        {
            bm = calcBranchMetric(input_set, getLabels(states[src], states[dst]));
            if (node[dst_id].pm == VOIDS || node[dst_id].pm > node[src_id].pm + bm)
            {
                node[dst_id].pm = node[src_id].pm + bm;
                node[dst_id].pre_id = src_id;
            }
        }
    }

    src_start_block = 1;
    src_end_block = src_start_block + STATES;
    dst_start_block = src_end_block;
    dst_end_block = dst_start_block + STATES;

    while (dst_end_block < NODE_SIZE)
    {
        /* 2bits from input */
        sprintf(input_set, "%d%d", encorded_signals[cnt], encorded_signals[cnt + 1]);
        cnt += 2;

        //   過去のdst(次のsrc)の中でpath metricが最小のnodeを探す.
        for (src_id = src_start_block; src_id < src_end_block; src_id++)
        {
            if (node[src_id].pm != VOIDS)
            {
                //---------- branch metricの計算 -------------
                for (dst_id = dst_start_block; dst_id < dst_end_block; dst_id++)
                {
                    src = (src_id - 1) % STATES;
                    dst = (dst_id - 1) % STATES;
                    // printf("%d --> %d\n", src, dst);

                    if (transitionsCheck(states[src], states[dst])) // その状態遷移がありえるか
                    {
                        bm = calcBranchMetric(input_set, getLabels(states[src], states[dst]));
                        // printf("%d : src( %d ) + bm( %d )\n", node[dst_id].pm, node[src_id].pm, bm);
                        //   printf("dst : %d\n", node[dst_id].pm);
                        // 行き先のパスメトリックがまだ計算されていない or より小さいパスメトリックを見つけた
                        if (node[dst_id].pm == VOIDS || ((node[dst_id].pm != VOIDS) && node[dst_id].pm > node[src_id].pm + bm))
                        {
                            node[dst_id].pm = node[src_id].pm + bm;
                            node[dst_id].pre_id = src_id;
                            // printf("dst : %d\n", dst_id);
                            //   printf("%d\n", src_id);
                            //   printf("%d", node[src_id].pm);
                        }
                    }
                }
            }
        }
        // printf("%d : %d\n", src_start_block, src_end_block);
        src_start_block = dst_start_block;
        src_end_block = src_start_block + STATES;

        dst_start_block = src_end_block;
        dst_end_block = dst_start_block + STATES;
    }

    // 最後に最小値が二つある場合のために
    dst_id = NODE_SIZE - 1;
    dst = 0;
    for (src_id = dst_id - 4; src_id < dst_id; src_id++)
    {
        if (node[src_id].pm != VOIDS)
        {
            src = (src_id - 1) % STATES;
            if (transitionsCheck(states[src], states[dst])) // その状態遷移がありえるか
            {
                bm = calcBranchMetric("00", getLabels(states[src], states[dst]));
                if (node[dst_id].pm == VOIDS || ((node[dst_id].pm != VOIDS) && node[dst_id].pm > node[src_id].pm + bm))
                {
                    node[dst_id].pm = node[src_id].pm + bm;
                    node[dst_id].pre_id = src_id;
                }
            }
        }
    }

    return node;
}

int *viterbi(Node *node)
{
    int *ans_nodes = (int *)malloc(BEFORE_DECORDE_LEN * sizeof(int));
    int min, min_id;

    string states[4] = {"00", "01", "10", "11"};

    int i;
    for (i = NODE_SIZE - 2; node[i].pm == VOIDS && i > NODE_SIZE - 2 - 4; i--)
        ;
    min = node[i].pm;
    min_id = i;
    int min_cnt = 0;

    for (int i = NODE_SIZE - 2; i > NODE_SIZE - 2 - 4; i--)
    {
        if (node[i].pm < min && node[i].pm != VOIDS && i != min_id)
        {
            min = node[i].pm;
            min_id = i;
        }
    }
    for (int i = NODE_SIZE - 2; i > NODE_SIZE - 2 - 4; i--)
    {
        if (node[i].pm == min)
            min_cnt++;
    }

    int cnt = BEFORE_DECORDE_LEN - 1;
    if (min_cnt > 1)
        ans_nodes[cnt] = node[NODE_SIZE - 1].pre_id;
    else
        ans_nodes[cnt] = min_id;
    int pre = ans_nodes[cnt];
    // printf("%d\n", pre);
    cnt--;

    for (int j = BEFORE_DECORDE_LEN - 2; j >= 0; j--)
    {
        ans_nodes[cnt] = node[pre].pre_id;
        pre = ans_nodes[cnt];
        // printf("%d\n", pre);
        cnt--;
    }
    ans_nodes[0] = 0;

    return ans_nodes;
}

int *getSurviversPath(int *ans_nodes)
{
    int *ans = (int *)malloc(sizeof(int) * AFTER_DECORDE_LEN);
    string label;
    string states[4] = {"00", "01", "10", "11"};
    string src_state, dst_state;
    int cnt = -1;
    for (int i = 0; i < BEFORE_DECORDE_LEN - 1; i++)
    {

        int src_node_id = ans_nodes[i];
        int dst_node_id = ans_nodes[i + 1];
        // printf("%d --> %d\n", src_node_id, dst_node_id);
        src_state = states[(src_node_id - 1) % 4];
        dst_state = states[(dst_node_id - 1) % 4];
        // printf("%s\n", src_state);
        // printf("%s\n", dst_state);
        label = getLabels(src_state, dst_state);
        // printf("%s\n", label);
        ans[++cnt] = label[0] - '0';
        ans[++cnt] = label[1] - '0';
    }
    return ans;
}

void printSurvivorPath(int *ans, signals *orign, string title)
{
    string states[4] = {"00", "01", "10", "11"};
    // printf("\n──────────── survivor path ────────────\n");

    printf("\e[m\n\n  \x1b[1m %s  | ", title);
    for (int i = 0; i < AFTER_DECORDE_LEN; i++)
    {
        printf("\x1b[1m%d", orign[i]);
        if ((i + 1) % 2 == 0)
            printf(" ");
    }
    printf("\e[m");
    printf("\n  \x1b[1m decorded  | ");
    printf("%s ", getLabels(states[0], states[(ans[1] - 1) % 4]));
    for (int i = 0; i < BEFORE_DECORDE_LEN; i++)
    {
        int src_node_id = ans[i];
        int dst_node_id = ans[i + 1];
        string src_state = states[(src_node_id - 1) % 4];
        string dst_state = states[(dst_node_id - 1) % 4];
        printf("\x1b[1m%s ", getLabels(src_state, dst_state));
    }
    printf("\e[m\n\n");
}
char *ascii2bin(char ascii_code)
{
    char *binary = (char *)malloc(BINARY_STRING_LENGTH * sizeof(char));
    if (binary == NULL)
    {
        perror("Failed to allocate memory");
        exit(EXIT_FAILURE);
    }

    binary[BINARY_STRING_LENGTH - 1] = '\0'; // Null terminator

    for (int i = BINARY_STRING_LENGTH - 2; i >= 0; --i)
    {
        binary[i] = (ascii_code & 1) ? '1' : '0';
        ascii_code >>= 1;
    }

    return binary;
}

int *asciiToBinSeries(string msg, int len)
{
    char *bin;
    int cnt = 0;
    int ans[(len - 1) * (BINARY_STRING_LENGTH - 1)];

    for (int j = 0; j < len - 1; j++)
    {
        printf("%c | \e[32m", msg[j]);
        bin = ascii2bin(msg[j]);
        for (int i = 0; i < BINARY_STRING_LENGTH - 1; i++)
        {
            printf("%d ", bin[i] - '0');
            ans[cnt++] = bin[i] - '0';
        }
        printf("\e[m\n");
    }
    return ans;
}

char binary2Ascii(const char *binary)
{
    if (strlen(binary) != 8)
    {
        fprintf(stderr, "Error: Binary string must be 8 bits long.\n");
        exit(EXIT_FAILURE);
    }

    char ascii_char = 0;
    for (int i = 0; i < 8; ++i)
    {
        if (binary[i] == '1')
        {
            ascii_char |= (1 << (7 - i));
        }
        else if (binary[i] != '0')
        {
            fprintf(stderr, "Error: Binary string must contain only '0' or '1'.\n");
            exit(EXIT_FAILURE);
        }
    }

    return ascii_char;
}
int *QPSK2bit(complex *qpsked, int len)
{
    int *ans, cnt = 0;
    ans = (int *)malloc(sizeof(int) * len * 2);

    for (int i = 0; i < len; i++)
    {
        if (qpsked[i].re == 1.0)
            ans[cnt] = 1;
        else if (qpsked[i].re == -1.0)
            ans[cnt] = 0;
        cnt++;

        if (qpsked[i].im == 1.0)
            ans[cnt] = 1;
        else if (qpsked[i].im == -1.0)
            ans[cnt] = 0;
        cnt++;
    }
    // printf("\n");
    return ans;
}

complex *convo2QPSK(signals *recived_signals, int len)
{
    complex *ans;
    ans = calloc_com(len / 2);
    int cnt = 0;
    for (int i = 0; i < len; i += 2)
    {

        if (recived_signals[i] == 0)
        {
            ans[cnt].re = -1.0;
        }
        else
        {
            ans[cnt].re = 1.0;
        }
        cnt++;
        // printf("%f\n", ans[cnt].re);

        // printf("%f\n", ans[cnt].re);
        if (recived_signals[i + 1] == 0)
        {
            ans[cnt].im = -1.0;
        }
        else
        {
            ans[cnt].im = 1.0;
        }
        cnt++;
        //  if (recived_signals[i] != 0 && recived_signals[i] != 1)
        //      printf("ERROR OCCUR!\n");
        //  if (recived_signals[i + 1] != 0 && recived_signals[i + 1] != 1)
        //      printf("ERROR OCCUR!\n");

        //    printf("%f\n\n", ans[cnt].im);

        // printf("%d\n", cnt);
    }
    return ans;
}