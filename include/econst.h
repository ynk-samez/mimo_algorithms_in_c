
#ifndef ECONST_H
#define ECONST_H
#include <limits.h>
#define VOIDS INT_MAX
#define TAIL_BITS 2
#define BEFORE_DECORDE_LEN 64
#define AFTER_DECORDE_LEN BEFORE_DECORDE_LEN * 2 + TAIL_BITS
#define NODE_SIZE BEFORE_DECORDE_LEN * 4 - 3 + 1
#define BINARY_STRING_LENGTH 9 // 8 bits + null terminator
#endif                         // ECONST_H
