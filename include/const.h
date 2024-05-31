#ifndef CONST_H
#define CONST_H

#define ti 0
#define ENMIN -10
#define ENMAX 20
#define LOOP 1e6

#define NT 4
#define NR 6
#define F_WAVES 8
// #define DOPPLER 0.0000005 // wfading 関数で使用
#define DOPPLER 0.0 // wfading 関数で使用
#define PATH_SIZE 512
#define Pilots 16
#define MODE 0
#define DEBUG(MODE) if (MODE)

#define VOIDS -999
#define BIT_LENGTH 4
#define STATE_00 0
#define STATE_01 1
#define STATE_10 2
#define STATE_11 3
#define STATES 4

#define BIT_LENGTH 4
#define PASS puts("\e[32m PASS \e[m\n")

#endif // CONST_H
