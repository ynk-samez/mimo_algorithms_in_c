#ifndef RFUNCTION_H
#define RFUNCTION_H

#include "complexArith.h"

void EstH(complex **y, complex **x_p, complex **ans, int symbol, int packl);
double VectEuclideanDis(complex *v1, complex *v2, int dim);
void MakeRep(complex **X_rep, int N_T);
void MLD(complex **x_rep, complex **h, complex *y, complex *ans, int N_T, int N_R);
void MMSE_w(complex **h, complex **w, double stv2, int nt, int nr);
void MMSE_SIC(complex **h, complex *y, complex *ans, double stv2);
void SQRD(int *p, complex **q, complex **r);
int ViterbiECount(int *a, int *out);
void ViterbiDecoding(complex y, int *t, int *out0, int *out1, double *eucdis, int KLEN);
void CountComparison(double **eucdis, double *sum, int *route, int STATE);
void FollowRoute(int **route, int *decoded, int LEN, int STATE);
void QPSK(int *code, complex *x, int len);
void SISO(complex *x, complex *y, int len, double stv, int ave);
void ConvertBit(complex *z, int **recv, int len);
void Convolution(int *t, int *out0, int *out1, int KLEN);
void PathComparison(double **eucdis, double *sum, int *route, int STATE);

#endif
