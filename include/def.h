#ifndef DEF_H
#define DEF_H
#include "complexArith.h"
void printInfo();
complex *genChannel(int width, int height, int stv, int mean);
double scalar(complex A);
void MatTrans(complex **A, complex **ans, int A_row, int A_column); // 共役転置行列の生成
void VectTrans(complex *A, complex **ans, int value);

complex wfading(complex h, long int t);

double ComAbsSqu(complex A);

double STV(double en); // 標準偏差の計算
// en = Eb/N0
double STV_2(double en);                                    // 雑音の標準偏差を算出
complex Hard(complex x);                                    // 硬判定
double count(complex x, complex y, double k);               // 誤りカウント
complex gaussianRand(double vd_stdev, double vd_meanValue); //  標準偏差σ  平均
complex rand_complex(complex x);                            // ±1をランダムで出力する
complex **calloc_com2d(int value1, int value2);             // 2次元の行列のメモリ領域確保
complex *calloc_com(int value);                             // 1次元の行列のメモリ領域確保
double *calloc_d(int value);
double **calloc_d2d(int value1, int value2);
int **calloc_i2d(int value1, int value2);
void InitComMat(complex **init, int row, int col);
//**：ポインタのポインタ
void InitComVect(complex *init, int value); // 複素数ベクトル行列の初期化

void InitMat(double **init, int col, int row);
void InitDMat(double **init, int row, int col);
void InitDVect(double *init, int value);
complex **recalloc_com2d(complex **c, int w, int h);                         // 2列以上の行列のメモリ領域の確保;
void MatVectMulCom(complex **A, complex *B, complex *ans, int row, int col); // ベクトル行列と行列の乗算
void VectAddCom(complex *A, complex *B, complex *ans, int value);            // ベクトル行列の加算
void VectSubCom(complex *A, complex *B, complex *ans, int value);            // ベクトル行列の減算
int *calloc_i(int value);
void Free_com2d(complex **c, int value);  // 2次元の行列のメモリ開放
void ShinFree2d(complex ***c, int value); // 2次元の行列のメモリ開放
void Free_d2d(double **d, int value);

void Free_i2d(int **i, int value);
void IdentityComMat(complex **imat, int n); // 単位行列の生成
void MatMatMulCom(complex **A, complex **B, complex **ans, int Arow, int Acol, int Bcol);
void Mat4MatMulCom(complex ****A, complex **B, complex **ans, int Arow, int Acol, int A3, int Bcol);
void MatMatAddCom(complex **A, complex **B, complex **ans, int row, int col);
void MatMatSubCom(complex **A, complex **B, complex **ans, int row, int col);
void InvComMat(complex **A, complex **inv_A, int n);       // 複素数逆行列の生成
void IdentitystvComMat(complex **imat, double stv, int n); // 単位複素数行列と雑音の標準偏差の乗算
void Identitystv2ComMat(complex **imat, double stv2, int n);
void EnseAve(complex **A, complex **ans, int row, int col, int n);

// 行列表示の関数
void PrintComVect(complex *A, int value);
void PrintComMat(complex **A, int row, int col);
void PrintComReMat(complex **A, int row, int col);
void Print_d_Mat(double **A, int row, int col);
void Print_i_Mat(int **A, int row, int col);

void excel_d_Vect(double *z, int n, int a);
void excel_d_Vect3(double *x, int x_n, double *y, int y_n, double *z, int z_n, int a);

void ChangeCol_i(int *A, int x, int y);
void ChangeCol_d(double *A, int x, int y);
void ChangeCol_com(complex *A, int x, int y);
void ChangeCol_MatCom(complex **A, int A_row, int x, int y);

void Cp_Mat(double **A, double **cp_A, int col, int row);
void Cp_Mat_i(int **A, int **cp_A, int col, int row);
void Cp_Vect(double *A, double *cp_A, int n);

#endif