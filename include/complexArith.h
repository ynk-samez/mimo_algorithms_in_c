#ifndef COMPLEXARRITH_H
#define COMPLEXARRITH_H
#define SWAP(type, a, b) \
    {                    \
        type temp = a;   \
        a = b;           \
        b = temp;        \
    }
typedef struct
{
    double re;
    double im;
} complex;

void printComplex(complex);
double norm(complex *vec, int h);
complex **exchangeColmn(complex **A, int w, int h, int col_a, int col_b);
int argmin(double *norms, int h);
complex **makeZeroMatrix(int w, int h);
int cmp(const void *x, const void *y);
void copyMatrix(complex **source, complex **destination, int rows, int cols);
complex **calcJmatrix(complex **W, complex **H, int w, int h);
complex **calcWeightMatrix(complex **H, int w, int h, double stv);
complex **hermitianMatrix(complex **Matrix, int w, int h);
complex **deleteCol(complex **H, int w, int h, int index);
complex *makeVecFromMat(complex **Mat, int h, int index);
complex *mulEach(complex *Vec, complex a);
void printPotinter(char *error_pointer, complex **A);
int minDiag(complex **A, int);

complex *complexVecSub(complex *A, complex *B, int Size);
complex **complexMatrixSub(complex **A, complex **B, int w, int h);
complex **complexMatrixAdd(complex **A, complex **B, int w, int h);
complex *conjugateVec(complex *, int);
complex **mulComplexVecVec(complex *A, complex *B, int sizeA, int sizeB);
complex **complexMatrixAve(complex **A, int w, int h);

complex **genE(int Size);

complex **scalarMat(complex **, int w, int h, double);
complex **addComplexMat(complex **, complex **, int w, int h);
complex **transposeMatrix(complex **A, int w, int h);
complex **conjugateMatrix(complex **, int w, int h);
complex **invMatrix(complex **, int Size);
void printAnyComplexMatrixColored(complex **A, int W, int H);
void printAnyComplexMatrix(complex **A, int W, int H);
void printAnyComplexVec(complex *A, int W);
complex *genRandomInput();
complex **genChannelComplex();
complex **mulComplexMat(complex **A, complex **B, int input1_w, int input1_h, int input2_w, int input2_h);
complex *mulComplexMatVec(complex **A, complex *B, int mat_w, int mat_h);

complex *addComplexVec(complex *, complex *);
complex *genNoiseComplex(double);
complex complex_add(complex x, complex y);
complex complex_sub(complex x, complex y);
complex complex_mul(complex x, complex y);
complex complex_div(complex x, complex y);

#endif