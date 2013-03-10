/* cblas.h */
/* cblas headers */

#ifndef __CBLAS_H__
#define __CBLAS_H__

#define CBLAS_INDEX size_t  /* this may vary between platforms */

typedef enum {
    CblasRowMajor=101, CblasColMajor=102
} CBLAS_ORDER;

typedef enum {
    CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113
} CBLAS_TRANSPOSE;

typedef enum {
    CblasUpper=121, CblasLower=122
} CBLAS_UPLO;

typedef enum {
    CblasNonUnit=131, CblasUnit=132
} CBLAS_DIAG;

typedef enum {
    CblasLeft=141, CblasRight=142
} CBLAS_SIDE;

void cblas_dgemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
                 const  CBLAS_TRANSPOSE TransB, 
                 const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

#endif /* __CBLAS_H__ */
