#ifndef _MISC_H_
#define _MISC_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     X               (0)
#define     Y               (1)
#define     Z               (2)
#define     CHN             (3)

#ifndef     PI
#define     PI              ((FLT)3.14159265358979323846)
#endif
#define     DEG2RAD         ((FLT)PI/(FLT)180.0)
#define     RAD2DEG         ((FLT)180.0/(FLT)PI)
#define     GRAVITY         ((FLT)9.80665)

#define     RE              (6378137.0)
#define     esqu            (0.00669437999013)
#define     RM(L)           (RE*(1-esqu)/pow((1-esqu*sin(L)*sin(L)),1.5))
#define     RN(L)           (RE/sqrt(1-esqu*sin(L)*sin(L)))

    DBL** mallocArray2D_DBL(U32 row, U32 col);
    U32 freeArray2D_DBL(DBL **pmatrix, U32 row, U32 col);
    U32 computeMeanStd(FLT* const mean, FLT* const std, const FLT array[][CHN], U32 count);
    void dcm2euler(const FLT cbn[3][3], FLT* const pyaw, FLT* const ppitch, FLT* const proll);
    void euler2dcm(FLT cbn[3][3], FLT fyaw, FLT fpitch, FLT froll);
    void euler2q(FLT q[4] ,FLT fyaw, FLT fpitch, FLT froll);
    void q2dcm(FLT q[4], FLT cbn[3][3]);
    void qNorm(FLT q[4]);
    void f3x3matrixTranspose(FLT matrix[3][3]);
    void fmatrixAeqI(FLT *A[], U32 rc);
    void f3x3matrixEqI(FLT matrix[3][3]);
    void f3x3matrixEqScalar(FLT matrix[3][3], FLT scalar);
    void f3x3matrixEqAxScalar(FLT A[][3], FLT scalar);
    void f3x3matrixEqMinusA(FLT A[][3]);
    U32 matrixMult(const DBL** const matrixA, const DBL** const matrixB, U32 rowA, U32 colA, U32 rowB, U32 colB, DBL** const matrixC);
    void fmatrixAeqInvA(FLT *A[], S32 iColInd[], S32 iRowInd[], S32 iPivot[], S32 isize, U32 *pierror);
    void eigencompute10(FLT A[][10], FLT eigval[], FLT eigvec[][10], S32 n);
    void udDecompose(DBL** const matrix, U32 n);
    void multPhimUp(const DBL** const phim, const DBL* const u, U32 n, DBL** const w);
    void storeUq(const U32 l, const U32 m, const DBL** const u, const U32 istart, const U32 jstart, DBL** const w);
    void udTimeUpdate(const U32 iw, const U32 jw, DBL** const w, const DBL* const dw, DBL* const u);
    void udMeasUpdate(DBL* const u, DBL* const x, const U32 n, const DBL r, const DBL* const a, const DBL z, DBL* const alpha, DBL* const res);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif