/**********************************************************************
* Descrition: kalman filter core
* Original Author: jpdeng
* Create Date: 2017/6/12
**********************************************************************/

#ifndef _KALMANLITE_H_
#define _KALMANLITE_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum kalmanStatus
{
    KF_CORE_INIT = 0,
    KF_CORE_RESET,
    KF_CORE_RESTART,
    KF_CORE_RUNNING,
} kalmanStatus_t;

typedef enum updateFlag
{
    UPDATE_NONE = 0,
    UPDATE_SAVE = 1<<0,
} updateFlag_t;

/* upper part of upper triangular matrix will be stored in vector type
   A = [a11 a12 a13
         0  a22 a23
         0   0  a33]
   is stored as [a11, a12, a13, a22, a23, a33]. About half memory size is saved.
*/
typedef struct kalmanInfo
{
    U8      status;                     // KF status flag
    U16     periodTms;                  // KF update period
    U32     msCnt;                      // The time KF update
    U16     stateNum;                   // the number of state
    U16     matrixNum;                  // the number of upper triangular matrix elements
    DBL     *X;                         // KF state variables
    DBL     *A;                         // KF state transfer matrix, upper triangular matrix
    DBL     *U_plus;                    // upper triangular matrix of UD decomposition of P(k) matrix
    DBL     *U_minus;                   // upper triangular matrix of UD decomposition of P(k|k-1) matrix
    DBL     *D_plus;                    // diagonal matrix of UD decomposition of P(k) matrix
    DBL     *D_minus;                   // diagonal matrix of UD decomposition of P(k|k-1) matrix
    DBL     *W1;                        // W1 = A*U, upper triangular matrix
    DBL     *W2;                        // W2 = I(N) when initialization; W = [W1 W2], W is needed in UD KF calculation
    DBL     **Q;                        // process noise covariance matrix
    // reserve space for ud update need
    DBL     *f;
    DBL     *v;
    DBL     *b;
} kalmanInfo_t;

int kalmanInit(kalmanInfo_t* const pKalmanInfo, const U16 uStateNum);
S16 uMatIdx(S16 i, S16 j, S16 n);
void uxuMat(DBL* U1, DBL* U2, DBL* U3, S16 n);
void uMatxVect(DBL* U, DBL* V1, DBL* V2, S16 n);
void uMatInit(DBL* U, S16 n);
//void uMat2Vector(DBL U1[][N_STATE], DBL* U2);
void udKfPredict(kalmanInfo_t* const pKalmanInfo);
FLT udKFUpdate(kalmanInfo_t* const pKalmanInfo, DBL* H, DBL* deltaX, DBL R, DBL ino, DBL varTh, updateFlag_t flag);
void getUdMatDiag(DBL* U, DBL* D, DBL* diag, S16 n);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif