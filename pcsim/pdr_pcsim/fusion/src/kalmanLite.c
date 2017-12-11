/**********************************************************************
* Descrition: kalman filter core
* Original Author: jpdeng
* Create Date: 2017/6/12
**********************************************************************/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "kalmanLite.h"

static DBL** mallocArray2D_DBL(U32 row, U32 col);
static U32 freeArray2D_DBL(DBL **pmatrix, U32 row, U32 col);

static DBL** mallocArray2D_DBL(U32 row, U32 col)
{
    DBL **pmatrix = NULL;
    U32 i = 0;
    U32 sizeofDBL = sizeof(DBL);
    U32 sizeofPtrDBL = sizeof(DBL *);

    pmatrix = (DBL **)malloc(sizeofPtrDBL * row);

    if (pmatrix == NULL)
    {
        return NULL;
    }

    for (i = 0; i < row; i++)
    {
        pmatrix[i] = (DBL *)malloc(sizeofDBL * col);

        if (pmatrix[i] == NULL)
        {
            U32 j;

            for (j = 0; j < i; j++)
            {
                free(pmatrix[j]);
            }
            free(pmatrix);

            return NULL;
        }
        memset(pmatrix[i], 0, sizeofDBL * col);
    }

    return pmatrix;
}

static U32 freeArray2D_DBL(DBL **pmatrix, U32 row, U32 col)
{
    U32 i = 0;
    PARAMETER_NOT_USED(col);

    for (i = 0; i < row; i++)
    {
        free(pmatrix[i]);
    }

    return 0;
}

int kalmanInit(kalmanInfo_t* const pKalmanInfo, const U16 uStateNum)
{
    U32 sizeofDBL = sizeof(DBL);

    pKalmanInfo->stateNum = uStateNum;
    pKalmanInfo->matrixNum = uStateNum * (uStateNum + 1) / 2 + 1; // start from suffix 1

    // malloc X array
    pKalmanInfo->X = (DBL *)malloc(sizeofDBL * uStateNum);
    if (pKalmanInfo->X == NULL)
    {
        return -1;
    }
    memset(pKalmanInfo->X, 0, sizeofDBL * uStateNum);

    // malloc A array
    pKalmanInfo->A = (DBL *)malloc(sizeofDBL * pKalmanInfo->matrixNum);
    if (pKalmanInfo->A == NULL)
    {
        free(pKalmanInfo->X);
        return -1;
    }
    memset(pKalmanInfo->A, 0, sizeofDBL * pKalmanInfo->matrixNum);

    // malloc U array
    pKalmanInfo->U_plus = (DBL *)malloc(sizeofDBL * pKalmanInfo->matrixNum);
    if (pKalmanInfo->U_plus == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        return -1;
    }
    memset(pKalmanInfo->U_plus, 0, sizeofDBL * pKalmanInfo->matrixNum);

    pKalmanInfo->U_minus = (DBL *)malloc(sizeofDBL * pKalmanInfo->matrixNum);
    if (pKalmanInfo->U_minus == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        return -1;
    }
    memset(pKalmanInfo->U_minus, 0, sizeofDBL * pKalmanInfo->matrixNum);

    // malloc D array
    pKalmanInfo->D_plus = (DBL *)malloc(sizeofDBL * (uStateNum + 1)); // start from suffix 1
    if (pKalmanInfo->D_plus == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        return -1;
    }
    memset(pKalmanInfo->D_plus, 0, sizeofDBL * (uStateNum + 1));

    pKalmanInfo->D_minus = (DBL *)malloc(sizeofDBL * (uStateNum + 1)); // start from suffix 1
    if (pKalmanInfo->D_minus == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        return -1;
    }
    memset(pKalmanInfo->D_minus, 0, sizeofDBL * (uStateNum + 1));

    // malloc W array
    pKalmanInfo->W1 = (DBL *)malloc(sizeofDBL * pKalmanInfo->matrixNum);
    if (pKalmanInfo->W1 == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        return -1;
    }
    memset(pKalmanInfo->W1, 0, sizeofDBL * pKalmanInfo->matrixNum);

    pKalmanInfo->W2 = (DBL *)malloc(sizeofDBL * pKalmanInfo->matrixNum);
    if (pKalmanInfo->W2 == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        free(pKalmanInfo->W1);
        return -1;
    }
    memset(pKalmanInfo->W2, 0, sizeofDBL * pKalmanInfo->matrixNum);

    // malloc reserve space
    pKalmanInfo->f = (DBL *)malloc(sizeofDBL * (uStateNum + 1)); // start from suffix 1
    if (pKalmanInfo->f == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        free(pKalmanInfo->W1);
        free(pKalmanInfo->W2);
        return -1;
    }
    memset(pKalmanInfo->f, 0, sizeofDBL * (uStateNum + 1));

    // malloc reserve space
    pKalmanInfo->v = (DBL *)malloc(sizeofDBL * (uStateNum + 1)); // start from suffix 1
    if (pKalmanInfo->v == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        free(pKalmanInfo->W1);
        free(pKalmanInfo->W2);
        free(pKalmanInfo->f);
        return -1;
    }
    memset(pKalmanInfo->v, 0, sizeofDBL * (uStateNum + 1));

    // malloc reserve space
    pKalmanInfo->b = (DBL *)malloc(sizeofDBL * (uStateNum + 1)); // start from suffix 1
    if (pKalmanInfo->b == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        free(pKalmanInfo->W1);
        free(pKalmanInfo->W2);
        free(pKalmanInfo->f);
        free(pKalmanInfo->v);
        return -1;
    }
    memset(pKalmanInfo->b, 0, sizeofDBL * (uStateNum + 1));

    // malloc Q array
    pKalmanInfo->Q = (DBL **)mallocArray2D_DBL(uStateNum, uStateNum);
    if (pKalmanInfo->Q == NULL)
    {
        free(pKalmanInfo->X);
        free(pKalmanInfo->A);
        free(pKalmanInfo->U_plus);
        free(pKalmanInfo->U_minus);
        free(pKalmanInfo->D_plus);
        free(pKalmanInfo->D_minus);
        free(pKalmanInfo->W1);
        free(pKalmanInfo->W2);
        free(pKalmanInfo->f);
        free(pKalmanInfo->v);
        free(pKalmanInfo->b);
        return -1;
    }
    uMatInit(pKalmanInfo->U_plus, uStateNum);
    uMatInit(pKalmanInfo->U_minus, uStateNum);

    return 0;
}

/**********************************************************************
* Function Name:    uMatIdx
*
* Description:
*    get the index of upper trianular matrix U(i,j) that stored in vector type
*
* Input:
*     i: row index of matrix
*     j: column index of matrix
*     N:  row and column length of matrix
*
* Return:
*     index of upper triangular that stored in vector type
*
* Author: jpdeng
**********************************************************************/
S16 uMatIdx(S16 i, S16 j, S16 N)
{
    //S16 i1, i2;
    //i1 = i - 1;
    //i2 = i - 2;

    //return i1 * N - (i1 * i2 >> 1) + j - i + 1;

    S16 i1;
    i1 = i - 1;

    return i1 * N - (i1 * i >> 1) + j;
}


/**********************************************************************
* Function Name:    uxuMat
*
* Description:
*    upper triangular matrix multiply upper triangular matrix
*
* Input:
*     U1: upper triangular matrix
*     U2: upper triangular matrix
*     U3: U3 = U1*U2
*     N:  row and column length of matrix
*
* Return:
*
* Author: jpdeng
**********************************************************************/
void uxuMat(DBL* U1, DBL* U2, DBL* U3, S16 N)
{
    S16 i, j, k, m, idx1, idx2;
    DBL tmp;

    m = 1;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            tmp = 0;

            for (k = i; k <= j; k++)
            {
                idx1 = uMatIdx(i, k, N);
                idx2 = uMatIdx(k, j, N);
                tmp += U1[idx1] * U2[idx2];
            }

            U3[m++] = tmp;
        }
    }
}

/**********************************************************************
* Function Name:    uMatInit
*
* Description:
*    initialize the diagonal value of upper triangular matrix to 1
*
* Input:
*     U: upper triangular matrix
*     N:  row and column length of matrix
*
* Return:
*
* Author: jpdeng
**********************************************************************/
void uMatInit(DBL* U, S16 N)
{
    S16 i, idx;

    for (i = 1; i <= N; i++)
    {
        idx = uMatIdx(i, i, N);
        U[idx] = 1;
    }
}

/**********************************************************************
* Function Name:    uMat2Vector
*
* Description:
*    store upper part of U1[N][N] to U2[N*(N+1)/2]
*
* Input:
*     U1: upper triangular matrix in 2 dimension
*     U2: upper triangular matrix in vector type
*
* Return:
*
* Author: jpdeng
**********************************************************************/
#if 0
void uMat2Vector(DBL U1[][N_STATE], DBL* U2)
{
    S16 i, j, idx, N;
    N = N_STATE;

    for (i = 1; i <= N; i++)
    {
        for (j = i; j <= N; j++)
        {
            idx = uMatIdx(i, j, N);
            U2[idx] = U1[i - 1][j - 1];
        }
    }
}
#endif

/**********************************************************************
* Function Name:    uMatxVect
*
* Description:
*    upper triangular matrix multiply vector, V2 = U*V1;
*
* Input:
*    U: upper triangular matrix, N*N
*    V1: vector array1, N*1, idx start from 0
*    V2: vector array2x, N*1, idx start from 0
*    N:  row and column length of matrix
*
* Return:
*
* Author: jpdeng
**********************************************************************/
void uMatxVect(DBL* U, DBL* V1, DBL* V2, S16 N)
{
    S16 i, j, idx;
    DBL tmp;

    for (i = 1; i <= N; i++)
    {
        idx = uMatIdx(i, i, N);
        tmp = 0;

        for (j = i; j <= N; j++)
        {
            tmp += U[idx] * V1[j - 1];
            idx++;
        }

        V2[i - 1] = tmp;
    }
}

/**********************************************************************
* Function Name:    getUdMatDiag
*
* Description:
*    get the diagonal components of U*D*D'
*
* Input:
*    U:     upper triangular matrix of UD decomposition, N*N
*    D:     diagonal matrix of UD decomposition, 1*N
*    diag:  diagonal components, 1*N
*    N:     row and column length of matrix
*
* Return:
*
* Author: jpdeng
**********************************************************************/
void getUdMatDiag(DBL* U, DBL* D, DBL* diag, S16 N)
{
    S16 i, k, idx;
    DBL tmp;

    for (i = 1; i <= N; i++)
    {
        tmp = 0;

        idx = uMatIdx(i, i, N);

        for (k = i; k <= N; k++)
        {
            tmp += U[idx] * D[k] * U[idx];
            idx = idx + 1;
        }

        diag[i - 1] = tmp;
    }
}

/**********************************************************************
* Function Name:    udKfPredict
*
* Description:
*    predict X(k|k-1) and P(k|k-1)
*
* Input:
*
* Return:
*
* Author: jpdeng
**********************************************************************/
void udKfPredict(kalmanInfo_t* const pKalmanInfo)
{
    S16 i, j, k, m, idx, idx1, idx_ij, idx_jj;
    DBL tmp, W1jk, WD;

    S16 N = pKalmanInfo->stateNum;
    S16 M = pKalmanInfo->matrixNum;
    /* KF state prediction, X(k|k-1) = A*X(k-1) */
    uMatxVect(pKalmanInfo->A, pKalmanInfo->X, pKalmanInfo->X, N);

    //memset(pKfCoreData->W1, 0, sizeof(DBL)*N_MAT);
    memset(pKalmanInfo->W2, 0, sizeof(DBL)*M);
    uMatInit(pKalmanInfo->W2, N);
    uxuMat(pKalmanInfo->A, pKalmanInfo->U_plus, pKalmanInfo->W1, N);

    /* P(k|k-1) = A*P(k-1)*A' + Q, calculate UD of P(k|k-1) */
    for (j = N; j >= 1; j--)
    {
        /* calculate D */
        tmp = 0;
        idx_jj = uMatIdx(j, j, N);
        idx = idx1 = idx_jj;

        for (k = j; k <= N; k++)
        {
            W1jk = pKalmanInfo->W1[idx];
            WD = W1jk * pKalmanInfo->D_plus[k];
            tmp += WD * W1jk;
            idx++;
        }

        for (k = j; k <= N; k++)
        {
            WD = 0;
            idx = idx_jj;

            for (m = j; m <= N; m++)
            {
                WD += pKalmanInfo->W2[idx] * pKalmanInfo->Q[m - 1][k - 1];
                idx++;
            }

            tmp += WD * pKalmanInfo->W2[idx1];
            idx1++;
        }

        pKalmanInfo->D_minus[j] = tmp;

        /* calculate U */
        for (i = 1; i <= j - 1; i++)
        {
            tmp = 0;
            idx_ij = uMatIdx(i, j, N);
            idx = idx_ij;
            idx1 = idx_jj;

            for (k = j; k <= N; k++)
            {
                WD = pKalmanInfo->W1[idx] * pKalmanInfo->D_plus[k];
                tmp += WD * pKalmanInfo->W1[idx1];
                idx++;
                idx1++;
            }

            idx1 = idx_jj;

            for (k = j; k <= N; k++)
            {
                WD = 0;
                idx = uMatIdx(i, i, N);

                for (m = i; m <= N; m++)
                {
                    WD += pKalmanInfo->W2[idx] * pKalmanInfo->Q[m - 1][k - 1];
                    idx++;
                }

                tmp += WD * pKalmanInfo->W2[idx1];
                idx1++;
            }

            idx = idx_ij;

            if (fabs(pKalmanInfo->D_minus[j]) > 1e-10)
            {
                pKalmanInfo->U_minus[idx] = tmp / pKalmanInfo->D_minus[j];
            }
            else
            {
                pKalmanInfo->U_minus[idx] = tmp * 1e10;
            }

            idx = idx_jj;
            idx1 = idx_ij;

            for (k = j; k <= N; k++)
            {
                pKalmanInfo->W1[idx1] -= pKalmanInfo->U_minus[idx_ij] * pKalmanInfo->W1[idx];
                pKalmanInfo->W2[idx1] -= pKalmanInfo->U_minus[idx_ij] * pKalmanInfo->W2[idx];
                idx++;
                idx1++;
            }
        }
    }//end for (j = N; j >= 1; j--)
}

/**********************************************************************
* Function Name:    udKFUpdate
*
* Description:
*    UD KF scalar update
*
* Input:
*     H:            measurement matrix, 1*N
*     R:            measurement noise variance, 1*1
*     ino:          innovation, 1*1
*     delatX:       state correction
*     varTh:        threshold of residual variance for update
*     flag:         1 -- save U,D results
*                   0 -- don't save U,D results
*
* Return:
*
*
* Author: jpdeng
**********************************************************************/
FLT udKFUpdate(kalmanInfo_t* const pKalmanInfo, DBL* H, DBL* deltaX, DBL R, DBL ino, DBL varTh, updateFlag_t flag)
{
    S16 i, j, k, idx;
    FLT normalizeRes;
    DBL p, a0, a1;
    
    S16 N = pKalmanInfo->stateNum;
    S16 M = pKalmanInfo->matrixNum;
    DBL *f = pKalmanInfo->f;
    DBL *v = pKalmanInfo->v;
    DBL *b = pKalmanInfo->b;


    /* Gain = P(k|k-1)*H'/[H*P(k|k-1)*H' + R], P(k) = (I - Gain*H)*P(k|k-1), calculate UD of P(k) */
    for (i = 1; i <= N; i++)
    {
        f[i] = 0;

        for (k = 1; k <= i; k++)
        {
            idx = uMatIdx(k, i, N);
            f[i] += H[k - 1] * pKalmanInfo->U_minus[idx];
        }
    }

    for (i = 1; i <= N; i++)
    {
        v[i] = pKalmanInfo->D_minus[i] * f[i];
    }

    a0 = R;

    for (i = 1; i <= N; i++)
    {
        a1 = a0 + f[i] * v[i];
        pKalmanInfo->D_plus[i] = pKalmanInfo->D_minus[i] * a0 / a1;
        b[i] = v[i];
        p = -f[i] / a0;
        a0 = a1;

        for (j = 1; j <= i - 1; j++)
        {
            idx = uMatIdx(j, i, N);
            pKalmanInfo->U_plus[idx] = pKalmanInfo->U_minus[idx] + b[j] * p;
            b[j] += pKalmanInfo->U_minus[idx] * v[i];
        }
    }

    /* b[] is gain */
    for (i = 1; i <= N; i++)
    {
        b[i] = b[i] / a1;
    }

    /* innovation error variance check */
    normalizeRes = (FLT)(fabs(ino) / sqrt(a1));

    if ((flag | UPDATE_SAVE) != 0)
    {
        if (normalizeRes > varTh)
        {
#ifndef ANDES
            printf("measurement cannot not pass KF normalizeRes check: %f\n", normalizeRes);
#endif
            return normalizeRes;
        }
    }

    /* calculate update */
    for (i = 1; i <= N; i++)
    {
        deltaX[i - 1] += (b[i] * ino);
    }

    if ((flag | UPDATE_SAVE) != 0)
    {
        /* prepare for next update */
        memcpy(pKalmanInfo->D_minus, pKalmanInfo->D_plus, sizeof(DBL) * (N + 1));
        memcpy(pKalmanInfo->U_minus, pKalmanInfo->U_plus, sizeof(DBL) * M);
    }

    return normalizeRes;
}