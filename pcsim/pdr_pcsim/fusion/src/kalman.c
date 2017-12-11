#include <stdlib.h>
#include <memory.h>
#include "types.h"
#include "kalman.h"
#include "misc.h"


static U32 UDInit(DBL *pud, U32 stateNum, const DBL rms[]);

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 kalmanInit(kalmanInfo_t *pKalmanInfo, U32 uStateNum, const DBL rms[])
{
    U32 sizeofDBL = sizeof(DBL);

    pKalmanInfo->uStateNum = uStateNum;
    pKalmanInfo->uUdNum = uStateNum * (uStateNum + 1) / 2;

    // malloc state array
    pKalmanInfo->pStateX = (DBL *)malloc(sizeofDBL * pKalmanInfo->uStateNum);
    if (pKalmanInfo->pStateX == NULL)
    {
        return -1;
    }
    memset(pKalmanInfo->pStateX, 0, sizeofDBL * pKalmanInfo->uStateNum);
    
    // malloc ud array
    pKalmanInfo->pUd = (DBL *)malloc(sizeofDBL * pKalmanInfo->uUdNum);
    if (pKalmanInfo->pUd == NULL)
    {
        free(pKalmanInfo->pStateX);
        return -1;
    }
    memset(pKalmanInfo->pUd, 0, sizeofDBL * pKalmanInfo->uUdNum);

    // malloc q array
    pKalmanInfo->pQd = (DBL **)mallocArray2D_DBL(pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
    if (pKalmanInfo->pQd == NULL)
    {
        free(pKalmanInfo->pStateX);
        free(pKalmanInfo->pUd);
        return -1;
    }

    // malloc phi array
    pKalmanInfo->pPhim = (DBL **)mallocArray2D_DBL(pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
    if (pKalmanInfo->pPhim == NULL)
    {
        free(pKalmanInfo->pStateX);
        free(pKalmanInfo->pUd);
        freeArray2D_DBL(pKalmanInfo->pQd, pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
        return -1;
    }

    // ud init
    UDInit(pKalmanInfo->pUd, pKalmanInfo->uStateNum, rms);

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void predict(kalmanInfo_t* const pkalmanInfo)
{
    U32 i = 0;
    U32 j = 0;
    U32 stateNum = pkalmanInfo->uStateNum;
    U32 stateNum2 = 2 * stateNum;
    DBL *pstateVector = NULL;
    DBL **pw = NULL;
    DBL *pddw = NULL;

    // predict x
    pstateVector = (DBL *)malloc(sizeof(DBL)*stateNum);
    for (i = 0; i < stateNum; i++)
    {
        pstateVector[i] = 0.0;
        for (j = 0; j < stateNum; j++)
        {
            pstateVector[i] += pkalmanInfo->pPhim[i][j] * pkalmanInfo->pStateX[j];
        }
    }

    for (i = 0; i < stateNum; i++)
    {
        pkalmanInfo->pStateX[i] = pstateVector[i];
    }

    // predict p
    // p(k/k-1)=phi*P(k-1)*phiT+Q(k)
    pw = mallocArray2D_DBL(stateNum, stateNum2);
    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum2; j++)
        {
            pw[i][j] = 0.0;
        }
    }

    // Prepare input arrays for subroutine WGS. They are:
    //
    //     w  = [ phi*u(p) | u(qd) ]
    //
    // and a vector
    //
    //     ddw = [ d(p) | d(qd) ]
    //
    // where u(p) is the upper triangular matrix in udu decomposition
    // of covariance p (with 1's on diagonal), and d(p) is the diagonal, i.e.
    //
    //     p = u(p)d(p)u(p)'
    //
    // Similarly
    //
    //     qd = u(qd)d(qd)u(qd)'

    // Compute  phi*u(p) and store in w(1..n,1..n)
    multPhimUp(pkalmanInfo->pPhim, pkalmanInfo->pUd, stateNum, pw);
    // Store u(qd) in the upper triangle of w(1..n,n+1..2*n)
    storeUq(0, stateNum, pkalmanInfo->pQd, 1, stateNum + 1, pw);

    // Store d(p) in ddw(1..l4)
    pddw = (DBL *)malloc(sizeof(DBL)*stateNum*2);
    for (i = 0; i < stateNum2; i++)
    {
        pddw[i] = 0.0;
    }

    j = 0;
    for (i = 1; i <= stateNum; i++)
    {
        j = j + i;
        pddw[i - 1] = pkalmanInfo->pUd[j - 1];
    }

    // Store d(qd) in ddw(n+1..2*n)
    for (i = 1; i <= stateNum; i++)
    {
        pddw[stateNum + i - 1] = pkalmanInfo->pQd[i - 1][i - 1];
    }

    // Compute udu' factors of extrapolated n by n block of covariance
    udTimeUpdate(stateNum, stateNum2, pw, pddw, pkalmanInfo->pUd);
    
    // release temporary arrays
    free(pstateVector);
    free(pddw);
    freeArray2D_DBL(pw, stateNum, stateNum2);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static U32 UDInit(DBL *pud, U32 stateNum, const DBL rms[])
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 1;

    U32 len1 = stateNum * (stateNum + 1) / 2;
    U32 len2 = stateNum;

    for (i = 0; i < len1; i++)
    {
        if (i == j)
        {
            if (k > len2)
            {
                return -1;
            }
            pud[i] = rms[k-1] * rms[k-1];
            k++;
            j = k * (k + 1) / 2 - 1;
        }
        else
        {
            pud[i] = 0;
        }
    }

    return 0;
}