#ifndef _KALMAN_H_
#define _KALMAN_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct kalmanInfo
    {
        U32 uStateNum;
        U32 uUdNum;
        DBL *pStateX;
        DBL *pUd;
        DBL **pQd;
        DBL **pPhim;
    } kalmanInfo_t;


    U32 kalmanInit(kalmanInfo_t *pkalmanInfo, U32 ustateNum);
    U32 UDInit(DBL *pud, U32 len1, const DBL rms[], U32 len2);
    void predict(kalmanInfo_t* const pkalmanInfo);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif