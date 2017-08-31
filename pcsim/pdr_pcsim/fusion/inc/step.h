#ifndef _STEP_H_
#define _STEP_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct stepInfo
    {
        U32 utime;
        //FLT curStepDet;
        FLT preStepDet;
        //U32 curSlop;
        S32 preSlop;
        U32 preStepTime;
        U32 stepDeltaTime;
        U32 stepCount;
        FLT stepLength;
        U32 moveCount;
        FLT stepDetAverage;
        FLT stepThreshold;
    } stepInfo_t;

    U32 stepInit(stepInfo_t* const pStepInfo);
    U32 stepDetection(U32 utime, FLT stepDet, stepInfo_t* const pStepInfo);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif