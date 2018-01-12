#include <memory.h>
#include "step.h"

#define STEP_DEFAULT    (0.8F)
#define AVERAGE_NUM     (10)
#define LOCKTIME        (300)
#define GRAVITY         ((FLT)9.80665)
/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 stepInit(stepInfo_t* const pStepInfo)
{
    memset(pStepInfo, 0, sizeof(stepInfo_t));
    pStepInfo->stepLength = STEP_DEFAULT;

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 stepDetection(U32 utime, FLT stepDet, stepInfo_t* const pStepInfo)
{
    U32 retval = 0;
    S32 slop = 0;
    U32 deltaStepTime = 0;
    U32 lastLoopTime = pStepInfo->utime;
    
    pStepInfo->utime = utime;

    // det move average calculate for threshold estimate and the window length is 10 (1s)
    if (pStepInfo->moveCount < AVERAGE_NUM)
    {
        pStepInfo->moveCount ++;
        pStepInfo->stepDetAverage = (pStepInfo->stepDetAverage * (pStepInfo->moveCount - 1) + stepDet) / pStepInfo->moveCount;
    }
    else
    {
        pStepInfo->moveCount = AVERAGE_NUM;
        pStepInfo->stepDetAverage = (pStepInfo->stepDetAverage * (AVERAGE_NUM - 1) + stepDet) / AVERAGE_NUM;
    }

    // determine the threshold
    if (pStepInfo->stepDetAverage - GRAVITY > 0.8)
    {
        pStepInfo->stepThreshold = (FLT)(0.15 * pStepInfo->stepDetAverage);
    }
    else
    {
        pStepInfo->stepDetAverage = GRAVITY;
        pStepInfo->stepThreshold = (FLT)(0.08 * GRAVITY);
    }

    // slop calculate
    if (stepDet - pStepInfo->preStepDet > 0)
    {
        slop = 1;
    }
    else
    {
        slop = -1;
    }

    // delta step time calculate
    if (utime < pStepInfo->preStepTime)
    {
        deltaStepTime = utime + 0xFFFFFFFF - pStepInfo->preStepTime;
    }
    else
    {
        deltaStepTime = utime - pStepInfo->preStepTime;
    }

    // step detect (wave summit & delta time & step det)
    if (pStepInfo->preSlop > 0 && slop < 0 && deltaStepTime > LOCKTIME && (pStepInfo->preStepDet - pStepInfo->stepDetAverage) > pStepInfo->stepThreshold)
    {
        pStepInfo->stepCount ++;
        pStepInfo->preStepTime = lastLoopTime;
        if (pStepInfo->stepCount == 1)
        {
            pStepInfo->stepDeltaTime = 500;
        }
        else
        {
            pStepInfo->stepDeltaTime = deltaStepTime;
        }
        

        retval = 1;
    }

    pStepInfo->preSlop = slop;
    pStepInfo->preStepDet = stepDet;

    return retval;
}