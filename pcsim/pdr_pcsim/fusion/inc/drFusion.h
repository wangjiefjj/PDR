#ifndef _DRFUSION_H_
#define _DRFUSION_H_

#include "types.h"
#include "kalmanLite.h"
#include "misc.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct drFusionData
    {
        U32 utime;              // dr fix time
        DBL fGnssLatitude;      // input measurement
        DBL fGnssLongitude;     // input measurement
        FLT fGnssHeading;       // input measurement
        DBL fPdrLatitude;
        DBL fPdrLongitude;
        FLT fPdrHeading;
        FLT fPdrFrequency;
        FLT fPdrStepLength;
    } drFusionData_t;

    typedef enum drFusionStatus
    {
        NoFix = 0,
        PosFix = 1<<1,
        HeadingFix = 1<<2,
    } drFusionStatus_t;

    U32 drKalmanInit(kalmanInfo_t* const pKalmanInfo);
    drFusionStatus_t drKalmanExec(U32 utime, kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData);


#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif