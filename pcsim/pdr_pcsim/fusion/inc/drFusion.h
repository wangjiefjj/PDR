#ifndef _DRFUSION_H_
#define _DRFUSION_H_

#include "kalman.h"
#include "types.h"
#include "misc.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct drFusionData
    {
        U32 utime;
        DBL fGnssLatitude;
        DBL fGnssLongitude;
        FLT fGnssHeading;
        DBL fPdrLatitude;
        DBL fPdrLongitude;
        FLT fPdrHeading;
    } drFusionData_t;

    U32 drKalmanInit(kalmanInfo_t* const pKalmanInfo);
    U32 drKalmanExec(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData);


#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif