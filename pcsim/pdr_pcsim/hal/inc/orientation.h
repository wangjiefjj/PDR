#ifndef _ORIENTATION_H_
#define _ORIENTATION_H_

#include "types.h"
#include "ahrs.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct pedestrianOrientation
{
    U32 uDeviceUpdateTime;
    FLT fDeviceHeadingRef;
    FLT fDeviceHeading;  
} pedestrianOrientation_t;

U32 getRelativeHeading(pedestrianOrientation_t* const pOrientation, FLT* const pHeading);
void updateReferenceOrientation(pedestrianOrientation_t* const pOrientation, const ahrsFixData_t* const pAhrsData);
void updateDeviceOrientation(pedestrianOrientation_t* const pOrientation, const ahrsFixData_t* const pAhrsData);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif

