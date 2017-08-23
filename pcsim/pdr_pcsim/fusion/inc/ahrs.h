#ifndef _AHRS_H_
#define _AHRS_H_

#include "types.h"
#include "misc.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     X               (0)
#define     Y               (1)
#define     Z               (2)
#define     CHN             (3)

    
    typedef struct quaternion
    {
        FLT  q0;	// scalar component
        FLT  q1;	// x vector component
        FLT  q2;	// y vector component
        FLT  q3;	// z vector component
    } quaternion_t;
	
    typedef struct ahrsFixData
    {
        U32   uTime;                    // time tag (ms)
        FLT   fPsiPl;					// yaw (deg)
        FLT   fThePl;					// pitch (deg)
        FLT   fPhiPl;					// roll (deg)    
        FLT   fCnb[3][3];               // a posteriori orientation matrix (Cnb)
        FLT   fCbn[3][3];               // a posteriori orientation matrix (Cbn)
        quaternion_t fqPl;              // a posteriori orientation quaternion
        FLT   fGyroBias[CHN];           // gyro bias (rad/s)
        FLT   fAccBias[CHN];            // acc bias (m/s2)
        FLT   fDelta;                   // inclination angle (rad)
    } ahrsFixData_t;

    void gyroCorrection(FLT gyro[], const ahrsFixData_t* const pAhrsFixData);
    void accCorrection(FLT acc[], const ahrsFixData_t* const pAhrsFixData);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif