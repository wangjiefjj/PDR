#ifndef _AHRS_H_
#define _AHRS_H_

#include "types.h"
#include "kalmanLite.h"
#include "misc.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     CHX               (0)
#define     CHY               (1)
#define     CHZ               (2)
#define     CHN               (3)

    
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
        FLT   fB;                       // current geomagnetic field magnitude (uT)
    } ahrsFixData_t;

    U32 ahrsInit(ahrsFixData_t* const pAhrsFixData);
    void gyroCorrection(FLT gyro[], const ahrsFixData_t* const pAhrsFixData);
    void accCorrection(FLT acc[], const ahrsFixData_t* const pAhrsFixData);
    U32 compassAlignment(const FLT acc[], const FLT mag[], ahrsFixData_t* const pAhrsFixData);
    U32 horizonAlignment(const FLT acc[], ahrsFixData_t* const pAhrsFixData);
    U32 headingAlignment(const FLT mag[], ahrsFixData_t* const pAhrsFixData);
    U32 quaternionIntegration (U32 utime, const FLT gyro[], ahrsFixData_t* const pAhrsFixData);
    U32 ahrsKalmanInit(kalmanInfo_t* const pKalmanInfo);
    U32 ahrsKalmanExec(U32 utime, const FLT acc[], const FLT mag[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif