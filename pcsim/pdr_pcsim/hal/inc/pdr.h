#ifndef _PDR_H_
#define _PDR_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     X       (0)
#define     Y       (1)
#define     Z       (2)
#define     CHN     (3)

    typedef enum dataType
    {
        GNSS_DATA = 0,
        SENSOR_DATA = 1,
    } dataType_t;

	typedef	struct sensorData
	{
        U32     uTime;      // ms
        FLT     fGyro[CHN];	// rad/s
        FLT     fAcc[CHN];	// m/s2
        FLT     fMag[CHN];	// uT
	} sensorData_t;

    typedef	struct gnssData
    {
        U32     uTime;      // ms
        DBL     fLatitude;	// rad
        DBL     fLongitude;	// rad
        FLT     fAltitude;	// m
        FLT     fVelE;      // m/s2
        FLT     fVelN;      // m/s2
        FLT     fVelU;      // m/s2
        FLT     fHeading;   // rad
    } gnssData_t;

    typedef struct pdrData
    {
        dataType_t   dataType;
        sensorData_t sensorData;
        gnssData_t   gnssData;
    } pdrData_t;

    U32 pdrNavInit();
    U32 pdrNavExec(const pdrData_t* const pdrData);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif