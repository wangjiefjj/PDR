#ifndef _PDR_H_
#define _PDR_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define     CHX             (0)
#define     CHY             (1)
#define     CHZ             (2)
#define     CHN             (3)


#define     MAG_SUPPORT     (0)
#ifndef     GPS_CH_NUM
#define     GPS_CH_NUM      (14)
#endif
    
    enum
    {
        GNSS_FIX_NONE = 0,
        GNSS_FIX_2D = 1,
        GNSS_FIX_3D = 2,
    };

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
        U32 	uGnssFix;   /* gnss fix type(NONE(0),2D(1),3D(2)) */
        DBL     fLatitude;	// rad
        DBL     fLongitude;	// rad
        FLT     fAltitude;	// m
        FLT     fVelE;      // m/s2
        FLT     fVelN;      // m/s2
        FLT     fVelU;      // m/s2
        FLT     fDOP[4];                        /* HDOP, VDOP, PDOP, TDOP */
        U8      uSvUsedNum[4];                  /* used satellite's number (GPS,GLONASS,Beidou,GLL) */
        U32     uSvUsed[4];                     /* used satellite's flag (GPS,GLONASS,Beidou,GLL) */
        U16     uSvCno[4][GPS_CH_NUM];     		/* satellite's C/No (GPS,GLONASS,Beidou,GLL) */
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