#include <math.h>
#include <memory.h>
#include "drFusion.h"

#define     STATE_NUM               (3)
#define     UD_NUM                  (STATE_NUM*(STATE_NUM+1)/2)
#define     MEAS_GNSS_NUM           (3)
#define     SIG_LAT                 (0.1)                       /* rms of pitch and latitude range (m) */
#define     SIG_LON                 (0.1)                       /* rms of pitch and longitude range (m) */
#define     SIG_HEADING             (10.0*PI/180)                /* rms of pitch and heading (rad) */
#define     SIGMA_LAT               (0.1 * 0.1)
#define     SIGMA_LON               (0.1 * 0.1)
#define     SIGMA_HEADING           (5.0*PI/180 * 5.0*PI/180)

static const DBL INIT_RMS[] = {SIG_LAT, SIG_LON, SIG_HEADING};

static void setPhimQd(kalmanInfo_t* const pKalmanInfo);
static U32 gnssMeasUpdate(kalmanInfo_t* const pKalmanInfo, const drFusionData_t* const pFusionData);
static void errCorrection(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData);

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 drKalmanInit(kalmanInfo_t* const pKalmanInfo)
{
    kalmanInit(pKalmanInfo, STATE_NUM, INIT_RMS);

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 drKalmanExec(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData)
{
    setPhimQd(pKalmanInfo);
    predict(pKalmanInfo);
    gnssMeasUpdate(pKalmanInfo, pFusionData);
    errCorrection(pKalmanInfo, pFusionData);

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void setPhimQd(kalmanInfo_t* const pKalmanInfo)
{
    pKalmanInfo->pPhim[0][0] = 1.0F;
    pKalmanInfo->pPhim[1][1] = 1.0F;
    pKalmanInfo->pPhim[2][2] = 1.0F;
    pKalmanInfo->pQd[0][0] = SIGMA_LAT;
    pKalmanInfo->pQd[1][1] = SIGMA_LON;
    pKalmanInfo->pQd[2][2] = SIGMA_HEADING;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static U32 gnssMeasUpdate(kalmanInfo_t* const pKalmanInfo, const drFusionData_t* const pFusionData)
{
    U32 i = 0;
    U32 j = 0;
    DBL zc = 0.0;
    DBL rc = 0.0;
    DBL hc[STATE_NUM] = {0.0};
    DBL z[MEAS_GNSS_NUM] = {0.0};
    DBL h[MEAS_GNSS_NUM][STATE_NUM] = {0.0};
    DBL r[MEAS_GNSS_NUM] = {0.0};
    DBL xSave[STATE_NUM];
    DBL udSave[UD_NUM];
    DBL ion = 0.0;
    DBL res = 0.0;
    DBL test = 0.0;

    DBL gnssLatitude = pFusionData->fGnssLatitude;
    DBL gnssLongitude = pFusionData->fGnssLongitude;
    FLT gnssHeading = pFusionData->fGnssHeading;
    DBL pdrLatitude = pFusionData->fPdrLatitude;
    DBL pdrLongitude = pFusionData->fPdrLongitude;
    FLT pdrHeading = pFusionData->fPdrHeading;

    z[0] = (gnssLatitude - pdrLatitude) * RM(gnssLatitude);
    z[1] = (gnssLongitude - pdrLongitude) * RN(gnssLongitude);

    if ((gnssHeading - pdrHeading) > PI)
    {
        gnssHeading -= 2*PI;
    }
    if ((gnssHeading - pdrHeading) < -PI)
    {
        gnssHeading += 2*PI;
    }
    z[2] = gnssHeading - pdrHeading;

    h[0][0] = 1.0F;
    h[1][1] = 1.0F;
    h[2][2] = 1.0F;

    r[0] = 20 * 20;
    r[1] = 20 * 20;
    r[2] = 10 * DEG2RAD * 10 * DEG2RAD;

    for (i = 0; i < MEAS_GNSS_NUM; i++)
    {
        zc = z[i];
        rc = r[i];

        for (j = 0; j < STATE_NUM; j++)
        {
            hc[j] = h[i][j];
        }

        // save x,p in case the measurement is rejected
        memcpy(xSave, pKalmanInfo->pStateX, sizeof(xSave));
        memcpy(udSave, pKalmanInfo->pUd, sizeof(udSave));

        // scalar measurement update
        udMeasUpdate(pKalmanInfo->pUd, pKalmanInfo->pStateX, pKalmanInfo->uStateNum, rc, hc, zc, &ion, &res);
        test = fabs(res) / sqrt(ion);

        // reject this measurement
        // 1. innovation test > 5, generally it is around 3.24
        if (test > 5)
        {
            memcpy(pKalmanInfo->pStateX, xSave, sizeof(xSave));
            memcpy(pKalmanInfo->pUd, udSave, sizeof(udSave));
        }
    }

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void errCorrection(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData)
{
    pFusionData->fPdrLatitude = pFusionData->fPdrLatitude + pKalmanInfo->pStateX[0]/RM(pFusionData->fGnssLatitude);
    pFusionData->fPdrLongitude = pFusionData->fPdrLongitude + pKalmanInfo->pStateX[1]/RN(pFusionData->fGnssLatitude);
    pKalmanInfo->pStateX[0] = 0.0F;
    pKalmanInfo->pStateX[1] = 0.0F;

    if (fabs(pKalmanInfo->pStateX[2]) < 10*DEG2RAD)
    {
        pFusionData->fPdrHeading = (FLT)(pFusionData->fPdrHeading + pKalmanInfo->pStateX[2]);

        if (pFusionData->fPdrHeading > PI)
        {
            pFusionData->fPdrHeading -= 2*PI;
        }
        if (pFusionData->fPdrHeading < -PI)
        {
            pFusionData->fPdrHeading += 2*PI;
        }
        pKalmanInfo->pStateX[2] = 0.0;
    }
}