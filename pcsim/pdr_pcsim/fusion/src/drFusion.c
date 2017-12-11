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
static drFusionStatus_t errCorrection(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData);

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 drKalmanInit(kalmanInfo_t* const pKalmanInfo)
{
    U8 i;

    kalmanInit(pKalmanInfo, STATE_NUM);
    for (i = 0; i < STATE_NUM; i++)
    {
        pKalmanInfo->D_plus[i+1] = INIT_RMS[i] * INIT_RMS[i];
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
drFusionStatus_t drKalmanExec(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData)
{
    drFusionStatus_t retvel;

    setPhimQd(pKalmanInfo);
    udKfPredict(pKalmanInfo);
    gnssMeasUpdate(pKalmanInfo, pFusionData);
    retvel = errCorrection(pKalmanInfo, pFusionData);

    return retvel;
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
    pKalmanInfo->A[uMatIdx(1, 1, pKalmanInfo->stateNum)] = 1.0F;
    pKalmanInfo->A[uMatIdx(2, 2, pKalmanInfo->stateNum)] = 1.0F;
    pKalmanInfo->A[uMatIdx(3, 3, pKalmanInfo->stateNum)] = 1.0F;
    pKalmanInfo->Q[0][0] = SIGMA_LAT;
    pKalmanInfo->Q[1][1] = SIGMA_LON;
    pKalmanInfo->Q[2][2] = SIGMA_HEADING;
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
    DBL test = 0.0;
    DBL deltaX[STATE_NUM] = {0.0};

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

        test = udKFUpdate(pKalmanInfo, hc, deltaX, rc, zc, 5, UPDATE_SAVE);
    }

    for (i = 0; i < STATE_NUM; i++)
    {
        pKalmanInfo->X[i] += deltaX[i];
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
static drFusionStatus_t errCorrection(kalmanInfo_t* const pKalmanInfo, drFusionData_t* const pFusionData)
{
    drFusionStatus_t retvel = NoFix;

    pFusionData->fPdrLatitude = pFusionData->fPdrLatitude + pKalmanInfo->X[0]/RM(pFusionData->fGnssLatitude);
    pFusionData->fPdrLongitude = pFusionData->fPdrLongitude + pKalmanInfo->X[1]/RN(pFusionData->fGnssLatitude);
    pKalmanInfo->X[0] = 0.0F;
    pKalmanInfo->X[1] = 0.0F;
    retvel |= PosFix;

    if (fabs(pKalmanInfo->X[2]) < 10*DEG2RAD)
    {
        pFusionData->fPdrHeading = (FLT)(pFusionData->fPdrHeading + pKalmanInfo->X[2]);

        if (pFusionData->fPdrHeading > PI)
        {
            pFusionData->fPdrHeading -= 2*PI;
        }
        if (pFusionData->fPdrHeading < -PI)
        {
            pFusionData->fPdrHeading += 2*PI;
        }
        pKalmanInfo->X[2] = 0.0;
        retvel |= HeadingFix;
    }

    return retvel;
}