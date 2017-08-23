#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "types.h"
#include "pdr.h"
#include "magcal.h"
#include "ahrs.h"


static magneticBuffer_t MagBuffer;
static magCalibration_t MagCalibration;
static ahrsFixData_t AhrsFixData;

static void seDataProc(const sensorData_t* const pSensorData);
static void gnssDataProc(const gnssData_t* const pGnssData);
static void sensorDataCorrection(FLT gyro[], FLT acc[], FLT mag[]);

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 pdrNavInit()
{
    if (magCalibrationInit(&MagCalibration, &MagBuffer))
    {
        printf("mag calibration init failed!\r\n");
        return -1;
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
U32 pdrNavExec(const pdrData_t* const pdrData)
{
    dataType_t dataType = pdrData->dataType;

    switch (dataType)
    {
        case SENSOR_DATA:
            seDataProc(&pdrData->sensorData);
            break;
        case GNSS_DATA:
            gnssDataProc(&pdrData->gnssData);
            break;
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
static void seDataProc(const sensorData_t* const pSensorData)
{
    U32 utime = 0;
    U32 i = 0;
    FLT fgyro[CHN] = {0};
    FLT facc[CHN] = {0};
    FLT fmag[CHN] = {0};
    static U32 LoopCounter = 0;

    utime = pSensorData->uTime;
    for (i = X; i <= Z; i++)
    {
        fgyro[i] = pSensorData->fGyro[i];
        facc[i] = pSensorData->fAcc[i];
        fmag[i] = pSensorData->fMag[i];
    }
    sensorDataCorrection(fgyro, facc, fmag);
    magBufferUpdate(&MagBuffer, pSensorData->fMag, fmag, LoopCounter);
    LoopCounter++;
    magCalibrationExec(&MagCalibration, &MagBuffer);
    if (MagCalibration.iValidMagCal != 0)
    {
#ifdef DEBUG
        // indicate mag calibration is valid
        printf("mag %dparameters calibration is completed.\r\n", MagCalibration.iValidMagCal);
        printf("fit error: %f%%\r\n", MagCalibration.fFitErrorpc);
        printf("geomagnetic field magnitude: %fuT\r\n", MagCalibration.fB);
        printf("hard iron offset: %fuT, %fuT, %fuT\r\n", MagCalibration.fV[0], MagCalibration.fV[1], MagCalibration.fV[2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[0][0], MagCalibration.finvW[0][1], MagCalibration.finvW[0][2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[1][0], MagCalibration.finvW[1][1], MagCalibration.finvW[1][2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[2][0], MagCalibration.finvW[2][1], MagCalibration.finvW[2][2]);
#endif
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void gnssDataProc(const gnssData_t* const pGnssData)
{
    ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return
  

 */
/*--------------------------------------------------------------------------*/
static void sensorDataCorrection(FLT gyro[], FLT acc[], FLT mag[])
{
    // gyro data correction
    gyroCorrection(gyro, &AhrsFixData);
    
    // acc data correction
    accCorrection(acc, &AhrsFixData);
    
    // mag data correction
    magCorrection(mag, &MagCalibration);
}