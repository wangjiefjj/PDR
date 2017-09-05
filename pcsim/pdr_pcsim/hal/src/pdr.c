#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "pdr.h"
#include "magcal.h"
#include "ahrs.h"
#include "step.h"

#define     GYRO_BUFFER_LEN     (50)
#define     AVE_NUM             (5)
#define     RE                  (6378137.0)
#define     esqu                (0.00669437999013)
#define     RM(L)               (RE*(1-esqu)/pow((1-esqu*sin(L)*sin(L)),1.5))
#define     RN(L)               (RE/sqrt(1-esqu*sin(L)*sin(L)))

typedef struct pdrCtrl
{
    U32   uStaticFlag;
    U32   uHorizonAlignFlag;
    U32   uHeadingAlignFlag;
    U32   uPdrNavFlag;
} pdrCtrl_t;

typedef struct pdrInfo
{
    U32 uTime;
    DBL fLatitude;
    DBL fLongitude;
    DBL fAltitude;
} pdrInfo_t;

static magneticBuffer_t MagBuffer;
static magCalibration_t MagCalibration;
static ahrsFixData_t AhrsFixData;
static kalmanInfo_t AhrsKalmanInfo;
static stepInfo_t StepInfo;
static pdrCtrl_t PdrCtrl;
static pdrInfo_t PdrInfo;
static FLT GyroSmoothBuffer[GYRO_BUFFER_LEN][CHN];
extern FILE* FpAhrs;
extern FILE* FpStep;
extern FILE* FpOutput;

static void seDataProc(const sensorData_t* const pSensorData);
static void gnssDataProc(const gnssData_t* const pGnssData);
static void sensorDataCorrection(FLT gyro[], FLT acc[], FLT mag[]);
static U32 initialAlignment(const FLT facc[], const FLT fmag[]);
static void gyroSmooth(FLT fgyro[], FLT gyroBuffer[][CHN]);
static U32 staticDetect(const FLT gyro[], const FLT acc[]);
static void gyroCalibration(FLT gyroBias[]);

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

    if (ahrsInit(&AhrsFixData))
    {
        printf("ahrs init failed!\r\n");
        return -1;
    }

    if (ahrsKalmanInit(&AhrsKalmanInfo))
    {
        printf("ahrs kalman init failed!\r\n");
        return -1;
    }

    if (stepInit(&StepInfo))
    {
        printf("step init failed!\r\n");
        return -1;
    }

    memset(&PdrCtrl, 0, sizeof(pdrCtrl_t));

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
    U32 retval;
    FLT fgyro[CHN] = {0};
    FLT facc[CHN] = {0};
    FLT fmag[CHN] = {0};
    static U32 LoopCounter = 0;
    static U32 AccLoopCounter = 0;
    static FLT AccDetSum = 0.0;

    /* sensor data correction */
    utime = pSensorData->uTime;
    for (i = X; i <= Z; i++)
    {
        fgyro[i] = pSensorData->fGyro[i];
        facc[i] = pSensorData->fAcc[i];
        fmag[i] = pSensorData->fMag[i];
    }
    sensorDataCorrection(fgyro, facc, fmag);

    /* mag calibration */
    magBufferUpdate(&MagBuffer, pSensorData->fMag, fmag, LoopCounter);
    LoopCounter++;
    magCalibrationExec(&MagCalibration, &MagBuffer);
    if (MagCalibration.iValidMagCal != 0)
    {
        AhrsFixData.fB = MagCalibration.fB;
#ifdef DEBUG
        // indicate mag calibration is valid
        /*printf("mag %dparameters calibration is completed.\r\n", MagCalibration.iValidMagCal);
        printf("fit error: %f%%\r\n", MagCalibration.fFitErrorpc);
        printf("geomagnetic field magnitude: %fuT\r\n", MagCalibration.fB);
        printf("hard iron offset: %fuT, %fuT, %fuT\r\n", MagCalibration.fV[0], MagCalibration.fV[1], MagCalibration.fV[2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[0][0], MagCalibration.finvW[0][1], MagCalibration.finvW[0][2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[1][0], MagCalibration.finvW[1][1], MagCalibration.finvW[1][2]);
        printf("inverse soft iron matrix: %5.3f, %5.3f, %5.3f\r\n", MagCalibration.finvW[2][0], MagCalibration.finvW[2][1], MagCalibration.finvW[2][2]);*/
#endif
    }

    /* static detect */
    PdrCtrl.uStaticFlag = staticDetect(fgyro, facc);

    /* initial alignment */
    // note: dip angle is computed in initial alignment which is executed only once.
    //       so if work for a long time, the dip angle need to be recomputed.
    if (retval = initialAlignment(facc, fmag))
    {
        AhrsFixData.uTime = utime;
#ifdef DEBUG
        printf("initial alignment is completed in case %d.\r\n", retval);
#endif
        return;
    }

    /* start AHRS loop */
    if (PdrCtrl.uHorizonAlignFlag == 1)
    {
        // gyro data smooth
        gyroSmooth(fgyro, GyroSmoothBuffer);
        for (i = X; i <= Z; i++)
        {
            if (fabsf(fgyro[i]) < 0.1)
            {
                fgyro[i] = 0;
            }
        }
        quaternionIntegration(utime, fgyro, &AhrsFixData);
        if (PdrCtrl.uHeadingAlignFlag == 1 && MagCalibration.iValidMagCal != 0)
        {
            FLT magVector[3] = {0.0};
            FLT magEstimate[3] = {0.0};
            FLT magResidual[3] = {0.0};
            // check mag vector residual
            magVector[0] = AhrsFixData.fB * cosf(AhrsFixData.fDelta);
            magVector[2] = AhrsFixData.fB * sinf(AhrsFixData.fDelta);
            for (i = X; i <= Z; i++)
            {
                magEstimate[i] = AhrsFixData.fCbn[i][X] * fmag[X] + AhrsFixData.fCbn[i][Y] * fmag[Y] + AhrsFixData.fCbn[i][Z] * fmag[Z];
            }

            magResidual[X] = magVector[X] - magEstimate[X];
            magResidual[Y] = magVector[Y] - magEstimate[Y];
            magResidual[Z] = magVector[Z] - magEstimate[Z];

            if (fabs(magResidual[X]) > 5)
            {
#ifdef DEBUG
                //printf("mag calibration invalid in %dms\r\n", utime);
#endif
                MagCalibration.iValidMagCal = 0;
                ahrsKalmanExec(utime, facc, NULL, &AhrsKalmanInfo, &AhrsFixData);
            }
            else
            {
                ahrsKalmanExec(utime, facc, fmag, &AhrsKalmanInfo, &AhrsFixData);
            }
        }
        else
        {
            ahrsKalmanExec(utime, facc, NULL, &AhrsKalmanInfo, &AhrsFixData);
        }
        AhrsFixData.uTime = utime;
#ifdef DEBUG
        fprintf(FpAhrs, "%f, %f, %f\n", AhrsFixData.fPsiPl*RAD2DEG, AhrsFixData.fThePl*RAD2DEG, AhrsFixData.fPhiPl*RAD2DEG);
#endif
    }
    
    /* start dead reckoning loop */
    if (PdrCtrl.uPdrNavFlag == 1)
    {
        FLT accDet = sqrtf(facc[0] * facc[0] + facc[1] * facc[1] + facc[2] * facc[2]);
        AccLoopCounter ++;

        if (AccLoopCounter % AVE_NUM != 0)
        {
            AccDetSum += accDet;
        }
        else
        {
            // 10Hz sample rate is enough for step detect
            accDet = (AccDetSum + accDet) / AVE_NUM;
            AccDetSum = 0.0F;
#ifdef DEBUG
            fprintf(FpStep, "%d, %f\n", utime, accDet);
#endif
            if (stepDetection(utime, accDet, &StepInfo))
            {
                DBL fLatitude = PdrInfo.fLatitude;
                DBL fLongitude = PdrInfo.fLongitude;
                DBL fAltitude = PdrInfo.fAltitude;

                fLatitude += StepInfo.stepLength * cosf(AhrsFixData.fPsiPl) / RM(fLatitude);
                fLongitude += StepInfo.stepLength * sinf(AhrsFixData.fPsiPl) / RN(fLatitude) / cosf(fLatitude);

                PdrInfo.uTime = utime;
                PdrInfo.fLatitude = fLatitude;
                PdrInfo.fLongitude = fLongitude;
#ifdef DEBUG
                printf("%d step occur in %dms.\r\n", StepInfo.stepCount, StepInfo.preStepTime);
                fprintf(FpOutput, "%d, %.8f, %.8f, %.8f\n", utime, fLatitude*RAD2DEG, fLongitude*RAD2DEG, fAltitude);
#endif
            }
        }
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
    if (PdrCtrl.uPdrNavFlag != 1)
    {
        if (PdrCtrl.uHeadingAlignFlag == 1 && PdrCtrl.uHorizonAlignFlag == 1)
        {
            PdrCtrl.uPdrNavFlag = 1;
            PdrInfo.uTime = pGnssData->uTime;
            PdrInfo.fLatitude = pGnssData->fLatitude;
            PdrInfo.fLongitude = pGnssData->fLongitude;
            PdrInfo.fAltitude = pGnssData->fAltitude;
        }
    }
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

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return
  

 */
/*--------------------------------------------------------------------------*/
enum
{
    None = 0,
    Case1 = 1,
    Case2 = 2,
    Case3 = 3,
};
static U32 initialAlignment(const FLT facc[], const FLT fmag[])
{
    if (MagCalibration.iValidMagCal != 0)
    {
        // initial alignment case 3: 
        // horizon alignment is completed before.
        // heading alignment can be executed when calibration is completed.
        if (PdrCtrl.uHeadingAlignFlag == 0 && PdrCtrl.uHorizonAlignFlag == 1)
        {
            headingAlignment(fmag, &AhrsFixData);
            PdrCtrl.uHeadingAlignFlag = 1;

            return Case3;
        }
    }

    if (PdrCtrl.uStaticFlag == 1)
    {
        gyroCalibration(AhrsFixData.fGyroBias);
        if (MagCalibration.iValidMagCal != 0)
        {
            // initial alignment case 1: 
            // device is static and mag calibration is completed.
            // horizon alignment and heading alignment can be executed simultaneously. 
            if (PdrCtrl.uHeadingAlignFlag == 0 && PdrCtrl.uHorizonAlignFlag == 0)
            {
                compassAlignment(facc, fmag, &AhrsFixData);
                PdrCtrl.uHeadingAlignFlag = 1;
                PdrCtrl.uHorizonAlignFlag = 1;

                return Case1;
            }
        }
        else
        {
            // initial alignment case 2: 
            // device is static and mag calibration is not completed.
            // only horizon alignment can be executed.
            if (PdrCtrl.uHorizonAlignFlag == 0)
            {
                horizonAlignment(facc, &AhrsFixData);
                PdrCtrl.uHorizonAlignFlag = 1;

                return Case2;
            }
        }
    }

    return None;
}

static void gyroSmooth(FLT fgyro[], FLT gyroBuffer[][CHN])
{
    U32 i = 0;
    U32 j = 0;
    FLT gyroSum[CHN] = {0.0};
    static U32 ucount = 0;

    if (ucount < GYRO_BUFFER_LEN)
    {
        for (i = X; i <= Z; i++)
        {
            gyroBuffer[ucount][i] = fgyro[i];
        }
    }
    else
    {
        ucount = GYRO_BUFFER_LEN;
        for (i = 0; i < GYRO_BUFFER_LEN - 1; i++)
        {
            for (j = X; j <= Z; j++)
            {
                gyroBuffer[i][j] = gyroBuffer[i+1][j];
            }
        }
        for (j = X; j <= Z; j++)
        {
            gyroBuffer[GYRO_BUFFER_LEN - 1][j] = fgyro[j];
        }
    }
    if (ucount > 1)
    {
        for (i = 0; i < ucount - 1; i++)
        {
            for (j = X; j <= Z; j++)
            {
                gyroSum[j] += gyroBuffer[i][j];
            }
        }
        for (i = X; i <= Z; i++)
        {
            fgyro[i] = gyroSum[i] / ucount;
        }
    }
    ucount ++;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
#define     ACC_STATIC      (0.1)
#define     GYRO_STATIC     (0.01)
#define     ALIGN_NUM       (100)
static FLT AlignGyroArray[ALIGN_NUM][CHN] = {0};
static FLT AlignAccArray[ALIGN_NUM][CHN] = {0};

static U32 staticDetect(const FLT gyro[], const FLT acc[])
{
    FLT gyro_mean = 0;
    FLT gyro_std = 0;
    FLT acc_mean = 0;
    FLT acc_std = 0;
    U32 i = 0;
    U32 j = 0;
    static U32 uCount = 0;

    if (uCount < ALIGN_NUM)
    {
        for (i = X; i <= Z; i++)
        {
            AlignGyroArray[uCount][i] = gyro[i];
            AlignAccArray[uCount][i] = acc[i];
        }
    }
    else
    {
        for(i = 0; i < ALIGN_NUM - 1; i++)
        {
            for (j = X; j <= Z; j++)
            {
                AlignGyroArray[i][j] = AlignGyroArray[i+1][j];
                AlignAccArray[i][j] = AlignAccArray[i+1][j];
            }
        }
        for (i = X; i <= Z; i++)
        {
            AlignGyroArray[ALIGN_NUM-1][i] = gyro[i];
            AlignAccArray[ALIGN_NUM-1][i] = acc[i];
        }
    }

    uCount++;
    if (uCount >= ALIGN_NUM)
    {
        uCount = ALIGN_NUM;
        if (computeMeanStd(&gyro_mean, &gyro_std, AlignGyroArray, ALIGN_NUM))
        {
            return -1;
        }

        if (computeMeanStd(&acc_mean, &acc_std, AlignAccArray, ALIGN_NUM))
        {
            return -1;
        }

        if (acc_std < ACC_STATIC && gyro_std < GYRO_STATIC)
        {
            // indicate static condition
            return 1;
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
static void gyroCalibration(FLT gyroBias[])
{
    U32 i = 0;
    FLT fgyroSum[CHN] = {0.0F};

    for (i = 0; i < ALIGN_NUM; i++)
    {
        fgyroSum[X] += AlignGyroArray[i][X];
        fgyroSum[Y] += AlignGyroArray[i][Y];
        fgyroSum[Z] += AlignGyroArray[i][Z];
    }
    gyroBias[X] += (FLT)(fgyroSum[X] * 1.0 / ALIGN_NUM);
    gyroBias[Y] += (FLT)(fgyroSum[Y] * 1.0 / ALIGN_NUM);
    gyroBias[Z] += (FLT)(fgyroSum[Z] * 1.0 / ALIGN_NUM);

    // clear gyro buffer
    for (i = 0; i < ALIGN_NUM; i++)
    {
        AlignGyroArray[i][X] = 0;
        AlignGyroArray[i][Y] = 0;
        AlignGyroArray[i][Z] = 0;
    }
}
