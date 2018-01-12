#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "pdr.h"
#include "magcal.h"
#include "ahrs.h"
#include "step.h"
#include "drFusion.h"
#include "orientation.h"

#define     GYRO_BUFFER_LEN     (50)
#define     AVE_NUM             (5)
#define     AHRS_INTERVAL       (10)

#define     ACC_STATIC          (0.1)
#define     GYRO_STATIC         (0.01)
#define     ACC_REST            (1.0)
#define     GYRO_REST           (0.1)
#define     ALIGN_NUM           (100)

enum
{
    UNSURE = 0,
    MOVE = 1,
    REST = 2,
    STATIC = 3,
};

enum
{
    None = 0,   // device alignment is not executed
    Case1 = 1,  // only device horizon alignment complete before pedestrian navigation start
    Case2 = 2,  // both device horizon/heading alignment complete before pedestrian navigation start
    Case3 = 3,  // device heading alignment complete after pedestrian navigation start
};

typedef struct pdrCtrl
{
    U32   uMotionFlag;
    U32   uDeviceHorizonAlignFlag;
    U32   uDeviceHeadingAlignFlag;
    U32   uPedestrianAlignFlag;
    U32   uPdrNavFlag;
} pdrCtrl_t;

typedef struct pdrInfo
{
    U32 uTime;
    DBL fLatitude;
    DBL fLongitude;
    FLT fAltitude;
    FLT fHeading;
} pdrInfo_t;

static magneticBuffer_t MagBuffer;
static magCalibration_t MagCalibration;
static ahrsFixData_t AhrsFixData;
static kalmanInfo_t AhrsKalmanInfo;
static kalmanInfo_t DrKalmanInfo;
static stepInfo_t StepInfo;
static pedestrianOrientation_t PdrOrientation;
static pdrCtrl_t PdrCtrl;
static pdrInfo_t PdrInfo;
static FLT GyroSmoothBuffer[GYRO_BUFFER_LEN][CHN];
static FLT AlignGyroArray[ALIGN_NUM][CHN] = {0};
static FLT AlignAccArray[ALIGN_NUM][CHN] = {0};

static void seDataProc(const sensorData_t* const pSensorData);
static void gnssDataProc(const gnssData_t* const pGnssData);
static void sensorDataCorrection(FLT gyro[], FLT acc[], FLT mag[]);
static U32 deviceAlignment(const FLT facc[], const FLT fmag[]);
static void gyroSmooth(FLT fgyro[], FLT gyroBuffer[][CHN]);
static U32 staticDetect(const FLT gyro[], const FLT acc[]);
static void gyroCalibration(FLT gyroBias[]);
static FLT magQualityControl(FLT mag[], const ahrsFixData_t* const pAhrsFixData);
static U32 pedestrianAlignment(const gnssData_t* const pGnssData, FLT *pHeading);

#ifndef ECOS

#ifdef DEBUG
    extern FILE* FpAhrs;
    extern FILE* FpStep;
    extern FILE* FpOutput;
#endif

#endif

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 pdrNavInit()
{
#if MAG_SUPPORT

    if (magCalibrationInit(&MagCalibration, &MagBuffer))
    {
        printf("mag calibration init failed!\r\n");
        return -1;
    }

#endif

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

    if (drKalmanInit(&DrKalmanInfo))
    {
        printf("dr kalman init failed!\r\n");
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
    FLT fgyro[CHN] = {0};  // Calibrated gyro data
    FLT facc[CHN] = {0};   // Calibrated acc data
    FLT fmag[CHN] = {0};   // Calibrated mag data
    static U32 LoopCounter = 0;
    static U32 AccLoopCounter = 0;
    static FLT AccDetSum = 0.0;
    static U8 AhrsLoop = 0;
    static FLT AccSum[CHN] = {0};
    static FLT MagSum[CHN] = {0};

    /* sensor data correction */
    utime = pSensorData->uTime;
    for (i = CHX; i <= CHZ; i++)
    {
        fgyro[i] = pSensorData->fGyro[i];
        facc[i] = pSensorData->fAcc[i];
        fmag[i] = pSensorData->fMag[i];
    }
    sensorDataCorrection(fgyro, facc, fmag);

#if MAG_SUPPORT

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

#endif

    /* static detect */
    PdrCtrl.uMotionFlag = staticDetect(pSensorData->fGyro, pSensorData->fAcc);

    /* Gyro calibration if device is strictly static */
    if (PdrCtrl.uMotionFlag == STATIC)
    {
        gyroCalibration(AhrsFixData.fGyroBias);
    }

    /* device initial alignment */
    // note: dip angle is computed in initial alignment which is executed only once.
    //       so if work for a long time, the dip angle need to be recomputed.
    if (retval = deviceAlignment(facc, fmag))
    {
        AhrsFixData.uTime = utime;
#ifdef DEBUG
        printf("initial alignment is completed in case %d.\r\n", retval);
#endif
        // Case3 happens when pedestrian has started,
        // so the following process need to be executed.
        if (retval != Case3)
        {
            return;
        }
    }

    /* start AHRS loop */
    if (PdrCtrl.uDeviceHorizonAlignFlag == 1)
    {
        // gyro data smooth
        gyroSmooth(fgyro, GyroSmoothBuffer);
        for (i = CHX; i <= CHZ; i++)
        {
            if (fabsf(fgyro[i]) < 0.1)
            {
                fgyro[i] = 0;
            }
        }

        // quaternion integration for AHRS
        quaternionIntegration(utime, fgyro, &AhrsFixData);

        // acc and mag aiding in 5 Hz
        if (AhrsLoop == AHRS_INTERVAL)
        {
            FLT fAccAverage[3];
            FLT fMagAverage[3];

            for (i = CHX; i <= CHZ; i++)
            {
                fAccAverage[i] = AccSum[i] / AHRS_INTERVAL;
                fMagAverage[i] = MagSum[i] / AHRS_INTERVAL;
            }

            if (PdrCtrl.uDeviceHeadingAlignFlag == 1 && MagCalibration.iValidMagCal != 0)
            {
                FLT ferror = 0.0F;

                // mag calibration quality check
                ferror = magQualityControl(fMagAverage, &AhrsFixData);
                if (fabs(ferror) > 5)
                {
#ifdef DEBUG
                    //printf("mag calibration invalid in %dms\r\n", utime);
#endif
                    MagCalibration.iValidMagCal = 0;
                    ahrsKalmanExec(utime, fAccAverage, NULL, &AhrsKalmanInfo, &AhrsFixData);
                }
                else
                {
                    ahrsKalmanExec(utime, fAccAverage, fMagAverage, &AhrsKalmanInfo, &AhrsFixData);
                }
            }
            else
            {
                ahrsKalmanExec(utime, fAccAverage, NULL, &AhrsKalmanInfo, &AhrsFixData);
            }

            // clear acc and mag sum value
            AhrsLoop = 0;
            for (i = CHX; i <= CHZ; i++)
            {
                AccSum[i] = 0;
                MagSum[i] = 0;
            } 
        }
        else
        {
            for (i = CHX; i <= CHZ; i++)
            {
                AccSum[i] += facc[i];
                MagSum[i] += fmag[i];
            }
            AhrsLoop++;
        }
        
        AhrsFixData.uTime = utime;
        updateDeviceOrientation(&PdrOrientation, &AhrsFixData);
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
                FLT fHeading = 0;
                FLT fRelativeHeading = 0;
                
                getRelativeHeading(&PdrOrientation, &fRelativeHeading);
                fHeading = PdrInfo.fHeading + fRelativeHeading;
                fHeading = fHeadingMod(fHeading);

                fLatitude += StepInfo.stepLength * cosf(fHeading) / RM(fLatitude);
                fLongitude += StepInfo.stepLength * sinf(fHeading) / RN(fLatitude) / cosf(fLatitude);

                PdrInfo.uTime = utime;
                PdrInfo.fLatitude = fLatitude;
                PdrInfo.fLongitude = fLongitude;
                PdrInfo.fHeading = fHeading;

                // update device reference heading
                updateReferenceOrientation(&PdrOrientation, &AhrsFixData);
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
    U32 utime = 0;
    drFusionData_t drFusionData;
    FLT gnssVel = 0.0F;

    if (pGnssData->uGnssFix != GNSS_FIX_NONE)
    {

    }
    else
    {
        // no gnss fix

        // return;
    }

    utime = pGnssData->uTime;
    memset(&drFusionData, 0, sizeof(drFusionData_t));

    if (PdrCtrl.uPdrNavFlag != 1)
    {
        if (PdrCtrl.uPedestrianAlignFlag != 1)
        {
            // heading alignment
            FLT heading = 0;

            if (pedestrianAlignment(pGnssData, &heading))
            {
                // pedestrian heading alignment succeed
                PdrCtrl.uPedestrianAlignFlag = 1;
                PdrInfo.fHeading = heading;
            }
        }

        // device alignment and pedestrian alignment succeed
        if (PdrCtrl.uDeviceHorizonAlignFlag == 1 &&
            PdrCtrl.uPedestrianAlignFlag == 1
            )
        {
            PdrCtrl.uPdrNavFlag = 1;
            PdrInfo.uTime = utime;
            PdrInfo.fLatitude = pGnssData->fLatitude;
            PdrInfo.fLongitude = pGnssData->fLongitude;
            PdrInfo.fAltitude = pGnssData->fAltitude;

            // update device heading reference
            updateReferenceOrientation(&PdrOrientation, &AhrsFixData);
        }

        return;
    }

    // gnss aiding process
    gnssVel = sqrtf(pGnssData->fVelE * pGnssData->fVelE + pGnssData->fVelN * pGnssData->fVelN + pGnssData->fVelU * pGnssData->fVelU);
    if (gnssVel > 1.5)
    {
        U32 status = NoFix;

        drFusionData.fGnssLatitude = pGnssData->fLatitude;
        drFusionData.fGnssLongitude = pGnssData->fLongitude;
        drFusionData.fGnssHeading = atan2f(pGnssData->fVelE, pGnssData->fVelN);
        drFusionData.fPdrLatitude = PdrInfo.fLatitude;
        drFusionData.fPdrLongitude = PdrInfo.fLongitude;
        drFusionData.fPdrHeading = PdrInfo.fHeading;
        drFusionData.fPdrFrequency = (FLT)(1000.0 / StepInfo.stepDeltaTime);

        status = drKalmanExec(utime, &DrKalmanInfo, &drFusionData);

        if ((status & PosFix) != 0)
        {
            drFusionData.utime = utime;
            // update pdr information
            PdrInfo.uTime = drFusionData.utime;
            PdrInfo.fLatitude = drFusionData.fPdrLatitude;
            PdrInfo.fLongitude = drFusionData.fPdrLongitude;
#ifdef DEBUG
            printf("gnss position aiding occur in %dms.\r\n", drFusionData.utime);
            fprintf(FpOutput, "%d, %.8f, %.8f, %.8f\n", PdrInfo.uTime, PdrInfo.fLatitude*RAD2DEG, PdrInfo.fLongitude*RAD2DEG, PdrInfo.fAltitude);
#endif
        }

        if ((status & HeadingFix) != 0)
        {
            drFusionData.utime = utime;
            PdrInfo.fHeading = drFusionData.fPdrHeading;
            // update device heading reference
            updateReferenceOrientation(&PdrOrientation, &AhrsFixData);
#ifdef DEBUG
            printf("gnss heading aiding occur in %dms.\r\n", drFusionData.utime);
#endif
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
static U32 deviceAlignment(const FLT facc[], const FLT fmag[])
{
    U32 retval = None;

    if (PdrCtrl.uMotionFlag >= REST)
    {
        /* device horizon alignment if device is rest */
        if (PdrCtrl.uDeviceHorizonAlignFlag == 0)
        {
            deviceHorizonAlignment(facc, &AhrsFixData);
            PdrCtrl.uDeviceHorizonAlignFlag = 1;
#if MAG_SUPPORT
            retval = Case1;
#else
            PdrCtrl.uDeviceHeadingAlignFlag = 1;
            retval = Case2;
#endif
        }
    }

    /* device heading alignment if:
       1. mag calibration is completed
       2. device horizon alignment is completed
    */
    if (PdrCtrl.uDeviceHeadingAlignFlag == 0)
    {
        if (PdrCtrl.uDeviceHorizonAlignFlag != 0 && MagCalibration.iValidMagCal != 0)
        {
            deviceHeadingAlignment(fmag, &AhrsFixData);
            PdrCtrl.uDeviceHeadingAlignFlag = 1;
            retval = Case2;

            /* Consider these 3 pedestrian alignment procedures:
                                    Proc1   Proc2   Proc3(*)
               Calibration:         1       2       3
               Device Horizon:      2       1       1
               Device Heading:      3       3       4
               Pedestrian Heading:  4       4       2
            */

            // process Procedure 3
            if (PdrCtrl.uPedestrianAlignFlag == 1)
            {
                // indicate pedestrian heading alignment is before device heading alignment
                // 1. feedback the old relative device heading to pedestrian heading
                // 2. update the new device heading reference
                FLT fRelativeHeading = 0;
                FLT fHeading = 0;

                getRelativeHeading(&PdrOrientation, &fRelativeHeading);
                fHeading = PdrInfo.fHeading + fRelativeHeading;
                fHeading = fHeadingMod(fHeading);
                PdrInfo.fHeading = fHeading;
                updateReferenceOrientation(&PdrOrientation, &AhrsFixData);

                retval = Case3;
            }
        }
    }

    return retval;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return
  

 */
/*--------------------------------------------------------------------------*/
#define HEADING_NUM 5

static U32 pedestrianAlignment(const gnssData_t* const pGnssData, FLT *pHeading)
{
    U8 i = 0;
    FLT heading = 0;
    FLT velocity = 0;
    FLT headingAverage = 0;
    FLT headingStd = 0;
    static FLT HeadingBuffer[HEADING_NUM];
    static U8 HeadingCount = 0;

    velocity = sqrtf(pGnssData->fVelE * pGnssData->fVelE + pGnssData->fVelN * pGnssData->fVelN);

    if (velocity > 1.0)
    {
        heading = atan2f(pGnssData->fVelE, pGnssData->fVelN);

        if (HeadingCount < HEADING_NUM)
        {
            HeadingBuffer[HeadingCount] = heading;
            HeadingCount++;
        }
        else
        {
            for (i = 0; i < HEADING_NUM - 1; i++)
            {
                HeadingBuffer[i] = HeadingBuffer[i+1];
            }
            HeadingBuffer[i] = heading;
        }

        if (HeadingCount == HEADING_NUM)
        {
            for (i = 0; i < HEADING_NUM; i++)
            {
                headingAverage += HeadingBuffer[i];
            }
            headingAverage = (FLT)(headingAverage * 1.0 / HEADING_NUM);

            for (i = 0; i < HEADING_NUM; i++)
            {
                headingStd += (HeadingBuffer[i] - headingAverage) * (HeadingBuffer[i] - headingAverage);
            }

            headingStd = (FLT)(headingStd * 1.0 / (HEADING_NUM - 1));

            if (headingStd < 0.1)
            {
                // pedestrian heading alignment complete
                *pHeading = headingAverage;

                return Case1;
            }
        }
    }

    return None;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return
  

 */
/*--------------------------------------------------------------------------*/
static void gyroSmooth(FLT fgyro[], FLT gyroBuffer[][CHN])
{
    U32 i = 0;
    U32 j = 0;
    FLT gyroSum[CHN] = {0.0};
    static U32 ucount = 0;

    if (ucount < GYRO_BUFFER_LEN)
    {
        for (i = CHX; i <= CHZ; i++)
        {
            gyroBuffer[ucount][i] = fgyro[i];
        }
    }
    else
    {
        ucount = GYRO_BUFFER_LEN;
        for (i = 0; i < GYRO_BUFFER_LEN - 1; i++)
        {
            for (j = CHX; j <= CHZ; j++)
            {
                gyroBuffer[i][j] = gyroBuffer[i+1][j];
            }
        }
        for (j = CHX; j <= CHZ; j++)
        {
            gyroBuffer[GYRO_BUFFER_LEN - 1][j] = fgyro[j];
        }
    }
    if (ucount > 1)
    {
        for (i = 0; i < ucount - 1; i++)
        {
            for (j = CHX; j <= CHZ; j++)
            {
                gyroSum[j] += gyroBuffer[i][j];
            }
        }
        for (i = CHX; i <= CHZ; i++)
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
        for (i = CHX; i <= CHZ; i++)
        {
            AlignGyroArray[uCount][i] = gyro[i];
            AlignAccArray[uCount][i] = acc[i];
        }
    }
    else
    {
        for(i = 0; i < ALIGN_NUM - 1; i++)
        {
            for (j = CHX; j <= CHZ; j++)
            {
                AlignGyroArray[i][j] = AlignGyroArray[i+1][j];
                AlignAccArray[i][j] = AlignAccArray[i+1][j];
            }
        }
        for (i = CHX; i <= CHZ; i++)
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
            return STATIC;
        }
        else if (acc_std < ACC_REST && gyro_std < GYRO_REST)
        {
            // indicate rest condition
            return REST;
        }
        else
        {
            // indicate move condition
            return MOVE;
        }
    }

    // indicate detecting
    return UNSURE;
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
        fgyroSum[CHX] += AlignGyroArray[i][CHX];
        fgyroSum[CHY] += AlignGyroArray[i][CHY];
        fgyroSum[CHZ] += AlignGyroArray[i][CHZ];
    }
    gyroBias[CHX] = (FLT)(fgyroSum[CHX] * 1.0 / ALIGN_NUM);
    gyroBias[CHY] = (FLT)(fgyroSum[CHY] * 1.0 / ALIGN_NUM);
    gyroBias[CHZ] = (FLT)(fgyroSum[CHZ] * 1.0 / ALIGN_NUM);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static FLT magQualityControl(FLT fmag[], const ahrsFixData_t* const pAhrsFixData)
{
    U32 i;
    FLT ferror = 0.0F;
    FLT magVector[3] = {0.0};
    FLT magEstimate[3] = {0.0};
    FLT magResidual[3] = {0.0};

    // check mag vector residual
    magVector[0] = AhrsFixData.fB * cosf(AhrsFixData.fDelta);
    magVector[2] = AhrsFixData.fB * sinf(AhrsFixData.fDelta);
    for (i = CHX; i <= CHZ; i++)
    {
        magEstimate[i] = AhrsFixData.fCbn[i][CHX] * fmag[CHX] + AhrsFixData.fCbn[i][CHY] * fmag[CHY] + AhrsFixData.fCbn[i][CHZ] * fmag[CHZ];
    }

    magResidual[CHX] = magVector[CHX] - magEstimate[CHX];
    magResidual[CHY] = magVector[CHY] - magEstimate[CHY];
    magResidual[CHZ] = magVector[CHZ] - magEstimate[CHZ];

    ferror = magResidual[CHX];

    return ferror;
}
