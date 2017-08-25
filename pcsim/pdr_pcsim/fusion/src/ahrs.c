#include <math.h>
#include <memory.h>
#include "ahrs.h"

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 ahrsInit(ahrsFixData_t* const pAhrsFixData)
{
    memset(pAhrsFixData, 0, sizeof(ahrsFixData_t));
    f3x3matrixEqI(pAhrsFixData->fCbn);
    f3x3matrixEqI(pAhrsFixData->fCnb);

    return 0;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void gyroCorrection(FLT gyro[], const ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;

    for (i = X; i <= Z; i++)
    {
        gyro[i] -= pAhrsFixData->fGyroBias[i];
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void accCorrection(FLT acc[], const ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;

    for (i = X; i <= Z; i++)
    {
        acc[i] -= pAhrsFixData->fAccBias[i];
    }
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

U32 staticDetect(const FLT gyro[], const FLT acc[])
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
void gyroCalibration(FLT gyroBias[])
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

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 compassAlignment(const FLT acc[], const FLT mag[], ahrsFixData_t* const pAhrsFixData)
{
    FLT fg[CHN] = {0.0F};
    FLT fmag[CHN] = {0.0F};
    FLT fR[3][3] = {0.0F};
    FLT fmod[3] = {0.0F};
    FLT ftmp = 0.0F;
    FLT fmodMag = 0.0F;
    FLT fGdotMag = 0.0F;
    FLT fyaw = 0.0F;
    FLT fpitch = 0.0F;
    FLT froll = 0.0F;
    FLT fq[4] = {0.0F};
    U32 i, j, retval;

    retval = 0;
    for (i = X; i <= Z; i++)
    {
        fg[i] = -acc[i];
        fmag[i] = mag[i];
        fR[i][Z] = fg[i];
        fR[i][X] = fmag[i];
    }

    // set y vector to vector product of z and x vectors
    fR[X][Y] = fR[Y][Z] * fR[Z][X] - fR[Z][Z] * fR[Y][X];
    fR[Y][Y] = fR[Z][Z] * fR[X][X] - fR[X][Z] * fR[Z][X];
    fR[Z][Y] = fR[X][Z] * fR[Y][X] - fR[Y][Z] * fR[X][X];

    // set x vector to vector product of y and z vectors
    fR[X][X] = fR[Y][Y] * fR[Z][Z] - fR[Z][Y] * fR[Y][Z];
    fR[Y][X] = fR[Z][Y] * fR[X][Z] - fR[X][Y] * fR[Z][Z];
    fR[Z][X] = fR[X][Y] * fR[Y][Z] - fR[Y][Y] * fR[X][Z];

    for (i = X; i <= Z; i++)
    {
        fmod[i] = sqrtf(fR[X][i] * fR[X][i] + fR[Y][i] * fR[Y][i] + fR[Z][i] * fR[Z][i]);
    }

    if (!((fmod[X] < 1e-5) || (fmod[Y] < 1e-5) || (fmod[Z] < 1e-5)))
    {
        for (j = X; j <= Z; j++)
        {
            ftmp = 1.0F / fmod[j];
            for (i = X; i <= Z; i++)
            {
                fR[i][j] *= ftmp;
            }
        }
    }
    else
    {
        // no solution is possible
        retval = -1;
    }
    memcpy(pAhrsFixData->fCnb, fR, sizeof(pAhrsFixData->fCnb));
    memcpy(pAhrsFixData->fCbn, fR, sizeof(pAhrsFixData->fCbn));
    f3x3matrixTranspose(pAhrsFixData->fCbn);

    fmodMag = sqrtf(fmag[X] * fmag[X] + fmag[Y] * fmag[Y] + fmag[Z] * fmag[Z]);
    fGdotMag = fg[X] * fmag[X] + fg[Y] * fmag[Y] + fg[Z] * fmag[Z];
    if (!((fmod[Z] == 0.0F) || (fmodMag == 0.0F)))
    {
        pAhrsFixData->fDelta = asinf(fGdotMag / (fmod[Z] * fmodMag));
    }
    dcm2euler(pAhrsFixData->fCbn, &fyaw, &fpitch, &froll);
    pAhrsFixData->fPsiPl = fyaw;
    pAhrsFixData->fThePl = fpitch;
    pAhrsFixData->fPhiPl = froll;

    euler2q(fq, fyaw, fpitch, froll);
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];

    horizonAlignment(acc, pAhrsFixData);

    return retval;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 horizonAlignment(const FLT acc[], ahrsFixData_t* const pAhrsFixData)
{
    FLT fmodG = 0.0F;
    FLT fq[4] = {0.0F};

    fmodG = sqrtf(acc[X] * acc[X] + acc[Y] * acc[Y] + acc[Z] * acc[Z]);
    if (fmodG < 1e-5)
    {
        // no solution is possible
        return -1;
    }
    pAhrsFixData->fPsiPl = 0.0F;
    pAhrsFixData->fThePl = -asinf(acc[X]/(-fmodG));
    pAhrsFixData->fPhiPl = atan2f(acc[Y]/(-fmodG), acc[Z]/(-fmodG));

    euler2q(fq, pAhrsFixData->fPsiPl, pAhrsFixData->fThePl, pAhrsFixData->fPhiPl);
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];

    q2dcm(fq, pAhrsFixData->fCbn);
    memcpy(pAhrsFixData->fCnb, pAhrsFixData->fCbn, sizeof(pAhrsFixData->fCnb));
    f3x3matrixTranspose(pAhrsFixData->fCnb);

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 headingAlignment(const FLT mag[], ahrsFixData_t* const pAhrsFixData)
{
    FLT mHorizon[CHN] = {0.0F};
    FLT mVector[CHN] = {0.0F};
    FLT fq[4] = {0.0F};

    mHorizon[X] = cosf(pAhrsFixData->fThePl) * mag[X] + sinf(pAhrsFixData->fPhiPl) * sinf(pAhrsFixData->fThePl) * mag[Y]
                + cosf(pAhrsFixData->fPhiPl) * sinf(pAhrsFixData->fThePl) * mag[Z];
    mHorizon[Y] = cosf(pAhrsFixData->fPhiPl) * mag[Y] - sinf(pAhrsFixData->fPhiPl) * mag[Z];
    mHorizon[Z] = -sinf(pAhrsFixData->fThePl) * mag[X] + sinf(pAhrsFixData->fPhiPl) * cosf(pAhrsFixData->fThePl) * mag[Y]
                + cosf(pAhrsFixData->fPhiPl) * cosf(pAhrsFixData->fThePl) * mag[Z];

    pAhrsFixData->fPsiPl = atan2f(-mHorizon[Y], mHorizon[X]);

    euler2q(fq, pAhrsFixData->fPsiPl, pAhrsFixData->fThePl, pAhrsFixData->fPhiPl);
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];

    q2dcm(fq, pAhrsFixData->fCbn);
    memcpy(pAhrsFixData->fCnb, pAhrsFixData->fCbn, sizeof(pAhrsFixData->fCnb));
    f3x3matrixTranspose(pAhrsFixData->fCnb);

    // estimate the dip angle
    mVector[X] = cosf(pAhrsFixData->fPsiPl) * mHorizon[X] - sinf(pAhrsFixData->fPsiPl) * mHorizon[Y];
    //mVector[Y] = sinf(pAhrsFixData->fPsiPl) * mHorizon[X] + cosf(pAhrsFixData->fPsiPl) * mHorizon[Y];
    mVector[Z] = mHorizon[Z];
    pAhrsFixData->fDelta = atanf(mVector[Z]/mVector[X]);
    
    return 0;
}
