#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "ahrs.h"

#define     STATE_NUM               (9)
#define     UD_NUM                  (STATE_NUM*(STATE_NUM+1)/2)
#define     MEAS_ACC_NUM            (3)
#define     MEAS_MAG_NUM            (3)

#define     SIG_PHI_E               (1.0*PI/180)                /* rms of pitch and roll */
#define     SIG_PHI_N               (1.0*PI/180)                /* rms of pitch and roll */
#define     SIG_PHI_U               (1.0*PI/180)                /* rms of heading */
#define     SIG_ACC                 (0.3)                       /* rms of acc error(m/(s.s)) */
#define     SIG_GYRO                (10.0*DEG2RAD/3600.0)       /* rms of gyro error  */

#define     GYRO_TIME_CONSTANT      (10.0F)
#define     ACC_TIME_CONSTANT       (10.0F)
#define     SIGMA_WIN               ((FLT)1.0e-6)
#define     SIGMA_ACC               ((FLT)((5.0e-4) * 9.78032667 * (5.0e-4) * 9.78032667))
#define     SIGMA_GYRO              ((FLT)(2.0 * PI / 180.0 / 3600 * 2.0 * PI / 180.0 / 3600))

static const DBL INIT_RMS[] = {SIG_PHI_E, SIG_PHI_N, SIG_PHI_U, SIG_GYRO, SIG_GYRO, SIG_GYRO, SIG_ACC, SIG_ACC, SIG_ACC};

static FLT dtCalculate(U32 timeNow, U32 timeLast);
static void setPhimQd(U32 utime, kalmanInfo_t* const pKalmanInfo, const ahrsFixData_t* const pAhrsFixData);
static U32 accMeasUpdate(const FLT acc[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData);
static U32 magMeasUpdate(const FLT mag[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData);
static void errCorrection(kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData);

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static FLT dtCalculate(U32 timeNow, U32 timeLast)
{
    if (timeLast > timeNow)
    {
        return (0xFFFFFFFF - timeLast + timeNow) / 1000.0F;
    }
    else
    {
        return (timeNow - timeLast) / 1000.0F;
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
FLT fHeadingMod(FLT fHeading)
{
    if (fHeading > PI)
    {
        fHeading -= 2*PI;
    }
    else if (fHeading < -PI)
    {
        fHeading += 2*PI;
    }
    else
    {

    }

    return fHeading;
}

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

    for (i = CHX; i <= CHZ; i++)
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

    for (i = CHX; i <= CHZ; i++)
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
    for (i = CHX; i <= CHZ; i++)
    {
        fg[i] = -acc[i];
        fmag[i] = mag[i];
        fR[i][CHZ] = fg[i];
        fR[i][CHX] = fmag[i];
    }

    // set y vector to vector product of z and x vectors
    fR[CHX][CHY] = fR[CHY][CHZ] * fR[CHZ][CHX] - fR[CHZ][CHZ] * fR[CHY][CHX];
    fR[CHY][CHY] = fR[CHZ][CHZ] * fR[CHX][CHX] - fR[CHX][CHZ] * fR[CHZ][CHX];
    fR[CHZ][CHY] = fR[CHX][CHZ] * fR[CHY][CHX] - fR[CHY][CHZ] * fR[CHX][CHX];

    // set x vector to vector product of y and z vectors
    fR[CHX][CHX] = fR[CHY][CHY] * fR[CHZ][CHZ] - fR[CHZ][CHY] * fR[CHY][CHZ];
    fR[CHY][CHX] = fR[CHZ][CHY] * fR[CHX][CHZ] - fR[CHX][CHY] * fR[CHZ][CHZ];
    fR[CHZ][CHX] = fR[CHX][CHY] * fR[CHY][CHZ] - fR[CHY][CHY] * fR[CHX][CHZ];

    for (i = CHX; i <= CHZ; i++)
    {
        fmod[i] = sqrtf(fR[CHX][i] * fR[CHX][i] + fR[CHY][i] * fR[CHY][i] + fR[CHZ][i] * fR[CHZ][i]);
    }

    if (!((fmod[CHX] < 1e-5) || (fmod[CHY] < 1e-5) || (fmod[CHZ] < 1e-5)))
    {
        for (j = CHX; j <= CHZ; j++)
        {
            ftmp = 1.0F / fmod[j];
            for (i = CHX; i <= CHZ; i++)
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

    fmodMag = sqrtf(fmag[CHX] * fmag[CHX] + fmag[CHY] * fmag[CHY] + fmag[CHZ] * fmag[CHZ]);
    fGdotMag = fg[CHX] * fmag[CHX] + fg[CHY] * fmag[CHY] + fg[CHZ] * fmag[CHZ];
    if (!((fmod[CHZ] == 0.0F) || (fmodMag == 0.0F)))
    {
        pAhrsFixData->fDelta = asinf(fGdotMag / (fmod[CHZ] * fmodMag));
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

    return retval;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 deviceHorizonAlignment(const FLT acc[], ahrsFixData_t* const pAhrsFixData)
{
    FLT fmodG = 0.0F;
    FLT fq[4] = {0.0F};

    fmodG = sqrtf(acc[CHX] * acc[CHX] + acc[CHY] * acc[CHY] + acc[CHZ] * acc[CHZ]);
    if (fmodG < 1e-5)
    {
        // no solution is possible
        return -1;
    }
    pAhrsFixData->fPsiPl = 0.0F;
    pAhrsFixData->fThePl = -asinf(acc[CHX]/(-fmodG));
    pAhrsFixData->fPhiPl = atan2f(acc[CHY]/(-fmodG), acc[CHZ]/(-fmodG));

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
U32 deviceHeadingAlignment(const FLT mag[], ahrsFixData_t* const pAhrsFixData)
{
    FLT mHorizon[CHN] = {0.0F};
    FLT mVector[CHN] = {0.0F};
    FLT fq[4] = {0.0F};

    mHorizon[CHX] = cosf(pAhrsFixData->fThePl) * mag[CHX] + sinf(pAhrsFixData->fPhiPl) * sinf(pAhrsFixData->fThePl) * mag[CHY]
                + cosf(pAhrsFixData->fPhiPl) * sinf(pAhrsFixData->fThePl) * mag[CHZ];
    mHorizon[CHY] = cosf(pAhrsFixData->fPhiPl) * mag[CHY] - sinf(pAhrsFixData->fPhiPl) * mag[CHZ];
    mHorizon[CHZ] = -sinf(pAhrsFixData->fThePl) * mag[CHX] + sinf(pAhrsFixData->fPhiPl) * cosf(pAhrsFixData->fThePl) * mag[CHY]
                + cosf(pAhrsFixData->fPhiPl) * cosf(pAhrsFixData->fThePl) * mag[CHZ];

    pAhrsFixData->fPsiPl = atan2f(-mHorizon[CHY], mHorizon[CHX]);

    euler2q(fq, pAhrsFixData->fPsiPl, pAhrsFixData->fThePl, pAhrsFixData->fPhiPl);
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];

    q2dcm(fq, pAhrsFixData->fCbn);
    memcpy(pAhrsFixData->fCnb, pAhrsFixData->fCbn, sizeof(pAhrsFixData->fCnb));
    f3x3matrixTranspose(pAhrsFixData->fCnb);

    // estimate the dip angle
    mVector[CHX] = cosf(pAhrsFixData->fPsiPl) * mHorizon[CHX] - sinf(pAhrsFixData->fPsiPl) * mHorizon[CHY];
    //mVector[CHY] = sinf(pAhrsFixData->fPsiPl) * mHorizon[CHX] + cosf(pAhrsFixData->fPsiPl) * mHorizon[CHY];
    mVector[CHZ] = mHorizon[CHZ];
    pAhrsFixData->fDelta = atanf(mVector[CHZ]/mVector[CHX]);
    
    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 quaternionIntegration (U32 utime, const FLT gyro[], ahrsFixData_t* const pAhrsFixData)
{
    FLT fq[4] = {0};
    FLT fdq[4] = {0};
    FLT fdt = 0;
    U32 i = 0;

    fq[0] = pAhrsFixData->fqPl.q0;
    fq[1] = pAhrsFixData->fqPl.q1;
    fq[2] = pAhrsFixData->fqPl.q2;
    fq[3] = pAhrsFixData->fqPl.q3;

    fdq[0] = -(gyro[0] * fq[1] + gyro[1] * fq[2] + gyro[2] * fq[3]) / 2.0F;
    fdq[1] =  (gyro[0] * fq[0] + gyro[2] * fq[2] - gyro[1] * fq[3]) / 2.0F;
    fdq[2] =  (gyro[1] * fq[0] - gyro[2] * fq[1] + gyro[0] * fq[3]) / 2.0F;
    fdq[3] =  (gyro[2] * fq[0] + gyro[1] * fq[1] - gyro[0] * fq[2]) / 2.0F;

    fdt = dtCalculate(utime, pAhrsFixData->uTime);
    for (i = 0; i < 4; i++)
    {
        fq[i] += fdq[i] * fdt;
    }
    qNorm(fq);

    //restore integration result
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];
    q2dcm(fq, pAhrsFixData->fCbn);
    dcm2euler(pAhrsFixData->fCbn, &pAhrsFixData->fPsiPl, &pAhrsFixData->fThePl, &pAhrsFixData->fPhiPl);
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
U32 ahrsKalmanInit(kalmanInfo_t* const pKalmanInfo)
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
U32 ahrsKalmanExec(U32 utime, const FLT acc[], const FLT mag[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    setPhimQd(utime, pKalmanInfo, pAhrsFixData);
    udKfPredict(pKalmanInfo);
    accMeasUpdate(acc, pKalmanInfo, pAhrsFixData);
    if (mag != NULL)
    {
        magMeasUpdate(mag, pKalmanInfo, pAhrsFixData);
    }
    errCorrection(pKalmanInfo, pAhrsFixData);

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void setPhimQd(U32 utime, kalmanInfo_t* const pKalmanInfo, const ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    U32 stateNum = STATE_NUM;
    DBL phim[STATE_NUM][STATE_NUM] = {0};
    DBL qdt[STATE_NUM][STATE_NUM] = {0};
    DBL G[STATE_NUM][STATE_NUM] = {0};         // the row and col of shaping matrix are related with model rather than fixed.
    DBL GT[STATE_NUM][STATE_NUM] = {0};        // the transpose of G matrix
    DBL M2[STATE_NUM][STATE_NUM] = {0};
    DBL temp[STATE_NUM][STATE_NUM] = {0};
    DBL *pRowA[STATE_NUM] = {0};
    DBL *pRowB[STATE_NUM] = {0};
    DBL *pRowC[STATE_NUM] = {0};
    FLT fdt = 0.0F;

    for (i = 0; i <stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)
        {
            phim[i][j] = 0.0F;
            qdt[i][j] = 0.0F;
        }
    }

    //set PHI matrix
    phim[0][3] = (DBL) -pAhrsFixData->fCbn[0][0];
    phim[0][4] = (DBL) -pAhrsFixData->fCbn[0][1];
    phim[0][5] = (DBL) -pAhrsFixData->fCbn[0][2];

    phim[1][3] = (DBL) -pAhrsFixData->fCbn[1][0];
    phim[1][4] = (DBL) -pAhrsFixData->fCbn[1][1];
    phim[1][5] = (DBL) -pAhrsFixData->fCbn[1][2];

    phim[2][3] = (DBL) -pAhrsFixData->fCbn[2][0];
    phim[2][4] = (DBL) -pAhrsFixData->fCbn[2][1];
    phim[2][5] = (DBL) -pAhrsFixData->fCbn[2][2];

    phim[3][3] = (DBL) - 1.0F / GYRO_TIME_CONSTANT;
    phim[4][4] = (DBL) - 1.0F / GYRO_TIME_CONSTANT;
    phim[5][5] = (DBL) - 1.0F / GYRO_TIME_CONSTANT;

    phim[6][6] = (DBL) - 1.0F / ACC_TIME_CONSTANT;
    phim[7][7] = (DBL) - 1.0F / ACC_TIME_CONSTANT;
    phim[8][8] = (DBL) - 1.0F / ACC_TIME_CONSTANT;

    //set Q matrix
    qdt[0][0] = (DBL)SIGMA_WIN;
    qdt[1][1] = (DBL)SIGMA_WIN;
    qdt[2][2] = (DBL)SIGMA_WIN;

    qdt[3][3] = (DBL)SIGMA_GYRO;
    qdt[4][4] = (DBL)SIGMA_GYRO;
    qdt[5][5] = (DBL)SIGMA_GYRO;

    qdt[6][6] = (DBL)SIGMA_ACC;
    qdt[7][7] = (DBL)SIGMA_ACC;
    qdt[8][8] = (DBL)SIGMA_ACC;

    // set G matrix
    G[0][0] = (DBL) -pAhrsFixData->fCbn[0][0];
    G[0][1] = (DBL) -pAhrsFixData->fCbn[0][1];
    G[0][2] = (DBL) -pAhrsFixData->fCbn[0][2];

    G[1][0] = (DBL) -pAhrsFixData->fCbn[1][0];
    G[1][1] = (DBL) -pAhrsFixData->fCbn[1][1];
    G[1][2] = (DBL) -pAhrsFixData->fCbn[1][2];

    G[2][0] = (DBL) -pAhrsFixData->fCbn[2][0];
    G[2][1] = (DBL) -pAhrsFixData->fCbn[2][1];
    G[2][2] = (DBL) -pAhrsFixData->fCbn[2][2];

    for (i = 3; i < stateNum; i++)
    {
        G[i][i] = 1.0;
    }

    // GT = G'
    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)
        {
            GT[i][j] = G[j][i];
        }
    }

    // qdt = G*w*G'
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = G[i];
        pRowB[i] = qdt[i];
        pRowC[i] = temp[i];
    }
    matrixMult(pRowA, pRowB, stateNum, stateNum, stateNum, stateNum, pRowC);
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = GT[i];
    }
    matrixMult(pRowC, pRowA, stateNum, stateNum, stateNum, stateNum, pRowB);

    // Q matrix discretization-2 order
    // M2=phi¡ÁM1£¬M1£½Q
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = phim[i];
        pRowB[i] = qdt[i];
        pRowC[i] = M2[i];
    }
    matrixMult(pRowA, pRowB, stateNum, stateNum, stateNum, stateNum, pRowC);

    fdt = dtCalculate(utime, pAhrsFixData->uTime);
    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)

        {
            qdt[i][j] = qdt[i][j] * fdt + (M2[i][j] + M2[j][i]) * fdt * fdt / 2.0;
        }
    }

    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)
        {
            pKalmanInfo->Q[i][j] = qdt[i][j];
        }
    }

    // phi matrix discretization-2 order
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = phim[i];
        pRowB[i] = temp[i];
    }
    matrixMult(pRowA, pRowA, stateNum, stateNum, stateNum, stateNum, pRowB);

    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)
        {
            phim[i][j] = phim[i][j] * fdt + temp[i][j] * fdt * fdt / 2.0; //second order phi matrix

            if (j == i)
            {
                phim[i][j] += 1.0;
            }
        }
    }

    pKalmanInfo->msCnt = utime;
    pKalmanInfo->periodTms = (U16)(fdt*1000 + 0.5);

    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)
        {
            if (j >= i)
            {
                pKalmanInfo->A[uMatIdx(i + 1, j + 1, stateNum)] = phim[i][j];
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
static U32 accMeasUpdate(const FLT acc[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    DBL zc = 0.0;
    DBL rc = 0.0;
    DBL hc[STATE_NUM] = {0.0};
    DBL z[MEAS_ACC_NUM] = {0.0};
    DBL h[MEAS_ACC_NUM][STATE_NUM] = {0.0};
    DBL r[MEAS_ACC_NUM] = {0.0};
    DBL gEstimate[3] = {0.0};
    DBL test = 0.0;
    DBL deltaX[STATE_NUM] = {0.0};

    h[0][1] = GRAVITY;
    h[1][0] = -GRAVITY;
    h[0][6] = pAhrsFixData->fCbn[0][0];
    h[0][7] = pAhrsFixData->fCbn[0][1];
    h[0][8] = pAhrsFixData->fCbn[0][2];
    h[1][6] = pAhrsFixData->fCbn[1][0];
    h[1][7] = pAhrsFixData->fCbn[1][1];
    h[1][8] = pAhrsFixData->fCbn[1][2];
    h[2][6] = pAhrsFixData->fCbn[2][0];
    h[2][7] = pAhrsFixData->fCbn[2][1];
    h[2][8] = pAhrsFixData->fCbn[2][2];

    r[0] = 5 * 5;
    r[1] = 5 * 5;
    r[2] = 5 * 5;

    for (i = CHX; i <= CHZ; i++)
    {
        gEstimate[i] = -(pAhrsFixData->fCbn[i][CHX] * acc[CHX] + pAhrsFixData->fCbn[i][CHY] * acc[CHY] + pAhrsFixData->fCbn[i][CHZ] * acc[CHZ]);
    }

    z[CHX] = 0 - gEstimate[CHX];
    z[CHY] = 0 - gEstimate[CHY];
    z[CHZ] = GRAVITY - gEstimate[CHZ];

    for (i = 0; i < MEAS_ACC_NUM; i++)
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
static U32 magMeasUpdate(const FLT mag[], kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    DBL zc = 0.0;
    DBL rc = 0.0;
    DBL hc[STATE_NUM] = {0.0};
    DBL z[MEAS_MAG_NUM] = {0.0};
    DBL h[MEAS_MAG_NUM][STATE_NUM] = {0.0};
    DBL r[MEAS_MAG_NUM] = {0.0};
    DBL magEstimate[3] = {0.0};
    DBL magVector[3] = {0.0};
    DBL test = 0.0;
    DBL deltaX[STATE_NUM] = {0};

    magVector[0] = pAhrsFixData->fB * cosf(pAhrsFixData->fDelta);
    magVector[2] = pAhrsFixData->fB * sinf(pAhrsFixData->fDelta);
    h[0][1] = magVector[2];
    h[1][0] = -magVector[2];
    h[1][2] = magVector[0];
    h[2][1] = -magVector[0];

    r[0] = 20 * 20;
    r[1] = 20 * 20;
    r[2] = 20 * 20;

    for (i = CHX; i <= CHZ; i++)
    {
        magEstimate[i] = pAhrsFixData->fCbn[i][CHX] * mag[CHX] + pAhrsFixData->fCbn[i][CHY] * mag[CHY] + pAhrsFixData->fCbn[i][CHZ] * mag[CHZ];
    }

    z[CHX] = magVector[CHX] - magEstimate[CHX];
    z[CHY] = magVector[CHY] - magEstimate[CHY];
    z[CHZ] = magVector[CHZ] - magEstimate[CHZ];

    for (i = 0; i < MEAS_MAG_NUM; i++)
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
static void errCorrection(kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    FLT fq[4] = {0.0};
    FLT deltaCbn[3][3] = {0.0};
    FLT temp = 0.0;
    FLT tempMatrix[3][3] = {0.0};

    euler2dcm(deltaCbn, (FLT)pKalmanInfo->X[2], (FLT)pKalmanInfo->X[1], (FLT)pKalmanInfo->X[0]);
    for (i = CHX; i <= CHZ; i++)
    {
        for (j = CHX; j <= CHZ; j++)
        {
            temp = 0.0F;
            for (k = CHX; k <= CHZ; k++)
            {
                temp += deltaCbn[i][k] * pAhrsFixData->fCbn[k][j];
            }
            tempMatrix[i][j] = temp;
        }
    }
    memcpy(pAhrsFixData->fCbn, tempMatrix, sizeof(pAhrsFixData->fCbn));
    memcpy(pAhrsFixData->fCnb, tempMatrix, sizeof(pAhrsFixData->fCnb));
    f3x3matrixTranspose(pAhrsFixData->fCnb);

    for (i = CHX; i <= CHZ; i++)
    {
        pAhrsFixData->fGyroBias[i] += (FLT)pKalmanInfo->X[i + 3];
        pAhrsFixData->fAccBias[i] += (FLT)pKalmanInfo->X[i + 6];
    }

    dcm2euler(pAhrsFixData->fCbn, &pAhrsFixData->fPsiPl, &pAhrsFixData->fThePl, &pAhrsFixData->fPhiPl);
    euler2q(fq, pAhrsFixData->fPsiPl, pAhrsFixData->fThePl, pAhrsFixData->fPhiPl);
    qNorm(fq);
    pAhrsFixData->fqPl.q0 = fq[0];
    pAhrsFixData->fqPl.q1 = fq[1];
    pAhrsFixData->fqPl.q2 = fq[2];
    pAhrsFixData->fqPl.q3 = fq[3];

    // clear x
    for (i = 0; i < STATE_NUM; i++)
    {
        pKalmanInfo->X[i] = 0.0;
    }
}
