#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "ahrs.h"

#define     ALIGN_NUM               (100)

#define     STATE_NUM               (9)
#define     UD_NUM                  (STATE_NUM*(STATE_NUM+1)/2)
#define     MEAS_ACC_NUM            (3)
#define     MEAS_MAG_NUM            (3)

#define     SIG_PHI_E               (1.0*PI/180)                /* rms of pitch and roll */
#define     SIG_PHI_N               (1.0*PI/180)                /* rms of pitch and roll */
#define     SIG_PHI_U               (1.0*PI/180)                /* (rad)0.001 rms of heading */
#define     SIG_ACC                 (0.3)                       /* rms of acc error(m/(s.s)) */
#define     SIG_GYRO                (1000.0*DEG2RAD/3600.0)     /* rms of gyro error  */

#define     GYRO_TIME_CONSTANT      (100.0F)
#define     ACC_TIME_CONSTANT       (100.0F)
#define     SIGMA_WIN               ((FLT)1.0e-6)
#define     SIGMA_ACC               ((FLT)((5.0e-4) * 9.78032667 * (5.0e-4) * 9.78032667))
#define     SIGMA_GYRO              ((FLT)(20.0 * PI / 180.0 / 3600 * 20.0 * PI / 180.0 / 3600))

static const DBL INIT_RMS[] = {SIG_PHI_E, SIG_PHI_N, SIG_PHI_U, SIG_GYRO, SIG_GYRO, SIG_GYRO, SIG_ACC, SIG_ACC, SIG_ACC};

static FLT dtCalculate(U32 timeNow, U32 timeLast);
static void setPhimQd(U32 utime, kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData);
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
    U32 sizeofDBL = sizeof(DBL);

    pKalmanInfo->uStateNum = STATE_NUM;
    pKalmanInfo->uUdNum = STATE_NUM * (STATE_NUM + 1) / 2;

    // malloc state array
    pKalmanInfo->pStateX = (DBL *)malloc(sizeofDBL * pKalmanInfo->uStateNum);
    if (pKalmanInfo->pStateX == NULL)
    {
        return -1;
    }
    memset(pKalmanInfo->pStateX, 0, sizeofDBL * pKalmanInfo->uStateNum);

    // malloc ud array
    pKalmanInfo->pUd = (DBL *)malloc(sizeofDBL * pKalmanInfo->uUdNum);
    if (pKalmanInfo->pUd == NULL)
    {
        free(pKalmanInfo->pStateX);
        return -1;
    }
    memset(pKalmanInfo->pUd, 0, sizeofDBL * pKalmanInfo->uUdNum);

    // malloc q array
    pKalmanInfo->pQd = (DBL **)mallocArray2D_DBL(pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
    if (pKalmanInfo->pQd == NULL)
    {
        free(pKalmanInfo->pStateX);
        free(pKalmanInfo->pUd);
        return -1;
    }

    // malloc phi array
    pKalmanInfo->pPhim = (DBL **)mallocArray2D_DBL(pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
    if (pKalmanInfo->pPhim == NULL)
    {
        free(pKalmanInfo->pStateX);
        free(pKalmanInfo->pUd);
        freeArray2D_DBL(pKalmanInfo->pQd, pKalmanInfo->uStateNum, pKalmanInfo->uStateNum);
        return -1;
    }

    UDInit(pKalmanInfo->pUd, pKalmanInfo->uUdNum, INIT_RMS, pKalmanInfo->uStateNum);

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
    predict(pKalmanInfo);
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
static void setPhimQd(U32 utime, kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    U32 stateNum = pKalmanInfo->uStateNum;
    DBL **phim = pKalmanInfo->pPhim;
    DBL **qdt = pKalmanInfo->pQd;
    FLT (*fCbn)[3] = pAhrsFixData->fCbn;
    DBL G[STATE_NUM][STATE_NUM] = {0};         // the row and col of shaping matrix are related with model rather than fixed.
    DBL GT[STATE_NUM][STATE_NUM] = {0};        // the transpose of G matrix
    DBL M2[STATE_NUM][STATE_NUM] = {0};
    DBL temp[STATE_NUM][STATE_NUM] = {0};
    DBL *pRowA[STATE_NUM] = {0};
    DBL *pRowB[STATE_NUM] = {0};
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
    phim[0][3] = (DBL) -fCbn[0][0];
    phim[0][4] = (DBL) -fCbn[0][1];
    phim[0][5] = (DBL) -fCbn[0][2];

    phim[1][3] = (DBL) -fCbn[1][0];
    phim[1][4] = (DBL) -fCbn[1][1];
    phim[1][5] = (DBL) -fCbn[1][2];

    phim[2][3] = (DBL) -fCbn[2][0];
    phim[2][4] = (DBL) -fCbn[2][1];
    phim[2][5] = (DBL) -fCbn[2][2];

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
    G[0][0] = (DBL) -fCbn[0][0];
    G[0][1] = (DBL) -fCbn[0][1];
    G[0][2] = (DBL) -fCbn[0][2];

    G[1][0] = (DBL) -fCbn[1][0];
    G[1][1] = (DBL) -fCbn[1][1];
    G[1][2] = (DBL) -fCbn[1][2];

    G[2][0] = (DBL) -fCbn[2][0];
    G[2][1] = (DBL) -fCbn[2][1];
    G[2][2] = (DBL) -fCbn[2][2];

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
    }
    for (i = 0; i < stateNum; i++)
    {
        pRowB[i] = temp[i];
    }
    matrixMult(pRowA, qdt, stateNum, stateNum, stateNum, stateNum, pRowB);
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = GT[i];
    }
    matrixMult(pRowB, pRowA, stateNum, stateNum, stateNum, stateNum, qdt);

    // Q matrix discretization-2 order
    // M2=phi¡ÁM1£¬M1£½Q
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = M2[i];
    }
    matrixMult(phim, qdt, stateNum, stateNum, stateNum, stateNum, pRowA);

    fdt = dtCalculate(utime, pAhrsFixData->uTime);
    for (i = 0; i < stateNum; i++)
    {
        for (j = 0; j < stateNum; j++)

        {
            qdt[i][j] = qdt[i][j] * fdt + (M2[i][j] + M2[j][i]) * fdt * fdt / 2.0;
        }
    }
    
    // ud decompose for Q matrix
    udDecompose(qdt, stateNum);

    // phi matrix discretization-2 order
    for (i = 0; i < stateNum; i++)
    {
        pRowA[i] = temp[i];
    }
    matrixMult(phim, phim, stateNum, stateNum, stateNum, stateNum, pRowA);

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
    DBL xSave[STATE_NUM];
    DBL udSave[UD_NUM];
    DBL ion = 0.0;
    DBL res = 0.0;
    DBL gEstimate[3] = {0.0};
    DBL test = 0.0;

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

    for (i = X; i <= Z; i++)
    {
        gEstimate[i] = -(pAhrsFixData->fCbn[i][X] * acc[X] + pAhrsFixData->fCbn[i][Y] * acc[Y] + pAhrsFixData->fCbn[i][Z] * acc[Z]);
    }

    z[X] = 0 - gEstimate[X];
    z[Y] = 0 - gEstimate[Y];
    z[Z] = GRAVITY - gEstimate[Z];

    for (i = 0; i < MEAS_ACC_NUM; i++)
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
    DBL xSave[STATE_NUM];
    DBL udSave[UD_NUM];
    DBL ion = 0.0;
    DBL res = 0.0;
    DBL magEstimate[3] = {0.0};
    DBL magVector[3] = {0.0};
    DBL test = 0.0;

    magVector[0] = pAhrsFixData->fB * cosf(pAhrsFixData->fDelta);
    magVector[2] = pAhrsFixData->fB * sinf(pAhrsFixData->fDelta);
    h[0][1] = magVector[2];
    h[1][0] = -magVector[2];
    h[1][2] = magVector[0];
    h[2][1] = -magVector[0];

    r[0] = 20 * 20;
    r[1] = 20 * 20;
    r[2] = 20 * 20;

    for (i = X; i <= Z; i++)
    {
        magEstimate[i] = pAhrsFixData->fCbn[i][X] * mag[X] + pAhrsFixData->fCbn[i][Y] * mag[Y] + pAhrsFixData->fCbn[i][Z] * mag[Z];
    }

    z[X] = magVector[X] - magEstimate[X];
    z[Y] = magVector[Y] - magEstimate[Y];
    z[Z] = magVector[Z] - magEstimate[Z];

    for (i = 0; i < MEAS_MAG_NUM; i++)
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
static void errCorrection(kalmanInfo_t* const pKalmanInfo, ahrsFixData_t* const pAhrsFixData)
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    FLT fq[4] = {0.0};
    FLT deltaCbn[3][3] = {0.0};
    FLT temp = 0.0;
    FLT tempMatrix[3][3] = {0.0};

    euler2dcm(deltaCbn, (FLT)pKalmanInfo->pStateX[2], (FLT)pKalmanInfo->pStateX[1], (FLT)pKalmanInfo->pStateX[0]);
    for (i = X; i <= Z; i++)
    {
        for (j = X; j <= Z; j++)
        {
            temp = 0.0F;
            for (k = X; k <= Z; k++)
            {
                temp += deltaCbn[i][k] * pAhrsFixData->fCbn[k][j];
            }
            tempMatrix[i][j] = temp;
        }
    }
    memcpy(pAhrsFixData->fCbn, tempMatrix, sizeof(pAhrsFixData->fCbn));
    memcpy(pAhrsFixData->fCnb, tempMatrix, sizeof(pAhrsFixData->fCnb));
    f3x3matrixTranspose(pAhrsFixData->fCnb);

    for (i = X; i <= Z; i++)
    {
        pAhrsFixData->fGyroBias[i] += (FLT)pKalmanInfo->pStateX[i + 3];
        pAhrsFixData->fAccBias[i] += (FLT)pKalmanInfo->pStateX[i + 6];
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
        pKalmanInfo->pStateX[i] = 0.0;
    }
}
