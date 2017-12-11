#include <math.h>
#include <stdio.h>
#include "misc.h"
#include "magcal.h"


static void updateCalibration4INV(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer);
static void updateCalibration7EIG(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer);


/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 magCalibrationInit(magCalibration_t* const pMagCalibration, magneticBuffer_t* const pMagBuffer)
{
    U32 i = 0;
    U32 j = 0;

    // initialize the calibration parameters
    f3x3matrixEqI(pMagCalibration->finvW);
    for (i = 0; i < 3; i++)
    {
        pMagCalibration->fV[i] = 0.0F;
    }
    pMagCalibration->fB = DEFAULTB;
    pMagCalibration->fFitErrorpc = 1000.0F;
    pMagCalibration->iValidMagCal = 0;

    // initialize the mag buffer
    pMagBuffer->iMagBufferCount = 0;
    for (i = 0; i < MAGBUFFSIZEX; i++)
    {
        for (j = 0; j < MAGBUFFSIZEY; j++)
        {
            pMagBuffer->index[i][j] = -1;
        }
    }
    pMagBuffer->fuTPerCount = MAGSENSITIVE;
    pMagBuffer->fCountsPeruT = (FLT)(1.0 / pMagBuffer->fuTPerCount);

    // initialize the array of (MAGBUFFSIZEX - 1) elements of 100 * tangents used for buffer indexing
    // entries cover the range 100 * tan(-PI/2 + PI/MAGBUFFSIZEX), 100 * tan(-PI/2 + 2*PI/MAGBUFFSIZEX) to
    // 100 * tan(-PI/2 + (MAGBUFFSIZEX - 1) * PI/MAGBUFFSIZEX).
    // for MAGBUFFSIZEX=14, the entries range in value from -438 to +438
    for (i = 0; i < (MAGBUFFSIZEX - 1); i++)
    {
        pMagBuffer->tanarray[i] = (S32) (100.0F * tanf(PI * (-0.5F + (i + 1) * 1.0F / MAGBUFFSIZEX)));
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
U32 magBufferUpdate(magneticBuffer_t * const pMagBuffer, const FLT* const pMagRaw, const FLT* const pMagCal, const U32 loopCounter)
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    U32 l = 0;
    U32 m = 0;
    S32 iMagRaw[CHN] = {0};
    S32 iMagCal[CHN] = {0};
    S32 itanj = 0;
    S32 itank = 0;
    U32 idelta = 0;
    U32 iclose = 0;

    // convert float mag data to int mag data to reduce multiplications
    for (i = CHX; i <= CHZ; i++)
    {
        iMagRaw[i] = (S32)(pMagRaw[i] * pMagBuffer->fCountsPeruT);
        iMagCal[i] = (S32)(pMagCal[i] * pMagBuffer->fCountsPeruT);
    }
    
    if (iMagCal[CHZ] == 0)
    {
        return -1;
    }
    itanj = 100 * iMagCal[CHX] / iMagCal[CHZ];
    itank = 100 * iMagCal[CHY] / iMagCal[CHZ];
    while ((j < (MAGBUFFSIZEX - 1) && (itanj >= pMagBuffer->tanarray[j]))) j++;
    while ((k < (MAGBUFFSIZEX - 1) && (itank >= pMagBuffer->tanarray[k]))) k++;
    if (iMagCal[CHX] < 0) k += MAGBUFFSIZEX;

    // case 1: buffer is full and this bin has a measurement: over-write without increasing number of measurements
    // this is the most common option at run time
    if ((pMagBuffer->iMagBufferCount == MAXMEASUREMENTS) && (pMagBuffer->index[j][k] != -1))
    {
        for (i = CHX; i <= CHZ; i++)
        {
            pMagBuffer->iMagRaw[i][j][k] = iMagRaw[i];
        }
        pMagBuffer->index[j][k] = loopCounter;

        return 0;
    }

    // case 2: the buffer is full and this bin does not have a measurement: store and retire the oldest
    // this is the second most common option at run time
    if ((pMagBuffer->iMagBufferCount == MAXMEASUREMENTS) && (pMagBuffer->index[j][k] == -1))
    {
        for (i = CHX; i <= CHZ; i++)
        {
            pMagBuffer->iMagRaw[i][j][k] = iMagRaw[i];
        }
        pMagBuffer->index[j][k] = loopCounter;

        // set l and m to the oldest active entry and disable it
        for (j = 0; j < MAGBUFFSIZEX; j++)
        {
            for (k = 0; k < MAGBUFFSIZEY; k++)
            {
                // check if the time stamp is older than the oldest found so far (normally fails this test)
                if (pMagBuffer->index[j][k] < (S32)loopCounter)
                {
                    // check if this bin is active (normally passes this test)
                    if (pMagBuffer->index[j][k] != -1)
                    {
                        // set l and m to the indices of the oldest entry found so far
                        l = j;
                        m = k;
                        // set i to the time stamp of the oldest entry found so far
                        i = pMagBuffer->index[l][m];
                    }
                }
            }
        }
        // deactivate the oldest measurement (no need to zero the measurement data)
        pMagBuffer->index[l][m] = -1;

        return 0;
    }

    // case 3: buffer is not full and this bin is empty: store and increment number of measurements
    if ((pMagBuffer->iMagBufferCount < MAXMEASUREMENTS) && (pMagBuffer->index[j][k] == -1))
    {
        for (i = CHX; i <= CHZ; i++)
        {
            pMagBuffer->iMagRaw[i][j][k] = iMagRaw[i];
        }
        pMagBuffer->index[j][k] = loopCounter;
        pMagBuffer->iMagBufferCount++;

        return 0;
    }

    // case 4: buffer is not full and this bin has a measurement: over-write if close or try to slot in
    // elsewhere if close to the other measurements so as to create a mesh at power up
    if ((pMagBuffer->iMagBufferCount < MAXMEASUREMENTS) && (pMagBuffer->index[j][k] != -1))
	{
		// calculate the vector difference between current measurement and the buffer entry
		idelta = 0;
		for (i = CHX; i <= CHZ; i++)
		{
			idelta += abs(iMagRaw[i] - pMagBuffer->iMagRaw[i][j][k]);
		}
		// check to see if the current reading is close to this existing magnetic buffer entry
		if (idelta < MESHDELTAUT * pMagBuffer->fCountsPeruT)
		{
			// simply over-write the measurement and return
			for (i = CHX; i <= CHZ; i++)
			{
				pMagBuffer->iMagRaw[i][j][k] = iMagRaw[i];
			}
			pMagBuffer->index[j][k] = loopCounter;
		}
		else
		{
			// reset the flag denoting that the current measurement is close to any measurement in the buffer
			iclose = 0;
			// to avoid compiler warning
			l = m = 0;
			// loop over the buffer j from 0 potentially up to MAGBUFFSIZEX - 1 
			j = 0;
			while (!iclose && (j < MAGBUFFSIZEX))
			{
				// loop over the buffer k from 0 potentially up to MAGBUFFSIZEY - 1 
				k = 0;
				while (!iclose && (k < MAGBUFFSIZEY))
				{
					// check whether this buffer entry already has a measurement or not
					if (pMagBuffer->index[j][k] != -1)
					{
						// calculate the vector difference between current measurement and the buffer entry
						idelta = 0;
						for (i = CHX; i <= CHZ; i++)
						{
							idelta += abs(iMagRaw[i] - pMagBuffer->iMagRaw[i][j][k]);
						}
						// check to see if the current reading is close to this existing magnetic buffer entry
						if (idelta < MESHDELTAUT)
						{
							// set the flag to abort the search
							iclose = 1;
						}
					}
					else
					{
						// store the location of this empty bin for future use
						l = j;
						m = k;
					} // end of test for valid measurement in this bin
					k++;
				} // end of k loop
				j++;
			} // end of j loop

			// if none too close, store the measurement in the last empty bin found and return
			// l and m are guaranteed to be set if no entries too close are detected
			if (!iclose)
			{
				for (i = CHX; i <= CHZ; i++)
				{
					pMagBuffer->iMagRaw[i][l][m] = iMagRaw[i];
				}
				pMagBuffer->index[l][m] = loopCounter;
				pMagBuffer->iMagBufferCount++;
			}
		}

		return 0;
    }

    return -1;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 magCalibrationExec(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer)
{
    U32 i = 0;
    U32 j = 0;
    U32 isolver = 0;
#ifdef DEBUG
    FILE *fp;
#endif

    // 4 element calibration case
    if (pMagBuffer->iMagBufferCount > MINMEASUREMENTS4CAL && pMagBuffer->iMagBufferCount < MINMEASUREMENTS7CAL)
    {
        if (pMagCalibration->iValidMagCal)
        {
            pMagCalibration->fFitErrorpc *= (FLT)(1.0F + INTERVAL4CAL * 1.0F / (SENSORFS * FITERRORAGINGSECS));	
        }

        isolver = 4;
#ifdef DEBUG
        // store the calibration data
        fopen_s(&fp, "../../data/magCalData.txt", "w");
        for (i = 0; i < MAGBUFFSIZEX; i++)
        {
            for (j = 0; j < MAGBUFFSIZEY; j++)
            {
                if (pMagBuffer->index[i][j] != -1)
                {
                    fprintf(fp, "%d %d %d\n", pMagBuffer->iMagRaw[CHX][i][j], pMagBuffer->iMagRaw[CHY][i][j], pMagBuffer->iMagRaw[CHZ][i][j]);
                }
            }
        }
        fclose(fp);
#endif
        updateCalibration4INV(pMagCalibration, pMagBuffer);
    }
    // 7 element calibration case
    else if (pMagBuffer->iMagBufferCount > MINMEASUREMENTS7CAL)
    {
        if (pMagCalibration->iValidMagCal)
        {
            pMagCalibration->fFitErrorpc *= (FLT)(1.0F + INTERVAL7CAL * 1.0F / (SENSORFS * FITERRORAGINGSECS));
        }
        isolver = 7;
#ifdef DEBUG
        // store the calibration data
        fopen_s(&fp, "../../data/magCalData.txt", "w");
        for (i = 0; i < MAGBUFFSIZEX; i++)
        {
            for (j = 0; j < MAGBUFFSIZEY; j++)
            {
                if (pMagBuffer->index[i][j] != -1)
                {
                    fprintf(fp, "%d %d %d\n", pMagBuffer->iMagRaw[CHX][i][j], pMagBuffer->iMagRaw[CHY][i][j], pMagBuffer->iMagRaw[CHZ][i][j]);
                }
            }
        }
        fclose(fp);
#endif
        updateCalibration7EIG(pMagCalibration, pMagBuffer);
    }

    // the trial geomagnetic field must be in range (earth is 22uT to 67uT)
    if ((pMagCalibration->ftrB >= MINBFITUT) && (pMagCalibration->ftrB <= MAXBFITUT))		
    {
        // always accept the calibration if i) no previous calibration exists ii) the calibration fit is reduced or
        // an improved solver was used giving a good trial calibration (4% or under)
        if ((pMagCalibration->iValidMagCal == 0) ||
            (pMagCalibration->ftrFitErrorpc <= pMagCalibration->fFitErrorpc) ||
            ((isolver > pMagCalibration->iValidMagCal) && (pMagCalibration->ftrFitErrorpc <= 4.0F)))
        {
            pMagCalibration->iValidMagCal = isolver;
            pMagCalibration->fFitErrorpc = pMagCalibration->ftrFitErrorpc;
            pMagCalibration->fB = pMagCalibration->ftrB;
            for (i = CHX; i <= CHZ; i++)
            {
                pMagCalibration->fV[i] = pMagCalibration->ftrV[i];
                for (j = CHX; j <= CHZ; j++)
                {
                    pMagCalibration->finvW[i][j] = pMagCalibration->ftrinvW[i][j];
                }
            }
        }
    }

    return isolver;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void magCorrection(FLT mag[], const magCalibration_t* const pMagCalibration)
{
    U32 i = 0;
    FLT ftemp[CHN] = {0.0};

    // remove the computed hard iron offsets (uT): ftmp[] = fmag[] - fV[]
    for (i = CHX; i <= CHZ; i++)
    {
        ftemp[i] = mag[i] - pMagCalibration->fV[i];
    }
    // remove the computed soft iron offsets (uT): fmag = inv(W)*(fmag[] - fV[])
    for (i = CHX; i <= CHZ; i++)
    {
        mag[i] = pMagCalibration->finvW[i][CHX]*ftemp[CHX] + pMagCalibration->finvW[i][CHY]*ftemp[CHY] + pMagCalibration->finvW[i][CHZ]*ftemp[CHZ];
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void updateCalibration4INV(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer)
{
    FLT fBs2;
    FLT fSumBs4;
    FLT fscaling;
    FLT fE;
    S32 iOffset[3];
    U32 iCount;
    U32 ierror;
    U32 i, j, k, l;

    FLT *pfRows[4];
    S32 iColInd[4];
    S32 iRowInd[4];
    S32 iPivot[4];

    fscaling = pMagBuffer->fuTPerCount / DEFAULTB;
    f3x3matrixEqI(pMagCalibration->ftrinvW);
    fSumBs4 = 0.0F;
    for (i = 0; i < 4; i++)
    {
        pMagCalibration->fvecB[i] = 0.0F;
        for (j = i; j < 4; j++)
        {
            pMagCalibration->fmatA[i][j] = 0.0F;
        }
    }
    iOffset[CHX] = iOffset[CHY] = iOffset[CHZ] = 0;
    iCount = 0;
    for (j = 0; j < MAGBUFFSIZEX; j++)
    {
        for (k = 0; k < MAGBUFFSIZEY; k++)
        {
            if (pMagBuffer->index[j][k] != -1)
            {
                if (iCount == 0)
                {
                    for (l = CHX; l <= CHZ; l++)
                    {
                        iOffset[l] = pMagBuffer->iMagRaw[l][j][k];
                    }
                }
                for (l = CHX; l <= CHZ; l++)
                {
                    pMagCalibration->fvecA[l] = (FLT)(pMagBuffer->iMagRaw[l][j][k] - iOffset[l]) * fscaling;
                    pMagCalibration->fvecA[l + 3] = pMagCalibration->fvecA[l] * pMagCalibration->fvecA[l];
                }
                fBs2 = pMagCalibration->fvecA[3] + pMagCalibration->fvecA[4] + pMagCalibration->fvecA[5];
                fSumBs4 += fBs2 * fBs2;
                for (l = CHX; l <= CHZ; l++)
                {
                    pMagCalibration->fvecB[l] += pMagCalibration->fvecA[l] * fBs2;
                }
                pMagCalibration->fvecB[3] += fBs2;
                pMagCalibration->fmatA[0][0] += pMagCalibration->fvecA[CHX + 3];
                pMagCalibration->fmatA[0][1] += pMagCalibration->fvecA[CHX] * pMagCalibration->fvecA[CHY];
                pMagCalibration->fmatA[0][2] += pMagCalibration->fvecA[CHX] * pMagCalibration->fvecA[CHZ];
                pMagCalibration->fmatA[0][3] += pMagCalibration->fvecA[CHX];
                pMagCalibration->fmatA[1][1] += pMagCalibration->fvecA[CHY + 3];
                pMagCalibration->fmatA[1][2] += pMagCalibration->fvecA[CHY] * pMagCalibration->fvecA[CHZ];
                pMagCalibration->fmatA[1][3] += pMagCalibration->fvecA[CHY];
                pMagCalibration->fmatA[2][2] += pMagCalibration->fvecA[CHZ + 3];
                pMagCalibration->fmatA[2][3] += pMagCalibration->fvecA[CHZ];
                iCount++;
            }
        }
    }
    pMagCalibration->fmatA[3][3] = (FLT)iCount;
    for (i = 0; i < 4; i++)
    {
        for (j = i; j < 4; j++)
        {
            pMagCalibration->fmatB[i][j] = pMagCalibration->fmatB[j][i] = pMagCalibration->fmatA[j][i] = pMagCalibration->fmatA[i][j];
        }
    }
    for (i = 0; i < 4; i++)
    {
        pfRows[i] = pMagCalibration->fmatB[i];
    }
    fmatrixAeqInvA(pfRows, iColInd, iRowInd, iPivot, 4, &ierror);
    for (i = 0; i < 4; i++)
    {
        pMagCalibration->fvecA[i] = 0.0F;
        for (k = 0; k < 4; k++)
        {
            pMagCalibration->fvecA[i] += pMagCalibration->fmatB[i][k] * pMagCalibration->fvecB[k];
        } 
    } 
    fE = 0.0F;
    for (i = 0; i < 4; i++)
    {
        fE += pMagCalibration->fvecA[i] * pMagCalibration->fvecB[i];
    }
    fE = fSumBs4 - 2.0F * fE;	
    for (i = 0; i < 4; i++)
    {
        pMagCalibration->fvecB[i] = 0.0F;
        for (k = 0; k < 4; k++)
        {
            pMagCalibration->fvecB[i] += pMagCalibration->fmatA[i][k] * pMagCalibration->fvecA[k];
        } 
    } 
    for (i = 0; i < 4; i++)
    {
        fE += pMagCalibration->fvecB[i] * pMagCalibration->fvecA[i];
    }
    for (l = CHX; l <= CHZ; l++)
    {
        pMagCalibration->ftrV[l] = 0.5F * pMagCalibration->fvecA[l];
    }
    pMagCalibration->ftrB = sqrtf(pMagCalibration->fvecA[3] + pMagCalibration->ftrV[CHX] * pMagCalibration->ftrV[CHX] +
        pMagCalibration->ftrV[CHY] * pMagCalibration->ftrV[CHY] + pMagCalibration->ftrV[CHZ] * pMagCalibration->ftrV[CHZ]);
    pMagCalibration->ftrFitErrorpc = sqrtf(fE * 1.0F / pMagBuffer->iMagBufferCount) * 100.0F /
        (2.0F * pMagCalibration->ftrB * pMagCalibration->ftrB);
    for (l = CHX; l <= CHZ; l++)
    {
        pMagCalibration->ftrV[l] = (FLT)(pMagCalibration->ftrV[l] * DEFAULTB + iOffset[l] * 1.0 * pMagBuffer->fuTPerCount);
    }
    pMagCalibration->ftrB *= DEFAULTB;

    return;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static void updateCalibration7EIG(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer)
{
    FLT det;
    FLT fscaling;
    FLT ftmp;
    S32 iOffset[3];
    U32 iCount;
    U32 i, j, k, l, m, n;

    fscaling = pMagBuffer->fuTPerCount / DEFAULTB;
    iOffset[CHX] = iOffset[CHY] = iOffset[CHZ] = 0;
    for (m = 0; m < 7; m++)
    {
        for (n = m; n < 7; n++)
        {
            pMagCalibration->fmatA[m][n] = 0.0F;
        }
    }
    iCount = 0;
    for (j = 0; j < MAGBUFFSIZEX; j++)
    {
        for (k = 0; k < MAGBUFFSIZEY; k++)
        {
            if (pMagBuffer->index[j][k] != -1)
            {
                if (iCount == 0)
                {
                    for (l = CHX; l <= CHZ; l++)
                    {
                        iOffset[l] = pMagBuffer->iMagRaw[l][j][k];
                    }
                }
                for (l = CHX; l <= CHZ; l++)
                {
                    pMagCalibration->fvecA[l + 3] = (FLT)(pMagBuffer->iMagRaw[l][j][k] - iOffset[l]) * fscaling;
                    pMagCalibration->fvecA[l] = pMagCalibration->fvecA[l + 3] * pMagCalibration->fvecA[l + 3];
                }
                for (m = 0; m < 6; m++)
                {
                    pMagCalibration->fmatA[m][6] += pMagCalibration->fvecA[m];
                }
                for (m = 0; m < 6; m++)
                {
                    for (n = m; n < 6; n++)
                    {
                        pMagCalibration->fmatA[m][n] += pMagCalibration->fvecA[m] * pMagCalibration->fvecA[n];
                    }
                }
                iCount++;
            }
        }
    }
    pMagCalibration->fmatA[6][6] = (FLT)iCount;
    for (m = 1; m < 7; m++)
    {
        for (n = 0; n < m; n++)
        {
            pMagCalibration->fmatA[m][n] = pMagCalibration->fmatA[n][m];
        }
    }
    eigencompute10(pMagCalibration->fmatA, pMagCalibration->fvecA, pMagCalibration->fmatB, 7);

    // find the smallest eigenvalue
    j = 0;
    for (i = 1; i < 7; i++)
    {
        if (pMagCalibration->fvecA[i] < pMagCalibration->fvecA[j])
        {
            j = i;
        }
    }
    f3x3matrixEqScalar(pMagCalibration->fA, 0.0F);
    det = 1.0F;
    for (l = CHX; l <= CHZ; l++)
    {
        pMagCalibration->fA[l][l] = pMagCalibration->fmatB[l][j];
        det *= pMagCalibration->fA[l][l];
        pMagCalibration->ftrV[l] = -0.5F * pMagCalibration->fmatB[l + 3][j] / pMagCalibration->fA[l][l];
    }
    if (det < 0.0F)
    {
        f3x3matrixEqMinusA(pMagCalibration->fA);
        pMagCalibration->fmatB[6][j] = -pMagCalibration->fmatB[6][j];
        det = -det;
    }
    ftmp = -pMagCalibration->fmatB[6][j];
    for (l = CHX; l <= CHZ; l++)
    {
        ftmp += pMagCalibration->fA[l][l] * pMagCalibration->ftrV[l] * pMagCalibration->ftrV[l];
    }
    pMagCalibration->ftrFitErrorpc = 50.0F * sqrtf(fabsf(pMagCalibration->fvecA[j]) / (FLT) pMagBuffer->iMagBufferCount) / fabsf(ftmp);
    f3x3matrixEqAxScalar(pMagCalibration->fA, powf(det, -(1.0F/3.0F)));
    pMagCalibration->ftrB = sqrtf(fabsf(ftmp)) * DEFAULTB * powf(det, -(1.0F/6.0F));
    f3x3matrixEqI(pMagCalibration->ftrinvW);
    for (l = CHX; l <= CHZ; l++)
    {
        pMagCalibration->ftrinvW[l][l] = sqrtf(fabsf(pMagCalibration->fA[l][l]));
        pMagCalibration->ftrV[l] = pMagCalibration->ftrV[l] * DEFAULTB + (FLT)iOffset[l] * pMagBuffer->fuTPerCount;
    }

    return;
}