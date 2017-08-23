#ifndef _MAGCAL_H_
#define _MAGCAL_H_

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// magnetic config
#define     MAGSENSITIVE                (1.0F)                      // MAGSENSITIVE (uT per counts)
#define     SENSORFS 			        (200.0F)                    // frequency of gyro sensor sampling process

// magnetic calibration constants
#define     MAGBUFFSIZEX                (14)					    // x dimension in magnetometer buffer (12x24 equals 288 elements)
#define     MAGBUFFSIZEY                (2 * MAGBUFFSIZEX)	        // y dimension in magnetometer buffer (12x24 equals 288 elements)
#define     MINMEASUREMENTS4CAL         (100)			            // minimum number of measurements for 4 element calibration
#define     MINMEASUREMENTS7CAL         (150)			            // minimum number of measurements for 7 element calibration
#define     MAXMEASUREMENTS             (240)				        // maximum number of measurements used for calibration
#define     MINBFITUT                   (10.0F)					    // minimum acceptable geomagnetic field B (uT) for valid calibration
#define     MAXBFITUT                   (90.0F)					    // maximum acceptable geomagnetic field B (uT) for valid calibration
#define     DEFAULTB                    (50.0F)					    // default geomagnetic field (uT)
#define     MESHDELTAUT                 (5)        				    // magnetic buffer mesh spacing in uT
#define     INTERVAL4CAL                (300)					    // 3s at 100Hz: 4 element interval (samples)
#define     INTERVAL7CAL                (1000)				        // 10s at 100Hz: 7 element interval (samples)
#define     MINBFITUT                   (10.0F)					    // minimum acceptable geomagnetic field B (uT) for valid calibration
#define     MAXBFITUT                   (90.0F)					    // maximum acceptable geomagnetic field B (uT) for valid calibration
#define     FITERRORAGINGSECS           (86400.0F)          		// 24 hours: time (s) for fit error to increase (age) by e=2.718


    typedef struct magneticBuffer
    {
        FLT fuTPerCount;						        // uT per count
        FLT fCountsPeruT;						        // counts per uT
        S32 iMagRaw[3][MAGBUFFSIZEX][MAGBUFFSIZEY];		// uncalibrated magnetometer readings
        S32 index[MAGBUFFSIZEX][MAGBUFFSIZEY];		    // array of time indices
        S32 tanarray[MAGBUFFSIZEX - 1];				    // array of tangents of (100 * angle)
        U32 iMagBufferCount;						    // number of magnetometer readings
    } magneticBuffer_t;

    typedef struct magCalibration
    {
        FLT fV[3];					// current hard iron offset x, y, z, (uT)
        FLT finvW[3][3];			// current inverse soft iron matrix
        FLT fB;						// current geomagnetic field magnitude (uT)
        FLT fFitErrorpc;			// current fit error %
        FLT ftrV[3];				// trial value of hard iron offset z, y, z (uT)
        FLT ftrinvW[3][3];			// trial inverse soft iron matrix size
        FLT ftrB;					// trial value of geomagnetic field magnitude in uT
        FLT ftrFitErrorpc;			// trial value of fit error %
        FLT fA[3][3];				// ellipsoid matrix A
        FLT finvA[3][3];			// inverse of ellipsoid matrix A
        FLT fmatA[10][10];			// scratch 10x10 matrix used by calibration algorithms
        FLT fmatB[10][10];			// scratch 10x10 matrix used by calibration algorithms
        FLT fvecA[10];				// scratch 10x1 vector used by calibration algorithms
        FLT fvecB[4];				// scratch 4x1 vector used by calibration algorithms
        U32 iValidMagCal;			// integer value 0, 4, 7 denoting both valid calibration and solver used
    } magCalibration_t;

    U32 magCalibrationInit(magCalibration_t* const pMagCalibration, magneticBuffer_t* const pMagBuffer);
    U32 magBufferUpdate(magneticBuffer_t * const pMagBuffer, const FLT* const pMagRaw, const FLT* const pMagCal, const U32 loopCounter);
    U32 magCalibrationExec(magCalibration_t* const pMagCalibration, const magneticBuffer_t* const pMagBuffer);
    void magCorrection(FLT mag[], const magCalibration_t* const pMagCalibration);

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif