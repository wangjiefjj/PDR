#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "types.h"
#include "pdr.h"


#define  MAX_BUFF_LEN	 (1024)
#define  MAX_SECTION     (20)


#ifdef DEBUG
FILE *FpOutput = NULL;
FILE* FpAhrs = NULL;
FILE* FpStep = NULL;
#endif

static U32 praseData(char *str, pdrData_t *pdrData);
/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
int main(int argc,char *argv[])
{
	FILE *fp;
	U8 line[MAX_BUFF_LEN];

	fopen_s(&fp, "../../data/pdr_data.log", "r");
#ifdef DEBUG
    fopen_s(&FpOutput, "../../data/output.txt", "w");
    fopen_s(&FpAhrs, "../../data/ahrsData.txt", "w");
    fopen_s(&FpStep, "../../data/stepData.txt", "w");
#endif
	if (fp == NULL)
	{
		printf("open file failed!\r\n");
		return -1;
	}

    if (pdrNavInit())
    {
        printf("pdr system init failed!\r\n");
        return -1;
    }

	while (fgets(line, MAX_BUFF_LEN, fp) != 0 )
	{
        pdrData_t pdrData;
#ifdef DEBUG
        //printf("\r\n--------------------------------------------------------------------\r\n");
        //printf("new data come:\r\n");
        //puts(line);
#endif
        memset(&pdrData, 0, sizeof(pdrData_t));
        if (praseData(line, &pdrData))
        {
            printf("parsing failed!\r\n");
            return -1;
        }
        pdrNavExec(&pdrData);
    }
#ifdef DEBUG
    fclose(fp);
    fclose(FpOutput);
    fclose(FpAhrs);
    fclose(FpStep);
#endif

    getchar();
	return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
static U32 praseData(char *str, pdrData_t *pdrData)
{
    U32  count = 0;
    U32  flag = 0;
    char *section[MAX_SECTION] = {NULL};
    char *buffer = str;
    char *ptr = NULL;

    while( (section[count] = strtok_s(buffer, " ", &ptr)) != NULL )
    {
        count ++;
        buffer = NULL;
    }

    // distinguish data types
    flag = (U32)atoi(section[0]);
    pdrData->dataType = flag;
    
    switch (flag)
    {
        case GNSS_DATA:
            pdrData->gnssData.uTime = atoi(section[1]);
            pdrData->gnssData.fLatitude = atof(section[2]);
            pdrData->gnssData.fLongitude = atof(section[3]);
            pdrData->gnssData.fAltitude = (FLT)atof(section[4]);
            pdrData->gnssData.fVelE = (FLT)atof(section[5]);
            pdrData->gnssData.fVelN = (FLT)atof(section[6]);
            pdrData->gnssData.fVelU = (FLT)atof(section[7]);
            //pdrData->gnssData.fHeading = (FLT)atof(section[8]);
            break;

        case SENSOR_DATA:
            pdrData->sensorData.uTime = atoi(section[1]);
            pdrData->sensorData.fAcc[CHX] = (FLT)atof(section[2]);
            pdrData->sensorData.fAcc[CHY] = (FLT)atof(section[3]);
            pdrData->sensorData.fAcc[CHZ] = (FLT)atof(section[4]);
            pdrData->sensorData.fGyro[CHX] = (FLT)atof(section[5]);
            pdrData->sensorData.fGyro[CHY] = (FLT)atof(section[6]);
            pdrData->sensorData.fGyro[CHZ] = (FLT)atof(section[7]);
            pdrData->sensorData.fMag[CHX] = (FLT)atof(section[8]);
            pdrData->sensorData.fMag[CHY] = (FLT)atof(section[9]);
            pdrData->sensorData.fMag[CHZ] = (FLT)atof(section[10]);
            break;
    }

    return 0;
}