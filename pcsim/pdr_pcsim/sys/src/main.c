#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "types.h"
#include "pdr.h"


#define  MAX_BUFF_LEN	 (1024)
#define  MAX_SECTION     (20)


#ifdef DEBUG
FILE *FpOutput;
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

	while ( fgets(line, MAX_BUFF_LEN, fp) != 0 )
	{
        pdrData_t pdrData;
#ifdef DEBUG
        printf("\r\n--------------------------------------------------------------------\r\n");
        printf("new data income:\r\n");
        puts(line);
#endif
        memset(&pdrData, 0, sizeof(pdrData_t));
        if ( praseData(line, &pdrData) )
        {
            printf("parsing failed!\r\n");
            return -1;
        }
        pdrNavExec(&pdrData);
    }
#ifdef DEBUG
    fclose(fp);
    fclose(FpOutput);
#endif

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
    pdrData->gnssData.uTime = atoi(section[1]);
    
    switch (flag)
    {
        case GNSS_DATA:
            pdrData->gnssData.fLatitude = atof(section[2]);
            pdrData->gnssData.fLongitude = atof(section[3]);
            pdrData->gnssData.fAltitude = (FLT)atof(section[4]);
            pdrData->gnssData.fVelE = (FLT)atof(section[5]);
            pdrData->gnssData.fVelN = (FLT)atof(section[6]);
            pdrData->gnssData.fVelU = (FLT)atof(section[7]);
            pdrData->gnssData.fHeading = (FLT)atof(section[8]);
            break;

        case SENSOR_DATA:
            pdrData->sensorData.fAcc[X] = (FLT)atof(section[2]);
            pdrData->sensorData.fAcc[Y] = (FLT)atof(section[3]);
            pdrData->sensorData.fAcc[Z] = (FLT)atof(section[4]);
            pdrData->sensorData.fGyro[X] = (FLT)atof(section[5]);
            pdrData->sensorData.fGyro[Y] = (FLT)atof(section[6]);
            pdrData->sensorData.fGyro[Z] = (FLT)atof(section[7]);
            pdrData->sensorData.fMag[X] = (FLT)atof(section[8]);
            pdrData->sensorData.fMag[Y] = (FLT)atof(section[9]);
            pdrData->sensorData.fMag[Z] = (FLT)atof(section[10]);
            break;
    }

    return 0;
}