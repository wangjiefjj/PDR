#include "ahrs.h"


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
