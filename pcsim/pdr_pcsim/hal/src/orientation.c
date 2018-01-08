#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "orientation.h"


/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 getRelativeHeading(pedestrianOrientation_t* const pOrientation, FLT* const pHeading)
{   
    *pHeading = pOrientation->fDeviceHeading - pOrientation->fDeviceHeadingRef;

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void updateReferenceOrientation(pedestrianOrientation_t* const pOrientation, const ahrsFixData_t* const pAhrsData)
{
    pOrientation->uDeviceUpdateTime = pAhrsData->uTime;
    pOrientation->fDeviceHeadingRef = pOrientation->fDeviceHeading = pAhrsData->fPsiPl;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void updateDeviceOrientation(pedestrianOrientation_t* const pOrientation, const ahrsFixData_t* const pAhrsData)
{
    pOrientation->uDeviceUpdateTime = pAhrsData->uTime;
    pOrientation->fDeviceHeading = pAhrsData->fPsiPl;
    
}






















