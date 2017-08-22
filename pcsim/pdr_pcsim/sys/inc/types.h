#ifndef _TYPES_H_
#define _TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

	typedef		signed char				S8;                 /* Signed 8 bits integer    */
	typedef		unsigned char			U8;                 /* Unsigned 8 bits integer  */
	typedef		signed short			S16;                /* Signed 16 bits integer   */
	typedef		unsigned short			U16;                /* Unsigned 16 bits integer */
	typedef		signed int				S32;                /* Signed 32 bits integer   */
	typedef		unsigned int			U32;                /* Unsigned 32 bits integer */
	typedef		signed long long		S64;                /* Signed 64 bits integer   */
	typedef		unsigned long long		U64;                /* Unsigned 64 bits integer */
	typedef		float					FLT;                /* 4 bytes floating point   */
	typedef		double					DBL;                /* 8 bytes floating point   */

#define PARAMETER_NOT_USED(p) ((void)(p))

#ifdef __cplusplus
}      /* extern "C" */
#endif

#endif