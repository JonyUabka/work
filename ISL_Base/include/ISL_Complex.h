/**
*	@file	ISL_Complex.h 
*	@brief	[Header file of Complex Functions], 复数计算函数；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-26
*	@refer	SU CWP
*/	

#ifndef PAI_FRAME_ISEISLIB_COMPLEX_H
#define PAI_FRAME_ISEISLIB_COMPLEX_H

#include "ISL_UserDefine.h"

namespace ISLib {

	/**
	* @brief	complex 结构类型定义
	* @param	r	实数	
	* @param	i	虚数
	*/
	typedef struct _complexStruct { /* complex number*/
		float r;
		float i;
	} complex;
	
	/**
	* @brief	利用complex结构进行复数计算
	*/
	   
	/**
	* @brief	构造一个 complex 结构
	* @param[in]	re	实数	
	* @param[in]	im	虚数
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cmplx(float re, float im);

	/**
	* @brief	计算 complex 的模
	* @param[in]	z	一个复数结构	
	* @return	返回一个 complex 模（浮点型）
	*/
	ISL_EXPORT float ISL_rcabs(complex z);

	/**
	* @brief	complex 的加法运算
	* @param[in]	a	复数a	
	* @param[in]	b	复数b
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cadd(complex a, complex b);

	/**
	* @brief	complex 的减法运算
	* @param[in]	a	复数a	
	* @param[in]	b	复数b
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_csub(complex a, complex b);

	/**
	* @brief	complex 的乘法运算
	* @param[in]	a	复数a	
	* @param[in]	b	复数b
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cmul(complex a, complex b);

	/**
	* @brief	complex 的除法运算
	* @param[in]	a	复数a	
	* @param[in]	b	复数b
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cdiv(complex a, complex b);

	/**
	* @brief	计算公扼的complex
	* @param[in]	z	复数z	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_conjg(complex z);

	/**
	* @brief	计算负的complex
	* @param[in]	z	复数z	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cneg(complex z);

	/**
	* @brief	计算归一化的complex
	* @param[in]	z	复数z	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cinv(complex z);

	/**
	* @brief	complex 的开方运算
	* @param[in]	z	复数z	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_csqrt(complex z);

	/**
	* @brief	complex 的指数运算
	* @param[in]	z	复数z	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cexp(complex z);

	/**
	* @brief	complex 乘常数运算
	* @param[in]	a	复数a	
	* @param[in]	x	浮点数x，常数	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_crmul(complex a, float x);

	/**
	* @brief	complex 的对数运算
	* @param[in]	b	复数b	
	* @return	返回一个 complex 结构
	*/
	ISL_EXPORT complex ISL_cln(complex b);


} /*End of ISLib*/

#endif
