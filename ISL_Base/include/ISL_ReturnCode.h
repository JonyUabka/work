/**
*	@file	ISL_ReturnCode.h 
*	@brief	[Header file of ISeisLib Return Codes], 算法库返回码头文件；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-26
*	@refer	SU CWP
*/	

#ifndef PAI_FRAME_ISEISLIB_RETURNCODE_H
#define PAI_FRAME_ISEISLIB_RETURNCODE_H

#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <assert.h>

using namespace std;

/**
* @brief	Linux 与 windows 都能够通用的库设置
* @param	No
* @return	No
*/

/*
#if defined(WIN32) || defined(WIN64) || defined(_WINDOWS)
    #ifdef  ISLBase_DLL
        #define ISLBase_EXPORT __declspec(dllexport)
    #else
        #define ISLBase_EXPORT __declspec(dllimport)
    #endif
	
	#ifdef  ISLAdv_DLL
        #define ISLAdv_EXPORT __declspec(dllexport)
    #else
        #define ISLAdv_EXPORT __declspec(dllimport)
    #endif
	
#else
    #define ISLBase_EXPORT
	#define ISLAdv_EXPORT
#endif
*/

#include "ISL_Export.h"

/**	
* @brief	[命名空间]：ISLib (ISeisLib算法库命名空间)
* 			[使用方式]：using namespace <命名空间名称>
* @param	No
* @return	No
*/
namespace ISLib {



} /*End of ISLib*/

#endif
