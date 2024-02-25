/**
*	@file	ISL_FileOperation.h
*	@brief	[Header file of File Operation Functions], 文件操作；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-02
*	@refer
*/

#ifndef PAI_FRAME_ISEISLIB_FILEOPERATION_H_
#define PAI_FRAME_ISEISLIB_FILEOPERATION_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
* @brief	从二进制文件中读取浮点型数组
*
* @param[in]	*fname				输入文件的文件名；
* @param[in]	num					输入数组的长度；
*
* @return	float * ，返回数据的指针，得到后请自己负责释放内存空间,数组的长度为num。
*/
ISL_EXPORT float * ISL_readFloatFromBinaryFile(char * fname, long num);


/**
* @brief	将一个浮点型数组保存为二进制文件
*
* @param[in]	*fname				输入文件的文件名；
* @param[in]	data[num]			输入的数组；
* @param[in]	num					输入数组的长度；
*
* @return	int，函数返回整型标志位，	若返回标志值小于0，则表示程序工作失败；
 *  							若返回标志值大于0；否则表示正常返回。
*/
ISL_EXPORT int ISL_writeFloatToBinaryFile(char * fname, float * data, long num);


/**
* @brief	清空一个二进制文件
*
* @param[in]	*fname				输入文件的文件名；
*
* @return	无
*/
ISL_EXPORT void ISL_clearBinaryFile(char * fname);



/**
* @brief	从一个文本文件逐行读取数据并放入浮点数组中
*
* @param[in]	*fname				输入文件的文件名；
* @param[in]	num					输入数组的长度；
*
* @return	float *，返回读取的数组，用完后，需自己释放内存空间，数组的长度为num。
*/
ISL_EXPORT float * ISL_readFloatFromTextFile(char * fname, long &num);


/**
* @brief	将一个浮点数组逐行写入一个文本文件中
*
* @param[in]	*fname				输入文件的文件名；
* @param[in]	data[num]			输入的数组；
* @param[in]	num					输入数组的长度；
*
* @return	无
*/
ISL_EXPORT void ISL_writeFloatToTextFile(char * fname, float * data, long num);



/**
* @brief	读取文件的大小
* @param[in] 	*f		FILE格式的文件指针；
* @return	无
*/
ISL_EXPORT long ISL_getFileSize(FILE * f);

}/* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_FILEOPERATION_H_ */
