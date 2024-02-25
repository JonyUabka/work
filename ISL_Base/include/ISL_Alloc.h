/**
 *	@file	ISL_Alloc.h
 *	@brief	[Header file of Memory Alloc Functions], 内存管理头文件（SU中所用，大量内存分配请不要用）；
 *	@see	ISeisLib Manual
 *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
 *	@date	2014-02-26
 *	@refer	SU CWP
 */


/**
 * C语言跟内存申请相关的函数主要有 alloca,calloc,malloc,free,realloc,sbrk等。其中alloca是向栈申请内存,因此无需释放。
 * malloc分配的内存是位于堆中的,并且没有初始化内存的内容,因此基本上malloc之后,调用函数memset来初始化这部分的内存空间。
 * calloc则将初始化这部分的内存,设置为0。而realloc则对malloc申请的内存进行大小的调整。
 * 申请的内存最终需要通过函数free来释放。而sbrk则是增加数据段的大小; malloc/calloc/free基本上都是C函数库实现的,跟OS无关。
 *
 * C函数库内部通过一定的结构来保存当前有多少可用内存。如果程序malloc的大小超出了库里所留存的空间,那么将首先调用brk系统调用来增加可用空间,然后再分配空间。
 * free时,释放的内存并不立即返回给OS,而是保留在内部结构中。
 *
 * 可以打个比方: brk类似于批发,一次性的向OS申请大的内存,而malloc等函数则类似于零售,满足程序运行时的要求。
 * 这套机制类似于缓冲。使用这套机制的原因: 系统调用不能支持任意大小的内存分配(有的系统调用只支持固定大小以及其倍数的内存申请,这样的话,对于小内存的分配会造成浪费）;
 * 系统调用申请内存代价昂贵,涉及到用户态和核心态的转换。
*/

#ifndef PAI_FRAME_ISEISLIB_ALLOC_H
#define PAI_FRAME_ISEISLIB_ALLOC_H

#include "ISL_UserDefine.h"
#include "ISL_Complex.h"

namespace ISLib
{

// =======================轻量级内存分配与删除======================================
/**
 * @brief	分配一维数组空间（轻量级）
 * @param[in]	n1	数组长度
 * @param[in]	size	单个数据类型字节数
 * @return	返回数组的指针
 */
ISL_EXPORT void * alloc1(int n1, int size);

/**
 * @brief	分配二维整型数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @param[in]	size	单个数据类型字节数
 * @return	返回数组的指针
 */
ISL_EXPORT void ** alloc2(int n1, int n2, int size);

/**
 * @brief	分配三维数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @param[in]	n3	二维数组长度
 * @param[in]	size	单个数据类型字节数
 * @return	返回数组的指针
 */
ISL_EXPORT void *** alloc3(int n1, int n2, int n3, int size);

/**
 * @brief	释放一维数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
ISL_EXPORT void free1(void *p);

/**
 * @brief	释放二维数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
ISL_EXPORT void free2(void **p);

/**
 * @brief	释放三维数组空间
 * @param[in]	***p	三维数组指针
 * @return	No
 */
ISL_EXPORT void free3(void ***p);

/**
 * @brief	重新分配数组空间（轻量级）
 * @param[in]	*v	数组的指针
 * @param[in]	n1	一维数组长度
 * @param[in]	size	单个数据类型字节数
 * @return	返回重新分配的数组指针
 */
ISL_EXPORT void * realloc1(void *v, int n1, int size);

/**
 * @brief	分配一维整型数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT int * alloc1int(int n1);

/**
 * @brief	重新分配一维整型数组空间（轻量级）
 * @param[in]	*v	数组的指针
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT int * realloc1int(int *v, int n1);

/**
 * @brief	分配二维整型数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT int ** alloc2int(int n1, int n2);

/**
 * @brief	分配三维整型数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @param[in]	n3	三维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT int *** alloc3int(int n1, int n2, int n3);

/**
 * @brief	释放一维整型数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
ISL_EXPORT void free1int(int *p);

/**
 * @brief	释放二维整型数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
ISL_EXPORT void free2int(int **p);

/**
 * @brief	释放三维整型数组空间
 * @param[in]	***p	三维数组指针
 * @return	No
 */
ISL_EXPORT void free3int(int ***p);

/**
 * @brief	重新分配一维浮点数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT float *realloc1float(float *v, int n1);

/**
 * @brief	分配一维浮点数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT float * alloc1float(int n1);

/**
 * @brief	分配二维浮点数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT float ** alloc2float(int n1, int n2);

/**
 * @brief	分配三维浮点数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @param[in]	n3	三维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT float *** alloc3float(int n1, int n2, int n3);

/**
 * @brief	释放一维浮点数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
ISL_EXPORT void free1float(float *p);

/**
 * @brief	释放二维浮点数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
ISL_EXPORT void free2float(float **p);

/**
 * @brief	释放三维浮点数组空间
 * @param[in]	***p	三维数组指针
 * @return	No
 */
ISL_EXPORT void free3float(float ***p);

/**
 * @brief	重新分配一维双精度数组空间（轻量级）
 * @param[in]	*v	数组的指针
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT double * realloc1double(double *v, int n1);

/**
 * @brief	分配一维双精度数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT double * alloc1double(int n1);

/**
 * @brief	分配二维双精度数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT double ** alloc2double(int n1, int n2);

/**
 * @brief	分配三维双精度数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @param[in]	n3	三维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT double *** alloc3double(int n1, int n2, int n3);

/**
 * @brief	释放一维双精度数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
ISL_EXPORT void free1double(double *p);

/**
 * @brief	释放二维双精度数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
ISL_EXPORT void free2double(double **p);

/**
 * @brief	释放三维双精度数组空间
 * @param[in]	***p	三维数组指针
 * @return	No
 */
ISL_EXPORT void free3double(double ***p);


// =======================大内存分配与删除======================================

/** @brief 分配整型二维数组 (new方式)
 *	@param[in]	n1	一维数组长度
 *	@param[in]	n2	二维数组长度
 *	@return	返回数组指针
 */
ISL_EXPORT int **new2dInt(int n1, int n2);


/** @brief 释放整型二维数组(delete方式)
 *  @param[in]	n1	一维数组长度
 *	@return	返回数组指针
 */
ISL_EXPORT void delete2dInt(int **&p, int n1);


/** @brief 分配浮点型二维数组 (new方式)
 *	@param[in]	n1	一维数组长度
 *	@param[in]	n2	二维数组长度
 *	@return	返回数组指针
 */
ISL_EXPORT float **new2dFloat(int n1, int n2) ;


/** @brief 释放浮点型二维数组(delete方式)
 *  @param[in]	n1	一维数组长度
 *	@return	返回数组指针
 */
ISL_EXPORT void delete2dFloat(float **&p, int n1);


// =======================复数分配与删除======================================

/**
 * @brief	分配一维复数数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT complex *alloc1complex( size_t n1 );


/**
 * @brief	释放一维复数数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
ISL_EXPORT void free1complex(complex *p);


/**
 * @brief	分配二维复数数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @return	返回数组指针
 */
ISL_EXPORT complex **alloc2complex (size_t n1, size_t n2);


/**
 * @brief	释放二维复数数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
ISL_EXPORT void free2complex (complex **p);

}/*End of ISLib*/

#endif
