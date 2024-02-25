/**
 *	@file	ISL_Alloc.h
 *	@brief	[Header file of Memory Alloc Functions], 内存管理头文件（SU中所用，大量内存分配请不要用）；
 *	@see	ISeisLib Manual
 *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
 *	@date	2014-02-26
 *	@refer	SU CWP
 */

#include "ISL_Alloc.h"

namespace ISLib
{

/* allocate a 1-d array */
void *alloc1(int n1, int size)
{
	void *p;

	if ((p = malloc(n1 * size)) == NULL)
		return NULL;

	return p;
}

/* allocate a 2-d array */
void **alloc2(int n1, int n2, int size)
{
	int i2;
	void **p;

	if ((p = (void**) malloc(n2 * sizeof(void*))) == NULL)
		return NULL;
	if ((p[0] = (void*) malloc(n2 * n1 * size)) == NULL)
	{
		free(p);
		return NULL;
	}
	for (i2 = 0; i2 < n2; i2++)
		p[i2] = (char*) p[0] + size * n1 * i2;

	return p;
}

/* allocate a 3-d array */
void ***alloc3(int n1, int n2, int n3, int size)
{
	int i3, i2;
	void ***p;

	if ((p = (void***) malloc(n3 * sizeof(void**))) == NULL)
		return NULL;
	if ((p[0] = (void**) malloc(n3 * n2 * sizeof(void*))) == NULL)
	{
		free(p);
		return NULL;
	}
	if ((p[0][0] = (void*) malloc(n3 * n2 * n1 * size)) == NULL)
	{
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3 = 0; i3 < n3; i3++)
	{
		p[i3] = p[0] + n2 * i3;
		for (i2 = 0; i2 < n2; i2++)
			p[i3][i2] = (char*) p[0][0] + size * n1 * (i2 + n2 * i3);
	}

	return p;
}

// === 根据数据类型 ===
/* allocate a 1-d array of ints */
int *alloc1int(int n1)
{
	return (int*) alloc1(n1, sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, int n1)
{
	return (int*) realloc1(v, n1, sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(int n1, int n2)
{
	return (int**) alloc2(n1, n2, sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**) p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(int n1, int n2, int n3)
{
	return (int***) alloc3(n1, n2, n3, sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***) p);
}

float *realloc1float(float *v, int n1)
{
	return (float*) realloc1(v, n1, sizeof(float));
}

/* allocate a 1-d array of floats */
float *alloc1float(int n1)
{
	return (float*) alloc1(n1, sizeof(float));
}

/* allocate a 2-d array of floats */
float **alloc2float(int n1, int n2)
{
	return (float**) alloc2(n1, n2, sizeof(float));
}

/* allocate a 3-d array of floats */
float ***alloc3float(int n1, int n2, int n3)
{
	return (float***) alloc3(n1, n2, n3, sizeof(float));
}

void *realloc1(void *v, int n1, int size)
{
	void *p;

	if ((p = realloc(v, n1 * size)) == NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1(void *p)
{
	free(p);
}

/* free a 2-d array */
void free2(void **p)
{
	free(p[0]);
	free(p);
}

/* free a 3-d array */
void free3(void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**) p);
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***) p);
}

double *alloc1double(int n1)
{
	return (double*) alloc1(n1, sizeof(double));
}

double *realloc1double(double *v, int n1)
{
	return (double*) realloc1(v, n1, sizeof(double));
}

void free1double(double *p)
{
	free1(p);
}

double **alloc2double(int n1, int n2)
{
	return (double**) alloc2(n1, n2, sizeof(double));
}

void free2double(double **p)
{
	free2((void**) p);
}

double ***alloc3double(int n1, int n2, int n3)
{
	return (double***) alloc3(n1, n2, n3, sizeof(double));
}

void free3double(double ***p)
{
	free3((void***) p);
}

// =======================大内存分配与删除======================================

/** @brief 分配整型二维数组 (new方式)
 *	@param[in]	n1	一维数组长度
 *	@param[in]	n2	二维数组长度
 *	@return	返回数组指针
 */
int **new2dInt(int n1, int n2)
{
	int **p = new int*[n1];
	for (int i = 0; i < n1; i++)
		p[i] = new int[n2];
	return p;
}

/** @brief 释放整型二维数组(delete方式)
 *  @param[in]	n1	一维数组长度
 *	@return	返回数组指针
 */
void delete2dInt(int **&p, int n1)
{
	for (int i = 0; i < n1; i++)
		delete[] p[i];
	delete[] p;
	p = NULL;
}

/** @brief 分配浮点型二维数组 (new方式)
 *	@param[in]	n1	一维数组长度
 *	@param[in]	n2	二维数组长度
 *	@return	返回数组指针
 */
float **new2dFloat(int n1, int n2)
{
	float **p = new float*[n1];
	for (int i = 0; i < n1; i++)
		p[i] = new float[n2];
	return p;
}

/** @brief 释放浮点型二维数组(delete方式)
 *  @param[in]	n1	一维数组长度
 *	@return	返回数组指针
 */
void delete2dFloat(float **&p, int n1)
{
	for (int i = 0; i < n1; i++)
		delete[] p[i];
	delete[] p;
	p = NULL;
}


/**
 * @brief	分配一维复数数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @return	返回数组指针
 */
complex *alloc1complex( size_t n1 )
{
	return (complex*) alloc1(n1, sizeof(complex));
}


/**
 * @brief	释放一维复数数组空间
 * @param[in]	*p	一维数组指针
 * @return	No
 */
void free1complex(complex *p)
{
	free1(p);
}

/**
 * @brief	分配二维复数数组空间（轻量级）
 * @param[in]	n1	一维数组长度
 * @param[in]	n2	二维数组长度
 * @return	返回数组指针
 */
complex **alloc2complex(size_t n1, size_t n2)
{
	return (complex**) alloc2(n1, n2, sizeof(complex));
}

/**
 * @brief	释放二维复数数组空间
 * @param[in]	**p	二维数组指针
 * @return	No
 */
void free2complex(complex **p)
{
	free(p[0]);
	free(p);
}

}/*End of ISLib*/

