/**
*	@file	ISL_SU_Matrix.h
*	@brief	[Header file of SU Matrix Functions], SU中的矩阵计算函数；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-10
*	@refer	SU CWP
*/	

#ifndef PAI_FRAME_ISEISLIB_SUMATRIX_H
#define PAI_FRAME_ISEISLIB_SUMATRIX_H

#include "ISL_MathBase.h"
#include "ISL_Alloc.h"

namespace ISLib {

#define TINY 0.000001

/* 三对角系统的方程 */
/**
* @brief	Solve a tridiagonal system of equations (float version)
*			Functions to solve tridiagonal systems of equations Tu=r for u.
*			
* @param[in]	n	dimension of system	
* @param[in]	a[]	array[n] of lower sub-diagonal of T (a[0] ignored)
* @param[in]	b[]	array[n] of diagonal of T
* @param[in]	c[]	array[n] of upper super-diagonal of T (c[n-1] ignored
* @param[in]	r[]	array[n] of right-hand-side column vector
* @param[out]	u[]	array[n] of solution (left-hand-side) column vector
* @return	no
*/
ISL_EXPORT void ISL_tridif (int n, float a[], float b[], float c[], float r[], float u[]);


/**
* @brief	Solve a tridiagonal system of equations (double version)
*			
* @param[in]	n	dimension of system	
* @param[in]	a[]	array[n] of lower sub-diagonal of T (a[0] ignored)
* @param[in]	b[]	array[n] of diagonal of T
* @param[in]	c[]	array[n] of upper super-diagonal of T (c[n-1] ignored
* @param[in]	r[]	array[n] of right-hand-side column vector
* @param[out]	u[]	array[n] of solution (left-hand-side) column vector
* @return	no
*/
ISL_EXPORT void ISL_tridid (int n, double a[], double b[], double c[], double r[], double u[]);


/**
* @brief	Solve a positive definite, symmetric tridiagonal system
*
* @param[in]	n	dimension of system			
* @param[in]	d	array[n], the diagonal of A	
* @param[in]	e	array[n], the superdiagonal of A
* @param[in]	b	array[n], the rhs column vector of Ax=b
* @param[out]	b	b is overwritten with the solution to Ax=b
* @return	no
*/
ISL_EXPORT void ISL_tripd(float *d, float *e, float *b, int n);


/**
* @brief	Solve an unsymmetric tridiagonal system that uses
*			Gaussian elimination with partial pivoting
*			
* @param[in]	n	dimension of system	
* @param[in]	d	diagonal vector of matrix
* @param[in]	e	upper-diagonal vector of matrix
* @param[in]	c	lower-diagonal vector of matrix
* @param[in]	b	right-hand vector
* @param[out]	b	solution vector
* @return	no
*/
ISL_EXPORT void ISL_tripp(int n, float *d, float *e, float *c, float *b);


/* 求解线性方程组Ax=b*/
/**
* @brief	Decompose a matrix (A) into a lower triangular (L)
*			and an upper triangular (U) such that A=LU
*			
* @param[in]	nrows	number of rows of matrix to invert	
* @param[in]	matrix	matrix of coefficients in linear system Ax=b
*
* @param[out]	matrix	matrix containing LU decomposition (original matrix destroyed)
* @param[out]	vector	recording the row permutations effected by partial pivoting
* @param[out]	d		+/- 1 depending on whether the number of row interchanges was even or odd
*
* @return	no
*/
ISL_EXPORT void ISL_LU_decomposition (int nrows, float **matrix, int *idx, float *d);


/**
* @brief	Apply backward substitution to an LU decomposed
*			matrix to solve the linear system of equations Ax=b
*			
* @param[in]	nrows	number of rows (and columns) of input matrix	
* @param[in]	matrix	matrix of coefficients (after LU decomposition)
* @param[in]	idx		permutation vector obtained from routine LU_decomposition
* @param[in]	b		right hand side vector in equation Ax=b
*
* @param[out]	b		vector with the solution
*
* @return	no
*/
ISL_EXPORT void ISL_backward_substitution (int nrows, float **matrix, int *idx, float *b);


/**
* @brief	compute the inverse of a square non-singular matrix
*			
* @param[in]	nrows	number of rows (and columns) of input matrix	
* @param[in]	matrix	matrix to invert
*
* @param[out]	matrix		inverse of input matrix
*
* @return	no
*/
ISL_EXPORT void ISL_inverse_matrix (int nrows, float **matrix);


/**
* @brief	computes the product A^(-1)*B without explicitely
*			computing the inverse matrix
*			
* @param[in]	nrows1		number of rows (and columns) of matrix to invert
* @param[in]	matrix1		square matrix to invert
* @param[in]	ncols2		number of coulmns of second matrix
* @param[in]	nrows2		number of rows of second matrix
* @param[in]	matrix		second matrix (multiplicator)
*
* @param[out]	out_matrix	matrix containing the product of the inverse of the first
*							matrix by the second one.
*
* @return	no
*/
ISL_EXPORT void ISL_inverse_matrix_multiply (	int nrows1, float **matrix1, int ncols2,
									int nrows2, float **matrix2, float **out_matrix);

									
/**
* @brief	Toeplitz矩阵的线性方程的对称性计算
*			Solve a symmetric Toeplitz linear system of equations Rf=g for f(Double version)
*			
* @param[in]	n			dimension of system
* @param[in]	r[]			array[n] of top row of Toeplitz matrix
* @param[in]	g[]			array[n] of right-hand-side column vector
*
* @param[out]	f[]			array[n] of solution (left-hand-side) column vector
* @param[out]	a[]			array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
*
* @return	no
*/
ISL_EXPORT void ISL_stoepd(int n, double r[], double g[], double f[], double a[]);


/**
* @brief	Toeplitz矩阵的线性方程的对称性计算
*			Solve a symmetric Toeplitz linear system of equations Rf=g for f(float version)
*			
* @param[in]	n			dimension of system
* @param[in]	r[]			array[n] of top row of Toeplitz matrix
* @param[in]	g[]			array[n] of right-hand-side column vector
*
* @param[out]	f[]			array[n] of solution (left-hand-side) column vector
* @param[out]	a[]			array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
*
* @return	no
*/
ISL_EXPORT void ISL_stoepf(int n, float *r, float *g, float *f, float *a);


// === determine index of x with respect to an array of x values ===
/**
* @brief	determine index of x with respect to an array of x values（判断某数值是否在数组内）
*			
* @param[in]	nx			number of x values in array ax
* @param[in]	ax[]		array[nx] of monotonically increasing or decreasing x values
* @param[in]	x			the value for which index is to be determined
* @param[in]	index		index determined previously (used to begin search)
*
* @param[out]	index		for monotonically increasing ax values, the largest index
*							for which ax[index]<=x, except index=0 if ax[0]>x;
*							for monotonically decreasing ax values, the largest index
*							for which ax[index]>=x, except index=0 if ax[0]<x
*							返回某个数组的下标，整型指针返回
*
* @return	no
*/
ISL_EXPORT void ISL_xindex (int nx, float ax[], float x, int *index);


/**
* @brief	计算n阶广义Hermite多项式
*			Compute n-th order generalized Hermite polynomial.
*			
* @param[in]	*h0			array of size nt holding H_{n-1}
* @param[in]	*h1			array of size nt holding H_{n}
* @param[in]	*t			array of size nt holding time vector
* @param[in]	nt			size of arrays, no. of samples
* @param[in]	n			order of polynomial
* @param[in]	sigma		variance
*
* @param[out]	*h			array of size nt holding H_{n+1}
*
* @return	no
*/
ISL_EXPORT void ISL_hermite_n_polynomial(double *h, double *h0, double *h1, double *t, int nt, int n, double sigma);


/**
* @brief	计算的最小二乘最佳正弦内插系数
*			The coefficients are a least-squares-best approximation to the ideal
*			sinc function for frequencies from zero up to a computed maximum
*			frequency.  For a given interpolator length, lsinc, mksinc computes
*			the maximum frequency, fmax (expressed as a fraction of the nyquist
*			frequency), using the following empirically derived relation (from
*			a Western Geophysical Technical Memorandum by Ken Larner):
*
*			fmax = min(0.066+0.265*log(lsinc),1.0)
*
*			Note that fmax increases as lsinc increases, up to a maximum of 1.0.
*			Use the coefficients to interpolate a uniformly-sampled function y(i)
*			as follows:
*
*				lsinc-1
*			y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
*				j=0
*
*			Interpolation error is greatest for d=0.5, but for frequencies less
*			than fmax, the error should be less than 1.0 percent.
*			
* @param[in]	d			fractional distance to interpolation point; 0.0<=d<=1.0
* @param[in]	lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20
*
* @param[out]	sinc[]		array[lsinc] containing interpolation coefficients
*
* @return	no
*/
ISL_EXPORT void ISL_mksinc (float d, int lsinc, float sinc[]);


/**
* @brief	一维数组的归一化
*			
* @param[in]	in[]		输入的数组array[n]
* @param[in]	n			样点个数
* @param[in]	normv		归一值
*
* @param[out]	out[]		归一值后的数组array[n]
*
* @return	no
*/
ISL_EXPORT void ISL_normalize_1d(float in[], int n, float normv, float out[]);


/**
* @brief	计算定期采样，单调递增函数x（y）的定期采样，单调增加的函数y（x）的逆线性插值。
*			compute a regularly sampled function x(y) from a regularly
*			sampled, monotonically increasing function y(x)
*
*			(1) dx>0.0 && dy>0.0
*			(2) y[0] < y[1] < ... < y[nx-1]
*			
* @param[in]	nx			number of samples of y(x)
* @param[in]	dx			x sampling interval; dx>0.0 is required
* @param[in]	fx			first x
* @param[in]	y			array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
* @param[in]	ny			number of samples of x(y)
* @param[in]	dy			y sampling interval; dy>0.0 is required
* @param[in]	fy			first y
* @param[in]	xylo		x value assigned to x(y) when y is less than smallest y(x)
* @param[in]	xyhi		x value assigned to x(y) when y is greater than largest y(x)
*
* @param[out]	x			array[ny] of x(y) values
*
* @return	no
*/
ISL_EXPORT void ISL_yxtoxy(	int nx, float dx, float fx,
					float y[], int ny, float dy, float fy,
					float xylo, float xyhi, float x[]);

} /*End of ISLib*/

#endif
