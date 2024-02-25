/**
*	@file	ISL_FFT.h 
*	@brief	[Header file of FFT], FFT；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-27
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_FFT_H
#define PAI_FRAME_ISEISLIB_FFT_H

//#include "ISL_Export.h"
#include "ISL_UserDefine.h"
#include "ISL_Complex.h"


namespace ISLib {


/**
* @brief	*** 非SU FFT ***			
*/

/**
* @brief	傅立叶级数逼近
*
* @param[in]	f[2n+1]			存放区间[0,2pi]内的2n+1个等距点处的函数值
* @param[in]	n				等距点数为2n+l
* @param[out]	a[n+1]			返回博里叶级数中的系数ak（k=0，l，…，n）
* @param[out]	b[n+1]			返到博里叶级数中的系数bk（k=0，l，…，n）
*
* @return	no
*/
ISL_EXPORT void ISL_four(double *f, int n, double *a, double *b);


/**
* @brief	快速傅里叶变换(单精度版本)
*			本函数用于计算地震数据的振幅谱、相位谱、幅角、模，同时也可用于时频分析、旋转相位等			
*
* @param[in]	pr			当l=0时，存放n个采样的实部，返回离散傅里叶变换的模；
*							当l=1时，存放傅里叶变换的n个实部，返回傅里叶变换的模
* @param[in]	pi			当l=0时，存放n个采样的虚部，返回离散傅里叶变换的幅角；
*							当l=1时，存放傅里叶变换的n个虚部，返回傅里叶变换的幅角。（其中幅角的单位为度）
* @param[in]	n			输入数据的采样点数
* @param[in]	k			满足n=2^Exp   由于Sample和Exp存在这样的关系，应用时受到一定的限制
* @param[in]	l			当l=0时，表示要求本函数计算傅里叶变换；
*							当l=1时，表示要求本函数计算逆傅里叶变换
* @param[in]	il			当il=0时，表示不要求本函数计算傅里叶变换或逆变换的模与幅角；
*							当il=1时，表示要求本函数计算傅里叶变换或逆变换的模与幅角
*
* @param[out]	pr			当l=0时，存放n个采样的实部，返回离散傅里叶变换的模；
*							当l=1时，存放傅里叶变换的n个实部，返回傅里叶变换的模
* @param[out]	pi			当l=0时，存放n个采样的虚部，返回离散傅里叶变换的幅角；
*							当l=1时，存放傅里叶变换的n个虚部，返回傅里叶变换的幅角。（其中幅角的单位为度）
* @param[out]	fr			当l=0时，返回傅里叶变换的实部；
*							当l=1时，返回傅里叶变换的实部
* @param[out]	fi			当l=0时，返回傅里叶变换的虚部；
*							当l=1时，返回傅里叶变换的虚部	
*
* @return	no
*/
ISL_EXPORT void ISL_kfft(float *pr, float *pi, int n, int k, float *fr, float *fi, int l, int il);

/**
* @brief	快速傅里叶变换(双精度版本)，参数设置请参考单精度版本
*/
ISL_EXPORT void ISL_kfft_d(double *pr, double *pi, int n, int k, double *fr, double *fi, int l, int il);



/**
* @brief	快速傅里叶变换（双精度，替换输入序列）
*
* @param[in]	x			存放n个采样的实部，返回傅里叶变换的实部；
* @param[in]	y			存放n个采样的虚部，返回傅里叶变换的虚部；
* @param[in]	n			输入数据的采样点数（满足n=2^Exp）
* @param[in]	sign		sign = -1，为反变换
*
* @return	no
*/
ISL_EXPORT void ISL_kfft_d(double *pr, double *pi, int n, int sign);
			
/**
* @brief	FFT辅助函数，由真实样点数求得FFT函数所需的样点数及EXP			
*
* @param[in]	ns			真实的样点数
*
* @param[out]	nPoints		输出的FFT所需样点数， nPoints
* @param[out]	nexp		输出的FFT，nExp	
*
* @return	no
*/
ISL_EXPORT void ISL_calcKFFTSampleAndExp(int ns, int &sample, int &exp);


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////		华丽的分割线		///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

class ISL_EXPORT su_fft
{
public:
	su_fft()
	{

	}

	~su_fft()
	{

	}
	/**
	* @brief	*** SU FFT ***
	*/

	/**
	* @brief	返回适用于素因子FFT所需的最小的n值，该值不小于输入参数nmin
	* @param[in]	nmin		返回值的下限
	* @return	适用于素因子FFT所需的n值
	*			返回值n是由素因子集{2,3,4,5,7,8,9,11,13,16}中的因子组成。
	*			n不会超过720720 = 5*7*9*11*13*16，因此如果输入的nmin超过了720720，则返回值就是720720。
	*/
	static int ISL_npfa (int nmin);


	/**
	* @brief	返回适用于素因子FFT所需的最优n值，该值在输入值nmin和nmax之间
	* @param[in]	nmin		返回值的下限
	* @param[in]	nmax		期望的返回值上限（但是不一定确保是）
	* @return	适用于素因子FFT所需的n值
	*			返回值n是由素因子集{2,3,4,5,7,8,9,11,13,16}中的因子组成。
	*			n不会超过720720 = 5*7*9*11*13*16，因此如果输入的nmin超过了720720，则返回值就是720720。
	*			选择最优值n的目的是降低计算成本，同时确保n值不超过nmax。
	*/
	static int ISL_npfao (int nmin, int nmax);


	/**
	* @brief	返回适用于实数到复数或者复数到实数素因子FFT所需的最小的n值，该值不小于输入参数nmin
	* @param[in]	nmin		返回值的下限
	* @return	返回适用于实数到复数或者复数到实数素因子FFT所需的n值
	*			目前实数到复数或者复数到实数素因子FFT要求变换长度n是偶数，并且n/2也是满足复数到复数的素因子FFT
	*/
	static int ISL_npfar (int nmin);


	/**
	* @brief	返回适用于实数到复数或者复数到实数素因子FFT所需的最小的n值，该值在输入值nmin和nmax之间
	* @param[in]	nmin		返回值的下限
	* @param[in]	nmax		期望的返回值上限（但是不一定确保是）
	* @return	返回适用于实数到复数或者复数到实数素因子FFT所需的n值。
	*			目前实数到复数或者复数到实数素因子FFT要求变换长度n是偶数，并且n/2也是满足复数到复数的素因子FFT。
	*/
	static int ISL_npfaro (int nmin, int nmax);


	/**
	* @brief	素因子FFT:复数到复数的FFT
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	n			n是变换的长度
	* @param[in]	z			z是需要变换的复数序列，长度是n
	* @param[out]	z			z变换后的复数序列
	*							n必须是由素因子集{2,3,4,5,7,8,9,11,13,16}中的素因子构成的合数，即n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u，
								并且0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1。这就意味着n的范围是1 <= n <= 720720 (= 5*7*9*11*13*16)。
	* @return	no
	*/
	static void ISL_pfacc(int isign, int n, complex *cz);


	/**
	* @brief	素因子FFT：复数到实数的变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	n			n是变换的长度
	* @param[in]	cz			cz是需要变换的复数序列，长度是[n/2+1]（也可能等于rz）
	* @param[out]	rz			变换后的实值序列（也可能等于cz）
	*							由于ISL_pfacr利用了ISL_pfacc进行计算，因此n必须是偶数，并且n/2对ISL_pfacc来说是合法的。
	*							一个简单的方法是利用ISL_npfar来获得n。
	* @return	no
	*/
	static void ISL_pfacr(int isign, int n, complex *cz, float *rz);


	/**
	* @brief	素因子FFT：实数到复数的变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	n			n是变换的长度，必须是偶数
	* @param[in]	rz			rz是需要变换的实数序列（可能等于cz）
	* @param[out]	cz			array[n/2+1] of complex values (may be equivalenced to rz)
	*							cz是变换后的复数序列（可能等于rz）
	*
	* @return		no
	*							由于ISL_pfacr利用了ISL_pfacc进行计算，因此n必须是偶数，并且n/2对ISL_pfacc来说是合法的。
	*							一个简单的方法是利用ISL_npfar来获得n。
	*/
	static void ISL_pfarc(int isign, int n, float *rz, complex *cz);


	/**
	* @brief	素因子FFT：复式复数到复数变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	n			n是每次变换中的复数的个数
	* @param[in]	nt			nt变换次数
	* @param[in]	k			k是一次变换中复数元素的步长
	* @param[in]	kt			kt是两次变换中复数元素的步长
	* @param[in]	cz			z是需要进行变换的复数元序列
	* @param[out]	cz			z是变换后的复数元序列
	*
	* @return		no
	*							n必须是由素因子集{2,3,4,5,7,8,9,11,13,16}中的素因子构成的合数，即n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u，
	*							并且0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1。这就意味着n的范围是1 <= n <= 720720 (= 5*7*9*11*13*16)。
	*							对一个n1*n2的复数矩阵进行二维变换（假定n1和n2都是合法的“n”），并且n1是快方向，n2是慢方向，
	*							ISL_pfamcc(isign,n1,n2,1,n1,z)；（进行第一维变换）
	*							ISL_pfamcc(isign,n2,n1,n1,1,z)；（进行第二维变换）
	*
	*/
	static void ISL_pfamcc(int isign, int n, int nt, int k, int kt, complex *cz);


	/**
	* @brief	素因子FFT：2维复数到复数的变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	idim		变换的维数，要么是1，要么是2
	* @param[in]	n1			n1需要变换的第一维（快）序列
	* @param[in]	n2			n2需要变换的第二维（慢）序列
	* @param[in]	z			z需要变换的复数矩阵[n1][n2]
	*
	* @param[out]	z			z变换后的复数矩阵[n1][n2]
	*
	* @return		no
	*							实际只对2D矩阵进行了一维（第一维或者第二维）变换。
	*							如果idim等于1，则对n1个复数元素进行了n2次变换，如果idim等于2，则对n2个复数元素进行了n1次变换。
	*							虽然z是一个1维数组的形式，但是它可以表示一个2维数组[n1][n2]
	*							变换长度n可以是n1，也可以是n2，n必须是由素因子集{2,3,4,5,7,8,9,11,13,16}中的素因子构成的合数，即n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u，
	*							并且0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1。这就意味着n的范围是1 <= n <= 720720 (= 5*7*9*11*13*16)。
	*							对一个n1*n2的复数矩阵进行二维变换（假定n1和n2都是合法的“n”），并且n1是快方向，n2是慢方向，
	*							ISL_pfa2cc(isign,1,n1,n2,z)
	*							ISL_pfa2cc(isign,2,n1,n2,z)
	*/
	static void ISL_pfa2cc(int isign, int idim, int n1, int n2, complex *cz);


	/**
	* @brief	素因子FFT：2D复数到实数的变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	idim		变换的维数，要么是1，要么是2
	* @param[in]	n1			n1需要变换的第一维（快）序列
	* @param[in]	n2			n2需要变换的第二维（慢）序列
	* @param[in]	cz			cz需要变换的复数序列（可能等于rz）
	*
	* @param[out]	rz			实数序列（可能等于cz）
	*
	* @return		no
	*							如果idim等于1，则对n1个实数进行n2次n1/2+1复数变换；如果idim等于2，则对n2个实数进行n1次n2/2+1复数变换。
	*							虽然rz是一个1维数组的形式，但是它可以表示一个2维数组[n1][n2]。
	*							同理，根据idim，cz可以表示一个2维复数数组[n1/2+1][n2]或者[n1][n2/2+1]
	*							根据idim,n既可以是n1也可以是n2。由于ISL_pfa2cr利用ISL_pfa2cc进行计算，n必须是偶数并且n/2对于ISL_pfa2cc是合法的。
	*							一个获得合法的n简单的方法是利用利用ISL_npfar。
	*/
	static void ISL_pfa2cr(int isign, int idim, int n1, int n2, complex *cz, float *rz);


	/**
	* @brief	素因子FFT：2维实数到复数的变换
	* @param[in]	isign		isign的符号是Fourier核中的指数的符号
	* @param[in]	idim		变换的维数，要么是1，要么是2
	* @param[in]	n1			n1需要变换的第一维（快）序列
	* @param[in]	n2			n2需要变换的第二维（慢）序列
	* @param[in]	rz			rz需要变换的实数序列（可能等于cz）
	*
	* @param[out]	cz			复数序列（可能等于rz）
	*
	* @return		no
	*							如果idim等于1，则对n1/2+1个复数进行n2次n1个实数变换；如果idim等于2，则对n2/2+1个复数进行n1次n2个实数变换。
	*							虽然rz是一个1维数组的形式，但是它可以表示一个2维数组[n1][n2]。
	*							同理，根据idim，cz可以表示一个2维复数数组[n1/2+1][n2]或者[n1][n2/2+1]。
	*							根据idim,n既可以是n1也可以是n2。由于ISL_pfa2cr利用ISL_pfa2cc进行计算，n必须是偶数并且n/2对于ISL_pfa2cc是合法的。
	*							一个获得合法的n简单的方法是利用利用ISL_npfar。。
	*/
	static void ISL_pfa2rc(int isign, int idim, int n1, int n2, float *rz, complex *cz);
};


}/*End of ISLib*/

#endif
