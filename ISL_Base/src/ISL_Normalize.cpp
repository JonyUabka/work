/**
*	@file	ISL_Normalize.cpp
*	@brief	[Header file of Normalize], 归一化；
*	@see	ISeisLib Manual
*	@author [Song ZhiXiang], 宋志翔；
*	@date	2014-06-11
*	@refer	SU CWP
*/


#include "ISL_Normalize.h"

namespace ISLib {


// === 对FFT进行归一化 ===
/**********************************************************************************************************************
 *
 *  功能：对FFT的计算结果进行归一化
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nw					In			输入的复数数组的长度
 *		complex  *			cpfft				In/Out		输入的复数数组
 *		int					type				In			计算类型
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_normalizeFFT(complex *cpfft, int nw, int type)
{
	int i;
	float hl, hv;

	if (type < 1 || type > 5)
		return;

	hl = 0.0;
	for (i = 0; i < nw; i++) {
		hv = fabs(ISL_rcabs(cpfft[i]));//表示cpfft[i]的摸的绝对值
		if (hl < hv)
			hl = hv;
	}

	for (i = 0; i < nw; i++) {
		switch (type) {
		case 1:
			cpfft[i].r = cpfft[i].r / hl;//对cpfft[i]进行单位化
			cpfft[i].i = cpfft[i].i / hl;
			break;

		case 2:
			hv = ISL_rcabs(cpfft[i]);//表示cpfft[i]的摸
			cpfft[i] = ISL_cmplx(hv / hl, 0.0);//cpfft[i]=hv+i0.0
			break;

		case 3:
			hv = exp(cpfft[i].r);
			cpfft[i].r = hv * cos(cpfft[i].i);
			cpfft[i].i = hv * sin(cpfft[i].i);
			break;
		}
	}
}

/**
* @brief	重排子波
*
* @param[in]	nw			输入的复数数组的长度
* @param[in]	*cpfft		输入的复数数组（个数为nw）
* @param[in]	ntfft
* @param[in]	*pfft		输入子波的数组（个数为ntfft）
*
* @param[out]	*pfft		输出子波的数组,输入时请分配好内存空间
* @param[out]	len			输出子波数组的长度
*
* @return		no
*/
void ISL_arrangeWavelet(int nw, complex *cpfft, int ntfft, float *pfft, int &len)
{
	int i, j;

	if (cpfft == NULL || ntfft <= 0 || nw <= 0)
		return;
	if (pfft == NULL || len <= 0)
		return;

	ISL_normalizeFFT(cpfft, nw, 1);

	su_fft::ISL_pfacr(-1, ntfft, cpfft, pfft);

	float val0 = pfft[len / 2];
	for (i = len - 1, j = len / 2; i >= len / 2; i--, j--) {
		pfft[i] = pfft[j];
		if (val0 < pfft[i])
			val0 = pfft[i];
	}

	for (i = 0, j = len - 1; i < len / 2; i++, j--) {
		pfft[i] = pfft[j];
	}

	pfft[len / 2] = val0 / 0.987;
}


/**********************************************************************************************************************
 *
 *  功能：对地震子波做归一化处理
 *
 *  说明：将地震子波做频率域归一化处理
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float*				inData[waveLen]		In/Out		输入的数组/输出排序后的数据
 *		int					waveLen				In			输入数据的采样点数
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_normalizeWaveletFre(float *inData, int waveLen)
{
	/////对子波做频率域归一化处理/////
	int powerNumbers = 1, nExp;
	for (int i = 1;; i++) {
		powerNumbers *= 2;
		nExp = i;
		if (powerNumbers >= waveLen)
			break;
	}

	int half = powerNumbers / 2 + 1;
	float *SpecAm = new float[half];
	float *tmpRe = new float[powerNumbers];
	float *tmpIm = new float[powerNumbers];
	float *fftRe = new float[powerNumbers];
	float *fftIm = new float[powerNumbers];
	memset(tmpRe, 0, sizeof(float) * powerNumbers);
	memset(tmpIm, 0, sizeof(float) * powerNumbers);
	memset(fftRe, 0, sizeof(float) * powerNumbers);
	memset(fftIm, 0, sizeof(float) * powerNumbers);
	memcpy(tmpRe, inData, sizeof(float) * waveLen);
	ISL_kfft(tmpRe, tmpIm, powerNumbers, nExp, fftRe, fftIm, 0, 0);

	for (int i = 0; i < half; i++)
		SpecAm[i] = sqrt(fftRe[i] * fftRe[i] + fftIm[i] * fftIm[i]);
	float Max = 0.;
	for (int i = 0; i < half; i++)
		if (Max < SpecAm[i])
			Max = SpecAm[i];

	for (int i = 0; i < powerNumbers; i++) {
		fftRe[i] /= Max;
		fftIm[i] /= Max;
	}
	memset(tmpRe, 0, sizeof(float) * powerNumbers);
	memset(tmpIm, 0, sizeof(float) * powerNumbers);
	memcpy(tmpRe, fftRe, sizeof(float) * powerNumbers);
	memcpy(tmpIm, fftIm, sizeof(float) * powerNumbers);
	ISL_kfft(tmpRe, tmpIm, powerNumbers, nExp, fftRe, fftIm, 1, 0);
	memcpy(inData, fftRe, sizeof(float) * waveLen);

	delete[] SpecAm;
	delete[] tmpRe;
	delete[] tmpIm;
	delete[] fftRe;
	delete[] fftIm;
}

} /*End of ISLib*/
