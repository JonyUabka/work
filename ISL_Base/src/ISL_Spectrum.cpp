///****************************************************************************
// **
// **				ISeis Lib [Base Algorithm] - 频谱分析
// **
// **				Writer By	Mouri Song
// **							Baihong Liu
// **							Qiang Yang
// **
// **				(Soft Center IGP)
// **
// **				DATA : 2014-03-19
// **
// ****************************************************************************/

///**
// *	@file	ISL_Spectrum.cpp
// *	@brief	[Header file of Spectrum Functions], 频谱分析；
// *	@see	ISeisLib Manual
// *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
// *			DATE : 2014-02-26
// */

//#include "ISL_Spectrum.h"
//#include "ISL_Transform.h"
//#include "ISL_Absorption.h"

//namespace ISLib
//{

//int Spectrum::ISL_powerSpectrum1(float *data, /* input seismic data */
//										int num, /* number of samples */
//										float dtime, /* sampling interval in milli-second */
//										float *&amps, float *&freqs, int flag)
//{
//	int i;
//	complex *ct = NULL;
//	int nfft, nf;
//	float df;

//	if (amps)
//	{
//		delete[] amps;
//		amps = NULL;
//	}

//	if (freqs)
//	{
//		delete[] freqs;
//		freqs = NULL;
//	}

//	nfft = su_fft::ISL_npfaro(num, num + num);
//	nf = nfft / 2 + 1;
//	amps = new float[nfft];
//	if (amps == NULL)
//	{
//		return 0;
//	}

//	freqs = new float[nf];
//	if (freqs == NULL)
//	{
//		delete[] amps;
//		amps = NULL;
//		return 0;
//	}

//	ct = new complex[nf];
//	if (ct == NULL)
//	{
//		delete[] amps;
//		amps = NULL;
//		delete[] freqs;
//		freqs = NULL;
//		return 0;
//	}

//	if (nfft > num)
//	{
//		for (i = 0; i < num; i++)
//			amps[i] = data[i];
//		for (i = num; i < nfft; i++)
//			amps[i] = 0.0;
//	}
//	else
//	{
//		for (i = 0; i < nfft; i++)
//			amps[i] = data[i];
//	}
//	su_fft::ISL_pfarc(1, nfft, amps, ct);

//	df = 1.0 / ((float) nfft * (dtime * 0.001));
//	float max = -99999.9;
//	for (i = 0; i < nf; i++)
//	{
//		amps[i] = ISL_rcabs(ct[i]);
//		freqs[i] = (float) i * df;
//		if (amps[i] > max)
//			max = amps[i];
//	}

//	delete[] ct;
//	ct = NULL;

//	if (!flag)
//		return nf;

//	for (i = 0; i < nf; i++)
//		amps[i] = amps[i] / max;

//	return nf;
//}

//int Spectrum::ISL_powerSpectrum2(float *data, /* input seismic data    */
//										int num, /* number of samples                 */
//										float dtime, /* sampling interval in milli-second */
//										float *&amps, /* amplitude values                  */
//										float *&freqs, /* frequency values                  */
//										float freq_step, /* frequency step                    */
//										int flag)
//{
//	if (amps)
//	{
//		delete[] amps;
//		amps = NULL;
//	}
//	if (freqs)
//	{
//		delete[] freqs;
//		freqs = NULL;
//	}

//	float *val1 = NULL, *val2 = NULL;

//	int nn = ISL_powerSpectrum1(data, num, dtime, val1, val2, flag);
//	if (val1 == NULL || val2 == NULL || nn <= 0)
//	{
//		if (val1)
//		{
//			delete[] val1;
//			val1 = NULL;
//		}
//		if (val2)
//		{
//			delete[] val2;
//			val2 = NULL;
//		}
//		return 0;
//	}

//	if (freq_step <= 0.0)
//	{
//		amps = val1;
//		freqs = val2;
//		return nn;
//	}

//	float fmin = val2[0];
//	float fmax = val2[nn - 1];
//	float fdelta = val2[1] - val2[0];
//	if (freq_step <= fdelta)
//	{
//		amps = val1;
//		freqs = val2;
//		return nn;
//	}

//	int count = (int) (fmin / freq_step);
//	fmin = (float) count * freq_step;
//	count = (int) ((fmax + freq_step) / freq_step);
//	fmax = (float) count * freq_step;

//	count = (fmax - fmin) / freq_step + 1;
//	float * xout = new float[count];

//	Interpolation::ISL_intlin_getx(nn, val2, count, xout);
//	amps = new float[count];
//	freqs = new float[count];
//	if (freqs == NULL || amps == NULL)
//	{
//		delete[] amps;
//		amps = NULL;
//		amps = val1;
//		freqs = val2;
//		return nn;
//	}

//	float max = 0, min = 0, top = 0;
//	ISL_findMaxValue(val1, nn, max, min, top);
//	Interpolation::ISL_intlin(nn, val2, val1, min, max, count, xout, amps, 1);

//	if(xout){ delete []xout; xout = NULL; }
//	delete[] val1;
//	val1 = NULL;
//	delete[] val2;
//	val2 = NULL;

//	for (int i = 0; i < count; i++)
//	{
//		freqs[i] = fmin + (float) i * freq_step;
//	}
//	return count;
//}




//}/* End Of namespace ISLIB */

