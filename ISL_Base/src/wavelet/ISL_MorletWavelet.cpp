/**
 *	@file	ISL_Wavelet.h
 *	@brief	[Header file of Morlet Wavelet], 莫莱子波；
 *	@see	ISeisLib Manual
 *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
 *	@date	2014-06-13
 *	@refer	SU CWP
 */

#include "ISL_Wavelet.h"
#include "ISL_Normalize.h"
#include "ISL_Absorption.h"

namespace ISLib {

/* Morlet小波, 本函数用于求取Morlet小波,pis2为Pi的根, pis4为pis的根*/
void ISL_morlet(float dt, float mainfreq, int nt, int level, float *wavelet)
{
	int half;
	float A = 1.;
	float PiS2, PiS4;
	float ss2, sdt, e0, tmp, tmp1, tmp2;

	PiS2 = sqrt((float)PI);
	PiS4 = 1 / sqrt(PiS2);

	half = nt / 2 + 1;
	ss2 = sqrt(mainfreq);
	sdt = mainfreq * dt;

	tmp = float(-level * level / 2.);
	e0 = exp(tmp);

	for (int i = 0; i < nt; i++) {
		tmp1 = sdt * (i - half);
		tmp2 = exp(-tmp1 * tmp1 / 2);
		wavelet[i] = PiS4 * ss2 * tmp2 * (cos(tmp1 * level) - e0);
	}

	float tmpMax = wavelet[0];
	for (int i = 0; i < nt; i++){
		if (tmpMax < wavelet[i]){
			tmpMax = wavelet[i];
		}
	}

	for (int i = 0; i < nt; i++){
		wavelet[i] = A * wavelet[i] / tmpMax;
	}
}

/*Morlet小波,利用时移、相位、尺度来计算Morlet小波*/
bool ISL_morlet2(float dt, float mainfreq, int nt, float timeDelay, float phase, float scale, float *wavelet)
{
	if ((mainfreq <= 0.) || (timeDelay <= 0.) || (scale <= 0.))
		return false;

	float time, item2, item3, angle;
	float item1 = -log(2.f) / pow(PI, 2.0f);
	float Tmp = float(2.0f * PI * mainfreq);
	float Tmp0 = Tmp * timeDelay;

	for (int i = 0; i < nt; i++) {
		time = i * dt;
		item2 = (Tmp * time - Tmp0) / scale;
		item3 = pow(item2, float(2.0f));

		angle = Tmp * time - Tmp0 + phase;

		wavelet[i] = exp(item1 * item3) * cos(angle);
	}

	return true;
}


/*计算地震子波的幅度*/
float ISL_morletAmp(float *pdata, int num, float dt, float meanF, float timeDelay, float phase, float scale)
{
	float *data = new float[num];
	memset(data, 0, num * sizeof(float));
	float amp = 0.;
	if ( ISL_morlet2(dt, meanF, num, timeDelay, phase, scale, data) == true) {
		float val1 = ISL_innerProduct(pdata, data, num);
		float val3 = ISL_innerProduct(data, data, num);
		amp = val1 / val3;
	}

	if (data) {
		delete data;
		data = NULL;
	}
	return amp;
}

/*修改Morlet子波相关参数,通过匹配原则，不断的搜索子波*/
void ISL_refineMorlet(float *pData, float dt, int num,
						float &mainfreq, float &timeDelay, float &phase, float &scale,
						float offsetT, float offsetF, float offsetP, float offsetS)
{
	float minF = mainfreq - offsetF;
	float maxF = mainfreq + offsetF;
	float deltaF = 1.;
	minF = (minF < 3.f) ? (3.0f) : minF;
	maxF = (maxF > 100.f) ? (100.f) : maxF;

	float minTimeDelay = timeDelay - offsetT * dt;
	float maxTimeDelay = timeDelay + offsetT * dt;

	minTimeDelay = (minTimeDelay < 5 * dt) ? (5 * dt) : minTimeDelay;
	maxTimeDelay = (maxTimeDelay > (num - 5) * dt) ? ((num - 5) * dt) : maxTimeDelay;

	float minPhase = phase - offsetP * PI / 180.f;
	float maxPhase = phase + offsetP * PI / 180.f;
	float deltaPhase = 5.f * PI / 180.f;
	minPhase = (minPhase < -PI) ? (-PI) : minPhase;
	maxPhase = (maxPhase > PI) ? (PI) : maxPhase;

	offsetS = 0;
	//float minScale = scale-offsetS;
	//float maxScale = scale+offsetS;
	//float deltaScale = 0.02f;

	float maxVal = 0.;
	float *data = new float[num];
	//for(float s = minScale;s <= maxScale;s+=deltaScale)
	{
		float s = scale;
		//for(float t = minTimeDelay;t <= maxTimeDelay;t+=deltaTimeDelay)
		float t = timeDelay;
		for (float p = minPhase; p <= maxPhase; p += deltaPhase) {
			for (float f = minF; f <= maxF; f += deltaF) {
				memset(data, 0, num * sizeof(float));
				if (ISL_morlet2(dt, f, num, t, p, s, data) == true) {
					float val1 = ISL_innerProduct(pData, data, num);
					float val2 = ISL_innerProduct(data, data, num);
					val2 = sqrt(val2);
					float val = val1 / val2;
					if (val > maxVal) {
						maxVal = val;
						mainfreq = f;
						timeDelay = t;
						phase = p;
						scale = s;
					}
				}
			}
		}
	}

	if (data) {
		delete[] data;
		data = NULL;
	}
}

} /* End Of namespace ISLIB */
