//*
//*	@file	ISL_TFAnalysis.cpp
//*	@brief	[Header file of time-frequency analysis Functions], 时频分析；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
//*	@date	2014-06-09
//*	@refer	SU CWP
//*/

//#include "ISL_TFAnalysis.h"
//#include "ISL_Absorption.h"
//#include "ISL_FFT.h"
//#include "ISL_Wavelet.h"
//#include "ISL_Normalize.h"
//#include "ISL_Transform.h"


//namespace ISLib {

///**********************************************************************************************************************
// *
// *  功能：广义S变换
// *
// *  说明：输入一个信号，输出信号的广义S变换时频谱
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		float*				inSignal[nPoints]	In			输入的地震信号
// *		int					nPoints				In			输入地震信号的长度
// *		float				fDT					In			输入地震信号的采样间隔(秒)
// *		float				fAlpha				In			广义S变换的参数(0.3f)
// *		int					nDelayTime			In			子波的时间延迟
// *		int					nFreqStrat			In			低通频
// *		int					nFreqEnd			In			高截频
// *		int					nFreqInterval		In			频率间隔
// *		float**				gstAmpSpec[nPoints][nF]		Out			输出广义S变换求得的时频谱,
// *																	nF = （起始频率-终止频率）/频率间隔+1
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_gsTransform(float * inSignal,
//						int nPoints,
//						float fDT,
//						float fAlpha,
//						int nDelayTime,
//						int nFreqStrat,
//						int nFreqEnd,
//						int nFreqInterval,
//						float ** outAmpSpec )
//{

//	int k, half;
//	int pNum = 1, nexp;
//	float Pi = (float)3.1415926f;
//	for (int i = 1;; i++) {
//		pNum *= 2;
//		nexp = i;
//		if (pNum >= nPoints)
//			break;
//	}

//	float *Gauss = new float[pNum];
//	float *fSigRe = new float[pNum];
//	float *fSigIm = new float[pNum];
//	float *fSpecRe = new float[pNum];
//	float *fSpecIm = new float[pNum];
//	float *stRe = new float[pNum];
//	float *stIm = new float[pNum];

//	memset(fSigRe, 0, sizeof(float) * pNum);
//	memset(fSigIm, 0, sizeof(float) * pNum);
//	memset(fSpecRe, 0, sizeof(float) * pNum);
//	memset(fSpecIm, 0, sizeof(float) * pNum);
//	memset(stRe, 0, sizeof(float) * pNum);
//	memset(stIm, 0, sizeof(float) * pNum);
//	memcpy(fSigRe, inSignal, sizeof(float) * nPoints);
//	ISL_kfft(fSigRe, fSigIm, pNum, nexp, fSpecRe, fSpecIm, 0, 0);

//	///// === 进行hilbert变换 === /////
//	half = pNum / 2;
//	for (int i = 1; i < pNum; i++) {
//		if (i <= half) {
//			fSpecRe[i] = 2 * fSpecRe[i];
//			fSpecIm[i] = 2 * fSpecIm[i];
//		} else {
//			fSpecRe[i] = 0.;
//			fSpecIm[i] = 0.;
//		}
//	}

//	int nLowPass = (int) (nFreqStrat * pNum * fDT);
//	float nfre = nLowPass;
//	if (nfre == 0)
//		nfre = 1.;

//	int nF = int((nFreqEnd - nFreqStrat) / nFreqInterval + 1);
//	float freStepLen = float(1.0 / pNum / fDT);
//	float fSep = float(nFreqInterval / freStepLen);

//	for (int ii = 0; ii < nF; ii++) {
//		Gauss[0] = 1;
//		for (int i = 1; i <= half; i++){
//			Gauss[i] = Gauss[pNum - i] = exp(-2.f * Pi * Pi * fAlpha * (i - nDelayTime) * (i - nDelayTime) / (nfre * nfre));
//		}

//		for (int i = 0; i < pNum; i++) {
//			k = int(nfre + i);
//			if (k >= pNum)
//				k = k - pNum;
//			fSigRe[i] = fSpecRe[k] * Gauss[i];
//			fSigIm[i] = fSpecIm[k] * Gauss[i];
//		}

//		ISL_kfft(fSigRe, fSigIm, pNum, nexp, stRe, stIm, 1, 0);

//		for (int j = 0; j < nPoints; j++){
//			outAmpSpec[j][ii] = sqrt(stRe[j] * stRe[j] + stIm[j] * stIm[j]);
//		}
//		nfre += fSep;
//	}

//	if(Gauss) { delete[] Gauss; Gauss = NULL; }
//	if(fSigRe) { delete[] fSigRe; fSigRe = NULL; }
//	if(fSigIm) { delete[] fSigIm; fSigIm = NULL; }
//	if(fSpecRe) { delete[] fSpecRe; fSpecRe = NULL; }
//	if(fSpecIm) { delete[] fSpecIm; fSpecIm = NULL; }
//	if(stRe) { delete[] stRe; stRe = NULL; }
//	if(stIm) { delete[] stIm; stIm = NULL; }
//}



///**********************************************************************************************************************
// *
// *  功能：三参数小波变换
// *
// *  说明：输入一个信号，输出信号的小波变换时频谱
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		float*				inSignal			In			地震信号
// *		int					nPoints				In			地震信号的长度
// *		float				DT					In			地震信号的采样间隔(秒)
// *		int					nScales				In			最大尺度
// *		int					freStart			In			起始频率
// *		int					nSep				In			尺度间隔
// *		float*				fWavRe				In			母小波的实部
// *		float*				fWavIm				In			母小波的虚部
// *		int					nWav				In			母小波的长度
// *		float				fFre				In			中心频率
// *		float**				outAmpSpec[nPoints][nF]			Out			输出三参数小波变换求得的时频谱,nF = （起始频率-终止频率）/频率间隔+1
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_triParWTransform(float *inSignal,
//							int nPoints,
//							float fDT,
//							int nFreqStrat,
//							int nFreqEnd,
//							int nFreqInterval,
//							float *fWavRe,
//							float *fWavIm,
//							int nWav,
//							float fFre,
//							float **outAmpSpec)
//{
//	/////小波的参数/////
//	int MaxTime = 8;    //母波的最大时间
//	int nWav1 = nWav - 1;
//	float df = 1.;    //频率域采样间隔
//	float Fre0 = 1;    //初始频率

//	int nn, nTmp, index;
//	float Sum1, Sum2, Coff;

//	/////计算最大尺度:int(Fre/初始频率/deltaT)+1/////
//	index = int((fFre / (Fre0 * fDT)) + 1) * MaxTime;
//	nTmp = index + nPoints;
//	float *fSWavRe = new float[index + 1];
//	float *fSWavIm = new float[index + 1];
//	float *fConvRe = new float[nTmp];
//	float *fConvIm = new float[nTmp];
//	float *fDeriRe = new float[nTmp - 1];
//	float *fDeriIm = new float[nTmp - 1];
//	memset(fSWavRe, 0, 4 * (index + 1));
//	memset(fSWavIm, 0, 4 * (index + 1));
//	memset(fConvRe, 0, 4 * (nTmp));
//	memset(fConvIm, 0, 4 * (nTmp));
//	memset(fDeriRe, 0, 4 * (nTmp - 1));
//	memset(fDeriIm, 0, 4 * (nTmp - 1));

//	int nF = int((nFreqEnd - nFreqStrat) / nFreqInterval + 1);

//	float TmpSum = 0.;
//	for (int fre = 0; fre < nF; fre++) {
//		TmpSum = float(fDT * ((fre + nFreqStrat) * df + Fre0));
//		Coff = float(fFre / TmpSum);
//		nTmp = int(MaxTime * Coff / nFreqInterval);

//		/////计算不同尺度下子波//////
//		for (int i = 0; i < nTmp + 1; i++) {
//			if (nTmp != 0) {
//				index = int(i * nWav1 / nTmp);
//				fSWavRe[nTmp - i] = fWavRe[index];
//				fSWavIm[nTmp - i] = fWavIm[index];
//			} else {
//				fSWavRe[0] = fWavRe[0];
//				fSWavIm[0] = fWavIm[0];
//			}
//		}

//		/////不同尺度下子波与输入信号卷积运算/////
//		for (int i = 0; i < nTmp + nPoints; i++) {
//			Sum1 = 0.;
//			Sum2 = 0.;
//			for (int j = 0; j <= i; j++) {
//				if (j < nPoints && i - j >= 0 && i - j < nTmp + 1) {
//					index = i - j;
//					Sum1 += inSignal[j] * fSWavRe[index];
//					Sum2 += inSignal[j] * fSWavIm[index];
//				}
//			}
//			fConvRe[i] = Sum1;
//			fConvIm[i] = Sum2;
//		}

//		/////求导/////
//		for (int i = 0; i < nTmp + nPoints - 1; i++) {
//			fDeriRe[i] = fConvRe[i + 1] - fConvRe[i];
//			fDeriIm[i] = fConvIm[i + 1] - fConvIm[i];
//		}

//		/////计算时频谱,把fDeriRe序列截成长度与输入信号一样的长度/////
//		nn = (int) ((nTmp - 1) / 2.);
//		for (int i = 0; i < nPoints; i++) {
//			index = nn + i;
//			outAmpSpec[i][fre] = sqrt(Coff)
//					* sqrt(fDeriRe[index] * fDeriRe[index] + fDeriIm[index] * fDeriIm[index]);
//		}
//	}

//	delete[] fSWavRe;
//	delete[] fSWavIm;
//	delete[] fConvRe;
//	delete[] fConvIm;
//	delete[] fDeriRe;
//	delete[] fDeriIm;
//}



///*短时傅里叶变换, 输入一个信号，输出信号的短时傅里叶变换的单频谱*/
//void ISL_uniFreSTFTransform(float *inSignal,
//							int nPoints,
//							float fDT,
//							float *inWinData,
//							int nWinLength,
//							int nFreqStrat,
//							int nFreqEnd,
//							int nFreqInterval,
//							float ** outAmpSpec)
//{
//	int Lh = int(0.5 * (nWinLength - 1));
//	int pNum = 1, nexp, index;
//	for (int i = 1;; i++) {
//		pNum *= 2;
//		nexp = i;
//		if (pNum >= nPoints)
//			break;
//	}

//	int nF = int((nFreqEnd - nFreqStrat) / nFreqInterval + 1);

//	float freStepLen = float(1.0 / pNum / fDT);
//	float fSep = float(nFreqInterval / freStepLen);

//	float *tmpRe = new float[pNum];
//	float *tmpIm = new float[pNum];
//	float *fftRe = new float[pNum];
//	float *fftIm = new float[pNum];

//	int L1, L2;
//	for (int j = 0; j < nPoints; j++) {
//		memset(tmpRe, 0, 4 * pNum);
//		memset(tmpIm, 0, 4 * pNum);
//		memset(fftRe, 0, 4 * pNum);
//		memset(fftIm, 0, 4 * pNum);

//		L1 = int(-min(min(pNum / 2 - 1, Lh), j));
//		L2 = int(min(min(pNum / 2 - 1, Lh), nPoints - j - 1));

//		for (int i = L1; i <= L2; i++) {
//			if (i < 0)
//				fftRe[pNum + i] = inSignal[j + i] * inWinData[i + Lh];
//			else
//				fftRe[i] = inSignal[j + i] * inWinData[i + Lh];
//		}

//		ISL_kfft(fftRe, fftIm, pNum, nexp, tmpRe, tmpIm, 0, 0);

//		for (int i = 0; i < nF; i++) {
//			index = int((i + nFreqStrat) * fSep);
//			outAmpSpec[j][i] = sqrt(tmpRe[index] * tmpRe[index] + tmpIm[index] * tmpIm[index]);
//		}
//		//cout<<AmpSpec[j][100]<<endl;
//	}

//	delete[] tmpRe;
//	delete[] tmpIm;
//	delete[] fftRe;
//	delete[] fftIm;
//}


///**********************************************************************************************************************
// *
// *  功能：卷积
// *
// *  说明：输入两个任意长度的离散点序列，返回一个卷积后的结果
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		float*				r					Int			输入的反射系数序列
// *		int					nr					Int			反射系数的维数
// *		float*				w					Int			输入的地震子波序列
// *		int					nw					Int			地震子波的维数
// *		float*				syn					Out			合成的地震记录
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_WTF_convolution(float *r, int nr, float *w, int nw, float *syn)
//{
//	int k;
//	float sum;
//	for (int i = 0; i < nr; i++) {
//		sum = 0.0;
//		for (int j = 0; j < nw; j++) {
//			k = i - j + nw / 2;
//			if (k >= 0 && k < nr)
//				sum += r[k] * w[j];
//		}
//		syn[i] = sum;
//	}
//}


///*小波变换, 输入一个信号，输出信号的小波变换的单频谱*/
//void ISL_uniFreWTransform(float *inSignal,
//							int nPoints,
//							float fDT,
//							int nFreqStrat,
//							int nFreqEnd,
//							int nFreqInterval,
//							int nWaveletLen,
//							float nWaveletDT,
//							float ** outAmpSpec,
//							int nGrade,
//							int nWFlag)
//{
//	float *seismic = new float[nPoints];
//	float *complexIm = new float[nPoints];
//	fDT = 0.002;

//	int Fre1 = nFreqStrat;
//	int Fre2 = nFreqEnd;
//	int fSep = nFreqInterval;
//	int flag = nWFlag;
//	int nFre = int((Fre2 - Fre1) / fSep + 1);

//	float waveLet[nWaveletLen];

//	int fre = Fre1 + fSep;
//	for (int j = 0; j < nFre; j++) {
//		if (flag == 1)
//			ISL_classicRicker(nWaveletLen, nWaveletDT, fre,  waveLet);
//		else if (flag == 2)
//			ISL_broadBandRicker(nWaveletLen, nWaveletDT, fre, fre + fSep, waveLet);
//		else if(flag == 3)
//			ISL_morlet(nWaveletDT, fre, nWaveletLen, nGrade, waveLet);

//		ISL_normalizeWaveletFre(waveLet, nWaveletLen);
//		ISL_WTF_convolution(inSignal, nPoints, waveLet, nWaveletLen, seismic);
//		ISL_hilbert(nPoints, seismic, complexIm);

//		for (int i = 0; i < nPoints; i++)
//			outAmpSpec[i][j] = sqrt(seismic[i] * seismic[i] + complexIm[i] * complexIm[i]);
//		fre += fSep;
//	}

//	delete[] seismic;
//	delete[] complexIm;
//}


///*匹配追踪, 对输入的信号进行时频分析分解，方法为匹配追踪方法，得到的时频分辨率较高。*/
///**
// * @brief	匹配追踪, 对输入的信号进行时频分析分解，方法为匹配追踪方法，得到的时频分辨率较高 ?
// *
// * @param[in]	inData[nPoints]				输入的地震信号；
// * @param[in]	nPoints						输入地震信号的长度；
// * @param[in]	fDT							地震信号的采样率(秒)
// * @param[in]	nIter						迭代次数
// * @param[in]	disS
// * @param[in]	disT
// * @param[in]	nLowFre						最小频率
// * @param[in]	nHighFre					最大频率
// * @param[in]	intFre						频率间隔
// * @param[in]	nBor						镶边长度
// * @param[out]	OutWav[((nHighFre - nLowFre) / intFre + 1)*nPoints]						输出的搜出来的子波
// * @param[out]	AmpWav[iter]				输出的子波的振幅
// * @param[out]	FreWav[iter]				输出的子波的频率
// * @param[out]	MPAmpSpec[((nHighFre - nLowFre) / intFre + 1)*nPoints]					输出的信号的时频谱
// *
// * @return	no
// */
///*
//void ISL_MPDTFS(float *inData, int nPoints, float fDT,
//				int nIter, float disS, int disT,
//				int disFre,	int disPh,
//				int nLowFre, int nHighFre, int intFre,
//				int nBor,
//				float *OutWav, float *AmpWav,
//				float *FreWav, float *MPAmpSpec)
//{
//	int loc;
//	float tmp, amp, scale, tShift, fPeak, phase;

//	int nFreNums = int(nHighFre - nLowFre) / intFre + 1;
//	float *insAmp = new float[nPoints];
//	float *insPh = new float[nPoints];
//	float *insFre = new float[nPoints];
//	float *tmpSei = new float[nPoints];
//	float *wav = new float[nPoints];
//	float *TFSpec = new float[nFreNums * nPoints];

//	memset(TFSpec, 0, sizeof(float) * nFreNums * nPoints);
//	memcpy(tmpSei, inData, sizeof(float) * nPoints);
//	float temp;
//	for (int i = 0; i < nBor; i++) {
//		temp = float(
//				exp(-0.5 * ((i - 0.5 * (nPoints - 1)) / (0.5 * 0.5 * (nPoints - 1)))
//					* ((i - 0.5 * (nPoints - 1)) / (0.5 * 0.5 * (nPoints - 1)))));    ////高斯窗
//		tmpSei[i] *= temp;
//		tmpSei[nPoints - i - 1] *= temp;
//	}

//	float max = 0;
//	long pos = 0;
//	for (int iter = 0; iter < nIter; iter++) {
//		ISL_insPhAmpFre(tmpSei, nPoints, fDT, 5, insAmp, insPh, insFre);
//		ISL_locateMax(insAmp, nPoints, max, pos);
//		loc = pos;

//		tShift = loc * fDT;
//		fPeak = insFre[loc];
//		phase = insPh[loc];
//		if (loc == 0 || loc == nPoints - 1)
//			tmpSei[loc] = 0.;

//		scale = 1.f;
//		//scale = getMorScale(tmpSei,dt,nPoints,fPeak,tShift,phase);
//		ISL_refineMorlet(tmpSei, fDT, nPoints, fPeak, tShift, phase, scale, disT, disFre, disPh, disS);
//		amp = ISL_morletAmp(tmpSei, fDT, nPoints, fPeak, tShift, phase, scale);
//		ISL_morlet2(fDT, fPeak, nPoints, tShift, phase, scale, wav);
//		ISL_timeFreqSpec(wav, fDT, nPoints, amp, fPeak, tShift, phase, scale, nLowFre, nHighFre,
//						intFre, TFSpec);

//		AmpWav[iter] = amp;
//		FreWav[iter] = fPeak;
//		for (int i = 0; i < nPoints; i++) {
//			tmp = amp * wav[i];
//			tmpSei[i] -= tmp;
//			OutWav[i + iter * nPoints] = tmp;
//		}
//	}

//	memcpy(inData, tmpSei, sizeof(float) * nPoints);
//	for (int j = 0; j < nFreNums; j++)
//		for (int i = 0; i < nPoints; i++)
//			MPAmpSpec[i + j * nPoints] = TFSpec[j + i * nFreNums];

//	delete[] insAmp;
//	delete[] insPh;
//	delete[] insFre;
//	delete[] tmpSei;
//	delete[] wav;
//	delete[] TFSpec;
//}
//*/

//} /* End Of namespace ISLIB
