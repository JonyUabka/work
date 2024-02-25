/**
*	@file	ISL_Absorption.h
*	@brief	[Header file of Absorption Functions], 吸收系数；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_ABSORPTION_H_
#define PAI_FRAME_ISEISLIB_ABSORPTION_H_

#include "ISL_Transform.h"

namespace ISLib {


/**
* @brief	计算相位
*
* @param[in]	pRe			输入谱的实部；
* @param[in]	pIm			输入谱的虚部；
* @param[out]	pPh			输出相位;
*
* @return	no
*/
ISL_EXPORT void ISL_calcPh(float pRe, float pIm, float& pPh);

/**
* @brief	计算频率
*
* @param[in]	ph1			输入的相位值1；
* @param[in]	ph2			输入的相位值2；
* @param[in]	dt			输入的采样间隔；（秒）
* @param[out]	fn			输出计算得到的频率值;
*
* @return	no
*/
ISL_EXPORT void ISL_calcFreq(float ph1, float ph2 ,float dt, float& fn);



/**
* @brief	计算两离散序列的内积
*
* @param[in]	inData1			离散序列1；
* @param[in]	inData2			离散序列2；
* @param[in]	num				离散序列长度;
* @param[out]	bABS			是否进行ABS，true=进行， false=不进行;
*
* @return	内积(float)
*/
ISL_EXPORT float ISL_innerProduct(float *inData1, float *inData2, int num, bool bABS = true);


/**
* @brief	计算振幅差值
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
*
* @return	计算得到的该点上的振幅差值
*/
ISL_EXPORT float ISL_ampDiff(float *inSpec, int num);


/**
* @brief	计算指定能量比
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
*
* @return	计算得到的该点上的指定能量比
*/
ISL_EXPORT float ISL_indEnRatio(float *inSpec, int num);


/**
* @brief	计算衰减频宽,计算每个时刻瞬时谱的衰减频宽
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	间隔的点数，间隔越小，衰减得越快
*/
ISL_EXPORT float ISL_attenGradFre(float *inSpec, int num, float deltaFre);


/**
* @brief	计算吸收系数,振幅谱的高频端进行曲线拟合
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	inFreq		输入的各个点上的瞬时频率
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	返回拟合后的指数系数
*/
ISL_EXPORT float ISL_absorbCoffAttenExp(float *inSpec, float *inFreq, int num, float deltaFre);


/**
* @brief	计算吸收系数,利用衰减梯度法计算吸收系数，衰减峰值振幅
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	返回按照一定的原则找振幅谱高频端两个点做直线的斜率
*/
ISL_EXPORT float ISL_absorbCoffAttenPeak(float *inSpec, int num, float deltaFre);


/**
* @brief	计算吸收系数——高频,利用衰减梯度法计算吸收系数，衰减总能量,反应高频端的衰减属性
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	计算得到的该点上的吸收系数的值
*/
ISL_EXPORT float ISL_absorbCoffAttenEnHigh(float *inSpec, int num, float deltaFre);


/**
* @brief	计算吸收系数——低频,利用衰减梯度法计算吸收系数，衰减总能量,反映低频端的低频增强特征
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	inFreq		输入的各个点上的瞬时频率
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	计算得到的该点上的吸收系数的值
*/
ISL_EXPORT float ISL_absorbCoffAttenEnLow(float *inSpec, int num, float deltaFre);


/**
* @brief	计算吸收系数, 利用谱比法计算吸收系数
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
*
* @return	计算得到的该点上的吸收系数的值
*/
ISL_EXPORT float ISL_absorbCoffSpec(float *inSpec, int num);


/**
* @brief	计算瞬时峰值频率
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	瞬时主极值频率
*/
ISL_EXPORT float ISL_insPeakFre(float *inSpec, int num, float deltaFre);


/**
* @brief	计算瞬时峰值振幅
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
*
* @return	瞬时峰值振幅
*/
ISL_EXPORT float ISL_insPeakAmp(float *inSpec, int num);


/**
* @brief	计算瞬时中心频率
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	瞬时中心频率
*/
ISL_EXPORT float ISL_insMidFre(float *inSpec, int num, float deltaFre);


/**
* @brief	计算瞬时谱带宽
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	瞬时谱带宽
*/
ISL_EXPORT float ISL_insSpecBandWidth(float *inSpec, int num, float deltaFre);


/**
* @brief	计算瞬时均方根频率
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	瞬时均方根频率
*/
ISL_EXPORT float ISL_insRMSFre(float *inSpec, int num, float deltaFre);


/**
* @brief	计算振幅截距,这个属性应该与高频衰减梯度、峰值频率等一起解释  同样的衰减梯度有可能会有不同的振幅截距  最终表现的是主频的宽度或频带范围特征
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	inFreq		输入的各个点上的瞬时频率
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	振幅截距
*/
ISL_EXPORT float ISL_ampCept(float *inSpec, int num, float deltaFre);


/**
* @brief	计算加权平均频率
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	加权平均频率
*/
ISL_EXPORT float ISL_weightedMeanFre(float *inSpec, int num, float deltaFre);


/**
* @brief	计算有效带宽能量
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
*
* @return	有效带宽能量
*/
ISL_EXPORT float ISL_validBWEnergy(float *inSpec, int num);


/**
* @brief	计算频谱上的两个点,本函数用于计算能量衰减到65%和85%时的两个点
*
* @param[in]	inFreq		输入点上的瞬时频
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[out]	Fre1		最大瞬时谱对应的频率
* @param[out]	Fre2		inFreq数组中最后的频率
* @param[out]	yCorrd		两个点的纵坐标
* @param[out]	nCorrd		纵坐标的个数
*
* @return	无
*/
ISL_EXPORT void ISL_attenDot(int *inFreq, float *inSpec, int num, int &Fre1, int &Fre2, float *&yCorrd, int &nCorrd);



/**
* @brief	高频端的频谱曲线最小二乘拟合
*
* @param[in]	inFreq		输入点上的瞬时频
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	num			瞬时谱的点数
* @param[out]	Fre1		最大瞬时谱对应的频率
* @param[out]	Fre2		inFreq数组中最后的频率
* @param[out]	yCorrd		两个点的纵坐标
* @param[out]	nCorrd		纵坐标的个数
*
* @return	无
*/
ISL_EXPORT void ISL_attenExp(int *inFreq, float *inSpec, int num, int &Fre1, int &Fre2, float *&yCorrd, int &nCorrd);



/**
* @brief	瞬时属性计算,计算瞬时属性，包括瞬时相位、瞬时频率
*
* @param[in]	InData		输入点上的瞬时频
* @param[in]	num			地震数据的样点数
* @param[in]	dt			地震数据的采样间隔
* @param[out]	nWins		nWins>1,则对瞬时频率采用中值滤波法
* @param[out]	insAmp		输出的瞬时属性
* @param[out]	insPh		输出的瞬时属性
* @param[out]	insFre		输出的瞬时属性
*
* @return	无
*/
ISL_EXPORT void ISL_insPhAmpFre(float *inData, int num,
						float dt, int nWins,
						float *insAmp, float *insPh, float *insFre);


} /* End Of namespace ISLIB */


#endif /* PAI_FRAME_ISEISLIB_ABSORPTION_H_ */
