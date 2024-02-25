/**
*	@file	ISL_Absorption.cpp
*	@brief	[Header file of time-frequency analysis Functions], 吸收系数；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#include "ISL_Absorption.h"
#include "ISL_Fitting.h"
#include "ISL_Filter.h"
#include "ISL_Statistics.h"
//namespace ROCKLIB {

/* 计算相位 */
void ISL_calcPh(float pRe, float pIm, float& pPh)
{
    float pai = (float) PI; //3.1415926f;
    float pai5 = (float) (0.5 * pai);
    float tmpRe, tmpIm, tmp;

    //if (pRe*pRe+pIm*pIm)
    //	pPh=atan2(pIm, pRe);
    //else
    //	pPh=0.;

    tmpRe = fabs(pRe);
    tmpIm = fabs(pIm);

    if (tmpRe >= tmpIm && pRe != 0) {
        if (pRe > 0)
            tmp = 0;
        if (pRe < 0 && pIm >= 0)
            tmp = pai;
        if (pRe < 0 && pIm < 0)
            tmp = -pai;
        pPh = atan(pIm / pRe) + tmp;
     }  else {
        if (pIm != 0) {
            if (pIm > 0)
                tmp = pai5;
            if (pIm < 0)
                tmp = -pai5;
            pPh = tmp - atan(pRe / pIm);
         }  else {
            if (pRe >= 0)
                pPh = 0;
            if (pRe < 0)
                pPh = pai;
         }
     }
 }

/* 计算频率 */
void ISL_calcFreq(float ph1, float ph2 ,float dt, float& fn)
{
    float pai, pai2, ph, f1, f2, f3, f4, f5;
    pai = (float) PI; //3.1415926f;
    pai2 = 2 * pai;
    ph = ph2 - ph1;

    f1 = fabs(ph + pai);
    f2 = fabs(ph - pai);
    f3 = fabs(ph + pai2);
    f4 = fabs(ph - pai2);
    f5 = fabs(ph);

    if (f4 < f5)
        f5 = f4;
    if (f3 < f5)
        f5 = f3;
    if (f2 < f5)
        f5 = f2;
    if (f1 < f5)
        f5 = f1;

    fn = f5 / (pai2 * dt);
 }




/*计算两离散序列的内积*/
float ISL_innerProduct(float *inData1, float *inData2, int num, bool bABS)
{
    float sum = 0.0;
    for (int i = 0; i < num; i++) {
        sum = sum + inData1[i] * inData2[i];
     }
    if(bABS == true)
        sum = fabs(sum);

    return sum;
 }


/*计算振幅差值*/
float ISL_ampDiff(float *inSpec, int num)
{
    float MaxValue;
    MaxValue = inSpec[0];
    /////计算主峰值/////
    for (int i = 1; i < num; i++)
        if (MaxValue < inSpec[i])
            MaxValue = inSpec[i];

    if (MaxValue == 0)
        return 0.;

    /////计算平均振幅/////
    float MeanVal = 0.;
    for (int i = 0; i < num; i++)
        MeanVal += inSpec[i];
    MeanVal /= num;

    return MaxValue - MeanVal;
 }


/* 计算指定能量比 */
float ISL_indEnRatio(float *inSpec,int num)
{
    int LocMax, LocRight;
    float MaxValue, TmpValue;

    LocRight = 0;
    /////判断最大值的位置/////
    MaxValue = inSpec[0];
    LocMax = 0;
    for (int i = 1; i < num; i++) {
        if (MaxValue < inSpec[i]) {
            MaxValue = inSpec[i];
            LocMax = i;
         }
     }

    if (MaxValue == 0)
        return 0.;

    /////计算总能量/////
    float TotalEn = 0.;
    for (int i = 0; i < num; i++)
        TotalEn += inSpec[i] * inSpec[i];

    if (TotalEn <= 0)
        return 0.;

    /////计算阀值点/////
    TmpValue = float(0.35 * MaxValue);
    if (TmpValue < inSpec[num - 1])
        LocRight = num - 1;
    else {
        for (int i = LocMax; i < num - 1; i++) {
            if (inSpec[i] >= TmpValue && inSpec[i + 1] < TmpValue) {
                LocRight = i;
                break;
             }
         }
     }

    if (LocRight == 0)
        return 0.;

    /////计算指定低频端能量/////
    float LowFreEn = 0.;
    for (int i = 0; i < LocRight; i++)
        LowFreEn += inSpec[i] * inSpec[i];

    return LowFreEn / TotalEn;
 }

/*计算衰减频宽,计算每个时刻瞬时谱的衰减频宽*/
float ISL_attenGradFre(float *inSpec, int num, float deltaFre)
{
    /////衰减总能量/////
    int LocLeft = 0, LocRight = 0;
    float TmpEnergy, TotalEnergy;
    float EnergyLeft, EnergyRight, valLeft, valRight;

    /////计算总能量/////
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    /////寻找衰减到65%时的频率/////
    EnergyLeft = 0.;
    TmpEnergy = float(0.65 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyLeft += inSpec[i] * inSpec[i];
        if (EnergyLeft > TmpEnergy) {
            LocLeft = i;
            valLeft = inSpec[i];
            break;
         }
     }

    /////寻找衰减到85%时的频率/////
    EnergyRight = 0.;
    TmpEnergy = float(0.85 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyRight += inSpec[i] * inSpec[i];
        if (EnergyRight > TmpEnergy) {
            LocRight = i;
            valRight = inSpec[i];
            break;
         }
     }

    return float(LocRight - LocLeft * 1.0) * deltaFre;
 }

/*计算吸收系数,振幅谱的高频端进行曲线拟合*/
float ISL_absorbCoffAttenExp(float *inSpec, float *inFreq, int num, float deltaFre)
{
    /////衰减最大振幅/////
    int MaxLoc;
    float MaxSpec;
    float eps = deltaFre / 2;

    /////判断最大值的位置/////
    MaxSpec = inSpec[0];
    MaxLoc = 0;
    for (int i = 0; i < num; i++) {
        if (MaxSpec < inSpec[i]) {
            MaxSpec = inSpec[i];
            MaxLoc = i;
         }
     }

    if (MaxSpec == 0)
        return 0.;

    int RightLoc = num - 1;
    int nterm = 3; //拟合项数=拟合次数+1
    int nNums = RightLoc - MaxLoc + 1;
    float *SpecX = new float[nNums];
    float *SpecY = new float[nNums];
    float *simu = new float[nNums];
    float *Coef = new float[nterm];
    for (int i = 0; i < nNums; i++) {
        SpecX[i] = inFreq[i + MaxLoc];
        if (SpecX[i] == 0)
            SpecX[i] = eps;
        SpecY[i] = inSpec[i + MaxLoc];
        if (SpecY[i] == 0)
            SpecY[i] = SpecY[i - 1];
        else
            SpecY[i] = log(SpecY[i] / pow(SpecX[i], 2));
     }
    /////曲线拟合—二次拟合/////
    ISL_lscf(SpecX, SpecY, nNums, Coef, nterm, simu);
    float Slope = Coef[nterm - 1]; //返回二次项的系数
    delete[] SpecX;
    delete[] SpecY;
    delete[] simu;
    delete[] Coef;

    return Slope;
 }

/*计算吸收系数,利用衰减梯度法计算吸收系数，衰减峰值振幅*/
float ISL_absorbCoffAttenPeak(float *inSpec, int num, float deltaFre)
{
    /////衰减最大振幅/////
    int LocMax, LocLeft = 0, LocRight = 0;
    float MaxValue, TmpValue;
    float valLeft = 0., valRight = 0.;

    /////判断最大值的位置/////
    MaxValue = inSpec[0];
    LocMax = 0;
    for (int i = 0; i < num; i++) {
        if (MaxValue < inSpec[i]) {
            MaxValue = inSpec[i];
            LocMax = i;
         }
     }

    if (MaxValue == 0)
        return 0.;

    /////寻找衰减到65%时的频率/////
    TmpValue = float(0.65 * MaxValue);
    if (inSpec[num - 1] > TmpValue) {
        return (inSpec[num - 1] - MaxValue) / (deltaFre * (num - 1 - LocMax));
     }  else {
        for (int i = LocMax; i < num - 1; i++) {
            if (inSpec[i] >= TmpValue && inSpec[i + 1] < TmpValue) {
                LocLeft = i;
                valLeft = inSpec[i];
                break;
             }
         }
     }

    if (LocLeft < LocMax) {
        LocLeft = LocMax;
        valLeft = MaxValue;
     }

    /////寻找衰减到35%时的频率/////
    TmpValue = float(0.35 * MaxValue);
    if (inSpec[num - 1] > TmpValue) {
        LocRight = num - 1;
        valRight = inSpec[num - 1];
     }  else {
        for (int i = LocMax; i < num - 1; i++) {
            if (inSpec[i] >= TmpValue && inSpec[i + 1] < TmpValue) {
                LocRight = i;
                valRight = inSpec[i];
                break;
             }
         }
     }

    if (LocRight == LocLeft)
        return 0.;
    float slope = (valRight - valLeft) / (deltaFre * (LocRight - LocLeft));
    if (slope >= 0)
        return 0.;
    else
        return slope;
 }


/* 计算吸收系数——高频,利用衰减梯度法计算吸收系数，衰减总能量 */
float ISL_absorbCoffAttenEnHigh(float *inSpec, int num, float deltaFre)
{
    /////衰减总能量/////
    int LocLeft = 0, LocRight = 0;
    float TmpEnergy, TotalEnergy;
    float EnergyLeft, EnergyRight, valLeft = 0., valRight = 0.;

    /////计算总能量/////
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    /////寻找衰减到65%时的频率/////
    EnergyLeft = 0.;
    TmpEnergy = float(0.65 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyLeft += inSpec[i] * inSpec[i];
        if (EnergyLeft > TmpEnergy) {
            LocLeft = i;
            valLeft = inSpec[i];
            break;
         }
     }

    /////寻找衰减到85%时的频率/////
    EnergyRight = 0.;
    TmpEnergy = float(0.85 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyRight += inSpec[i] * inSpec[i];
        if (EnergyRight > TmpEnergy) {
            LocRight = i;
            valRight = inSpec[i];
            break;
         }
     }

    if (LocRight - LocLeft == 0)
        return 0.;
    else {
        float slope = (inSpec[LocRight] - inSpec[LocLeft]) / (deltaFre * (LocRight - LocLeft));
        return slope;
     }
 }


/*计算吸收系数——低频,利用衰减梯度法计算吸收系数，衰减总能量*/
float ISL_absorbCoffAttenEnLow(float *inSpec, int num, float deltaFre)
{
    /////衰减总能量/////
    int LocLeft = 0, LocRight = 0;
    float TmpEnergy, TotalEnergy;
    float EnergyLeft, EnergyRight;

    /////计算总能量/////
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    /////寻找衰减到15%时的频率/////
    EnergyLeft = 0.;
    TmpEnergy = float(0.15 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyLeft += inSpec[i] * inSpec[i];
        if (EnergyLeft > TmpEnergy) {
            LocLeft = i;
            break;
         }
     }

    /////寻找衰减到35%时的频率/////
    EnergyRight = 0.;
    TmpEnergy = float(0.35 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyRight += inSpec[i] * inSpec[i];
        if (EnergyRight > TmpEnergy) {
            LocRight = i;
            break;
         }
     }

    if (LocRight - LocLeft == 0)
        return 0.;
    else {
        float slope = (inSpec[LocRight] - inSpec[LocLeft]) / (deltaFre * (LocRight - LocLeft));
        return slope;
     }
 }


/*计算吸收系数, 利用谱比法计算吸收系数*/
float ISL_absorbCoffSpec(float *inSpec, int num)
{
    int LocMax, LocLeft = 0, LocRight = 0;
    float Max, Tmp;
    float AreaLeft, AreaRight;

    /////判断最大值的位置/////
    Max = inSpec[0];
    LocMax = 0;
    for (int i = 0; i < num; i++) {
        if (Max < inSpec[i]) {
            Max = inSpec[i];
            LocMax = i;
         }
     }

    if (Max == 0)
        return 0.;

    Tmp = float(0.2 * Max);
    /////判断左边临界值的位置/////
    if (inSpec[0] >= Tmp)
        LocLeft = 0;
    else {
        for (int i = 0; i < num; i++) {
            if (inSpec[i] <= Tmp && inSpec[i + 1] > Tmp) {
                LocLeft = i;
                break;
             }
         }
     }

    /////判断右边临界值的位置/////
    if (inSpec[num - 1] > Tmp)
        LocRight = num - 1;
    else {
        for (int i = num - 1; i >= 0; i--) {
            if (inSpec[i] <= Tmp && inSpec[i - 1] > Tmp) {
                LocRight = i;
                break;
             }
         }
     }

    AreaLeft = 0.;
    for (int i = LocLeft; i < LocMax; i++)
        AreaLeft += inSpec[i] * inSpec[i];

    AreaRight = 0.;
    for (int i = LocMax; i < LocRight; i++)
        AreaRight += inSpec[i] * inSpec[i];

    if (AreaRight < 0.01)
        return 0.;
    else
        return AreaLeft / AreaRight;
 }

/*计算瞬时峰值频率*/
float ISL_insPeakFre(float *inSpec, int num, float deltaFre)
{
    int MaxLoc;
    float MaxSpec;
    MaxLoc = 0;
    MaxSpec = inSpec[0];
    for (int i = 1; i < num; i++) {
        if (MaxSpec < inSpec[i]) {
            MaxLoc = i;
            MaxSpec = inSpec[i];
         }
     }

    if (MaxSpec == 0)
        return 0.;

    return float(MaxLoc * deltaFre);
 }

/*计算瞬时峰值振幅*/
float ISL_insPeakAmp(float *inSpec, int num)
{
    int MaxLoc;
    float MaxSpec;

    MaxLoc = 0;
    MaxSpec = inSpec[0];
    for (int i = 1; i < num; i++) {
        if (MaxSpec < inSpec[i]) {
            MaxLoc = i;
            MaxSpec = inSpec[i];
         }
     }

    if (MaxSpec == 0)
        return 0.;

    return MaxSpec;
 }

/*计算瞬时中心频率*/
float ISL_insMidFre(float *inSpec, int num, float deltaFre)
{
    float TmpEnergy, TotalEnergy;
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    TmpEnergy = 0.;
    for (int i = 0; i < num; i++)
        TmpEnergy += deltaFre * i * inSpec[i] * inSpec[i];

    return float(TmpEnergy / TotalEnergy);
 }

/*计算瞬时谱带宽*/
float ISL_insSpecBandWidth(float *inSpec, int num, float deltaFre)
{
    float TmpBW, TmpEnergy, TotalEnergy;
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    TmpEnergy = 0.;
    for (int i = 0; i < num; i++)
        TmpEnergy += deltaFre * i * inSpec[i] * inSpec[i];

    TmpBW = float(TmpEnergy / TotalEnergy);
    TmpEnergy = 0.;
    for (int i = 0; i < num; i++)
        TmpEnergy += (deltaFre * i - TmpBW) * (deltaFre * i - TmpBW) * inSpec[i] * inSpec[i];

    return sqrt(float(TmpEnergy / TotalEnergy));
 }


/*计算瞬时均方根频率*/
float ISL_insRMSFre(float *inSpec, int num, float deltaFre)
{
    float TmpEnergy = 0, TotalEnergy = 0;
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    TmpEnergy = 0.;
    for (int i = 0; i < num; i++)
        TmpEnergy += deltaFre * i * deltaFre * i * inSpec[i] * inSpec[i];

    return sqrt(float(TmpEnergy / TotalEnergy));
 }


/*计算振幅截距*/
float ISL_ampCept(float *inSpec, int num, float deltaFre)
{
    /////衰减总能量/////
    int LocLeft = 0, LocRight = 0;
    float TmpEnergy, TotalEnergy;
    float EnergyLeft, EnergyRight, valLeft = 0., valRight = 0.;

    /////计算总能量/////
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i] * inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    /////寻找衰减到65%时的频率/////
    EnergyLeft = 0.;
    TmpEnergy = float(0.65 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyLeft += inSpec[i] * inSpec[i];
        if (EnergyLeft > TmpEnergy) {
            LocLeft = i;
            valLeft = inSpec[i];
            break;
         }
     }

    /////寻找衰减到85%时的频率/////
    EnergyRight = 0.;
    TmpEnergy = float(0.85 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyRight += inSpec[i] * inSpec[i];
        if (EnergyRight > TmpEnergy) {
            LocRight = i;
            valRight = inSpec[i];
            break;
         }
     }
    if (LocRight == LocLeft)
        return 0.;
    float slope = (inSpec[LocRight] - inSpec[LocLeft]) / (deltaFre * (LocRight - LocLeft));
    float cept = inSpec[LocLeft] - slope * LocLeft * deltaFre;
    return cept;
 }


/*计算加权平均频率*/
float ISL_weightedMeanFre(float *inSpec, int num, float deltaFre)
{
    float TmpEnergy, TotalEnergy = 0;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i];

    if (TotalEnergy == 0)
        return 0.;

    TmpEnergy = 0.;
    for (int i = 0; i < num; i++)
        TmpEnergy += deltaFre * i * inSpec[i];

    return float(TmpEnergy / TotalEnergy);
 }


/*计算有效带宽能量*/
float ISL_validBWEnergy(float *inSpec, int num)
{
    float TotalEnergy = 0;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i];

    return TotalEnergy;
 }


/**********************************************************************************************************************
 *
 *  功能：计算频谱上的两个点
 *
 *  说明：本函数用于计算能量衰减到65%和85%时的两个点
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float*				inFreq				In			输入点上的瞬时频率
 *		float*				inSpec				In			输入点上的瞬时谱
 *		int					num					In			瞬时谱的点数
 *		int					Fre1				Out			最大瞬时谱对应的频率
 *		int					Fre2				Out			inFreq数组中最后的频率
 *		float*				yCorrd				Out			两个点的纵坐标
 *		int					nCorrd				Out			纵坐标的个数
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_attenDot(int *inFreq, float *inSpec, int num, int &Fre1, int &Fre2, float *&yCorrd, int &nCorrd)
{
    float x1, x2, y1, y2;
    float TmpEnergy, TotalEnergy;
    float EnergyLeft, EnergyRight;
    x1 = inFreq[0];
    x2 = inFreq[num - 1];
    y1 = inSpec[0];
    y2 = inSpec[num - 1];

    /////判断最大值的位置/////
    float MaxSpec = inSpec[0];
    int MaxLoc = 0;
    for (int i = 0; i < num; i++) {
        if (MaxSpec < inSpec[i]) {
            MaxSpec = inSpec[i];
            MaxLoc = i;
         }
     }
    Fre1 = inFreq[MaxLoc];
    Fre2 = inFreq[num - 1];

    /////计算总能量/////
    TotalEnergy = 0.;
    for (int i = 0; i < num; i++)
        TotalEnergy += inSpec[i];

    /////寻找衰减到65%时的频率/////
    EnergyLeft = 0.;
    TmpEnergy = float(0.65 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyLeft += inSpec[i];
        if (EnergyLeft > TmpEnergy) {
            x1 = inFreq[i];
            y1 = inSpec[i];
            break;
         }
     }

    /////寻找衰减到85%时的频率/////
    EnergyRight = 0.;
    TmpEnergy = float(0.85 * TotalEnergy);
    for (int i = 0; i < num; i++) {
        EnergyRight += inSpec[i];
        if (EnergyRight > TmpEnergy) {
            x2 = inFreq[i];
            y2 = inSpec[i];
            break;
         }
     }

    if (x1 < MaxLoc) {
        x1 = inFreq[MaxLoc];
        y1 = inSpec[MaxLoc];
     }

    float Tmp;
    if (y1 < y2) {
        Tmp = y1;
        y1 = y2;
        y2 = Tmp;
     }

    nCorrd = num - MaxLoc;
    if(yCorrd){  delete []yCorrd; yCorrd = NULL;  }
    yCorrd = new float[nCorrd];

    float CoeffA = float((y2 - y1) / (x2 - x1));
    float CoeffB = y1 - CoeffA * x1;
    for (int i = 0; i < nCorrd; i++) {
        yCorrd[i] = CoeffA * inFreq[i + MaxLoc] + CoeffB;
        if (yCorrd[i] < 0) {
            Fre2 = inFreq[i + MaxLoc];
            break;
         }
     }
 }


/*高频端的频谱曲线最小二乘拟合*/
void ISL_attenExp(int *inFreq, float *inSpec, int num, int &Fre1, int &Fre2, float *&yCorrd, int &nCorrd)
{
    int MaxLoc;
    float MaxSpec;
    float eps = 0.5f;

    /////判断最大值的位置/////
    MaxSpec = inSpec[0];
    MaxLoc = 0;
    for (int i = 0; i < num; i++) {
        if (MaxSpec < inSpec[i]) {
            MaxSpec = inSpec[i];
            MaxLoc = i;
         }
     }
    Fre1 = inFreq[MaxLoc];

    int RightLoc = num - 1;
    float RightValue = float(MaxSpec * 0.01);
    if (RightValue < inSpec[num - 1])
        RightLoc = num - 1;
    else {
        for (int i = num - 1; i > MaxLoc; i--) {
            if (inSpec[i] < RightValue && inSpec[i - 1] > RightValue) {
                RightLoc = i;
                break;
             }
         }
     }
    Fre2 = inFreq[RightLoc];

    int nterm = 3;
    nCorrd = RightLoc - MaxLoc + 1;

    if(yCorrd){  delete []yCorrd; yCorrd = NULL;  }
    yCorrd = new float[nCorrd];
    float *xCoord = new float[nCorrd];
    float *simu = new float[nCorrd];
    float *Coef = new float[nterm];
    for (int i = 0; i < nCorrd; i++) {
        xCoord[i] = inFreq[i + MaxLoc]; //取对应的频率值
        if (xCoord[i] == 0)
            xCoord[i] = eps;
        yCorrd[i] = inSpec[i + MaxLoc]; //取对应的能量值
        if (yCorrd[i] == 0)
            yCorrd[i] = yCorrd[i - 1];
        else
            yCorrd[i] = log(yCorrd[i] / pow(xCoord[i], 2)); //k=2,N=1
     }

    ISL_lscf(xCoord, yCorrd, nCorrd, Coef, nterm, simu);
    for (int i = 0; i < nCorrd; i++)
        yCorrd[i] = xCoord[i] * xCoord[i] * exp(simu[i]);

    delete[] xCoord;
    delete[] simu;
    delete[] Coef;
 }



/*瞬时属性计算,计算瞬时属性，包括瞬时相位、瞬时频率*/
void ISL_insPhAmpFre(float *inData, int num,
                        float dt, int nWins,
                        float *insAmp, float *insPh, float *insFre)
{
    float * comIm = new float[num];
    float * nInstFre = new float[num];
    memset(comIm, 0, sizeof(float) * num);
    memset(nInstFre, 0, sizeof(float) * num);
    ISL_hilbert(num, inData, comIm);

    float phase = 0;
    for (int i = 0; i < num; i++) {
        insAmp[i] = sqrt(inData[i] * inData[i] + comIm[i] * comIm[i]);
        ISL_calcPh(inData[i], comIm[i], insPh[i]);
        ISL_calcFreq(phase, insPh[i], dt, nInstFre[i]);
        phase = insPh[i];
     }
    nInstFre[0] = nInstFre[1];
    if (nWins <= 1)
        memcpy(insFre, nInstFre, sizeof(float) * num);
    else
        ISL_windowMidValue(nInstFre, num, nWins, insFre);

    delete[] comIm;
    delete[] nInstFre;
 }
// } /* End Of namespace ISLIB */
