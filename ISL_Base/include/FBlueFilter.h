#ifndef _FBLUEFILTER_H
#define _FBLUEFILTER_H
#include <math.h>
#include <stdlib.h>

/*
* @brief    	蓝色滤波的参数输入
*
* @param[in]		nw		时窗个数; 
* @param[in]            llf           	滤波窗口长度;
* @param[in]            time_point	时窗值;
* @param[in]            coe		加强算子系数;
* @param[in]            f_wide		加强频率时窗宽度;
* @param[in]            f_expand	加强频率;
* @param[in]	        f_high_cut	高通频率;
* @param[in]		slope_h		高通频率斜坡;	
* @param[in]		f_low_cut	低通频率;
* @param[in]		slope_l		低通频率斜坡;
*/

extern "C" void bluefilter_initialize_(int *nw
		, int *llf
		, float * time_point
		, float * coe
		, float * f_wide
		, float * f_expand
		, float * f_high_cut
		, float * slope_h
		, float * f_low_cut
		, float * slope_l);

/*
 * @brief        对地震道进行蓝色滤波
 *
 * @param[in/out]            trace		道数据;
 * @param[in]            trace_lt	道数据长度;
 * @param[in]            opt
 *
*/

extern "C" void bluefilter_process_(float *trace, int * trace_lt,float * opt);
extern "C" void bluefilter_finish_();


#endif

