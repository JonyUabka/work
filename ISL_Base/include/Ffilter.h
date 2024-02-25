
// 引用此 文件  #include "Ffilter.h"

#ifndef _FFILTER_H
#define _FFILTER_H
#include <math.h>
#include <stdlib.h>


/*
*+------------------------------------------------------------
*+    程序功能：计算滤波因子
*+    参数：
*+        LMS--     滤波算子长度（ms）
*+        ISI--     采样间隔（ms）
*+        JOB--     滤波因子类型
*+                  =1  低通滤波
*+                  =2  高通滤波
*+                  =3  陷频滤波
*+                  =4  带通滤波
*+        RIT1--    带通滤波参数1
*+        RIT2--    带通滤波参数2
*+        RIT3--    带通滤波参数3
*+        RIT4--    带通滤波参数4
*+        FF--      存放滤波算子
*+------------------------------------------------------------
      SUBROUTINE CALFIL(LMS,ISI,JOB,RIT1,RIT2,RIT3,RIT4,FF)
      REAL  FF(*)
*/

/**
 * @brief      计算滤波算子 
 *
 * @param[in]   length          滤波算子长度 (毫秒) 输出数据滤波算子filt的个数由length和dt计算得出；
 * @param[in]   dt              采样间隔  (毫秒)；
 * @param[in]   icode           滤波因子类型 
 * 				=1  低通滤波 
 * 				=2  高通滤波 
 * 				=3  陷频滤波 ；
 * 				=4  带通滤波
 * @param[in]   f1              滤波档; 
 * @param[in]   f2              滤波档;
 * @param[in]   f3              滤波档;
 * @param[in]   f4              滤波档;
 * @param[out]   filt		输出数据，滤波算子;
 *
 */
extern "C" void calfil_(int* length, int* dt, int* icode, float* f1, float* f2, float* f3, float* f4, float* filt);

/*
*+------------------------------------------------------------
*+    程序功能：对地震道进行滤波
*+    参数：
*+        filt --  滤波因子
*+        nt   --  滤波因子长度
*+        LSAMP--  地震道长
*+        trace--  地震道
*+------------------------------------------------------------
      subroutine filterexec( FILT,NT,LSAMP, TRACE)
      real TRACE(*), FILT(*), WKPADI(LSAMP)
*/

/**
 * @brief      对地震道进行滤波 
 *
 * @param[in]   	filt         	 滤波算子；
 * @param[in]   	lfilt         	 滤波算子长度；
 * @param[in]   	numOfSamples     滤波数据的长度; 
 * @param[in/out]	trace_data	 需要计算的滤波数据，返回时改变为滤波结果;
 *
 */
extern "C" void filter_exec_(float* filt, int* lfilt, int* numOfSamples, float* trace_data);



/*
*
* *+------------------------------------------------------------
* *+    程序功能：相关或褶积运算
* *+    参数：
* *+        Y--      输出数组
* *+        LY--     输出数组长度
* *+        X--      输入数组
* *+        LX--     输入数组长度
* *+        F--      滤波算子
* *+        LF--     滤波算子长度
* *+        IFLAG--  相关或褶积
* *+                 >=0，相关
* *+                 其它，褶积
* *+        ISHIFT-- 時移量
* *+------------------------------------------------------------
	SUBROUTINE CONVOX(Y,LY,X,LX,F,LF,IFLAG,ISHIFT)
	REAL  Y(*),X(*),F(*)
*/

/**
 * @brief     相关或褶积运算 
 *
 * @param[out]	y		输出数组; 
 * @param[in]   ly		输出数组长度;  
 * @param[in]   x		输入数组;
 * @param[in]   lx              输入数组长度; 
 * @param[in]   f		滤波算子;
 * @param[in]   lf		滤波算子长度;
 * @param[in]   iflag		相关或褶积	>=0，相关;其它，褶积; 
 * @param[in]   ishift		時移量;
 *
 */
extern "C" void convox_(float* y, int *ly, float* x, int *lx, float* f, int *lf, int *iflag, int *ishift);



#endif
