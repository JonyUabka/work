// 引用此 文件  #include "prdecon.h"

#ifndef _PRDECON_H
#define _PRDECON_H
#include <math.h>
#include <stdlib.h>

/**
* @brief 完成预测反褶积功能
* @brief 完成预测反褶积模块的参数分析，工作缓冲区长度计算等功能。
*
* @param[in]  address1   实例地址
* @param[in]  tmtable    参数表
* @param[in]  params     界面参数
* @param[in]  dt  	 采样间隔
* @param[in]  nt   	 道长
* @param[out] lbuf1   	 工作缓冲区长度
*
* SUBROUTINE PRDECON_AM (address1,params,IPAR, DT, NT, LC )
*/
extern "C" void prdecon_am_(long* address1, int* tmtable, int* params, int* dt, int* nt, int* lbuf1);



/**
* @brief 完成预测反褶积模块的核心算法。
*
* @param[in] 	  buf1	  	工作缓冲区
* @param[in] 	  libm 	  	切除数据表
* @param[in/out]  ih  	  	道头
* @param[in/out]  trace_data   	地震道
*
* SUBROUTINE PRDECON_PM (KBUF, LIBM, IH, TRACE )
*/
extern "C" void prdecon_pm_(int* buf1, int* libm, int* ih, float* trace_data);



#endif
