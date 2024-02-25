/**
*	@file	ISL_UserDefine.h 
*	@brief	[Header file of ISeisLib User-defined macros], 算法库用户定义头文件；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-10
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_USERDEFINE_H
#define PAI_FRAME_ISEISLIB_USERDEFINE_H

#include "ISL_Export.h"
#include "ISL_ReturnCode.h"


namespace ISLib {

	#ifndef NFAX
	#define NFAX 10
	#endif

	#ifndef NTAB
	#define NTAB 240
	#endif

	#ifndef NullValue
	#define NullValue -9999.99
	#endif

	#ifndef MinDataError
	#define MinDataError 0.0001
	#endif
	
	#ifndef MinDistanceError
	#define MinDistanceError 5.0
	#endif

	#ifndef MaxDistanceError
	#define MaxDistanceError 20.0
	#endif

	#ifndef MaxPermValue
	#define MaxPermValue 1e4
	#endif

	#ifndef MinDeltaError
	#define MinDeltaError 3
	#endif

	#ifndef INF
	#define INF 1E200
	#endif

	#ifndef EP
	#define EP 1E-10
	#endif

	#ifndef MAXV
	#define MAXV 300
	#endif

	#ifndef INITIAL_T
	#define INITIAL_T 999999
	#endif

	/**	@brief look factor			  */
	#ifndef LOOKFAC
	#define LOOKFAC 2
	#endif

	/**	@brief Largest allowed nfft		  */
	#ifndef PFA_MAX
	#define PFA_MAX 720720
	#endif

	/* default pnoise value		*/
	#ifndef PNOISE
	#define PNOISE 0.001
	#endif

	#ifndef	F1
	#define	F1	1.125
	#endif

	#ifndef	F2
	#define	F2	-0.04166667
	#endif

	/**	@brief	[PI] : 圆周率 */
	#ifndef PI
	#define PI 3.141592653589793
	#endif

	/**	@brief	[ABS(x)] : 求绝对值	*/
	#ifndef ABS
	#define ABS(x) ((x) < 0 ? -(x) : (x))
	#endif
	
	/**	@brief	[MAX(x,y)] : 两数最大值	*/
	#ifndef MAX
	#define	MAX(x,y) ((x) > (y) ? (x) : (y))
	#endif

	/**	@brief	[MIN(x,y)] : 两数最小值	*/
	#ifndef MIN
	#define	MIN(x,y) ((x) < (y) ? (x) : (y))
	#endif

	#ifndef NINT
	#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
	#endif

	/**	@brief	[ISODD(n)] : 奇数判断 */
	#define ISODD(n) ((n) & 01)
	
	/**	@brief	[ISIZE] : 整型数据空间字节数 */
	#define ISIZE sizeof(int)
	
	/**	@brief	[FSIZE] : 浮点型数据空间字节数 */
	#define FSIZE sizeof(float)
	
	/**	@brief	[DSIZE] : 双精度型数据空间字节数	*/
	#define DSIZE sizeof(double)

	/**	@brief	字节的定义	*/
	typedef unsigned char ISL_byte;
	typedef short ISL_Int16;
	typedef int ISL_Int32;
	typedef long long ISL_Int64;
	
	/**	@brief	颜色的定义 */
	struct ISL_color
	{
		short r;
		short g;
		short b;
		short a;
	};
	
	#define P120 0.120536680
	#define P142 0.142314838
	#define P173 0.173648178
	#define P222 0.222520934
	#define P239 0.239315664
	#define P281 0.281732557
	#define P342 0.342020143
	#define P354 0.354604887
	#define P382 0.382683432
	#define P415 0.415415013
	#define P433 0.433883739
	#define P464 0.464723172
	#define P540 0.540640817
	#define P559 0.559016994
	#define P568 0.568064747
	#define P587 0.587785252
	#define P623 0.623489802
	#define P642 0.642787610
	#define P654 0.654860734
	#define P663 0.663122658
	#define P707 0.707106781
	#define P748 0.748510748
	#define P755 0.755749574
	#define P766 0.766044443
	#define P781 0.781831482
	#define P822 0.822983866
	#define P841 0.841253533
	#define P866 0.866025404
	#define P885 0.885456026
	#define P900 0.900968868
	#define P909 0.909631995
	#define P923 0.923879533
	#define P935 0.935016243
	#define P939 0.939692621
	#define P951 0.951056516
	#define P959 0.959492974
	#define P970 0.970941817
	#define P974 0.974927912
	#define P984 0.984807753
	#define P989 0.989821442
	#define P992 0.992708874

}/*End of ISLib*/

#endif
