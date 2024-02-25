/****************************************************************************
**
**				ISeis Lib [Base Algorithm] - 图像处理
**
**				Writer By	Mouri Song
**							Baihong Liu
**							Qiang Yang
**
**				(Soft Center IGP)
**
**				DATA : 2014-08-20
**
****************************************************************************/

/**
*	@file	ISL_ImageProcessing
*	@brief	[Header file of ImageProcessing Functions], 图像处理；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-08-20
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_IMAGEPROCESSING_H_
#define PAI_FRAME_ISEISLIB_IMAGEPROCESSING_H_

#include "ISL_UserDefine.h"

namespace ISLib
{


// BMP 头
/**
 * @brief	BMP图像格式文件头
 */
typedef struct
{
	char id[2];						/** bmp 文件标志 "BM" */
	unsigned long fileSize;			/** 文件大小 */
	unsigned long reserved0;
	unsigned long bitmapDataOffset;
	unsigned long bitmapHeaderSize;
	unsigned long width;			/** 图像宽度 */
	unsigned long height;			/** 图像高度 */
	unsigned short planes;
	unsigned short bitsPerPixel;	/** 每个像素站多少位 */
	unsigned long compression;		/** 是否压缩 */
	unsigned long bitmapDataSize;
	unsigned long hRes;
	unsigned long vRes;
	unsigned long colors;
	unsigned long importantColors;
	unsigned char palette[256][4];	/** 调色板数据,24位及以上像素没有该数据 */
} BMPHeader_t;


/**
 * @brief	读取BMP格式图片
 *
 * @brief	注意bmp文件是由左到右, 由下到上的循序存储的像素.每个像素颜色顺序为B G R (A)
 * 			BMP的行象素是4字节对齐的，不足的补0，比如有个图像每行宽度是63象素，
 * 			BMP文件存储时会每行存储64个字节，最后一个字节用0补齐
 *
 * 			pic为返回值，使用完毕后，请用 free()来释放空间
 * 			pic是一维数组，每四个数据存储一个像素点，例如：
 * 			pic[0]、 pic[1]、 pic[2]、 pic[3] 分别存储的是第1个像素点的 R G B （A） ，以后的数据依次类推
 *
 * @param[in]	name	输入文件的文件名
 * @param[in]	pic		返回 R G B(A) 类型像素数据, 数组的大小为 width*height*4
 * @param[in]	width	返回图像修正后的宽度
 * @param[in]	height	返回图像修正后的高度
 *
 * @return	true	读取成功
 * 			false	读取失败
 */
ISL_EXPORT bool ISL_loadBMP(const char *name, ISL_byte *& pic, int &width, int &height);



/**
 * @brief	得到灰度值(单个像素)
 *
 * @param[in]	r	输入红色值
 * @param[in]	g	输入绿色值
 * @param[in]	b	输入蓝色值
 *
 * @return	灰度值
 */
ISL_EXPORT ISL_byte ISL_IMG_grayPixel(ISL_byte r, ISL_byte g, ISL_byte b);



/**
 * @brief	得到黑白值, mid = 黑白判断阈值（0 ~ 255）(单个像素)
 *
 * @param[in]	in	输入的原始颜色值
 * @param[in]	mid	黑白判断阈值（0 ~ 255）
 *
 * @return	得到黑白值
 */
ISL_EXPORT ISL_byte ISL_IMG_blackPixel(ISL_byte in, int mid = 127);



/**
 * @brief	亮度设置(单个像素)
 *
 * @param[in]	in	输入的原始颜色值
 * @param[in]	light	输入的亮度值（0~255）
 *
 * @return	修正后的值
 */
ISL_EXPORT ISL_byte ISL_IMG_lightPixel(ISL_byte in, int light);




/**
 * @brief	亮度反转(底片色)(单个像素)
 *
 * @param[in]	in	输入的原始颜色值
 *
 * @return	修正后的值
 */
ISL_EXPORT ISL_byte ISL_IMG_lightReversePixel(ISL_byte in);



/**
 * @brief	图像曝光处理(单个像素)
 *
 * @param[in]	in	输入的原始颜色值
 *
 * @return	修正后的值
 */
ISL_EXPORT ISL_byte ISL_IMG_exposalPixel(ISL_byte in);



/**
 * @brief	用户值域范围调整(单个像素)
 *
 * @param[in]	ri	输入三原色，红色分量
 * @param[in]	gi	输入三原色，绿色分量
 * @param[in]	bi	输入三原色，蓝色分量
 * @param[in]	incement	调整值
 * @param[in]	ro	输出三原色，红色分量
 * @param[in]	go	输出三原色，绿色分量
 * @param[in]	bo	输出三原色，蓝色分量
 *
 * @return	无
 */
ISL_EXPORT void ISL_IMG_contrastAlterPixel(	ISL_byte ri,
							ISL_byte gi,
							ISL_byte bi,
							int incement,
							ISL_byte &ro,
							ISL_byte &go,
							ISL_byte &bo );



/**
 * @brief	对图像使用阈值法进行着色处理(单个像素)
 *
 * @param[in]	ri	输入三原色，红色分量
 * @param[in]	gi	输入三原色，绿色分量
 * @param[in]	bi	输入三原色，蓝色分量
 * @param[in]	ru	用户定义三原色，红色分量
 * @param[in]	gu	用户定义三原色，绿色分量
 * @param[in]	bu	用户定义三原色，蓝色分量
 * @param[in]	ro	输出三原色，红色分量
 * @param[in]	go	输出三原色，绿色分量
 * @param[in]	bo	输出三原色，蓝色分量
 *
 * @return	无
 */
ISL_EXPORT void ISL_IMG_colorAlterPixel(	ISL_byte ri, ISL_byte gi, ISL_byte bi, /* 原始rgb值 */
							ISL_byte ru, ISL_byte gu, ISL_byte bu, /* 用户定义 rgb值 */
							ISL_byte &ro, ISL_byte &go, ISL_byte &bo /* 返回的rgb值 */ );


/** 以下函数是针对完整图像 */


/**
 * @brief	用户值域范围调整(图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[in]	incement	调整值
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_contrastAlter( ISL_byte *& pic, int width, int height, int incement );



/**
 * @brief	对图像使用阈值法进行着色处理 (图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[in]	ru			用户定义三原色，红色分量
 * @param[in]	gu			用户定义三原色，绿色分量
 * @param[in]	bu			用户定义三原色，蓝色分量
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_colorAlter(	ISL_byte *& pic, int width, int height,
							ISL_byte ru, ISL_byte gu, ISL_byte bu/* 用户定义 rgb值 */ );



/**
 * @brief	使图像产生霓虹处理效果 (图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_neonLight(ISL_byte *& pic, int width, int height);



/**
 * @brief	使图像灰化 (图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_gray(ISL_byte *& pic, int width, int height);



/**
 * @brief	使图像黑白  (图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[in]	mid			黑白判断阈值（0 ~ 255）
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_black(ISL_byte *& pic, int width, int height, int mid = 127);



/**
 * @brief	使灰度图像黑白  (图片)
 *
 * @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
 * @param[in]	width		图像宽度
 * @param[in]	height		图像高度
 * @param[in]	mid			黑白判断阈值（0 ~ 255）
 * @param[out]	pic			返回修改后的图像
 *
 * @return	一切正常返回 true ，否则返回 false
 */
ISL_EXPORT bool ISL_IMG_grap2black(ISL_byte *& pic, int width, int height, int mid = 127);




/** @brief	使图像亮化  (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[in]	light		输入的亮度值（0~255）
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_light(ISL_byte *& pic, int width, int height, int light);



/** @brief	图像亮化反转(底片色)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_lightReverse(ISL_byte *& pic, int width, int height);



/** @brief	图像曝光处理 (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_exposal(ISL_byte *& pic, int width, int height);




/** @brief	使图像平滑处理  (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_smoothness(ISL_byte *& pic, int width, int height);



/** @brief	产生图像浮雕处理效果  (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_embossment(ISL_byte *& pic, int width, int height);



/** @brief	图像扩散处理  (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_spread(ISL_byte *& pic, int width, int height);



/** @brief	图像锐化处理  (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_sharp(ISL_byte *& pic, int width, int height);



/** @brief	对图像使用阈值法进行高通滤波(3X3)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[in]	mode		1 =  基本高通, 2 =  中等高通, 3 =  过量高通
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_highLVBO(ISL_byte *& pic, int width, int height, int mode = 1);



/** @brief	实现图像低通滤波(3X3)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_lowLVBO3X3(ISL_byte *& pic, int width, int height);



/** @brief	实现图像低通滤波(5X5)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_lowLVBO5X5(ISL_byte *& pic, int width, int height);



/** @brief	使图像水平增强(3X1)(可用于层位增强)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_horizontalGROW(ISL_byte *& pic, int width, int height);



/** @brief	使图像垂直增强(3X1)(可用于断层增强)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_verticalGROW(ISL_byte *& pic, int width, int height);\



/** @brief	使图像双向增强(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_doubleGROW(ISL_byte *& pic, int width, int height);



/** @brief	使图像产生马赛克效果(5x5)(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_mosaic(ISL_byte *& pic, int width, int height);



/** @brief	使图像产生轮廓 (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[in]	mid 		黑白判断阈值（0 ~ 255）
* @param[out]	pic			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_outLine(ISL_byte *& pic, int width, int height, int mid = 127);



/** @brief	图像数据二值化， 仅对黑白色有效 (图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[out]	out			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_binaryzation(ISL_byte *pic, int width, int height, float *out);



/** @brief	Sobel 边缘检测(图片)
*
* @param[in]	pic			输入的图像值/输出的处理后的图像值,数组大小为 width*height*4
* @param[in]	width		图像宽度
* @param[in]	height		图像高度
* @param[in]	type		当type为true时，差分结果取水平和垂直方向差分中较大者，否则取平均值
* @param[in]	scale		扫描尺度（0 ~ 1）
* @param[out]	out			返回修改后的图像
*
* @return	一切正常返回 true ，否则返回 false
*/
ISL_EXPORT bool ISL_IMG_SideSobel(ISL_byte *& pic, int width, int height, bool type = true, double scale = 0.5);

}/*End of ISLib*/

#endif /* PAI_FRAME_ISEISLIB_IMAGEPROCESSING_H_ */
