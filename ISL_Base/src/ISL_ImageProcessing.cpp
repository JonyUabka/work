//***************************************************************************
//**
//**				ISeis Lib [Base Algorithm] - 图像处理
//**
//**				Writer By	Mouri Song
//**							Baihong Liu
//**							Qiang Yang
//**
//**				(Soft Center IGP)
//**
//**				DATA : 2014-08-20
//**
//****************************************************************************/

///**
//*	@file	ISL_ImageProcessing
//*	@brief	[Source file of ImageProcessing Functions], 图像处理；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
//*	@date	2014-08-20
//*	@refer	SU CWP
//*/


//#include "ISL_ImageProcessing.h"
//#include "ISL_FileOperation.h"


//namespace ISLib
//{

///**  读取BMP格式图片  */
//bool ISL_loadBMP(const char *name, ISL_byte *& pic, int &width, int &height)
//{
//	int columns, rows, numPixels;
//	ISL_byte * pixbuf = NULL;
//	int row, column;
//	ISL_byte * buf_p = NULL;
//	ISL_byte * buffer = NULL;
//	unsigned long length;
//	BMPHeader_t bmpHeader;
//	ISL_byte *bmpRGBA = NULL;
//	pic = NULL;

//	// load the file
//	FILE* pfile = fopen(name, "rb");
//	if (pfile == NULL) {
//		printf("LoadBMP: Load BMP files failed (%s)\n", name);
//		return false;
//	}
//	length = ISL_getFileSize(pfile);
//	buffer = (ISL_byte*) malloc(length * sizeof(ISL_byte));
//	fread(buffer, 1, length, pfile);
//	if (!buffer) {
//		printf("LoadBMP: Memory alloc failed (%s)\n", name);
//		return false;
//	}

//	buf_p = buffer;

//	bmpHeader.id[0] = *buf_p++;
//	bmpHeader.id[1] = *buf_p++;
//	bmpHeader.fileSize = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.reserved0 = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.bitmapDataOffset = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.bitmapHeaderSize = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.width = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.height = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.planes = (*(short *) buf_p);
//	buf_p += 2;
//	bmpHeader.bitsPerPixel = (*(short *) buf_p);
//	buf_p += 2;
//	bmpHeader.compression = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.bitmapDataSize = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.hRes = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.vRes = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.colors = (*(long *) buf_p);
//	buf_p += 4;
//	bmpHeader.importantColors = (*(long *) buf_p);
//	buf_p += 4;

//	memcpy(bmpHeader.palette, buf_p, sizeof(bmpHeader.palette));

//	cout << "bmpHeader.id[0] = " << bmpHeader.id[0] << endl;
//	cout << "bmpHeader.id[1] = " << bmpHeader.id[1] << endl;
//	cout << "bmpHeader.fileSize = " << bmpHeader.fileSize << endl;
//	cout << "bmpHeader.reserved0 = " << bmpHeader.reserved0 << endl;
//	cout << "bmpHeader.bitmapDataOffset = " << bmpHeader.bitmapDataOffset << endl;
//	cout << "bmpHeader.bitmapHeaderSize = " << bmpHeader.bitmapHeaderSize << endl;
//	cout << "bmpHeader.width = " << bmpHeader.width << endl;
//	cout << "bmpHeader.height = " << bmpHeader.height << endl;
//	cout << "bmpHeader.planes = " << bmpHeader.planes << endl;
//	cout << "bmpHeader.bitsPerPixel = " << bmpHeader.bitsPerPixel << endl;
//	cout << "bmpHeader.compression = " << bmpHeader.compression << endl;
//	cout << "bmpHeader.bitmapDataSize = " << bmpHeader.bitmapDataSize << endl;
//	cout << "bmpHeader.hRes = " << bmpHeader.hRes << endl;
//	cout << "bmpHeader.vRes = " << bmpHeader.vRes << endl;
//	cout << "bmpHeader.colors = " << bmpHeader.colors << endl;
//	cout << "bmpHeader.importantColors = " << bmpHeader.importantColors << endl;

//	if (bmpHeader.bitsPerPixel == 8)
//		buf_p += 1024;

//	if (bmpHeader.id[0] != 'B' && bmpHeader.id[1] != 'M') {
//		printf("LoadBMP: only Windows-style BMP files supported (%s)\n", name);
//		return false;
//	}
//	if (bmpHeader.fileSize != length) {
//		printf("LoadBMP: header size does not match file size (%d vs. %d) (%s)\n",
//				(int) bmpHeader.fileSize, (int) length, name);
//		return false;
//	}
//	if (bmpHeader.compression != 0) {
//		printf("LoadBMP: only uncompressed BMP files supported (%s)\n", name);
//		return false;
//	}
//	if (bmpHeader.bitsPerPixel < 8) {
//		printf("LoadBMP: monochrome and 4-bit BMP files not supported (%s)\n", name);
//		return false;
//	}

//	columns = bmpHeader.width;
//	rows = bmpHeader.height;

//	if (bmpHeader.bitsPerPixel == 24) {
//		if ((columns & 3) != 0) //检查宽度是否为4倍数
//			columns = (columns & ~3) + 4; //修正位图宽度值,对齐到4的倍数,不然图像会变形
//	}

//	if (rows < 0)
//		rows = -rows;
//	numPixels = columns * rows;

//	width = columns;
//	height = rows;

//	bmpRGBA = (ISL_byte*) malloc(numPixels * 4);
//	pic = bmpRGBA;
//	buf_p = buffer + 54;

//	for (row = rows - 1; row >= 0; row--) {
//		pixbuf = bmpRGBA + row * columns * 4;

//		for (column = 0; column < columns; column++) {
//			unsigned char red, green, blue, alpha;
//			int palIndex;
//			unsigned short shortPixel;

//			switch (bmpHeader.bitsPerPixel)
//			{
//				case 8: {
//					palIndex = *buf_p++;
//					*pixbuf++ = bmpHeader.palette[palIndex][2];
//					*pixbuf++ = bmpHeader.palette[palIndex][1];
//					*pixbuf++ = bmpHeader.palette[palIndex][0];
//					*pixbuf++ = 0xff;
//				}
//				break;

//				case 16: {
//					shortPixel = *(unsigned short *) pixbuf;
//					pixbuf += 2;
//					*pixbuf++ = (shortPixel & (31 << 10)) >> 7;
//					*pixbuf++ = (shortPixel & (31 << 5)) >> 2;
//					*pixbuf++ = (shortPixel & (31)) << 3;
//					*pixbuf++ = 0xff;
//				}
//				break;

//				case 24: {
//					blue = *buf_p++;
//					green = *buf_p++;
//					red = *buf_p++;
//					*pixbuf++ = red;
//					*pixbuf++ = green;
//					*pixbuf++ = blue;
//					*pixbuf++ = 255;
//				}
//				break;

//				case 32: {
//					blue = *buf_p++;
//					green = *buf_p++;
//					red = *buf_p++;
//					alpha = *buf_p++;
//					*pixbuf++ = red;
//					*pixbuf++ = green;
//					*pixbuf++ = blue;
//					*pixbuf++ = alpha;
//				}
//				break;

//				default: {
//					printf("LoadBMP: illegal pixel_size '%d' in file '%s'\n",
//							bmpHeader.bitsPerPixel, name);
//					return false;
//				}
//			}
//		}
//	}

//	free(buffer);
//	fclose(pfile);
//	return true;
//}


///** 得到灰度值 */
//ISL_byte ISL_IMG_grayPixel(ISL_byte r, ISL_byte g, ISL_byte b)
//{
//	int gray = 0;
//	if (r > g)
//		gray = r;
//	else
//		gray = g;

//	if (gray < b)
//		gray = b;

//	return gray;
//}


///** 得到黑白值, mid = 黑白判断阈值（0 ~ 255）*/
//ISL_byte ISL_IMG_blackPixel(ISL_byte in, int mid)
//{
//	int gray = 0;
//	if ((int) in > mid)
//		gray = 255;
//	else
//		gray = 0;

//	return gray;
//}


///** 亮度设置 */
//ISL_byte ISL_IMG_lightPixel(ISL_byte in, int light)
//{
//	int a = 0;
//	a = (int) (in * light / 100);

//	if (a < 0)
//		a = 0;
//	if (a > 255)
//		a = 255;

//	return a;
//}


///** 亮度反转(底片色) */
//ISL_byte ISL_IMG_lightReversePixel(ISL_byte in)
//{
//	int a = 255 - in;
//	return a;
//}


///** 图像曝光处理 */
//ISL_byte ISL_IMG_exposalPixel(ISL_byte in)
//{
//	int a = (int) in;
//	a = (a > 128) ? a : (255 - a); //调整

//	return a;
//}


///** 用户值域范围调整 */
//void ISL_IMG_contrastAlterPixel(	ISL_byte ri,
//							ISL_byte gi,
//							ISL_byte bi,
//							int incement,
//							ISL_byte &ro,
//							ISL_byte &go,
//							ISL_byte &bo )
//{
//	int nHigh = 255 - incement;

//	//对于极端情况加以处理
//	if (nHigh < incement) {
//		nHigh = 127;
//		incement = 120;
//	}
//	if (incement < -127)
//		incement = -120;

//	//扩展或压缩区间的长度
//	int nStretch = 255;
//	if (incement >= 0)
//		nStretch = 255 - 2 * incement;
//	else
//		nStretch = 255 + 2 * incement;

//	if (incement >= 0) /* m_Increment>=0时 */
//	{
//		//取得当前点（蓝色）的值，调整
//		if (bi <= incement)
//			bo = 0;
//		else if (bi > nHigh)
//			bo = 255;
//		else
//			bo = (ISL_byte)((((int) bi - incement) * 255) / nStretch);

//		//取得当前点（绿色）的值，调整
//		if (gi <= incement)
//			go = 0;
//		else if (gi > nHigh)
//			go = 255;
//		else
//			go = (ISL_byte)((((int) gi - incement) * 255) / nStretch);

//		//取得当前点（红色）的值，调整
//		if (ri <= incement)
//			ro = 0;
//		else if (ri > nHigh)
//			ro = 255;
//		else
//			ro = (ISL_byte)((((int) ri - incement) * 255) / nStretch);

//	} else { /* m_Increment < 0 时 */
//		bo = (ISL_byte)((((int) (bi) * nStretch) / 255) - incement);
//		go = (ISL_byte)((((int) (gi) * nStretch) / 255) - incement);
//		ro = (ISL_byte)((((int) (ri) * nStretch) / 255) - incement);
//	}
//}



///** 对图像使用阈值法进行着色处理 */
//void ISL_IMG_colorAlterPixel(	ISL_byte ri, ISL_byte gi, ISL_byte bi, /* 原始rgb值 */
//							ISL_byte ru, ISL_byte gu, ISL_byte bu, /* 用户定义 rgb值 */
//							ISL_byte &ro, ISL_byte &go, ISL_byte &bo /* 返回的rgb值 */ )
//{
//	ISL_byte gray = (ISL_byte)(((short) ri * 59 + (short) gi * 30 + (short) bi * 11) / 100);

//	bo = (ISL_byte)((bu * gray) / 255);
//	go = (ISL_byte)((gu * gray) / 255);
//	ro = (ISL_byte)((ru * gray) / 255);
//}



///** 用户值域范围调整(图片) */
//bool ISL_IMG_contrastAlter( ISL_byte *& pic, int width, int height, int incement )
//{
//	if (pic == NULL || width < 0 || height < 0)
//		return false;

//	long idx = 0;
//	ISL_byte r = 0, g = 0, b = 0;

//	for (int j = 0; j < height; j++) {
//		for (int i = 0; i < width; i++) {
//			ISL_IMG_contrastAlterPixel(pic[idx], pic[idx + 1], pic[idx + 2], incement, r, g, b);
//			pic[idx] = r;
//			pic[idx + 1] = g;
//			pic[idx + 2] = b;

//			idx = idx + 4;
//		}
//	}
//	return true;
//}



///** 对图像使用阈值法进行着色处理 (图片) */
//bool ISL_IMG_colorAlter( 	ISL_byte *& pic, int width, int height,
//							ISL_byte ru, ISL_byte gu, ISL_byte bu/* 用户定义 rgb值 */ )
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	long idx = 0;
//	ISL_byte r=0, g=0, b=0;

//	for (int j = 0; j < height; j++) {
//		for (int i = 0; i < width; i++) {
//			ISL_IMG_colorAlterPixel(	pic[idx], pic[idx+1], pic[idx+2],
//								ru, gu, bu,
//								r, g, b);
//			pic[idx] = r;
//			pic[idx + 1] = g;
//			pic[idx + 2] = b;

//			idx = idx + 4;
//		}
//	}

//	return true;
//}


///** 使图像产生霓虹处理效果(图片) */
//bool ISL_IMG_neonLight(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0 ) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	for (int j = 0; j < height - 4; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 1; i++) // 每列
//		{

//			int pby_pt = 0;
//			//对像素执行算法
//			pby_pt = (*(pic + (height - j - 1) * DibWidth + i) - *(pic + (height - j - 1) * DibWidth + i + 3))
//					* (*(pic + (height - j - 1) * DibWidth + i) - *(pic + (height - j - 1) * DibWidth + i + 3))
//					+ (*(pic + (height - j - 1) * DibWidth + i) - *(pic + (height - j - 2) * DibWidth + i))
//					* (*(pic + (height - j - 1) * DibWidth + i) - *(pic + (height - j - 2) * DibWidth + i));

//			*(p_temp + (height - j - 1) * DibWidth + i) = 2 * int(sqrt(pby_pt));

//			/** 判断合法性 */
//			if ( *(p_temp + (height - j - 1) * DibWidth + i) < 0)
//				 *(p_temp + (height - j - 1) * DibWidth + i) = 0;
//			if ( *(p_temp + (height - j - 1) * DibWidth + i) > 255)
//				 *(p_temp + (height - j - 1) * DibWidth + i) = 255;
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像灰化 (图片) */
//bool ISL_IMG_gray(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0)
//		return false;

//	long idx = 0;
//	ISL_byte gray = 0;
//	for (int j = 0; j < height; j++) {
//		for (int i = 0; i < width; i++) {
//			gray = ISL_IMG_grayPixel(pic[idx], pic[idx + 1], pic[idx + 2]);
//			pic[idx] = gray;
//			pic[idx + 1] = gray;
//			pic[idx + 2] = gray;

//			idx = idx + 4;
//		}
//	}

//	return true;
//}


///** 使图像黑白 (图片) */
//bool ISL_IMG_black(ISL_byte *& pic, int width, int height, int mid)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	if(ISL_IMG_gray(pic, width, height) == false)
//		return false;

//	for (long j = 0; j < height*width*4; j++) {
//		pic[j] = ISL_IMG_blackPixel(pic[j], mid);
//	}

//	return true;
//}


///** 使灰度图像黑白(图片) , mid = 黑白判断阈值（0 ~ 255）*/
//bool ISL_IMG_grap2black(ISL_byte *& pic, int width, int height, int mid)
//{
//	if (pic == NULL || width < 0 || height < 0)
//		return false;

//	for (long j = 0; j < height * width * 4; j++) {
//		pic[j] = ISL_IMG_blackPixel(pic[j], mid);
//	}

//	return true;
//}


///** 使图像亮化 (图片) */
//bool ISL_IMG_light(ISL_byte *& pic, int width, int height, int light)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	for (long j = 0; j < height*width*4; j++) {
//		pic[j] = ISL_IMG_lightPixel(pic[j], light);
//	}

//	return true;
//}


///** 图像亮化反转(底片色)(图片) */
//bool ISL_IMG_lightReverse(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0)
//		return false;

//	for (long j = 0; j < height * width * 4; j++) {
//		pic[j] = ISL_IMG_lightReversePixel(pic[j]);
//	}

//	return true;
//}


///** 图像曝光处理 (图片) */
//bool ISL_IMG_exposal(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0)
//		return false;

//	for (long j = 0; j < height * width * 4; j++) {
//		pic[j] = ISL_IMG_exposalPixel(pic[j]);
//	}

//	return true;
//}


///** 使图像平滑处理  (图片) */
//bool ISL_IMG_smoothness(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[3][3];////定义(3x3)矩阵
//	h[0][0] = 1;	h[0][1] = 1;	h[0][2] = 1;
//	h[1][0] = 1;	h[1][1] = 1;	h[1][2] = 1;
//	h[2][0] = 1;	h[2][1] = 1;	h[2][2] = 1;

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			double pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6))
//			//对应的第1行的值乘以矩阵对应值，再相加
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[1][1] * (*(pic + (height - j - 2) * DibWidth + i + 3))
//					+ h[1][2] * (*(pic + (height - j - 2) * DibWidth + i + 6))
//			//对应的第2行的值乘以矩阵对应值，再相加
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i))
//					+ h[2][1] * (*(pic + (height - j - 3) * DibWidth + i + 3))
//					+ h[2][2] * (*(pic + (height - j - 3) * DibWidth + i + 6));

//			*(p_temp + (height - j - 2) * DibWidth + i + 3) = abs(int(pby_pt / 9));//取总和的的平均值
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 产生图像浮雕处理效果  (图片) */
//bool ISL_IMG_embossment(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	for (int j = 0; j < height; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 4; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对像素得每个分量执行算法
//			pby_pt = *(pic + (height - j - 1) * DibWidth + i) -
//					 *(pic + (height - j - 1) * DibWidth + i + 3) + 128;

//			*(p_temp + (height - j - 1) * DibWidth + i + 3) = pby_pt;

//			//检验合法性
//			if (*(p_temp + (height - j - 1) * DibWidth + i + 3) < 0)
//				*(p_temp + (height - j - 1) * DibWidth + i + 3) = 0;

//			else if (*(p_temp + (height - j - 1) * DibWidth + i + 3) > 255)
//				*(p_temp + (height - j - 1) * DibWidth + i + 3) = 255;
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}

///** 图像扩散处理  (图片) */
//bool ISL_IMG_spread(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	for (int j = 0; j < height - 4; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 14; i++) // 每列
//		{
//			int m = 0, n = 0;
//			m = rand() % 5; //取得行随机数
//			n = rand() % 5; //取得列随机数
//			int pby_pt = 0;
//			pby_pt = *(pic + (height - j - 1 - m) * DibWidth + i + 3 * n);//得到对应随机像素值
//			*(p_temp + (height - j - 3) * DibWidth + i + 6) = pby_pt;
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 图像锐化处理  (图片) */
//bool ISL_IMG_sharp(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	for (int j = 0; j < height - 4; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 14; i++) // 每列
//		{
//			int pby_pt = 0;
//			pby_pt = *(pic + (height - j - 2) * DibWidth + i + 3)
//					- *(pic + (height - j - 1) * DibWidth + i);

//			*(p_temp + (height - j - 2) * DibWidth + i + 3) = *(pic
//					+ (height - j - 2) * DibWidth + i + 3) + abs( (int)(pby_pt / 4) );

//			if (*(p_temp + (height - j - 2) * DibWidth + i + 3) > 255)
//				*(p_temp + (height - j - 2) * DibWidth + i + 3) = 255;
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}



///** 对图像使用阈值法进行高通滤波(3X3)(图片) */
//bool ISL_IMG_highLVBO(ISL_byte *& pic, int width, int height, int mode)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[3][3];  ////定义(3x3)矩阵
//	if (mode == 1) { //矩阵1（基本高通）
//		h[0][0] = 1;	h[0][1] = -2;	h[0][2] = 1;
//		h[1][0] = -2;	h[1][1] = 5;	h[1][2] = -2;
//		h[2][0] = 1;	h[2][1] = -2;	h[2][2] = 1;
//	} else if (mode == 2) { //矩阵2（中等高通）
//		h[0][0] = 0;	h[0][1] = -1;	h[0][2] = 0;
//		h[1][0] = -1;	h[1][1] = 5;	h[1][2] = -1;
//		h[2][0] = 0;	h[2][1] = -1;	h[2][2] = 0;
//	} else { //矩阵3（过量高通）
//		h[0][0] = -1;	h[0][1] = -1;	h[0][2] = -1;
//		h[1][0] = -1;	h[1][1] = 9;	h[1][2] = -1;
//		h[2][0] = -1;	h[2][1] = -1;	h[2][2] = -1;
//	}

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6))
//			//对应的第1行的值乘以矩阵对应值，再相加
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[1][1] * (*(pic + (height - j - 2) * DibWidth + i + 3))
//					+ h[1][2] * (*(pic + (height - j - 2) * DibWidth + i + 6))
//			//对应的第2行的值乘以矩阵对应值，再相加
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i))
//					+ h[2][1] * (*(pic + (height - j - 3) * DibWidth + i + 3))
//					+ h[2][2] * (*(pic + (height - j - 3) * DibWidth + i + 6));

//			*(p_temp + (height - j - 2) * DibWidth + i + 3) = abs(pby_pt);

//			if (pby_pt > 255) //判断是否越界
//				*(p_temp + (height - j - 2) * DibWidth + i + 3) = 255;
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}

///** 实现图像低通滤波(3X3)(图片) */
//bool ISL_IMG_lowLVBO3X3(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	double h[3][3];////定义(3x3)矩阵
//	h[0][0] = 0.1;	h[0][1] = 0.1;	h[0][2] = 0.1;
//	h[1][0] = 0.1;	h[1][1] = 0.2;	h[1][2] = 0.1;
//	h[2][0] = 0.1;	h[2][1] = 0.1;	h[2][2] = 0.1;

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			double pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6))
//			//对应的第0行的值乘以矩阵对应值，再相加
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[1][1] * (*(pic + (height - j - 2) * DibWidth + i + 3))
//					+ h[1][2] * (*(pic + (height - j - 2) * DibWidth + i + 6))
//			//对应的第0行的值乘以矩阵对应值，再相加
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i))
//					+ h[2][1] * (*(pic + (height - j - 3) * DibWidth + i + 3))
//					+ h[2][2] * (*(pic + (height - j - 3) * DibWidth + i + 6));

//			*(p_temp + (height - j - 2) * DibWidth + i + 3) = abs(int(pby_pt));
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 实现图像低通滤波(5X5)(图片) */
//bool ISL_IMG_lowLVBO5X5(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[5][5];//定义(5x5)矩阵
//	h[0][0] = 1;  h[0][1] = 1; h[0][2] = 1; h[0][3] = 1; h[0][4] = 1;
//	h[1][0] = 1;  h[1][1] = 2; h[1][2] = 2; h[1][3] = 2; h[1][4] = 1;
//	h[2][0] = 1;  h[2][1] = 2; h[2][2] = 3; h[2][3] = 2; h[2][4] = 1;
//	h[3][0] = 1;  h[3][1] = 2; h[3][2] = 2; h[3][3] = 2; h[3][4] = 1;
//	h[4][0] = 1;  h[4][1] = 1; h[4][2] = 1; h[4][3] = 1; h[4][4] = 1;

//	for (int j = 0; j < height - 4; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 14; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6))
//					+ h[0][3] * (*(pic + (height - j - 1) * DibWidth + i + 9))
//					+ h[0][4] * (*(pic + (height - j - 1) * DibWidth + i + 12))
//			//对应的第1行的值乘以矩阵对应值，再相加
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[1][1] * (*(pic + (height - j - 2) * DibWidth + i + 3))
//					+ h[1][2] * (*(pic + (height - j - 2) * DibWidth + i + 6))
//					+ h[1][3] * (*(pic + (height - j - 2) * DibWidth + i + 9))
//					+ h[1][4] * (*(pic + (height - j - 2) * DibWidth + i + 12))
//			//对应的第2行的值乘以矩阵对应值，再相加
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i))
//					+ h[2][1] * (*(pic + (height - j - 3) * DibWidth + i + 3))
//					+ h[2][2] * (*(pic + (height - j - 3) * DibWidth + i + 6))
//					+ h[2][3] * (*(pic + (height - j - 3) * DibWidth + i + 9))
//					+ h[2][4] * (*(pic + (height - j - 3) * DibWidth + i + 12))
//			//对应的第3行的值乘以矩阵对应值，再相加
//					+ h[3][0] * (*(pic + (height - j - 4) * DibWidth + i))
//					+ h[3][1] * (*(pic + (height - j - 4) * DibWidth + i + 3))
//					+ h[3][2] * (*(pic + (height - j - 4) * DibWidth + i + 6))
//					+ h[3][3] * (*(pic + (height - j - 4) * DibWidth + i + 9))
//					+ h[3][4] * (*(pic + (height - j - 4) * DibWidth + i + 12))
//			//对应的第4行的值乘以矩阵对应值，再相加
//					+ h[4][0] * (*(pic + (height - j - 5) * DibWidth + i))
//					+ h[4][1] * (*(pic + (height - j - 5) * DibWidth + i + 3))
//					+ h[4][2] * (*(pic + (height - j - 5) * DibWidth + i + 6))
//					+ h[4][3] * (*(pic + (height - j - 5) * DibWidth + i + 9))
//					+ h[4][4] * (*(pic + (height - j - 5) * DibWidth + i + 12));

//			//为了计算方便我们把除以35（矩阵权和）放在求总和之后
//			*(p_temp + (height - j - 3) * DibWidth + i + 6) = abs(int(pby_pt / 35));
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像水平增强(3X1)(可用于层位增强)(图片) */
//bool ISL_IMG_horizontalGROW(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[3][1];//定义(3x1)矩阵
//	h[0][0] = -1;	h[1][0] = 2;	h[2][0] = -1;

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对应的3行的值乘分别以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i));

//			if (pby_pt > 20)
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt) + 100;
//			else
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt);
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像垂直增强(3X1)(可用于断层增强)(图片) */
//bool ISL_IMG_verticalGROW(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[1][3];//定义(1x3)矩阵
//	h[0][0] = -1;	h[0][1] = 2;	h[0][2] = -1;

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6));

//			if (pby_pt > 20)
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt) + 100;
//			else
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt);
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像双向增强(图片) */
//bool ISL_IMG_doubleGROW(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	int h[3][3];//定义(3x3)矩阵
//	h[0][0] = -1;	h[0][1] = -1;	h[0][2] = -1;
//	h[1][0] = -1;	h[1][1] = 8;	h[1][2] = -1;
//	h[2][0] = -1;	h[2][1] = -1;	h[2][2] = -1;

//	for (int j = 0; j < height - 2; j++) // 每行
//	{
//		for (int i = 0; i < DibWidth - 8; i++) // 每列
//		{
//			int pby_pt = 0;
//			//对应的第0行的值乘以矩阵对应值，再相加
//			pby_pt = h[0][0] * (*(pic + (height - j - 1) * DibWidth + i))
//					+ h[0][1] * (*(pic + (height - j - 1) * DibWidth + i + 3))
//					+ h[0][2] * (*(pic + (height - j - 1) * DibWidth + i + 6))
//			//对应的第1行的值乘以矩阵对应值，再相加
//					+ h[1][0] * (*(pic + (height - j - 2) * DibWidth + i))
//					+ h[1][1] * (*(pic + (height - j - 2) * DibWidth + i + 3))
//					+ h[1][2] * (*(pic + (height - j - 2) * DibWidth + i + 6))
//			//对应的第2行的值乘以矩阵对应值，再相加
//					+ h[2][0] * (*(pic + (height - j - 3) * DibWidth + i))
//					+ h[2][1] * (*(pic + (height - j - 3) * DibWidth + i + 3))
//					+ h[2][2] * (*(pic + (height - j - 3) * DibWidth + i + 6));

//			if (pby_pt > 20)
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt) + 100;
//			else
//				*(p_temp + (height - j - 2) * DibWidth + i) = abs(pby_pt);
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像产生马赛克效果(5x5)(图片) */
//bool ISL_IMG_mosaic(ISL_byte *& pic, int width, int height)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];

//	for (int j = 0; j < height - 4; j += 5) // 每行
//	{
//		for (int i = 0; i < DibWidth - 14; i += 15) // 每列
//		{ //对应周围(5x5)矩阵蓝色值求和平均
//			int pby_pt = 0;
//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 15; n += 3) {
//					pby_pt += *(pic + (height - j - 1 - m) * DibWidth + i + n);
//				}

//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 14; n += 3) {
//					*(p_temp + (height - j - 1 - m) * DibWidth + i + n) = int( pby_pt / 25);
//				}
//			//对应周围(5x5)矩阵绿色值求和平均
//			pby_pt = 0;
//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 15; n += 3) {
//					pby_pt += *(pic + (height - j - 1 - m) * DibWidth + i + n + 1);
//				}
//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 14; n += 3) {
//					*(p_temp + (height - j - 1 - m) * DibWidth + i + n + 1) = int(pby_pt / 25);
//				}
//			//对应周围(5x5)矩阵红色值求和平均
//			pby_pt = 0;
//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 15; n += 3) {
//					pby_pt += *(pic + (height - j - 1 - m) * DibWidth + i + n + 2);
//				}
//			for (int m = 0; m < 5; m++)
//				for (int n = 0; n < 14; n += 3) {
//					*(p_temp + (height - j - 1 - m) * DibWidth + i + n + 2) = int(pby_pt / 25);
//				}
//		}
//	}

//	memcpy(pic, p_temp, height * DibWidth); // 复制处理后的图像
//	delete[] p_temp; //删除暂时分配内存
//	p_temp = NULL;

//	return true;
//}


///** 使图像产生轮廓 (图片) */
//bool ISL_IMG_outLine(ISL_byte *& pic, int width, int height, int mid)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int DibWidth = width * 4;
//	ISL_byte *p_temp = new ISL_byte[height * DibWidth];
//	memcpy(p_temp, pic, height * DibWidth);

//	if(ISL_IMG_smoothness(p_temp, width, height) == false)
//		return false;

//	for (long j = 0; j < height*width*4; j++) {
//		pic[j] = pic[j] - p_temp[j];
//	}

//	if(ISL_IMG_black(pic, width, height, mid) == false)
//		return false;

//	return true;
//}


///** 图像数据二值化， 仅对黑白色有效 (图片) */
//bool ISL_IMG_binaryzation(ISL_byte *pic, int width, int height, float *out)
//{
//	if (pic == NULL || width < 0 || height < 0 || out == NULL)
//		return false;

//	long idx = 0;

//	for (int j = 0; j < height * width; j++) {
//		if (pic[idx] == 0)
//			out[j] = 1;
//		else
//			out[j] = 0;

//		idx = idx + 4;
//	}

//	return true;
//}


///** Sobel 边缘检测 ,当type为true时，差分结果取水平和垂直方向差分中较大者，否则取平均值(图片) */
//bool ISL_IMG_SideSobel(ISL_byte *& pic, int width, int height, bool type, double scale)
//{
//	if (pic == NULL || width < 0 || height < 0) return false;

//	int x, y, a, aHr, aHg, aHb, aVr, aVg, aVb, aH, aV;
//	long n;

//	int w = width;
//	int h = height;

//	//依次处理每个像素
//	for (y = 1; y < h - 1; y++){
//		for (x = 1; x < w - 1; x++) {
//			//计算像素的偏移位置
//			n = (y * w + x) * 4;
//			//计算红色分量水平灰度差
//			aHr = abs((pic[n - w * 4 - 4] + pic[n - 4] * 2 + pic[n + w * 4 - 4])
//					- (pic[n - w * 4 + 4] + pic[n + 4] * 2 + pic[n + w * 4 + 4]));
//			//计算红色分量垂直灰度差
//			aVr = abs((pic[n - w * 4 - 4] + pic[n - w * 4] * 2 + pic[n - w * 4 + 4])
//					- (pic[n + w * 4 - 4] + pic[n + w * 4] * 2 + pic[n + w * 4 + 4]));
//			//计算绿色分量水平灰度差
//			aHg = abs((pic[n - w * 4 - 4 + 1] + pic[n - 4 + 1] * 2 + pic[n + w * 4 - 4 + 1])
//					- (pic[n - w * 4 + 4 + 1] + pic[n + 4 + 1] * 2 + pic[n + w * 4 + 4 + 1]));
//			//计算绿色分量垂直灰度差
//			aVg = abs((pic[n - w * 4 - 4 + 1] + pic[n - w * 4 + 1] * 2 + pic[n - w * 4 + 4 + 1])
//					- (pic[n + w * 4 - 4 + 1] + pic[n + w * 4 + 1] * 2 + pic[n + w * 4 + 4 + 1]));
//			//计算蓝色分量水平灰度差
//			aHb = abs((pic[n - w * 4 - 4 + 2] + pic[n - 4 + 2] * 2 + pic[n + w * 4 - 4 + 2])
//					- (pic[n - w * 4 + 4 + 2] + pic[n + 4 + 2] * 2 + pic[n + w * 4 + 4 + 2]));
//			//计算蓝色分量垂直灰度差
//			aVb = abs((pic[n - w * 4 - 4 + 2] + pic[n - w * 4 + 2] * 2 + pic[n - w * 4 + 4 + 2])
//					- (pic[n + w * 4 - 4 + 2] + pic[n + w * 4 + 2] * 2 + pic[n + w * 4 + 4 + 2]));

//			//计算水平综合灰度差
//			aH = aHr + aHg + aHb;
//			//计算垂直综合灰度差
//			aV = aVr + aVg + aVb;

//			if (type) {
//				//取水平和垂直方向差分中较大者
//				if (aH > aV)
//					a = aH;
//				else
//					a = aV;
//			} else {
//				//取水平和垂直方向差分的平均值
//				a = (aH + aV) / 2;
//			}

//			a = a * scale;

//			a = a > 255 ? 255 : a;
//			//生成边缘扫描结果
//		//	SetPixel(image1, n, a);

//			pic[n - w * 4 - 4] = a;
//			pic[n - w * 4 - 4 + 1] = a;
//			pic[n - w * 4 - 4 + 2] = a;
//		}
//	}
//	return true;
//}

//}/*End of ISLib
