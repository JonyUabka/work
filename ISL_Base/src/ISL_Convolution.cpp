///**
// *	@file	ISL_Convolution.h
// *	@brief	[Header file of Convolution Functions], 卷积函数；
// *	@see	ISeisLib Manual
// *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
// *	@date	2014-02-26
// *	@refer	SU CWP
// */

//#include "ISL_Convolution.h"
//#include "ISL_FFT.h"

//namespace ISLib{
//static int _conv_proc_type = 0;
//static int _conv_scl_flag = 1;

///* internal function optimized for short x */
//void ISL_convs(int lx, int ifx, float *x, int ly, int ify, float *y, int lz,
//               int ifz, float *z)
//{
//    int ilx = ifx + lx - 1, ily = ify + ly - 1, ilz = ifz + lz - 1, i, j, ilow,
//            ihigh, jlow, jhigh;
//    float x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15,
//            x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28,
//            x29, ya, yb, z0, z1, sum;

//    x -= ifx;
//    y -= ify;
//    z -= ifz;

//    /* OFF LEFT:  i < ifx+ify */
//    ilow = ifz;
//    ihigh = ify + ifx - 1;
//    if (ihigh > ilz)
//        ihigh = ilz;
//    for (i = ilow; i <= ihigh; ++i)
//        z[i] = 0.0;

//    /* ROLLING ON:  ify+ifx <= i < ify+ilx */
//    ilow = ify + ifx;
//    if (ilow < ifz)
//        ilow = ifz;
//    ihigh = ify + ilx - 1;
//    if (ihigh > ilz)
//        ihigh = ilz;
//    jlow = ifx;
//    jhigh = ilow - ify;
//    for (i = ilow; i <= ihigh; ++i, ++jhigh)
//    {
//        for (j = jlow, sum = 0.0; j <= jhigh; ++j)
//            sum += x[j] * y[i - j];
//        z[i] = sum;
//    }

//    /* MIDDLE:  ify+ilx <= i <= ily+ifx */
//    ilow = ify + ilx;
//    if (ilow < ifz)
//        ilow = ifz;
//    ihigh = ily + ifx;
//    if (ihigh > ilz)
//        ihigh = ilz;
//    if (lx == 1)
//    {
//        x0 = x[ifx];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 2)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 3)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 4)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 5)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 6)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 7)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 8)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 9)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 10)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 11)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 12)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 13)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 14)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 15)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 16)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 17)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 18)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 19)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 20)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 21)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 22)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 23)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 24)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 25)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 26)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        x25 = x[ifx + 25];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z1 += x25 * yb;
//            ya = y[i - ifx - 25];
//            z0 += x25 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 27)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        x25 = x[ifx + 25];
//        x26 = x[ifx + 26];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z1 += x25 * yb;
//            ya = y[i - ifx - 25];
//            z0 += x25 * ya;
//            z1 += x26 * ya;
//            yb = y[i - ifx - 26];
//            z0 += x26 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 28)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        x25 = x[ifx + 25];
//        x26 = x[ifx + 26];
//        x27 = x[ifx + 27];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z1 += x25 * yb;
//            ya = y[i - ifx - 25];
//            z0 += x25 * ya;
//            z1 += x26 * ya;
//            yb = y[i - ifx - 26];
//            z0 += x26 * yb;
//            z1 += x27 * yb;
//            ya = y[i - ifx - 27];
//            z0 += x27 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 29)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        x25 = x[ifx + 25];
//        x26 = x[ifx + 26];
//        x27 = x[ifx + 27];
//        x28 = x[ifx + 28];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z1 += x25 * yb;
//            ya = y[i - ifx - 25];
//            z0 += x25 * ya;
//            z1 += x26 * ya;
//            yb = y[i - ifx - 26];
//            z0 += x26 * yb;
//            z1 += x27 * yb;
//            ya = y[i - ifx - 27];
//            z0 += x27 * ya;
//            z1 += x28 * ya;
//            yb = y[i - ifx - 28];
//            z0 += x28 * yb;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    else if (lx == 30)
//    {
//        x0 = x[ifx];
//        x1 = x[ifx + 1];
//        x2 = x[ifx + 2];
//        x3 = x[ifx + 3];
//        x4 = x[ifx + 4];
//        x5 = x[ifx + 5];
//        x6 = x[ifx + 6];
//        x7 = x[ifx + 7];
//        x8 = x[ifx + 8];
//        x9 = x[ifx + 9];
//        x10 = x[ifx + 10];
//        x11 = x[ifx + 11];
//        x12 = x[ifx + 12];
//        x13 = x[ifx + 13];
//        x14 = x[ifx + 14];
//        x15 = x[ifx + 15];
//        x16 = x[ifx + 16];
//        x17 = x[ifx + 17];
//        x18 = x[ifx + 18];
//        x19 = x[ifx + 19];
//        x20 = x[ifx + 20];
//        x21 = x[ifx + 21];
//        x22 = x[ifx + 22];
//        x23 = x[ifx + 23];
//        x24 = x[ifx + 24];
//        x25 = x[ifx + 25];
//        x26 = x[ifx + 26];
//        x27 = x[ifx + 27];
//        x28 = x[ifx + 28];
//        x29 = x[ifx + 29];
//        for (i = ilow; i <= ihigh - 1; i += 2)
//        {
//            ya = y[i + 1 - ifx];
//            z1 = x0 * ya;
//            yb = y[i - ifx];
//            z0 = x0 * yb;
//            z1 += x1 * yb;
//            ya = y[i - ifx - 1];
//            z0 += x1 * ya;
//            z1 += x2 * ya;
//            yb = y[i - ifx - 2];
//            z0 += x2 * yb;
//            z1 += x3 * yb;
//            ya = y[i - ifx - 3];
//            z0 += x3 * ya;
//            z1 += x4 * ya;
//            yb = y[i - ifx - 4];
//            z0 += x4 * yb;
//            z1 += x5 * yb;
//            ya = y[i - ifx - 5];
//            z0 += x5 * ya;
//            z1 += x6 * ya;
//            yb = y[i - ifx - 6];
//            z0 += x6 * yb;
//            z1 += x7 * yb;
//            ya = y[i - ifx - 7];
//            z0 += x7 * ya;
//            z1 += x8 * ya;
//            yb = y[i - ifx - 8];
//            z0 += x8 * yb;
//            z1 += x9 * yb;
//            ya = y[i - ifx - 9];
//            z0 += x9 * ya;
//            z1 += x10 * ya;
//            yb = y[i - ifx - 10];
//            z0 += x10 * yb;
//            z1 += x11 * yb;
//            ya = y[i - ifx - 11];
//            z0 += x11 * ya;
//            z1 += x12 * ya;
//            yb = y[i - ifx - 12];
//            z0 += x12 * yb;
//            z1 += x13 * yb;
//            ya = y[i - ifx - 13];
//            z0 += x13 * ya;
//            z1 += x14 * ya;
//            yb = y[i - ifx - 14];
//            z0 += x14 * yb;
//            z1 += x15 * yb;
//            ya = y[i - ifx - 15];
//            z0 += x15 * ya;
//            z1 += x16 * ya;
//            yb = y[i - ifx - 16];
//            z0 += x16 * yb;
//            z1 += x17 * yb;
//            ya = y[i - ifx - 17];
//            z0 += x17 * ya;
//            z1 += x18 * ya;
//            yb = y[i - ifx - 18];
//            z0 += x18 * yb;
//            z1 += x19 * yb;
//            ya = y[i - ifx - 19];
//            z0 += x19 * ya;
//            z1 += x20 * ya;
//            yb = y[i - ifx - 20];
//            z0 += x20 * yb;
//            z1 += x21 * yb;
//            ya = y[i - ifx - 21];
//            z0 += x21 * ya;
//            z1 += x22 * ya;
//            yb = y[i - ifx - 22];
//            z0 += x22 * yb;
//            z1 += x23 * yb;
//            ya = y[i - ifx - 23];
//            z0 += x23 * ya;
//            z1 += x24 * ya;
//            yb = y[i - ifx - 24];
//            z0 += x24 * yb;
//            z1 += x25 * yb;
//            ya = y[i - ifx - 25];
//            z0 += x25 * ya;
//            z1 += x26 * ya;
//            yb = y[i - ifx - 26];
//            z0 += x26 * yb;
//            z1 += x27 * yb;
//            ya = y[i - ifx - 27];
//            z0 += x27 * ya;
//            z1 += x28 * ya;
//            yb = y[i - ifx - 28];
//            z0 += x28 * yb;
//            z1 += x29 * yb;
//            ya = y[i - ifx - 29];
//            z0 += x29 * ya;
//            z[i + 1] = z1;
//            z[i] = z0;
//        }
//    }
//    if (ihigh >= ilow && (ihigh - ilow) % 2 == 0)
//    {
//        ilow = ihigh;
//        jlow = ifx;
//        jhigh = ilx;
//        for (i = ilow; i <= ihigh; ++i)
//        {
//            for (j = jlow, sum = 0.0; j <= jhigh; ++j)
//                sum += x[j] * y[i - j];
//            z[i] = sum;
//        }
//    }

//    /* ROLLING OFF:  ily+ifx < i <= ily+ilx */
//    ilow = ily + ifx + 1;
//    if (ilow < ifz)
//        ilow = ifz;
//    ihigh = ily + ilx;
//    if (ihigh > ilz)
//        ihigh = ilz;
//    jlow = ilow - ily;
//    jhigh = ilx;
//    for (i = ilow; i <= ihigh; ++i, ++jlow)
//    {
//        for (j = jlow, sum = 0.0; j <= jhigh; ++j)
//            sum += x[j] * y[i - j];
//        z[i] = sum;
//    }

//    /* OFF RIGHT:  ily+ilx < i */
//    ilow = ily + ilx + 1;
//    if (ilow < ifz)
//        ilow = ifz;
//    ihigh = ilz;
//    for (i = ilow; i <= ihigh; ++i)
//        z[i] = 0.0;
//}

///**********************************************************************************************************************
// *
// *  功能：互相关函数
// *
// *  说明：Compute z = x 与 y 的相关;
// *			ifx+lx-1
// *	z[i] =   sum    x[j]*y[i+j]  ;  i = ifz,...,ifz+lz-1
// *			j=ifx
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					lx					In			length of x array
// *		int					ifx					In			sample index of first x
// *		float *				x					In			array[lx] to be cross-correlated with y
// *		int					ly					In			length of y array
// *		int					ify					In			sample index of first y
// *		float *				y					In			array[ly] with which x is to be cross-correlated
// *		int					lz					In			length of z array
// *		int					ifz					In			sample index of first z
// *		float *				z					Out			array[lz] containing x cross-correlated with y
// *
// *  返回：无
// *
// **********************************************************************************************************************/
//void ISL_xcor(int lx, int ifx, float *x, int ly, int ify, float *y, int lz,
//              int ifz, float *z)
//{
//    int i, j;

//    float *xr = new float[lx];
//    if (xr == NULL)
//        return;

//    for (i = 0, j = lx - 1; i < lx; ++i, --j)
//        xr[i] = x[j];

//    _conv_proc_type = 0;
//    _conv_scl_flag = 1;

//    ISL_conv(lx, 1 - ifx - lx, xr, ly, ify, y, lz, ifz, z);

//    delete[] xr;
//}

///**********************************************************************************************************************
// *
// *  功能：卷积函数
// *
// *  说明：Compute z = x 与 y 的卷积;
// *			ifx+lx-1
// *	z[i] =   sum    x[j]*y[i+j]  ;  i = ifz,...,ifz+lz-1
// *			j=ifx
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					lx					In			length of x array
// *		int					lfx					In			sample index of first x
// *		float *				x					In			array[lx] to be convolved with y
// *		int					ly					In			length of y array
// *		int					lfy					In			sample index of first y
// *		float *				y					In			array[ly] with which x is to be convolved
// *		int					lz					In			length of z array
// *		int					lfz					In			sample index of first z
// *		float *				z					Out			array[lz] containing x convolved with y
// *
// *  返回：无
// *
// **********************************************************************************************************************/
//void ISL_conv(int lx, //length of x array
//              int ifx, //sample index of first x
//              float *x, //array[lx] to be convolved with y
//              int ly, //length of y array
//              int ify, //sample index of first y
//              float *y, //array[ly] with which x is to be convolved
//              int lz, //length of z array
//              int ifz, //sample index of first z
//              float *z //array[lz] containing x convolved with y【Out】
//              )
//{
//    if (_conv_proc_type == 0)
//    {
//        int ilx = ifx + lx - 1, ily = ify + ly - 1, ilz = ifz + lz - 1, i, j,
//                jlow, jhigh;
//        float sum = 0;

//        x -= ifx;
//        y -= ify;
//        z -= ifz;
//        for (i = ifz; i <= ilz; ++i)
//        {
//            jlow = i - ily;
//            if (jlow < ifx)
//                jlow = ifx;
//            jhigh = i - ify;
//            if (jhigh > ilx)
//                jhigh = ilx;
//            for (j = jlow, sum = 0.0; j <= jhigh; ++j)
//                sum += x[j] * y[i - j];
//            z[i] = sum;
//        }

//        if (_conv_scl_flag == 0)
//            return;

//        float xamp = 1.0;
//        for (i = ifz; i <= ilz; ++i)
//            if (sum < fabs(z[i]))
//                sum = fabs(z[i]);
//        xamp = xamp / sum;
//        for (i = ifz; i <= ilz; ++i)
//            z[i] *= xamp;

//    }
//    else
//    {
//        int ilx = ifx + lx - 1, ily = ify + ly - 1, ilz = ifz + lz - 1, i, j,
//                ilow, ihigh, jlow, jhigh;
//        float sa, sb, xa, xb, ya, yb, *t;

//        /* if x is longer than y, swap x and y */
//        if (lx > ly)
//        {
//            i = ifx;
//            ifx = ify;
//            ify = i;
//            i = ilx;
//            ilx = ily;
//            ily = i;
//            i = lx;
//            lx = ly;
//            ly = i;
//            t = x;
//            x = y;
//            y = t;
//        }

//        /* handle short x with special code */
//        if (lx >= 1 && lx <= 30)
//        {
//            ISL_convs(lx, ifx, x, ly, ify, y, lz, ifz, z);
//            return;
//        }

//        /* adjust pointers for indices of first samples */
//        x -= ifx;
//        y -= ify;
//        z -= ifz;

//        /* OFF LEFT:  i < ify+ifx */

//        /* zero output for all i */
//        ilow = ifz;
//        ihigh = ify + ifx - 1;
//        if (ihigh > ilz)
//            ihigh = ilz;
//        for (i = ilow; i <= ihigh; ++i)
//            z[i] = 0.0;

//        /* ROLLING ON:  ify+ifx <= i < ify+ilx */

//        /* if necessary, do one i so that number of j in overlap is odd */
//        if (i < ify + ilx && i <= ilz)
//        {
//            jlow = ifx;
//            jhigh = i - ify;
//            if ((jhigh - jlow) % 2)
//            {
//                sa = 0.0;
//                for (j = jlow; j <= jhigh; ++j)
//                    sa += x[j] * y[i - j];
//                z[i++] = sa;
//            }
//        }

//        /* loop over pairs of i and j */
//        ilow = i;
//        ihigh = ilx + ify - 1;
//        if (ihigh > ilz)
//            ihigh = ilz;
//        jlow = ifx;
//        jhigh = ilow - ify;
//        for (i = ilow; i < ihigh; i += 2, jhigh += 2)
//        {
//            sa = sb = 0.0;
//            xb = x[jhigh + 1];
//            yb = 0.0;
//            for (j = jhigh; j >= jlow; j -= 2)
//            {
//                sa += xb * yb;
//                ya = y[i - j];
//                sb += xb * ya;
//                xa = x[j];
//                sa += xa * ya;
//                yb = y[i + 1 - j];
//                sb += xa * yb;
//                xb = x[j - 1];
//            }
//            z[i] = sa;
//            z[i + 1] = sb;
//        }

//        /* if number of i is odd */
//        if (i == ihigh)
//        {
//            jlow = ifx;
//            jhigh = i - ify;
//            sa = 0.0;
//            for (j = jlow; j <= jhigh; ++j)
//                sa += x[j] * y[i - j];
//            z[i++] = sa;
//        }

//        /* MIDDLE:  ify+ilx <= i <= ily+ifx */

//        /* determine limits for i and j */
//        ilow = i;
//        ihigh = ily + ifx;
//        if (ihigh > ilz)
//            ihigh = ilz;
//        jlow = ifx;
//        jhigh = ilx;

//        /* if number of j is even, do j in pairs with no leftover */
//        if ((jhigh - jlow) % 2)
//        {
//            for (i = ilow; i < ihigh; i += 2)
//            {
//                sa = sb = 0.0;
//                yb = y[i + 1 - jlow];
//                xa = x[jlow];
//                for (j = jlow; j < jhigh; j += 2)
//                {
//                    sb += xa * yb;
//                    ya = y[i - j];
//                    sa += xa * ya;
//                    xb = x[j + 1];
//                    sb += xb * ya;
//                    yb = y[i - 1 - j];
//                    sa += xb * yb;
//                    xa = x[j + 2];
//                }
//                z[i] = sa;
//                z[i + 1] = sb;
//            }
//            /* else, number of j is odd, so do j in pairs with leftover */
//        }
//        else
//        {
//            for (i = ilow; i < ihigh; i += 2)
//            {
//                sa = sb = 0.0;
//                yb = y[i + 1 - jlow];
//                xa = x[jlow];
//                for (j = jlow; j < jhigh; j += 2)
//                {
//                    sb += xa * yb;
//                    ya = y[i - j];
//                    sa += xa * ya;
//                    xb = x[j + 1];
//                    sb += xb * ya;
//                    yb = y[i - 1 - j];
//                    sa += xb * yb;
//                    xa = x[j + 2];
//                }
//                z[i] = sa + x[jhigh] * y[i - jhigh];
//                z[i + 1] = sb + x[jhigh] * y[i + 1 - jhigh];
//            }
//        }

//        /* if number of i is odd */
//        if (i == ihigh)
//        {
//            sa = 0.0;
//            for (j = jlow; j <= jhigh; ++j)
//                sa += x[j] * y[i - j];
//            z[i++] = sa;
//        }

//        /* ROLLING OFF:  ily+ifx < i <= ily+ilx */

//        /* if necessary, do one i so that number of j in overlap is even */
//        if (i <= ily + ilx && i <= ilz)
//        {
//            jlow = i - ily;
//            jhigh = ilx;
//            if (!((jhigh - jlow) % 2))
//            {
//                sa = 0.0;
//                for (j = jlow; j <= jhigh; ++j)
//                    sa += x[j] * y[i - j];
//                z[i++] = sa;
//            }
//        }

//        /* number of j is now even, so loop over both i and j in pairs */
//        ilow = i;
//        ihigh = ily + ilx;
//        if (ihigh > ilz)
//            ihigh = ilz;
//        jlow = ilow - ily;
//        jhigh = ilx - 2; /* Dave's new patch */
//        for (i = ilow; i < ihigh; i += 2, jlow += 2)
//        {
//            sa = sb = 0.0;
//            xa = x[jlow];
//            yb = 0.0;
//            for (j = jlow; j < jhigh; j += 2)
//            {
//                sb += xa * yb;
//                ya = y[i - j];
//                sa += xa * ya;
//                xb = x[j + 1];
//                sb += xb * ya;
//                yb = y[i - 1 - j];
//                sa += xb * yb;
//                xa = x[j + 2];
//            }
//            sb += xa * yb;
//            ya = y[i - j];
//            sa += xa * ya;
//            xb = x[j + 1];
//            sb += xb * ya;
//            yb = y[i - 1 - j];
//            sa += xb * yb;
//            z[i] = sa;
//            z[i + 1] = sb;
//        }

//        /* if number of i is odd */
//        if (i == ihigh)
//        {
//            jlow = i - ily;
//            jhigh = ilx;
//            sa = 0.0;
//            for (j = jlow; j <= jhigh; ++j)
//                sa += x[j] * y[i - j];
//            z[i++] = sa;
//        }
//        /* OFF RIGHT:  ily+ilx < i */

//        /* zero output for all i */
//        ilow = i;
//        ihigh = ilz;
//        for (i = ilow; i <= ihigh; ++i)
//            z[i] = 0.0;
//    }
//}


///*用FFT计算两个有限长序列的线性卷积*/
//void ISL_conv(double *x, double *y, int m, int n, int len)
//{
//    int i, len2;
//    double t;
//    for (i = m; i < len; i++)
//        x[i] = 0.0;
//    for (i = n; i < len; i++)
//        y[i] = 0.0;

//    double pi[len];
//    for (i = 0; i < len; i++)
//        pi[i] = 0.0;

//    ISL_kfft_d(x, pi, len, 1);
//    ISL_kfft_d(y, pi, len, 1);

//    len2 = len / 2;
//    x[0] = x[0] * y[0];
//    x[len2] = x[len2] * y[len2];
//    for (i = 1; i < len2; i++) {
//        t = x[i] * y[i] - x[len - i] * y[len - i];
//        x[len - i] = x[i] * y[len - i] + x[len - i] * y[i];
//        x[i] = t;
//    }
//    ISL_kfft_d(x, pi, len, -1);
//}


///*计算合成记录的简单卷积， 输入两个任意长度的离散点序列，返回一个卷积后的结果*/
//void ISL_synConv(float *r, int nr, float *w, int nw, float *syn)
//{
//    int k;
//    float sum;
//    for (int i = 0; i < nr; i++) {
//        sum = 0.0;
//        for (int j = 0; j < nw; j++) {
//            k = i - j + nw / 2;
//            if (k >= 0 && k < nr)
//                sum += r[k] * w[j];
//        }
//        syn[i] = sum;
//    }
//}
//}
