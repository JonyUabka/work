/**
 *	@file	ISL_Complex.h
 *	@brief	[Header file of Complex Functions], 复数计算函数；
 *	@see	ISeisLib Manual
 *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
 *	@date	2014-02-26
 *	@refer	SU CWP
 */

#include "ISL_Complex.h"

namespace ISLib
{

/******************************************************************
 *******			complex 相关函数的用户借口说明 				*******
 *******************************************************************/

// 计算complex 的摸
float ISL_rcabs(complex z)
{
	float x, y, ans, temp;
	x = fabs(z.r);
	y = fabs(z.i);

	if (x == 0.0)
		ans = y;
	else if (y == 0.0)
		ans = x;
	else if (x > y)
	{
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	}
	else
	{
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}

// complex 的加法运算
complex ISL_cadd(complex a, complex b)
{
	complex c;
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return c;
}

// complex 的减法运算
complex ISL_csub(complex a, complex b)
{
	complex c;
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return c;
}

// complex 的乘法运算
complex ISL_cmul(complex a, complex b)
{
	complex c;
	c.r = a.r * b.r - a.i * b.i;
	c.i = a.i * b.r + a.r * b.i;
	return c;
}

// complex 的除法运算
complex ISL_cdiv(complex a, complex b)
{
	complex c;
	float r, den;

	if (fabs(b.r) >= fabs(b.i))
	{
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else
	{
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}

	return c;
}

// 构造一个complex 结构
complex ISL_cmplx(float re, float im)
{
	complex c;
	c.r = re;
	c.i = im;
	return c;
}

// 计算公扼的complex
complex ISL_conjg(complex z)
{
	complex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

// 计算负的complex
complex ISL_cneg(complex z)
{
	complex c;
	c.r = -z.r;
	c.i = -z.i;
	return c;
}

//计算归一化的complex
complex ISL_cinv(complex z)
{
	complex c;
	float s;
    s = 1.0 / (z.r * z.r + z.i * z.i);
    c.r = z.r * s;
	c.i = -z.i * s;
	return c;
}
//complex 的开方运算
complex ISL_csqrt(complex z)
{
	float x, y, w, r;
	complex a;

	if (z.r == 0.0 && z.i == 0.0)
	{
		a.r = 0.0;
		a.i = 0.0;
		return a;
	}
	else
	{
		x = fabs(z.r);
		y = fabs(z.i);
		if (x >= y)
		{
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		}
		else
		{
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}
		if (z.r >= 0.0)
		{
			a.r = w;
			a.i = z.i / (2.0 * w);
		}
		else
		{
			a.i = (z.i >= 0.0) ? w : -w;
			a.r = z.i / (2.0 * a.i);
		}
		return a;
	}
}

//complex 的指数运算
complex ISL_cexp(complex z)
{
	float a;
	complex c;
	a = exp(z.r);
	c.r = a * cos(z.i);
	c.i = a * sin(z.i);
	return c;
}

//complex 乘常数运算
complex ISL_crmul(complex a, float x)
{
	complex c;
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

//complex 的对数运算
complex ISL_cln(complex b)
{
	float a;
	complex c;

	if (fabs(b.r) <= 0.00001)
	{
		c.r = 0.0;
		c.i = 0.0;
	}
	else
	{
		a = ISL_rcabs(b);
		c.r = log(a);
		c.i = (float) atan2(b.i, b.r);
	}
	return c;
}

}/*End of ISLib*/
