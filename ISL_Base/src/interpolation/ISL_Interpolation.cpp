/****************************************************************************
 **
 **				ISeis Lib [Base Algorithm] - 插值与拟合
 **
 **				Writer By	Mouri Song
 **							Baihong Liu
 **							Qiang Yang
 **
 **				(Soft Center IGP)
 **
 **				DATA : 2014-03-19
 **
 ****************************************************************************/

/**
 *	@file	ISL_Interpolation.h
 *	@brief	[Header file of Interpolation Functions], 插值与拟合；
 *	@see	ISeisLib Manual
 *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
 *	@date	2014-03-20
 *	@refer	SU CWP
 */

#include "ISL_Interpolation.h"

namespace ISLib
{

/* 等间隔一维线性插值得输出X位置(辅助函数) */
void Interpolation::ISL_intlin_getx(int nin, float *xin, int nout, float *xout)
{
	if (nin == nout)
	{
		for (int i = 0; i < nout; i++)
		{
			xout[i] = xin[i];
		}
		return;
	}

	float diff, dx;
	float first, end;

	first = xin[0];
	end = xin[nin - 1];
	diff = xin[nin - 1] - xin[0];
	dx = diff / (float) nout;

	for (int i = 0; i < nout; i++)
	{
		xout[i] = xin[0] + (float) i * dx;
	}
}

/* 8系数一维插值(浮点版本)  */
void Interpolation::ISL_intt8r(int ntable, float table[][8], int nxin,
		float dxin, float fxin, float yin[], float yinl, float yinr, int nxout,
		float xout[], float yout[])
{
	int ioutb, nxinm8, ixout, ixoutn, kyin, ktable, itable;
	float xoutb, xoutf, xouts, xoutn, frac, fntablem1, yini, sum, *yin0,
			*table00, *pyin, *ptable;

	/* compute constants */
	ioutb = -3 - 8;
	xoutf = fxin;
	xouts = 1.0 / dxin;
	xoutb = 8.0 - xoutf * xouts;
	fntablem1 = (float) (ntable - 1);
	nxinm8 = nxin - 8;
	yin0 = &yin[0];
	table00 = &table[0][0];

	/* loop over output samples */
	for (ixout = 0; ixout < nxout; ixout++)
	{

		/* determine pointers into table and yin */
		xoutn = xoutb + xout[ixout] * xouts;
		ixoutn = (int) xoutn;
		kyin = ioutb + ixoutn;
		pyin = yin0 + kyin;
		frac = xoutn - (float) ixoutn;
		ktable = frac >= 0.0 ? frac * fntablem1 + 0.5 : (frac + 1.0)
				* fntablem1 - 0.5;
		ptable = table00 + ktable * 8;

		/* if totally within input array, use fast method */
		if (kyin >= 0 && kyin <= nxinm8)
		{
			yout[ixout] = pyin[0] * ptable[0] + pyin[1] * ptable[1] + pyin[2]
					* ptable[2] + pyin[3] * ptable[3] + pyin[4] * ptable[4]
					+ pyin[5] * ptable[5] + pyin[6] * ptable[6] + pyin[7]
					* ptable[7];

			/* else handle end effects with care */
		}
		else
		{

			/* sum over 8 tabulated coefficients */
			for (itable = 0, sum = 0.0; itable < 8; itable++, kyin++)
			{
				if (kyin < 0)
					yini = yinl;
				else if (kyin >= nxin)
					yini = yinr;
				else
					yini = yin[kyin];
				sum += yini * (*ptable++);
			}
			yout[ixout] = sum;
		}
	}
}


/* 计算连续的二阶导数的三次样条插值系数, 	和ISL_intcub配合使用 */
void Interpolation::ISL_csplin ( int n, float x[], float y[], float yd[][4] )
{
	int i;
	float h1,h2,del1,del2,dmax,hsum,w1,w2,divdf3,sleft,sright,alpha,t;

	/* if n=1, then use constant interpolation */
	if (n==1) {
		yd[0][0] = y[0];
		yd[0][1] = 0.0;
		yd[0][2] = 0.0;
		yd[0][3] = 0.0;
		return;

	/* else, if n=2, then use linear interpolation */
	} else if (n==2) {
		yd[0][0] = y[0];  yd[1][0] = y[1];
		yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
		yd[0][2] = yd[1][2] = 0.0;
		yd[0][3] = yd[1][3] = 0.0;
		return;
	}

	/* set left end derivative via shape-preserving 3-point formula */
	h1 = x[1]-x[0];
	h2 = x[2]-x[1];
	hsum = h1+h2;
	del1 = (y[1]-y[0])/h1;
	del2 = (y[2]-y[1])/h2;
	w1 = (h1+hsum)/hsum;
	w2 = -h1/hsum;
	sleft = w1*del1+w2*del2;
	if (sleft*del1<=0.0)
		sleft = 0.0;
	else if (del1*del2<0.0) {
		dmax = 3.0*del1;
		if (ABS(sleft)>ABS(dmax)) sleft = dmax;
	}

	/* set right end derivative via shape-preserving 3-point formula */
	h1 = x[n-2]-x[n-3];
	h2 = x[n-1]-x[n-2];
	hsum = h1+h2;
	del1 = (y[n-2]-y[n-3])/h1;
	del2 = (y[n-1]-y[n-2])/h2;
	w1 = -h2/hsum;
	w2 = (h2+hsum)/hsum;
	sright = w1*del1+w2*del2;
	if (sright*del2<=0.0)
		sright = 0.0;
	else if (del1*del2<0.0) {
		dmax = 3.0*del2;
		if (ABS(sright)>ABS(dmax)) sright = dmax;
	}

	/* compute tridiagonal system coefficients and right-hand-side */
	yd[0][0] = 1.0;
	yd[0][2] = 2.0*sleft;
	for (i=1; i<n-1; i++) {
		h1 = x[i]-x[i-1];
		h2 = x[i+1]-x[i];
		del1 = (y[i]-y[i-1])/h1;
		del2 = (y[i+1]-y[i])/h2;
		alpha = h2/(h1+h2);
		yd[i][0] = alpha;
		yd[i][2] = 3.0*(alpha*del1+(1.0-alpha)*del2);
	}
	yd[n-1][0] = 0.0;
	yd[n-1][2] = 2.0*sright;

	/* solve tridiagonal system for slopes */
	t = 2.0;
	yd[0][1] = yd[0][2]/t;
	for (i=1; i<n; i++) {
		yd[i][3] = (1.0-yd[i-1][0])/t;
		t = 2.0-yd[i][0]*yd[i][3];
		yd[i][1] = (yd[i][2]-yd[i][0]*yd[i-1][1])/t;
	}
	for (i=n-2; i>=0; i--)
		yd[i][1] -= yd[i+1][3]*yd[i+1][1];

	/* copy ordinates into output array */
	for (i=0; i<n; i++)
		yd[i][0] = y[i];

	/* compute 2nd and 3rd derivatives of cubic polynomials */
	for (i=0; i<n-1; i++) {
		h2 = x[i+1]-x[i];
		del2 = (y[i+1]-y[i])/h2;
		divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
		yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
		yd[i][3] = (divdf3/h2)*(6.0/h2);
	}
	yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
	yd[n-1][3] = yd[n-2][3];
}


/* 单点线性插值  */
float Interpolation::ISL_line(float *x, float *y, int n, float t)
{
	static int idx;

	ISL_xindex(n, x, t, &idx);
	return y[idx] + (t - x[idx]) * (y[idx + 1] - y[idx])
			/ (x[idx + 1] - x[idx]);

}

/* 一元三点插值  */
float Interpolation::ISL_lg3(float *x, float *y, int n, float t)
{
	int i, j, k, m;
	float z, s;
	z = 0.0;
	if (n < 1)
		return (z);
	if (n == 1)
	{
		z = y[0];
		return (z);
	}
	if (n == 2)
	{
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return (z);
	}
	if (t <= x[1])
	{
		k = 0;
		m = 2;
	}
	else if (t >= x[n - 2])
	{
		k = n - 3;
		m = n - 1;
	}
	else
	{
		k = 1;
		m = n;
		while (m - k != 1)
		{
			i = (k + m) / 2;
			if (t < x[i - 1])
				m = i;
			else
				k = i;
		}
		k = k - 1;
		m = m - 1;
		if (fabs(t - x[k]) < fabs(t - x[m]))
			k = k - 1;
		else
			m = m + 1;
	}
	z = 0.0;
	for (i = k; i <= m; i++)
	{
		s = 1.0;
		for (j = k; j <= m; j++)
			if (j != i)
				s = s * (t - x[j]) / (x[i] - x[j]);
		z = z + s * y[i];
	}
	return (z);
}

/* 拉格朗日（Lagrange）插值  */
float Interpolation::ISL_lgr(float *x, float *y, int n, float t)
{
	int i, j, k, m;
	float z, s;

	z = 0.0;
	if (n < 1)
		return (z);

	if (n == 1)
	{
		z = y[0];
		return (z);
	}

	if (n == 2)
	{
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return (z);
	}

	i = 0;
	while ((x[i] < t) && (i < n))
		i = i + 1;
	k = i - 4;

	if (k < 0)
		k = 0;
	m = i + 3;

	if (m > n - 1)
		m = n - 1;

	for (i = k; i <= m; i++)
	{
		s = 1.0;
		for (j = k; j <= m; j++)
			if (j != i)
				s = s * (t - x[j]) / (x[i] - x[j]);
		z = z + s * y[i];
	}

	return (z);
}

/* 连分式插值  */
float Interpolation::ISL_pqs(float *x, float *y, int n, float t)
{
	int i, j, k, m, l;
	float z, h, b[8];
	z = 0.0;

	if (n < 1)
		return (z);
	if (n == 1)
	{
		z = y[0];
		return (z);
	}

	if (n <= 8)
	{
		k = 0;
		m = n;
	}
	else if (t < x[4])
	{
		k = 0;
		m = 8;
	}
	else if (t > x[n - 5])
	{
		k = n - 8;
		m = 8;
	}
	else
	{
		k = 1;
		j = n;
		while (j - k != 1)
		{
			i = (k + j) / 2;
			if (t < x[i - 1])
				j = i;
			else
				k = i;
		}
		k = k - 4;
		m = 8;
	}

	b[0] = y[k];
	for (i = 2; i <= m; i++)
	{
		h = y[i + k - 1];
		l = 0;
		j = 1;
		while ((l == 0) && (j <= i - 1))
		{
			if (fabs(h - b[j - 1]) + 1.0 == 1.0)
				l = 1;
			else
				h = (x[i + k - 1] - x[j + k - 1]) / (h - b[j - 1]);
			j = j + 1;
		}
		b[i - 1] = h;
		if (l != 0)
			b[i - 1] = 1.0e+35;
	}

	z = b[m - 1];
	for (i = m - 1; i >= 1; i--)
		z = b[i - 1] + (t - x[i + k - 1]) / z;

	return (z);
}

/*埃特金插值*/
float Interpolation::ISL_atk(float *x, float *y, int n, float t, float eps)
{
	int i, j, k, m, l = 0;
	float z, xx[10], yy[10];
	z = 0.0;

	if (n < 1)
		return (z);

	if (n == 1)
	{
		z = y[0];
		return (z);
	}

	m = 10;
	if (m > n)
		m = n;

	if (t <= x[0])
		k = 1;
	else if (t >= x[n - 1])
		k = n;
	else
	{
		k = 1;
		j = n;
		while ((k - j != 1) && (k - j != -1))
		{
			l = (k + j) / 2;
			if (t < x[l - 1])
				j = l;
			else
				k = l;
		}
		if (fabs(t - x[l - 1]) > fabs(t - x[j - 1]))
			k = j;
	}

	j = 1;
	l = 0;
	for (i = 1; i <= m; i++)
	{
		k = k + j * l;
		if ((k < 1) || (k > n))
		{
			l = l + 1;
			j = -j;
			k = k + j * l;
		}
		xx[i - 1] = x[k - 1];
		yy[i - 1] = y[k - 1];
		l = l + 1;
		j = -j;
	}

	i = 0;
	do
	{
		i = i + 1;
		z = yy[i];
		for (j = 0; j <= i - 1; j++)
			z = yy[j] + (t - xx[j]) * (yy[j] - z) / (xx[j] - xx[i]);
		yy[i] = z;
	} while ((i != m - 1) && (fabs(yy[i] - yy[i - 1]) > eps));

	return (z);
}

/* 埃尔米特插值  */
float Interpolation::ISL_hmt(float *x, float *y, float *dy, int n, float t)
{
	int i, j;
	float z, p, q, s;
	z = 0.0;

	for (i = 1; i <= n; i++)
	{
		s = 1.0;
		for (j = 1; j <= n; j++)
			if (j != i)
				s = s * (t - x[j - 1]) / (x[i - 1] - x[j - 1]);

		s = s * s;
		p = 0.0;

		for (j = 1; j <= n; j++)
			if (j != i)
				p = p + 1.0 / (x[i - 1] - x[j - 1]);

		q = y[i - 1] + (t - x[i - 1]) * (dy[i - 1] - 2.0 * y[i - 1] * p);

		z = z + q * s;
	}
	return (z);
}

/* 一维插值( 对序列y(x[0]), y(x[1])进行插值 ) */
void Interpolation::ISL_intlin(int nin, float xin[], float yin[], float yinl,
		float yinr, int nout, float xout[], float yout[], int flag, float eps)
{
	//	static int idx;
	int jout;
	float x;

	/* if input x values are monotonically increasing, then */
	if (xin[0] <= xin[nin - 1])
	{
		for (jout = 0; jout < nout; jout++)
		{
			x = xout[jout];
			if (x < xin[0])
				yout[jout] = yinl;
			else if (x > xin[nin - 1])
				yout[jout] = yinr;
			else if (x == xin[nin - 1] || nin == 1)
				yout[jout] = yin[nin - 1];
			else
			{

				switch (flag)
				{
				case 1:
					yout[jout] = Interpolation::ISL_line(xin, yin, nin, x);
					break;

				case 2:
					yout[jout] = Interpolation::ISL_lg3(xin, yin, nin, x);
					break;

				case 3:
					yout[jout] = Interpolation::ISL_lgr(xin, yin, nin, x);
					break;

				case 4:
					yout[jout] = Interpolation::ISL_pqs(xin, yin, nin, x);
					break;

				case 5:
					yout[jout] = Interpolation::ISL_atk(xin, yin, nin, x, eps);
					break;

				default:
					yout[jout] = Interpolation::ISL_line(xin, yin, nin, x);
					break;
				}

				//				ISL_xindex(nin,xin,x,&idx);
				//				yout[jout] = yin[idx]+(x-xin[idx])
				//					*(yin[idx+1]-yin[idx])
				//					/(xin[idx+1]-xin[idx]);
			}
		}

		/* else, if input x values are monotonically decreasing, then */
	}
	else
	{
		for (jout = 0; jout < nout; jout++)
		{
			x = xout[jout];
			if (x > xin[0])
				yout[jout] = yinl;
			else if (x < xin[nin - 1])
				yout[jout] = yinr;
			else if (x == xin[nin - 1] || nin == 1)
				yout[jout] = yin[nin - 1];
			else
			{

				switch (flag)
				{
				case 1:
					yout[jout] = Interpolation::ISL_line(xin, yin, nin, x);
					break;

				case 2:
					yout[jout] = Interpolation::ISL_lg3(xin, yin, nin, x);
					break;

				case 3:
					yout[jout] = Interpolation::ISL_lgr(xin, yin, nin, x);
					break;

				case 4:
					yout[jout] = Interpolation::ISL_pqs(xin, yin, nin, x);
					break;

				case 5:
					yout[jout] = Interpolation::ISL_atk(xin, yin, nin, x, eps);
					break;

				default:
					yout[jout] = Interpolation::ISL_line(xin, yin, nin, x);
					break;
				}

				//				ISL_xindex(nin,xin,x,&idx);
				//				yout[jout] = yin[idx]+(x-xin[idx])
				//					*(yin[idx+1]-yin[idx])
				//					/(xin[idx+1]-xin[idx]);
			}
		}
	}
}

/* 正弦一维插值(浮点版本)  */
void Interpolation::ISL_ints8r(int nxin, float dxin, float fxin, float yin[],
		float yinl, float yinr, int nxout, float xout[], float yout[])
{
	static float table[NTABLE2][LTABLE];
	static int tabled = 0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled)
	{
		for (jtable = 1; jtable < NTABLE2 - 1; jtable++)
		{
			frac = (float) jtable / (float) (NTABLE2 - 1);
			ISL_mksinc(frac, LTABLE, &table[jtable][0]);
		}
		for (jtable = 0; jtable < LTABLE; jtable++)
		{
			table[0][jtable] = 0.0;
			table[NTABLE2 - 1][jtable] = 0.0;
		}
		table[0][LTABLE / 2 - 1] = 1.0;
		table[NTABLE2 - 1][LTABLE / 2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	ISL_intt8r(NTABLE2, table, nxin, dxin, fxin, yin, yinl, yinr, nxout, xout,
			yout);
}

/* y(x), y'(x), y''(x) , 分段插值 */
void Interpolation::ISL_intcub (	int ideriv,
					int nin, float xin[], float ydin[][4],
					int nout, float xout[], float yout[] )
{
	static int idx;
	int iout;
	float delx;

	/* y(x) is desired, then */
	if (ideriv==0) {
		for (iout=0; iout<nout; iout++) {
			ISL_xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][0]+delx*
				(ydin[idx][1]+delx*
				(ydin[idx][2]*O2+delx*
				(ydin[idx][3]*O6))));
		}

	/* else, if y'(x) is desired, then */
	} else if (ideriv==1) {
		for (iout=0; iout<nout; iout++) {
			ISL_xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][1]+delx*
				(ydin[idx][2]+delx*
				(ydin[idx][3]*O2)));
		}

	/* else, if y''(x) is desired, then */
	} else if (ideriv==2) {
		for (iout=0; iout<nout; iout++) {
			ISL_xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][2]+delx*
				(ydin[idx][3]));
		}

	/* else, if y'''(x) is desired, then */
	} else if (ideriv==3) {
		for (iout=0; iout<nout; iout++) {
			ISL_xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][3]);
		}

	/* else, if any other derivative is desired, then */
	} else {
		for (iout=0; iout<nout; iout++)
			yout[iout] = 0.0;
	}
}

}

