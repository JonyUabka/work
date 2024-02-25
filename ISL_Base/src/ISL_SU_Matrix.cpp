//***************************************************************************
// **
// **				ISeis Lib [Base Algorithm] - SU中的矩阵计算
// **
// **				Writer By	Mouri Song
// **							Baihong Liu
// **							Qiang Yang
// **
// **				(Soft Center IGP)
// **
// **				DATA : 2014-02-10
// **
// ****************************************************************************/

//#include "ISL_SU_Matrix.h"

//namespace ISLib
//{
///*****************************************************************************
// TRIDIAGONAL - Functions to solve tridiagonal systems of equations Tu=r for u.

// tridif		Solve a tridiagonal system of equations (float version)
// tridid		Solve a tridiagonal system of equations (double version)
// tripd		Solve a positive definite, symmetric tridiagonal system
// tripp		Solve an unsymmetric tridiagonal system that uses
// Gaussian elimination with partial pivoting

// ******************************************************************************
// Function Prototypes:
// void tridif (int n, float a[], float b[], float c[], float r[], float u[]);
// void tridid (int n, double a[], double b[], double c[], double r[], double u[]);
// void tripd (float *d, float *e, float *b, int n);
// void tripp(int n, float *d, float *e, float *c, float *b);

// ******************************************************************************
// tridif, tridid:
// Input:
// n		dimension of system
// a		array[n] of lower sub-diagonal of T (a[0] ignored)
// b		array[n] of diagonal of T
// c		array[n] of upper super-diagonal of T (c[n-1] ignored)
// r		array[n] of right-hand-side column vector

// Output:
// u		array[n] of solution (left-hand-side) column vector

// *****************************************************************************
// tripd:
// Input:
// d	array[n], the diagonal of A
// e	array[n], the superdiagonal of A
// b	array[n], the rhs column vector of Ax=b

// Output:
// b	b is overwritten with the solution to Ax=b
// *****************************************************************************
// tripp:
// Input:
// d	diagonal vector of matrix
// e       upper-diagonal vector of matrix
// c       lower-diagonal vector of matrix
// b       right-hand vector
// n       dimension of matrix

// Output:
// b       solution vector

// ******************************************************************************
// Notes:
// For example, a tridiagonal system of dimension 4 is specified as:

// |b[0]    c[0]     0       0  | |u[0]|     |r[0]|
// |a[1]    b[1]    c[1]     0  | |u[1]|  =  |r[1]|
// | 0      a[2]    b[2]    c[2]| |u[2]|     |r[2]|
// | 0       0      a[3]    b[3]| |u[3]|     |r[3]|

// The tridiagonal matrix is assumed to be non-singular.

// tripd:
// Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
// n-vector b, following algorithm overwrites b with the solution to Ax = b

// ******************************************************************************
// Authors:  tridif, tridid: Dave Hale, Colorado School of Mines, 10/03/89
// tripd, tripp: Zhenyue Liu, Colorado School of Mines,  1992-1993
// *****************************************************************************/
//void ISL_tridif(int n, float a[], float b[], float c[], float r[], float u[])
//{
//	int j;
//	float t, *w;

//	w = (float*) malloc(n * sizeof(float));
//	t = b[0];
//	u[0] = r[0] / t;
//	for (j = 1; j < n; j++)
//	{
//		w[j] = c[j - 1] / t;
//		t = b[j] - a[j] * w[j];
//		u[j] = (r[j] - a[j] * u[j - 1]) / t;
//	}
//	for (j = n - 2; j >= 0; j--)
//		u[j] -= w[j + 1] * u[j + 1];
//	free(w);
//}

//void ISL_tridid(int n, double a[], double b[], double c[], double r[],
//		double u[])
//{
//	int j;
//	double t, *w;

//	w = (double*) malloc(n * sizeof(double));
//	t = b[0];
//	u[0] = r[0] / t;
//	for (j = 1; j < n; j++)
//	{
//		w[j] = c[j - 1] / t;
//		t = b[j] - a[j] * w[j];
//		u[j] = (r[j] - a[j] * u[j - 1]) / t;
//	}
//	for (j = n - 2; j >= 0; j--)
//		u[j] -= w[j + 1] * u[j + 1];
//	free(w);
//}

//void ISL_tripd(float *d, float *e, float *b, int n)
//{
//	int k;
//	float temp;

//	/* decomposition */
//	for (k = 1; k < n; ++k)
//	{
//		temp = e[k - 1];
//		e[k - 1] = temp / d[k - 1];
//		d[k] -= temp * e[k - 1];
//	}

//	/* substitution	*/
//	for (k = 1; k < n; ++k)
//		b[k] -= e[k - 1] * b[k - 1];

//	b[n - 1] /= d[n - 1];
//	for (k = n - 1; k > 0; --k)
//		b[k - 1] = b[k - 1] / d[k - 1] - e[k - 1] * b[k];
//}

//static void exch(float x, float y)
//{
//	float t;
//	t = x;
//	x = y;
//	y = t;
//}

//void ISL_tripp(int n, float *d, float *e, float *c, float *b)
//{
//	int k;
//	float temp;

//	/*      elimination   */
//	for (k = 0; k < n - 1; ++k)
//	{
//		c[k] = 0;
//		if (ABS(d[k]) < ABS(c[k+1]))
//		{
//			exch(d[k], c[k + 1]);
//			exch(e[k], d[k + 1]);
//			exch(c[k], e[k + 1]);
//			exch(b[k], b[k + 1]);
//		}

//		if (d[k] == 0)
//			fprintf(stderr, "coefficient matrix is singular!\n");
//		temp = c[k + 1] / d[k];
//		d[k + 1] -= temp * e[k];
//		e[k + 1] -= temp * c[k];
//		b[k + 1] -= temp * b[k];
//	}

//	/*      substitution      */
//	if (d[n - 1] == 0)
//		fprintf(stderr, "coefficient matrix is singular!\n");
//	b[n - 1] = b[n - 1] / d[n - 1];
//	b[n - 2] = (b[n - 2] - b[n - 1] * e[n - 2]) / d[n - 2];
//	for (k = n - 3; k >= 0; --k)
//		b[k] = (b[k] - b[k + 1] * e[k] - b[k + 2] * c[k]) / d[k];

//}

///*************************************************************************
// AXB - Functions to solve a linear system of equations Ax=b by LU
// decomposition, invert a square matrix or directly multiply an
// inverse matrix by another matrix (without explicitely computing
// the inverse).

// LU_decomposition	Decompose a matrix (A) into a lower triangular (L)
// and an upper triangular (U) such that A=LU

// backward_substitution	Apply backward substitution to an LU decomposed
// matrix to solve the linear system of equations Ax=b

// inverse_matrix		compute the inverse of a square non-singular matrix

// inverse_matrix_multiply	computes the product A^(-1)*B without explicitely
// computing the inverse matrix

// ***************************************************************************
// Function prototypes:
// void LU_decomposition (int nrows, float **matrix, int *idx, float *d);
// void backward_substitution (int nrows, float **matrix, int *idx, float *b);
// void inverse_matrix (int nrows, float **matrix);
// void inverse_matrix_multiply (int nrows1, float **matrix1, int ncols2,
// int nrows2, float **matrix2, float **out_matrix);

// ***************************************************************************
// LU_decomposition:
// Input:
// nrows		number of rows of matrix to invert
// matrix		matrix of coefficients in linear system Ax=b

// Output:
// matrix		matrix containing LU decomposition (original matrix destroyed)
// idx		vector recording the row permutations effected by partial
// pivoting
// d		+/- 1 depending on whether the number of row interchanges
// was even or odd
// *************************************************************************
// backward_substitution
// Input:
// nrows		number of rows (and columns) of input matrix
// matrix		matrix of coefficients (after LU decomposition)
// idx		permutation vector obtained from routine LU_decomposition
// b		right hand side vector in equation Ax=b

// Output:
// b		vector with the solution
// **************************************************************************
// inverse_matrix
// Input:
// nrows		number of rows (and columns) of input matrix
// matrix		matrix to invert

// Output:
// matrix		inverse of input matrix
// **************************************************************************
// inverse_matrix_multiply
// nrows1          number of rows (and columns) of matrix to invert
// matrix1         square matrix to invert
// ncols2          number of coulmns of second matrix
// nrows2			number of rows of second matrix
// matrix          second matrix (multiplicator)

// Output Parameters:
// out_matrix      matrix containing the product of the inverse of the first
// matrix by the second one.
// Note:
// matrix1 and matrix2 are not destroyed, (not clobbered)
// *************************************************************************
// Notes:
// To solve the set of linear equations Ax=b, first do the LU decomposition of
// A (which will clobber A with its LU decomposition) and then do the backward
// substitution with this new matrix and the right-hand side vector b. The vector
// b will be clobbered with the solution. Both, the original matrix and vector B,
// will have been destroyed.

// The LU decomposition is carried out with the Crout's method with implicit
// partial pivoting that guaratees that the maximum pivot is used in every
// step of the algorithm.

// The operation count to solve a linear system of equations via LU decomposition
// is 1/3N^3 and is a factor of 3 better than the standard Gauss-Jordan algorithm
// To invert a matrix the count is the same with both algorithms: N^3.

// Once a linear system Ax=b has been solved, to solve another linear system
// with the same matrix A but with different vetor b, ONLY the back substitution
// has to be repeated with the new b (remember that the matrix in backsubstitution
// is not the original matrix but its LU decomposition)

// If you want to compute A^(-1)*B from matrices A and B, it is better to
// use the subroutine inverse_matrix_multiply rather than explicitely computing
// the inverse. This saves a whole martix multiplication and is also more accurate.

// ***************************************************************************
// Refferences:
// Press, Teukolsky, Vettering and Flannery, Numerical Recipes in C:
// The art of scientific computing. Cambridge University Press.
// second edition. (1992).
// Golub and Van Loan, Matrix Computations. John Hopkins University Press.
// Second Edition. (1989).
// Horn and Johnson, Matrix Analysis. Cambridge University Press. (1985).
// *************************************************************************
// Credits:
// Adapted from discussions in Numerical Recipes, by Gabriel Alvarez (1995)
// *************************************************************************/
//void ISL_LU_decomposition(int nrows, float **matrix, int *idx, float *d)
//{
//	int i, j, k; /* loop counters for rows, columns and summs */
//	int imax = 0; /* index of maximum pivot */
//	float big; /* largest number in input matrix */
//	float dum; /* pivot scale factor */
//	float sum; /* auxiliary variable for summ  */
//	float temp; /* auxiliary variable */
//	float *vv = NULL; /* vector to store implicit scaling for each row */

//	/* allocate working space */
//	vv = alloc1float(nrows);

//	/* initialize interchanges counter */
//	*d = 1.0;

//	/* loop over rows to get implicit scaling information */
//	for (i = 0; i < nrows; i++)
//	{
//		big = 0.0;
//		for (j = 0; j < nrows; j++)
//			if ((temp = ABS(matrix[i][j])) > big)
//				big = temp;
//		if (big == 0.0)
//			fprintf(stderr, "error, singular matrix in LU decomposition\n");

//		/* save the scaling */
//		vv[i] = 1.0 / big;
//	}

//	/* loop over columns (Crout's method) */
//	for (j = 0; j < nrows; j++)
//	{
//		for (i = 0; i < j; i++)
//		{
//			sum = matrix[i][j];
//			for (k = 0; k < i; k++)
//				sum -= matrix[i][k] * matrix[k][j];
//			matrix[i][j] = sum;
//		}

//		/* initialize for the search for largest pivot element */
//		big = 0.0;
//		for (i = j; i < nrows; i++)
//		{
//			sum = matrix[i][j];
//			for (k = 0; k < j; k++)
//				sum -= matrix[i][k] * matrix[k][j];
//			matrix[i][j] = sum;

//			/*  Is new pivot better than best so far? */
//			if ((dum = vv[i] * ABS(sum)) >= big)
//			{
//				big = dum;
//				imax = i;
//			}
//		}

//		/* Do we need to interchange rows */
//		if (j != imax)
//		{
//			for (k = 0; k < nrows; k++)
//			{
//				dum = matrix[imax][k];
//				matrix[imax][k] = matrix[j][k];
//				matrix[j][k] = dum;
//			}

//			/* change the parity of d */
//			*d = -(*d);

//			/* interchange the scale factor */
//			vv[imax] = vv[j];
//		}
//		idx[j] = imax;

//		/* if matrix becomes singular don't use pivot=0 */
//		if (matrix[j][j] == 0.0)
//			matrix[j][j] = TINY;

//		if (j != nrows)
//		{
//			/* divide by the pivot element */
//			dum = 1.0 / matrix[j][j];
//			if (j < nrows)
//			{
//				for (i = j + 1; i < nrows; i++)
//					matrix[i][j] *= dum;
//			}
//		}
//	}

//	/* free workspace */
//	free1float(vv);
//}

//void ISL_backward_substitution(int nrows, float **matrix, int *idx, float *b)
//{
//	int i, ii = -1, j; /* loop counters */
//	int ip; /* index of first nonvanishing element of b */
//	float sum; /* auxiliary variable for partial sums */

//	for (i = 0; i < nrows; i++)
//	{
//		ip = idx[i];

//		/* do forward substitution */
//		sum = b[ip];
//		b[ip] = b[i];
//		if (ii != -1)
//			for (j = ii; j < i; j++)
//				sum -= matrix[i][j] * b[j];
//		else if (sum != 0.0)
//			ii = i;
//		b[i] = sum;
//	}

//	/* now, do the backward substitution */
//	for (i = nrows - 1; i >= 0; i--)
//	{
//		sum = b[i];
//		for (j = i + 1; j < nrows; j++)
//			sum -= matrix[i][j] * b[j];

//		/* store results in output vector */
//		b[i] = sum / matrix[i][i];
//	}
//}

//void ISL_inverse_matrix(int nrows, float **matrix)
//{
//	int i, j; /* loop counters */
//	float d; /* +/-1 depending on row interchanges even/odd*/
//	int *idx = NULL; /* vector of row permutations */
//	float *column = NULL; /* unit vector for backward substitution*/
//	float **inverse = NULL; /* array to hold the inverse matrix */

//	/* allocate working space */
//	idx = alloc1int(nrows);
//	column = alloc1float(nrows);
//	inverse = alloc2float(nrows, nrows);

//	/* first, do the LU decomposition of input matrix */
//	ISL_LU_decomposition(nrows, matrix, idx, &d);

//	/* find inverse by columns */
//	for (j = 0; j < nrows; j++)
//	{

//		/* unit vector corresponding to current column */
//		for (i = 0; i < nrows; i++)
//			column[i] = 0.0;
//		column[j] = 1.0;

//		/* backward substitution column by column */
//		ISL_backward_substitution(nrows, matrix, idx, column);

//		/* compute inverse matrix column by column */
//		for (i = 0; i < nrows; i++)
//			inverse[i][j] = column[i];
//	}

//	/* clobber original matrix with its inverse */
//	for (i = 0; i < nrows; i++)
//		for (j = 0; j < nrows; j++)
//			matrix[i][j] = inverse[i][j];

//	/* free allocated space */
//	free1int(idx);
//	free1float(column);
//	free2float(inverse);
//}

//void ISL_inverse_matrix_multiply(int nrows1, float **matrix1, int ncols2,
//		int nrows2, float **matrix2, float **out_matrix)
//{
//	int i, j; /* loop counters for rows and coulmns */
//	float d; /* to use in LU decomposition */
//	int *idx = NULL; /* to record permutations by partial pivoting*/
//	float **scratch1 = NULL; /* array to hold input matrix1 */
//	float *scratch2 = NULL; /* vector to hold column of input matrix2 */

//	/* allocate working space */
//	idx = alloc1int(nrows1);
//	scratch1 = alloc2float(nrows1, nrows1);
//	scratch2 = alloc1float(nrows2);

//	/* copy input matrix1 to scratch to avoid clobbering */
//	for (i = 0; i < nrows1; i++)
//		for (j = 0; j < nrows1; j++)
//			scratch1[i][j] = matrix1[i][j];

//	/* do the LU decomposition */
//	ISL_LU_decomposition(nrows1, scratch1, idx, &d);

//	/* find inverse by columns */
//	for (j = 0; j < ncols2; j++)
//	{

//		/* copy column of second input matrix to scratch vector */
//		for (i = 0; i < nrows2; i++)
//			scratch2[i] = matrix2[i][j];

//		/* do backward substitution */
//		ISL_backward_substitution(nrows1, scratch1, idx, scratch2);

//		/* copy results to output matrix */
//		for (i = 0; i < nrows1; i++)
//			out_matrix[i][j] = scratch2[i];
//	}

//	/* free allocated space */
//	free2float(scratch1);
//	free1float(scratch2);
//}

///**********************************************************************************************************************
// *
// *  功能：Toeplitz矩阵的线性方程的对称性计算
// *
// *  说明：Solve a symmetric Toeplitz linear system of equations Rf=g for f(Double version)
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					n					In			dimension of system
// *		float				r[]					In			array[n] of top row of Toeplitz matrix
// *		float				g[]					In			array[n] of right-hand-side column vector
// *
// *		float				f[]					Out			array[n] of solution (left-hand-side) column vector
// *		float				a[]					Out			array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
// *
// *  返回：无
// *  Notes:
// *  	This routine does NOT solve the case when the main diagonal is zero, it
// *  	just silently returns.
// *
// *  	The left column of the Toeplitz matrix is assumed to be equal to the top
// *  	row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
// *
// **********************************************************************************************************************/
//void ISL_stoepd(int n, double r[], double g[], double f[], double a[])
//{
//	int i, j;
//	double v, e, c, w, bot;

//	if (r[0] == 0.0)
//		return;

//	a[0] = 1.0;
//	v = r[0];
//	f[0] = g[0] / r[0];

//	for (j = 1; j < n; j++)
//	{

//		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
//		a[j] = 0.0;
//		f[j] = 0.0;
//		for (i = 0, e = 0.0; i < j; i++)
//			e += a[i] * r[j - i];
//		c = e / v;
//		v -= c * e;
//		for (i = 0; i <= j / 2; i++)
//		{
//			bot = a[j - i] - c * a[i];
//			a[i] -= c * a[j - i];
//			a[j - i] = bot;
//		}

//		/* use a and v above to get f[i], i = 0,1,2,...,j */
//		for (i = 0, w = 0.0; i < j; i++)
//			w += f[i] * r[j - i];
//		c = (w - g[j]) / v;
//		for (i = 0; i <= j; i++)
//			f[i] -= c * a[j - i];
//	}
//}

///**********************************************************************************************************************
// *
// *  功能：Toeplitz矩阵的线性方程的对称性计算
// *
// *  说明：Solve a symmetric Toeplitz linear system of equations Rf=g for f(float version)
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					n					In			dimension of system
// *		float				r[]					In			array[n] of top row of Toeplitz matrix
// *		float				g[]					In			array[n] of right-hand-side column vector
// *
// *		float				f[]					Out			array[n] of solution (left-hand-side) column vector
// *		float				a[]					Out			array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
// *
// *  返回：无
// *  Notes:
// *  	This routine does NOT solve the case when the main diagonal is zero, it
// *  	just silently returns.
// *
// *  	The left column of the Toeplitz matrix is assumed to be equal to the top
// *  	row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
// *
// **********************************************************************************************************************/
//void ISL_stoepf(int n, float *r, float *g, float *f, float *a)
//{
//	int i, j;
//	float v, e, c, w, bot;

//	if (r[0] == 0.0)
//		return;

//	a[0] = 1.0;
//	v = r[0];
//	f[0] = g[0] / r[0];

//	for (j = 1; j < n; j++)
//	{

//		a[j] = 0.0;
//		f[j] = 0.0;
//		for (i = 0, e = 0.0; i < j; i++)
//			e += a[i] * r[j - i];
//		c = e / v;
//		v -= c * e;
//		for (i = 0; i <= j / 2; i++)
//		{
//			bot = a[j - i] - c * a[i];
//			a[i] -= c * a[j - i];
//			a[j - i] = bot;
//		}

//		for (i = 0, w = 0.0; i < j; i++)
//			w += f[i] * r[j - i];
//		c = (w - g[j]) / v;
//		for (i = 0; i <= j; i++)
//			f[i] -= c * a[j - i];
//	}
//}

///**********************************************************************************************************************
// *
// *  功能：计算n阶广义Hermite多项式
// *
// *  说明：Compute n-th order generalized Hermite polynomial.
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		double				*h0					In			array of size nt holding H_{n-1}
// *		double				*h1					In			array of size nt holding H_{n}
// *		double 				*t					In			array of size nt holding time vector
// *		double				nt					In			size of arrays, no. of samples
// *		double				n					In			order of polynomial
// *		double				sigma				In			variance

// *		double				*h					Out			array of size nt holding H_{n+1}
// *
// *  返回：无
// *  Note:	Note that n in the function call is the order of the derivative and
// *  		j in the code below is the n in the recurrence relation
// *
// **********************************************************************************************************************/
//void ISL_hermite_n_polynomial(double *h, double *h0, double *h1, double *t,
//		int nt, int n, double sigma)
//{
//	int i; /* loop variable */
//	int j = 1; /* recurrence counter */

//	/* as long as necessary use recurrence relation */
//	do
//	{
//		/* current instance of recurrence relation */
//		for (i = 0; i < nt; ++i)
//			h[i] = (t[i] * h1[i] - j * h0[i]) / sigma;

//		/* update inputs to recurrence relation */
//		memcpy((void *) h0, (const void *) h1, DSIZE * nt);
//		memcpy((void *) h1, (const void *) h, DSIZE * nt);

//		/* update counter */
//		++j;

//	} while (j < n);
//}

///**********************************************************************************************************************
// *
// *  功能：
// *
// *  说明：determine index of x with respect to an array of x values
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					nx					In			number of x values in array ax
// *		int					ax[]				In			array[nx] of monotonically increasing or decreasing x values
// *		float 				x					In			the value for which index is to be determined
// *		int					index				in			index determined previously (used to begin search)

// *		int					index				Out			for monotonically increasing ax values, the largest index
// *															for which ax[index]<=x, except index=0 if ax[0]>x;
// *															for monotonically decreasing ax values, the largest index
// *															for which ax[index]>=x, except index=0 if ax[0]<x
// *
// *  返回：无
// *  	This function is designed to be particularly efficient when called
// *		repeatedly for slightly changing x values; in such cases, the index
// *		returned from one call should be used in the next.
// *
// **********************************************************************************************************************/
//void ISL_xindex(int nx, float ax[], float x, int *index)
//{
//	int lower, upper, middle, step;

//	/* initialize lower and upper indices and step */
//	lower = *index;
//	if (lower < 0)
//		lower = 0;
//	if (lower >= nx)
//		lower = nx - 1;
//	upper = lower + 1;
//	step = 1;

//	/* if x values increasing */
//	if (ax[nx - 1] > ax[0])
//	{

//		/* find indices such that ax[lower] <= x < ax[upper] */
//		while (lower > 0 && ax[lower] > x)
//		{
//			upper = lower;
//			lower -= step;
//			step += step;
//		}
//		if (lower < 0)
//			lower = 0;
//		while (upper < nx && ax[upper] <= x)
//		{
//			lower = upper;
//			upper += step;
//			step += step;
//		}
//		if (upper > nx)
//			upper = nx;

//		/* find index via bisection */
//		while ((middle = (lower + upper) >> 1) != lower)
//		{
//			if (x >= ax[middle])
//				lower = middle;
//			else
//				upper = middle;
//		}

//		/* else, if not increasing */
//	}
//	else
//	{

//		/* find indices such that ax[lower] >= x > ax[upper] */
//		while (lower > 0 && ax[lower] < x)
//		{
//			upper = lower;
//			lower -= step;
//			step += step;
//		}
//		if (lower < 0)
//			lower = 0;
//		while (upper < nx && ax[upper] >= x)
//		{
//			lower = upper;
//			upper += step;
//			step += step;
//		}
//		if (upper > nx)
//			upper = nx;

//		/* find index via bisection */
//		while ((middle = (lower + upper) >> 1) != lower)
//		{
//			if (x <= ax[middle])
//				lower = middle;
//			else
//				upper = middle;
//		}
//	}
//	/* return lower index */
//	*index = lower;
//}

///**********************************************************************************************************************
// *
// *  功能：计算的最小二乘最佳正弦内插系数
// *
// *  说明：The coefficients are a least-squares-best approximation to the ideal
// *  sinc function for frequencies from zero up to a computed maximum
// *  frequency.  For a given interpolator length, lsinc, mksinc computes
// *  the maximum frequency, fmax (expressed as a fraction of the nyquist
// *  frequency), using the following empirically derived relation (from
// *  a Western Geophysical Technical Memorandum by Ken Larner):
// *
// *  	fmax = min(0.066+0.265*log(lsinc),1.0)
// *
// *  Note that fmax increases as lsinc increases, up to a maximum of 1.0.
// *  Use the coefficients to interpolate a uniformly-sampled function y(i)
// *  as follows:
// *
// *			lsinc-1
// *	y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
// *			  j=0
// *
// *	Interpolation error is greatest for d=0.5, but for frequencies less
// *	than fmax, the error should be less than 1.0 percent.
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		float				d					In			fractional distance to interpolation point; 0.0<=d<=1.0
// *		int					lsinc				In			length of sinc approximation; lsinc%2==0 and lsinc<=20
// *
// *		float				sinc[]				Out			array[lsinc] containing interpolation coefficients
// *
// *  返回：无
// *
// **********************************************************************************************************************/
//void ISL_mksinc(float d, int lsinc, float sinc[])
//{
//	int j;
//	double s[20], a[20], c[20], work[20], fmax;

//	/* compute auto-correlation and cross-correlation arrays */
//	fmax = 0.066 + 0.265 * log((double) lsinc);
//	fmax = (fmax < 1.0) ? fmax : 1.0;
//	for (j = 0; j < lsinc; j++)
//	{
//		a[j] = ISL_dsinc(fmax * j);
//		c[j] = ISL_dsinc(fmax * (lsinc / 2 - j - 1 + d));
//	}

//	/* solve symmetric Toeplitz system for the sinc approximation */
//	ISL_stoepd(lsinc, a, c, s, work);
//	for (j = 0; j < lsinc; j++)
//		sinc[j] = s[j];
//}

///**********************************************************************************************************************
// *
// *  功能：一维数组的归一化
// *
// *  说明：
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int 				in[]				In			输入的数组array[n]
// *		int 				n					In			样点个数
// *		float				normv				In			归一值
// *
// *		float				data[]				Out			归一值后的数组array[n]
// *
// *  返回：无
// *
// **********************************************************************************************************************/
//void ISL_normalize_1d(float in[], int n, float normv, float out[])
//{
//	float max = fabs(in[0]);

//	for (int i = 1; i < n; i++)
//	{
//		if (max < fabs(in[i]))
//			max = fabs(in[i]);
//	}

//	for (int i = 0; i < n; i++)
//	{
//		out[i] = in[i] / normv;
//	}
//}

///**********************************************************************************************************************
// *
// *  功能：compute a regularly sampled function x(y) from a regularly
// sampled, monotonically increasing function y(x)
// *
// *  说明：	(1) dx>0.0 && dy>0.0
// *  	(2) y[0] < y[1] < ... < y[nx-1]
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int 				nx					In			number of samples of y(x)
// *		float 				dx					In			x sampling interval; dx>0.0 is required
// *		float				fx					In			first x
// *		float				y					In			array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
// *		float				ny					In			number of samples of x(y)
// *		float				dy					In			y sampling interval; dy>0.0 is required
// *		float				fy					In			first y
// *		float				xylo				In			x value assigned to x(y) when y is less than smallest y(x)
// *		float				xyhi				In			x value assigned to x(y) when y is greater than largest y(x)
// *
// *		float				x					Out			array[ny] of x(y) values
// *
// *  返回：无
// *
// **********************************************************************************************************************/
//void ISL_yxtoxy(int nx, float dx, float fx, float y[], int ny, float dy,
//		float fy, float xylo, float xyhi, float x[])
//{
//	int nxi, nyo, jxi1, jxi2, jyo;
//	float dxi, fxi, dyo, fyo, fyi, yo, xi1, yi1, yi2;

//	nxi = nx;
//	dxi = dx;
//	fxi = fx;
//	nyo = ny;
//	dyo = dy;
//	fyo = fy;
//	fyi = y[0];

//	/* loop over output y less than smallest input y */
//	for (jyo = 0, yo = fyo; jyo < nyo; jyo++, yo += dyo)
//	{
//		if (yo >= fyi)
//			break;
//		x[jyo] = xylo;
//	}

//	/* loop over output y between smallest and largest input y */
//	if (jyo == nyo - 1 && yo == fyi)
//	{
//		x[jyo++] = fxi;
//		yo += dyo;
//	}
//	jxi1 = 0;
//	jxi2 = 1;
//	xi1 = fxi;
//	while (jxi2 < nxi && jyo < nyo)
//	{
//		yi1 = y[jxi1];
//		yi2 = y[jxi2];
//		if (yi1 <= yo && yo <= yi2)
//		{
//			x[jyo++] = xi1 + dxi * (yo - yi1) / (yi2 - yi1);
//			yo += dyo;
//		}
//		else
//		{
//			jxi1++;
//			jxi2++;
//			xi1 += dxi;
//		}
//	}

//	/* loop over output y greater than largest input y */
//	while (jyo < nyo)
//		x[jyo++] = xyhi;
//}

///***************************************
// *
// *
// ****************************************/
//}/*End of ISLib
