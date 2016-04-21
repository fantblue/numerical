#nclude "inc/atlas_enum.h"
#include "inc/atlas_blas.h"

void ATL_sger
(
   const int                  M,
   const int                  N,
   const float                ALPHA,
   const float                * X,
   const int                  INCX,
   const float                * Y,
   const int                  INCY,
   float                      * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * ATL_srefger performs the rank 1 operation
 *
 *    A := alpha * x * y' + A,
 *
 * where alpha is a scalar,  x is an m-element vector, y is an n-element
 * vector and A is an m by n matrix.
 *
 * Arguments
 * =========
 *
 * M       (input)                       const int
 *         On entry,  M  specifies the number of rows of  the matrix  A.
 *         M must be at least zero. Unchanged on exit.
 *
 * N       (input)                       const int
 *         On entry, N  specifies the number of columns of the matrix A.
 *         N  must be at least zero. Unchanged on exit.
 *
 * ALPHA   (input)                       const float
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied as zero then the arrays X and Y need not be set on
 *         input. Unchanged on exit.
 *
 * X       (input)                       const float *
 *         On entry,  X  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( m - 1 ) * abs( INCX ) ) * sizeof(   float   ),
 *         that contains the vector x. Unchanged on exit.
 *
 * INCX    (input)                       const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero. Unchanged on exit.
 *
 * Y       (input)                       const float *
 *         On entry,  Y  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCY ) ) * sizeof(   float   ),
 *         that contains the vector y. Unchanged on exit.
 *
 * INCY    (input)                       const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero. Unchanged on exit.
 *
 * A       (input/output)                float *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than   LDA * n * sizeof(   float   ).  Before entry, the lea-
 *         ding  m by n  part of the array  A  must  contain the  matrix
 *         coefficients.  On exit,  A  is overwritten by the updated ma-
 *         trix.
 *
 * LDA     (input)                       const int
 *         On entry, LDA  specifies the first dimension of A as declared
 *         in the calling (sub) program. LDA must be at least  max(1,m).
 *         Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
	register float t0;
	int            i, iaij, ix, j, jaj, jy;

	if((M==0) || (N==0) || (ALPHA == ATL_sZERO) ) return;

	for( j = 0, jaj = 0, jy = 0; j < N; j++, jaj += LDA, jy+= INCY )
	{
		t0 = ALPHA * Y[jy];
		for( i = 0, iaij  = jaj, ix = 0; i < M; i++, iaij += 1, ix += INCX )
		{ A[iaij] += X[ix] * t0; }
	}
	
}

