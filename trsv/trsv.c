#include <stdio.h>
#include <stdlib.h>
#include "inc/nc.h"



int main(int argc, char *argv[]) 
{
	int n;
	float *A;
	int i, j;


	return 0;
}


/* Function: trsvLNN
 * UPLO:  U
 * TRANS: N
 * DIAG:  N
 **********************/
	static void my_strsvLNN
(
 const int                  N,
 const float               * A,
 const int                  LDA,
 float                     * X,
 const int                  INCX
 )
{
	int                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
	register float            t0;

	for( j = 0, jaj = 0, jx  = 0; j < N; j++, jaj += ldap1, jx += INCX )
	{
		X[jx] /= A[jaj]; t0 = X[jx];
		for( i = j+1,    iaij  = jaj+1, ix  = jx + INCX;
				i < N; i++, iaij += 1,     ix += INCX ) {
			X[ix] -= t0 * A[iaij];
		}
	}
	return;
}


/* Function: trsvLNU
 * UPLO:  U
 * TRANS: N
 * DIAG:  U
 **********************/
	static void my_strsvLNU
(
 const int                  N,
 const float               * A,
 const int                  LDA,
 float                     * X,
 const int                  INCX
 )
{
	int                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
	register float            t0;

	for( j = 0, jaj = 0, jx  = 0; j < N; j++, jaj += ldap1, jx += INCX )
	{
		t0 = X[jx];
		for( i = j+1,    iaij  = jaj+1, ix  = jx + INCX;
				i < N; i++, iaij += 1,     ix += INCX ) {
			X[ix] -= t0 * A[iaij];
		}
	}
	return;
}

/* Function: trsvLTU
 * UPLO:  L
 * TRANS: T
 * DIAG:  N
 **********************/
static void my_strsvLTN
(
   const int                  N,
   const float               * A,
   const int                  LDA,
   float                     * X,
   const int                  INCX
)
{
   int                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register float            t0;

   for( j = N-1,     jaj  = (N-1)*(ldap1), jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= ldap1,         jx -= INCX )
   {
      t0 = X[jx];
      for( i = j+1,    iaij  = 1+jaj, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { t0 -= A[iaij] * X[ix]; }
      t0 /= A[jaj]; X[jx] = t0;
   }
}

/* Function: trsvLTU
 * UPLO:  L
 * TRANS: T
 * DIAG:  U
 **********************/
static void my_strsvLTU
(
   const int                  N,
   const float               * A,
   const int                  LDA,
   float                     * X,
   const int                  INCX
)
{
   int                        i, iaij, ix, j, jaj, jx, ldap1 = LDA + 1;
   register float            t0;

   for( j = N-1,     jaj  = (N-1)*(ldap1), jx  = (N-1)*INCX;
        j >= 0; j--, jaj -= ldap1,         jx -= INCX )
   {
      t0 = X[jx];
      for( i = j+1,    iaij  = 1+jaj, ix  = jx + INCX;
           i < N; i++, iaij += 1,     ix += INCX ) { t0 -= A[iaij] * X[ix]; }
      X[jx] = t0;
   }
}

/* Function: trsvUNN
 * UPLO:  U
 * TRANS: N
 * DIAG:  N
 **********************/
	static void my_strsvUNN
(
 const int                  N,
 const float               * A,
 const int                  LDA,
 float                     * X,
 const int                  INCX
 )
{
	int                        i, iaij, ix, j, jaj, jx;
	register float            t0;

	for( j = N-1,     jaj  = (N-1)*LDA, jx  = (N-1)*INCX;
			j >= 0; j--, jaj -= LDA,       jx -= INCX )
	{
		t0 = X[jx];
		for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
		{ X[ix] -= t0 * A[iaij]; }
	}

	return;
}


/* Function: trsvUNU
 * UPLO:  U
 * TRANS: N
 * DIAG:  U
 **********************/
	static void my_strsvUNU
(
 const int                  N,
 const float               * A,
 const int                  LDA,
 float                     * X,
 const int                  INCX
 )
{
	int                        i, iaij, ix, j, jaj, jx;
	register float            t0;

	for( j = N-1,     jaj  = (N-1)*LDA, jx  = (N-1)*INCX;
			j >= 0; j--, jaj -= LDA,       jx -= INCX )
	{
		X[jx] /= A[j+jaj]; t0 = X[jx];
		for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
		{ X[ix] -= t0 * A[iaij]; }
	}

	return;
}

/* Function: trsvUTN
 * UPLO:  U
 * TRANS: T
 * DIAG:  N
 **********************/
static void my_strsvUTN
(
   const int                  N,
   const float               * A,
   const int                  LDA,
   float                     * X,
   const int                  INCX
)
{
   int                        i, iaij, ix, j, jaj, jx;
   register float            t0;

   for( j = 0, jaj = 0,jx = 0; j < N; j++, jaj += LDA, jx += INCX )
   {
      t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { t0 -= A[iaij] * X[ix]; }
      t0 /= A[iaij]; X[jx] = t0;
   }
}

/* Function: trsvUTU
 * UPLO:  U
 * TRANS: T
 * DIAG:  U
 **********************/
static void my_strsvUTU
(
   const int                  N,
   const float               * A,
   const int                  LDA,
   float                     * X,
   const int                  INCX
)
{
   int                        i, iaij, ix, j, jaj, jx;
   register float            t0;

   for( j = 0, jaj = 0, jx = 0; j < N; j++, jaj += LDA, jx += INCX )
   {
      t0 = X[jx];
      for( i = 0, iaij = jaj, ix = 0; i < j; i++, iaij += 1, ix += INCX )
      { t0 -= A[iaij] * X[ix]; }
      X[jx] = t0;
   }
}


	static void my_strsv0
(
 const enum MY_UPLO        UPLO,
 const enum MY_TRANS       TRANS,
 const enum MY_DIAG        DIAG,
 const int                  N,
 const float               * A,
 const int                  LDA,
 float                     * X,
 const int                  INCX
 ) 
{
	if( N == 0 ) return;

	if( UPLO == MyUpper )
	{
		if( TRANS == MyNoTrans )
		{
			if( DIAG == MyNonUnit ) { my_strsvUNN( N,    A, LDA, X, INCX ); }
			else                     { my_strsvUNU( N,    A, LDA, X, INCX ); }
		}
		else
		{
			if( DIAG == MyNonUnit ) { my_strsvUTN( N,    A, LDA, X, INCX ); }
			else                     { my_strsvUTU( N,    A, LDA, X, INCX ); }
		}
	}
	else
	{
		if( TRANS == MyNoTrans )
		{
			if( DIAG == MyNonUnit ) { my_strsvLNN( N,    A, LDA, X, INCX ); }
			else                     { my_strsvLNU( N,    A, LDA, X, INCX ); }
		}
		else
		{
			if( DIAG == MyNonUnit ) { my_strsvLTN( N,    A, LDA, X, INCX ); }
			else                     { my_strsvLTU( N,    A, LDA, X, INCX ); }
		}
	}
}


/* 
 * Purpose
 * =======
 *
 * my_dtrsv solves one of the systems of equations
 *  
 *     A * x = b,   or   A^T * x = b,
 *  
 * where b and x are n-element vectors and  A  is an n by n non-unit, or
 * unit, upper or lower triangular matrix.
 *  
 * No test for  singularity  or  near-singularity  is included  in  this
 * routine. Such tests must be performed before calling this routine.
 *
 * Arguments
 * =========
 *
 * ORDER   (local input)                 const enum MY_ORDER
 *         On entry, ORDER  specifies the storage format of the operands
 *         as follows:                                                  
 *            ORDER = MyRowMajor,                                      
 *            ORDER = MyColumnMajor.                                   
 *
 * UPLO    (local input)                 const enum MY_UPLO
 *         On  entry,   UPLO   specifies  whether  the  upper  or  lower
 *         triangular  part  of the array  A  is to be referenced.  When
 *         UPLO==MyUpper, only  the upper triangular part of A is to be
 *         referenced, otherwise only the lower triangular part of A is 
 *         to be referenced. 
 *
 * TRANS   (local input)                 const enum MY_TRANS
 *         On entry,  TRANS  specifies  the equations  to  be  solved as
 *         follows:
 *            TRANS==MyNoTrans     A   * x = b,
 *            TRANS==MyTrans       A^T * x = b.
 *
 * DIAG    (local input)                 const enum MY_DIAG
 *         On entry,  DIAG  specifies  whether  A  is unit triangular or
 *         not. When DIAG==MyUnit,  A is assumed to be unit triangular,
 *         and otherwise, A is not assumed to be unit triangular.
 *
 * N       (local input)                 const int
 *         On entry, N specifies the order of the matrix A. N must be at
 *         least zero.
 *
 * A       (local input)                 const float *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than LDA * n. Before entry with  UPLO==MyUpper,  the leading
 *         n by n upper triangular  part of the array A must contain the
 *         upper triangular  matrix and the  strictly  lower  triangular
 *         part of A is not referenced.  When  UPLO==MyLower  on entry,
 *         the  leading n by n lower triangular part of the array A must
 *         contain the lower triangular matrix  and  the  strictly upper
 *         triangular part of A is not referenced.
 *          
 *         Note  that  when  DIAG==MyUnit,  the diagonal elements of  A
 *         not referenced  either,  but are assumed to be unity.
 *
 * LDA     (local input)                 const int
 *         On entry,  LDA  specifies  the  leading  dimension  of  A  as
 *         declared  in  the  calling  (sub) program.  LDA  must  be  at
 *         least MAX(1,n).
 *
 * X       (local input/output)          float *
 *         On entry,  X  is an incremented array of dimension  at  least
 *         ( 1 + ( n - 1 ) * abs( INCX ) )  that  contains the vector x.
 *         Before entry,  the  incremented array  X  must contain  the n
 *         element right-hand side vector b. On exit,  X  is overwritten
 *         with the solution vector x.
 *
 * INCX    (local input)                 const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero.
 *
 * ---------------------------------------------------------------------
 */ 
void my_strsv
(
   const enum MY_ORDER             ORDER,
   const enum MY_UPLO              UPLO,
   const enum MY_TRANS             TRANS,
   const enum MY_DIAG              DIAG,
   const int                        N,
   const float *                   A,
   const int                        LDA,
   float *                         X,
   const int                        INCX
)
{
   if( ORDER == MyColumnMajor )
   {
      my_strsv0( UPLO, TRANS, DIAG, N, A, LDA, X, INCX );
   }
   else
   {
      my_strsv0( ( UPLO  == MyUpper   ? MyLower : MyUpper   ),
                  ( TRANS == MyNoTrans ? MyTrans : MyNoTrans ),
                  DIAG, N, A, LDA, X, INCX );
   }

}

