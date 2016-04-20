#ifndef __NC_H_
#define __NC_H_
/*
 * ---------------------------------------------------------------------
 * typedef definitions
 * ---------------------------------------------------------------------
 */
enum MY_ORDER
{  MyRowMajor = 101,  MyColumnMajor  = 102 };
enum MY_TRANS
{  MyNoTrans  = 111,  MyTrans        = 112,  MyConjTrans    = 113 };
enum MY_UPLO
{  MyUpper    = 121,  MyLower        = 122 };
enum MY_DIAG
{  MyNonUnit  = 131,  MyUnit         = 132 };
enum MY_SIDE
{  MyLeft     = 141,  MyRight        = 142 }; 

#endif //__NC_H_
