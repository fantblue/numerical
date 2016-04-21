#ifndef ATLAS_BLAS_H
#define ATLAS_BLAS_H
/*
 * =====================================================================
 * Prototypes for Level 1 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefrotg
(
   float *,
   float *,
   float *,
   float *
);

void       ATL_srefrotmg
(
   float *,
   float *,
   float *,
   const float,
   float *
);

float      ATL_srefnrm2
(
   const int,
   const float *,          const int
);

float      ATL_srefasum
(
   const int,
   const float *,          const int
);

int        ATL_isrefamax
(
   const int,
   const float *,          const int
);

void       ATL_srefscal
(
   const int,
   const float,
   float *,                const int
);

void       ATL_srefswap
(
   const int,
   float *,                const int,
   float *,                const int
);

void       ATL_srefcopy
(
   const int,
   const float *,          const int,
   float *,                const int
);

void       ATL_srefaxpy
(
   const int,
   const float,
   const float *,          const int,
   float *,                const int
);

void       ATL_srefrot
(
   const int,
   float *,                const int,
   float *,                const int,
   const float,
   const float
);

void       ATL_srefrotm
(
   const int,
   float *,                const int,
   float *,                const int,
   const float *
);

float      ATL_srefdot
(
   const int,
   const float *,          const int,
   const float *,          const int
);

float      ATL_sdsrefdot
(
   const int,
   const float,
   const float *,          const int,
   const float *,          const int
);

double     ATL_dsrefdot
(
   const int,
   const float *,          const int,
   const float *,          const int
);

void       ATL_drefrotg
(
   double *,
   double *,
   double *,
   double *
);

void       ATL_drefrotmg
(
   double *,
   double *,
   double *,
   const double,
   double *
);

double     ATL_drefnrm2
(
   const int,
   const double *,         const int
);

double     ATL_drefasum
(
   const int,
   const double *,         const int
);

int        ATL_idrefamax
(
   const int,
   const double *,         const int
);

void       ATL_drefscal
(
   const int,
   const double,
   double *,               const int
);

void       ATL_drefswap
(
   const int,
   double *,               const int,
   double *,               const int
);

void       ATL_drefcopy
(
   const int,
   const double *,         const int,
   double *,               const int
);

void       ATL_drefaxpy
(
   const int,
   const double,
   const double *,         const int,
   double *,               const int
);

void       ATL_drefrot
(
   const int,
   double *,               const int,
   double *,               const int,
   const double,
   const double
);

void       ATL_drefrotm
(
   const int,
   double *,               const int,
   double *,               const int,
   const double *
);

double     ATL_drefdot
(
   const int,
   const double *,         const int,
   const double *,         const int
);

void       ATL_crefrotg
(
   float *,
   const float *,
   float *,
   float *
);

float      ATL_screfnrm2
(
   const int,
   const float *,          const int
);

float      ATL_screfasum
(
   const int,
   const float *,          const int
);

int        ATL_icrefamax
(
   const int,
   const float *,          const int
);

void       ATL_crefscal
(
   const int,
   const float *,
   float *,                const int
);

void       ATL_csrefscal
(
   const int,
   const float,
   float *,                const int
);

void       ATL_crefswap
(
   const int,
   float *,                const int,
   float *,                const int
);

void       ATL_crefcopy
(
   const int,
   const float *,          const int,
   float *,                const int
);

void       ATL_crefaxpy
(
   const int,
   const float *,
   const float *,          const int,
   float *,                const int
);

void       ATL_csrefrot
(
   const int,
   float *,                const int,
   float *,                const int,
   const float,
   const float
);

void       ATL_crefdotc_sub
(
   const int,
   const float *,          const int,
   const float *,          const int,
   float *
);

void       ATL_crefdotu_sub
(
   const int,
   const float *,          const int,
   const float *,          const int,
   float *
);

void       ATL_zrefrotg
(
   double *,
   const double *,
   double *,
   double *
);

double     ATL_dzrefnrm2
(
   const int,
   const double *,         const int
);

double     ATL_dzrefasum
(
   const int,
   const double *,         const int
);

int        ATL_izrefamax
(
   const int,
   const double *,         const int
);

void       ATL_zrefscal
(
   const int,
   const double *,
   double *,               const int
);

void       ATL_zdrefscal
(
   const int,
   const double,
   double *,               const int
);

void       ATL_zrefswap
(
   const int,
   double *,               const int,
   double *,               const int
);

void       ATL_zrefcopy
(
   const int,
   const double *,         const int,
   double *,               const int
);

void       ATL_zrefaxpy
(
   const int,
   const double *,
   const double *,         const int,
   double *,               const int
);

void       ATL_zdrefrot
(
   const int,
   double *,               const int,
   double *,               const int,
   const double,
   const double
);

void       ATL_zrefdotc_sub
(
   const int,
   const double *,         const int,
   const double *,         const int,
   double *
);

void       ATL_zrefdotu_sub
(
   const int,
   const double *,         const int,
   const double *,         const int,
   double *
);

/*
 * =====================================================================
 * Prototypes for Level 2 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgbmv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmv
(
  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpr
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefger
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsbmv
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspmv
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspr
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *
);

void       ATL_srefspr2
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *
);

void       ATL_srefsymv
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsyr2
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_sreftpsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_sreftrmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgbmv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmv
(
  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpr
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefger
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsbmv
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspmv
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspr
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *
);

void       ATL_drefspr2
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *
);

void       ATL_drefsymv
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsyr2
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_dreftpsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_dreftrmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgbmv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmv
(
  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgprc
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgpru
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgerc
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgeru
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhbmv
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpmv
(
  const enum ATLAS_UPLO,
  const int,
  const float *,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpr
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *
);

void       ATL_crefhpr2
(
  const enum ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *
);

void       ATL_crefhemv
(
  const enum ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefher
(
  const enum ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefher2
(
  const enum ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_creftpsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_creftrmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgbmv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmv
(
  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemv
(
  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgprc
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgpru
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgerc
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgeru
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhbmv
(
  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpmv
(
  const enum ATLAS_UPLO,
  const int,
  const double *,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpr
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *
);

void       ATL_zrefhpr2
(
  const enum ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *
);

void       ATL_zrefhemv
(
  const enum ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefher
(
  const enum ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefher2
(
  const enum ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_zreftpsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_zreftrmv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsv
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);
/*
 * =====================================================================
 * Prototypes for Level 3 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_sreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_dreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefherk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_creftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefherk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

#endif /*ATLAS_blas.h*/
