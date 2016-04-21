#ifndef ATLAS_ENUM_H
#define ATLAS_ENUM_H
/*
 * ---------------------------------------------------------------------
 * typedef definitions
 * ---------------------------------------------------------------------
 */
#ifndef ATLAS_ENUM_DEFINED_H
#define ATLAS_ENUM_DEFINED_H
enum ATLAS_ORDER
{  AtlasRowMajor = 101,  AtlasColumnMajor  = 102 };
enum ATLAS_TRANS
{  AtlasNoTrans  = 111,  AtlasTrans        = 112,  AtlasConjTrans    = 113 };
enum ATLAS_UPLO
{  AtlasUpper    = 121,  AtlasLower        = 122 };
enum ATLAS_DIAG
{  AtlasNonUnit  = 131,  AtlasUnit         = 132 };
enum ATLAS_SIDE
{  AtlasLeft     = 141,  AtlasRight        = 142 }; 
#endif

#endif //ATLAS_ENUM_H
