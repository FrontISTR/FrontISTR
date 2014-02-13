/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/TypeDef.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef TYPE_DEF_HH_5B7AEEB4_B0BA_491b_83F7_F5AEFDCFD675
#define TYPE_DEF_HH_5B7AEEB4_B0BA_491b_83F7_F5AEFDCFD675
#include "CommonStd.h"

//--
// uiint, iint 定義 (32bit,64bit : プリプロセッサによりサイズ変更)
//--
#ifdef MSVC
#ifdef INT64
typedef unsigned __int64 uiint;
typedef __int64 iint;
#else
typedef unsigned __int32 uiint;
typedef __int32 iint;
#endif
#else
#ifdef INT64
typedef uint64_t uiint;
typedef  int64_t iint;
#else
typedef uint32_t uiint;
typedef  int32_t iint;
#endif
#endif

//--
// int8, int16, int32, int64 定義
//--
#ifdef MSVC
typedef __int8  int8;
typedef __int16 int16;
typedef __int32 int32;
typedef __int64 int64;
typedef unsigned __int8  uint8;
typedef unsigned __int16 uint16;
typedef unsigned __int32 uint32;
typedef unsigned __int64 uint64;
#else
typedef  int8_t    int8;
typedef  int16_t   int16;
typedef  int32_t   int32;
typedef  int64_t   int64;
typedef uint8_t   uint8;
typedef uint16_t  uint16;
typedef uint32_t  uint32;
typedef uint64_t  uint64;
#endif

#define IINT32_MAX  2147483647
#define IINT32_MIN -2147483647
#define UIINT32_MAX  4294967295
#define UIINT32_MIN  0
#define IINT64_MAX  9223372036854775807
#define IINT64_MIN -9223372036854775807
#define UIINT64_MAX  18446744073709551615
#define UIINT64_MIN  0

//--
// IINT_MAX,IINT_MIN,  UIINT_MAX,UIINT_MIN 定義
//--
#ifdef INT64
#define IINT_MAX IINT64_MAX
#define IINT_MIN IINT64_MIN
#define UIINT_MAX UIINT64_MAX
#define UIINT_MIN UIINT64_MIN
#else
#define IINT_MAX IINT32_MAX
#define IINT_MIN IINT32_MIN
#define UIINT_MAX UIINT32_MAX
#define UIINT_MIN UIINT32_MIN
#endif

//--
// vector利用変数 定義
//--
typedef vector<iint> vint;
typedef vector<vector<iint> > vvint;

typedef vector<uiint> vuint;
typedef vector<vector<uiint> > vvuint;
typedef vector<vector<vector<uiint> > > vvvuint;

typedef vector<double> vdouble;
typedef vector<vector<double> > vvdouble;
typedef vector<vector<vector<double> > > vvvdouble;
typedef vector<vector<vector<vector<double> > > > v4double;
typedef vector<vector<vector<vector<vector<double> > > > > v5double;
typedef vector<string> vstring;

//--
// 定数
//--
#define ASCII  0
#define BINARY 1

#define MW_ERROR 0
#define MW_SUCCESS 1

#endif
