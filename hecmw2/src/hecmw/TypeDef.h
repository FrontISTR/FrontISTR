//
//  TypeDef.h
//
//              2011.03.24
//              2008.12.01
//              k.Takeda
#ifndef TYPE_DEF_HH_5B7AEEB4_B0BA_491b_83F7_F5AEFDCFD675
#define TYPE_DEF_HH_5B7AEEB4_B0BA_491b_83F7_F5AEFDCFD675

#include "CommonStd.h"


///*--- MW3 ---*/
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

///*--- File ---*/
#ifdef MSVC
 typedef __int32 int32;
 typedef __int64 int64;
 typedef unsigned __int32 uint32;
 typedef unsigned __int64 uint64;
#else
 typedef  int32_t   int32;
 typedef  int64_t   int64;
 typedef uint32_t  uint32;
 typedef uint64_t  uint64;
#endif

///*--- MAX MIN ---*/
#define IINT32_MAX  2147483647
#define IINT32_MIN -2147483647
#define UINT32_MAX  4294967295
#define UINT32_MIN  0
#define IINT64_MAX  9223372036854775807
#define IINT64_MIN -9223372036854775807
#define UINT64_MAX  18446744073709551615
#define UINT64_MIN  0


#ifdef INT64
 #define IINT_MAX IINT64_MAX
 #define IINT_MIN IINT64_MIN
 #define UIINT_MAX UINT64_MAX
 #define UIINT_MIN UINT64_MIN
#else
 #define IINT_MAX IINT32_MAX
 #define IINT_MIN IINT32_MIN
 #define UIINT_MAX UINT32_MAX
 #define UIINT_MIN UINT32_MIN
#endif

  


///*--- vector ---*/
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

// File Format
#define ASCII  0
#define BINARY 1

#endif

