/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.1 beta
|
|   ./src/CommonStd.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef COMMON_STD_HH_3AF0ABDA_F467_454f_B9DB_A487484A060F
#define COMMON_STD_HH_3AF0ABDA_F467_454f_B9DB_A487484A060F

#define BUFFERLENGTH 256

#ifndef MSVC
#include <stdint.h>
#endif

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "DeleteObject.h"
#ifdef  MSVC
 #define _USE_MATH_DEFINES
#endif
#include <cmath>
using namespace std;
#endif 
