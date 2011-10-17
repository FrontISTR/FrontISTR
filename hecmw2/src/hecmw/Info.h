/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Info.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _INFO_H_62c11211_b4a0_4c3a_97c9_ec8f15c9330c
#define	_INFO_H_62c11211_b4a0_4c3a_97c9_ec8f15c9330c
namespace Utility{
struct Info{
    static const char* Header(){return
    " --------------------------------------------------------------------------- \n"
    "|                                                                           |\n"
    "|                     HEC_MW ver4.0   MW3 ver0.2 Beta                       |\n"
    "|                                                                           |\n"
    "|        High-end fundamental middleware for finite element analysis        |\n"
    "|                                                                           |\n"
    "|               Contact address : IIS,The University of Tokyo, CISS         |\n"
    "|                                                                           |\n"
    " --------------------------------------------------------------------------- \n"
    ;}
};
}
#endif	/* _INFO_H */
