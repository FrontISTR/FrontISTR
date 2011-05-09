/* 
 * File:   Info.h
 * Author: ktakeda
 *
 * Created on 2009/05/12, 15:05
 */

#ifndef _INFO_H_62c11211_b4a0_4c3a_97c9_ec8f15c9330c
#define	_INFO_H_62c11211_b4a0_4c3a_97c9_ec8f15c9330c

namespace Utility{
union Info{
    static const char* Header(){return
    " --------------------------------------------------------------------------- \n"
    "|                                                                           |\n"
    "|                              HEC_MW3                                      |\n"
    "|                                                                           |\n"
    "|                                                                           |\n"
    "|                          FEM Static Library                               |\n"
    "|                                               Version 0.0                 |\n"
    "|                                                                           |\n"
    "|                                                                           |\n"
    "|                                                               2009.05.12  |\n"
    "|                                                               2009.05.12  |\n"
    " --------------------------------------------------------------------------- ";}
};
}

#endif	/* _INFO_H */

