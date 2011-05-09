/* 
 * File:   SolutionType.h
 * Author: ktakeda
 *
 * Created on 2010/08/27, 14:29
 */

namespace pmw{
#ifndef SOLUTIONTYPE_H
#define	SOLUTIONTYPE_H
struct SolutionType{
    enum {
        FEM,
        FVM,
        Limit
    };
};
#endif	/* SOLUTIONTYPE_H */
}


