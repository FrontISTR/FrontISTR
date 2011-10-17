/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/SolutionType.h
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
