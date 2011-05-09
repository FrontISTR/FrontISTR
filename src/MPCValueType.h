/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MPCValueType.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
namespace pmw{
#ifndef _MPCVALUETYPE_H
#define	_MPCVALUETYPE_H    
struct MPCValueType{
    enum{
        Displacement,
        Scalar
    };
};
#endif	/* _MPCVALUETYPE_H */
}
