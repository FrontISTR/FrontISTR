/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   NodeType.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _NODETYPE_H_0ca01a33_f12a_4ae4_9dde_3da45f62b934
#define	_NODETYPE_H_0ca01a33_f12a_4ae4_9dde_3da45f62b934
namespace pmw{
struct NodeType{
    enum{
        Scalar,
        Vector,
        ScalarVector
    };
};
}
#endif	/* _NODETYPE_H */
