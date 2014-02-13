/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CodeType.h
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
namespace pmw
{
#ifndef CODE_TYPE_H
#define	CODE_TYPE_H
struct CodeType {
    enum {
        Solid=0,
        Flow,
        Limit
    };
};
#endif	/* CODE_TYPE_H */
}

