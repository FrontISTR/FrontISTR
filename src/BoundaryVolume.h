/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryVolume.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _BOUNDARYELEMENT_H_af5368fd_f5a0_4058_bd57_36e9b590b973
#define	_BOUNDARYELEMENT_H_af5368fd_f5a0_4058_bd57_36e9b590b973
#include "GeneralBoundary.h"
#include "Logger.h"
namespace pmw{
class CBoundaryVolume:public CGeneralBoundary{
public:
    CBoundaryVolume();
    virtual ~CBoundaryVolume();
private:
    uint mnID;
public:
    virtual void initialize(const uint& bndType, const uint& dof);
    void setID(const uint& id){ mnID = id;}
    uint& getID(){ return mnID;}
};
}
#endif	/* _BOUNDARYELEMENT_H */
