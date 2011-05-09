/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryFace.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _BOUNDARYFACE_H_c5c754e6_c2be_4199_854f_2929f2afa847
#define	_BOUNDARYFACE_H_c5c754e6_c2be_4199_854f_2929f2afa847
#include "GeneralBoundary.h"
#include "Logger.h"
namespace pmw{
class CBoundaryFace:public CGeneralBoundary{
public:
    CBoundaryFace();
    virtual ~CBoundaryFace();
private:
    uint mnID;
    uint mnFaceID;
public:
    void initialize(const uint& bndType,const uint& dof);
    void setElementID(const uint& id){ mnID = id;}
    uint& getElementID(){ return mnID;}
    void setFaceID(const uint& id){ mnFaceID= id;}
    uint& getFaceID(){ return mnFaceID;}
};
}
#endif	/* _BOUNDARYFACE_H */
