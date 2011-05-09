/* 
 * File:   BoundaryFace.h
 * Author: ktakeda
 *
 * Created on 2009/05/18, 15:11
 */

#ifndef _BOUNDARYFACE_H_c5c754e6_c2be_4199_854f_2929f2afa847
#define	_BOUNDARYFACE_H_c5c754e6_c2be_4199_854f_2929f2afa847

//#include "Polyhedral.h"
#include "GeneralBoundary.h"

#include "Logger.h"

namespace pmw{
//class CBoundaryFace:public CPolyhedral,public CGeneralBoundary{
class CBoundaryFace:public CGeneralBoundary{
public:
    CBoundaryFace();
    virtual ~CBoundaryFace();

private:
    uint mnID;// Meshの要素Index番号
    uint mnFaceID;//要素のローカル面番号

    //圧力(垂直荷重)
    //面方向の面力(粘性:Viscous Force)
    //熱流束(面)
public:
    void initialize(const uint& bndType,const uint& dof);

    // Element ID
    void setElementID(const uint& id){ mnID = id;}
    uint& getElementID(){ return mnID;}

    // Face ID
    void setFaceID(const uint& id){ mnFaceID= id;}
    uint& getFaceID(){ return mnFaceID;}
};
}

#endif	/* _BOUNDARYFACE_H */

