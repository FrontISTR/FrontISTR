/* 
 * File:   BoundaryVolume.h
 * Author: ktakeda
 *
 * Modify     2009/05/18
 * Created on 2009/05/13, 16:08
 */

#ifndef _BOUNDARYELEMENT_H_af5368fd_f5a0_4058_bd57_36e9b590b973
#define	_BOUNDARYELEMENT_H_af5368fd_f5a0_4058_bd57_36e9b590b973

//#include "Polyhedral.h"
#include "GeneralBoundary.h"

#include "Logger.h"

namespace pmw{
//class CBoundaryVolume:public CPolyhedral,public CGeneralBoundary{
class CBoundaryVolume:public CGeneralBoundary{
public:
    CBoundaryVolume();
    virtual ~CBoundaryVolume();

private:
    uint mnID;// Meshの要素Index番号
//    // fix parameters
//    //重力
//    //加速度
//    //遠心力 回転軸(点1、点2)、角速度 ：パラメータ==7
//    //発熱

public:
    virtual void initialize(const uint& bndType, const uint& dof);

    // Element ID
    void setID(const uint& id){ mnID = id;}
    uint& getID(){ return mnID;}
};
}


#endif	/* _BOUNDARYELEMENT_H */

