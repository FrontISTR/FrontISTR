/* 
 * File:   BoundaryNode.h
 * Author: ktakeda
 *
 * Modify     2009/05/25
 * Created on 2009/05/13, 15:04
 */

#ifndef _BOUNDARYNODE_H_37946911_224e_49ea_90c9_472937aab367
#define	_BOUNDARYNODE_H_37946911_224e_49ea_90c9_472937aab367

#include "Vertex.h"
#include "GeneralBoundary.h"

#include "Logger.h"

namespace pmw{
class CBoundaryNode:public CVertex, public CGeneralBoundary{
public:
    CBoundaryNode();
    virtual ~CBoundaryNode();

private:
    //荷重
    //変位
    //速度
    //加速度
    //固定温度(熱応力に利用)
    //集中熱流束
public:
    virtual void initialize(const uint& bndType,const uint& dof);
};
}

#endif	/* _BOUNDARYNODE_H */

