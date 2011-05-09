/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   NodeConnectNodeTree.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _NODECONNECTNODETREE_H_209a9301_802f_414d_b5c4_c879ec87e1d1
#define	_NODECONNECTNODETREE_H_209a9301_802f_414d_b5c4_c879ec87e1d1
#include "CommonStd.h"
#include "TypeDef.h"
#include <iostream>
namespace pmw{
class CNodeConnectNodeTree{
private:
    CNodeConnectNodeTree();
public:
    static CNodeConnectNodeTree* Instance(){
        static CNodeConnectNodeTree moNodeConnect;
        return &moNodeConnect;
    }
    virtual ~CNodeConnectNodeTree();
private:
    vvuint mvHexaConnectNode;
    vvuint mvTetraConnectNode;
    vvuint mvPrismConnectNode;
    vvuint mvPyramidConnectNode;
    vvuint mvQuadConnectNode;
    vvuint mvTriangleConnectNode;
    vvuint mvBeamConnectNode;
public:
    vuint& getHexaConnectNode(const uint& localID){ return mvHexaConnectNode[localID];}
    vuint& getTetraConnectNode(const uint& localID){ return mvTetraConnectNode[localID];}
    vuint& getPrismConnectNode(const uint& localID){ return mvPrismConnectNode[localID];}
    vuint& getPyramidConnectNode(const uint& localID){ return mvPyramidConnectNode[localID];}
    vuint& getQuadConnectNode(const uint& localID){ return mvQuadConnectNode[localID];}
    vuint& getTriangleConnectNode(const uint& localID){ return mvTriangleConnectNode[localID];}
    vuint& getBeamConnectNode(const uint& localID){ return mvBeamConnectNode[localID];}
};
}
#endif	/* _NODECONNECTNODETREE_H */
