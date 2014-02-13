/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/NodeConnectNodeTree.h
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
#ifndef _NODECONNECTNODETREE_H_209a9301_802f_414d_b5c4_c879ec87e1d1
#define	_NODECONNECTNODETREE_H_209a9301_802f_414d_b5c4_c879ec87e1d1
#include "CommonStd.h"
#include "TypeDef.h"
#include <iostream>
namespace pmw
{
class CNodeConnectNodeTree
{
private:
    CNodeConnectNodeTree();
public:
    static CNodeConnectNodeTree* Instance() {
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
    vuint& getHexaConnectNode(const uiint& localID) {
        return mvHexaConnectNode[localID];
    }
    vuint& getTetraConnectNode(const uiint& localID) {
        return mvTetraConnectNode[localID];
    }
    vuint& getPrismConnectNode(const uiint& localID) {
        return mvPrismConnectNode[localID];
    }
    vuint& getPyramidConnectNode(const uiint& localID) {
        return mvPyramidConnectNode[localID];
    }
    vuint& getQuadConnectNode(const uiint& localID) {
        return mvQuadConnectNode[localID];
    }
    vuint& getTriangleConnectNode(const uiint& localID) {
        return mvTriangleConnectNode[localID];
    }
    vuint& getBeamConnectNode(const uiint& localID) {
        return mvBeamConnectNode[localID];
    }
};
}
#endif	/* _NODECONNECTNODETREE_H */
