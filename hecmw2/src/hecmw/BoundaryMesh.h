/* 
 * File:   BoundaryMesh.h
 * Author: ktakeda
 *
 * Created on 2010/05/06, 13:57
 */
#include "TypeDef.h"
#include "BoundaryNode.h"
#include <map>

namespace pmw{
#ifndef _BOUNDARYMESH_H
#define	_BOUNDARYMESH_H
class CBoundaryMesh{
public:
    CBoundaryMesh();
    virtual ~CBoundaryMesh();

protected:
    uiint mnID;     //境界ID(BoundaryID)
    uiint mnBndType;//境界種類(Dirichlet | Neumann) :=> BNodeへの配分が異なる.

    string msName; //境界名称

    uiint mMGLevel;   //自身の階層Level
    uiint mMaxMGLevel;//最大Level(BNode境界条件の領域確保に使用)

    vector<CBoundaryNode*> mvBNode;
    map<uiint, uiint, less<uiint> > mmBNodeID2Index;

    vuint  mvDOF;                            //DOF番号の配列
    map<uiint, uiint, less<uiint> > mmDOF2Index;//DOF番号 => DOFインデックス

    ////bool mbSecondOrder;//2次要素の境界メッシュBOOLEAN(mvBNodeのカウントで利用)
    uiint mnEdgeNodeCount;//2次要素が混じった際に、progBMeshのmvBNodeサイズが二重カウントにならないように途中でカウントアップする

public:
    //境界ID(BoundaryID)
    void setID(const uiint& id){ mnID= id;}
    uiint& getID(){ return mnID;}

    //境界名称
    void setName(const string& name){ msName = name;}
    string& getName(){ return msName;}

    //DOFの追加・代入・配列領域
    void addDOF(const uiint& dof);
    void setDOF(const uiint& index, const uiint& dof);
    void resizeDOF(const uiint& res_size);
    //DOFの提供・DOFインデックス・DOF数
    uiint& getDOF(const uiint& index);
    uiint& getDOF_Index(const uiint& dof);
    uiint getNumOfDOF();

    //境界種類(Dirichlet | Neumann)
    void setBndType(const uiint& bndType){ mnBndType= bndType;}
    uiint& getBndType(){ return mnBndType;}

    //階層Level
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}

    //最大Level(BNode境界条件の領域確保に使用)
    void setMaxMGLevel(const uiint& maxLevel){ mMaxMGLevel= maxLevel;}
    uiint& getMaxMGLevel(){ return mMaxMGLevel;}

    //CMWのRefineでコースグリッドのBNode境界条件の領域確保に使用 2011.04.22
    void resizeCGrid_BNodeValue(const uiint& maxLevel);


    //境界Node :=> BNode
    void resizeBNode(const uiint& res_size){ mvBNode.resize(res_size);}
    void setBNode(const uiint& index, CBoundaryNode *pBNode);
    void addBNode(CBoundaryNode* pBNode);
    uiint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNodeIX(const uiint& index){ return mvBNode[index];}
    CBoundaryNode* getBNodeID(const uiint& id){
        uiint index= mmBNodeID2Index[id];
        return mvBNode[index];
    }
    uiint& getBNodeIndex(const uiint& id){ return mmBNodeID2Index[id];}

protected:
    virtual void distNeumannValue()=0;
    virtual void distDirichletValue()=0;
//    virtual void distDirichletValue_at_CGrid()=0;//Level==0(coarse grid)
//    virtual void distDirichletValue_at_FGrid()=0;//Level>=1(fine grid)

public:
    void distValueBNode();//distNeumannValue, distDirichletValue呼び出し
};
#endif	/* _BOUNDARYMESH_H */
}


