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
    uint mnID;     //境界ID(BoundaryID)
    uint mnBndType;//境界種類(Dirichlet | Neumann) :=> BNodeへの配分が異なる.

    uint mMGLevel;   //自身の階層Level
    uint mMaxMGLevel;//最大Level(BNode境界条件の領域確保に使用)

    vector<CBoundaryNode*> mvBNode;
    map<uint, uint, less<uint> > mmBNodeID2Index;

    vuint  mvDOF;                            //DOF番号の配列
    map<uint, uint, less<uint> > mmDOF2Index;//DOF番号 => DOFインデックス

public:
    //境界ID(BoundaryID)
    void setID(const uint& id){ mnID= id;}
    uint& getID(){ return mnID;}

    //DOFの追加・代入・配列領域
    void addDOF(const uint& dof);
    void setDOF(const uint& index, const uint& dof);
    void resizeDOF(const uint& res_size);
    //DOFの提供・DOFインデックス・DOF数
    uint& getDOF(const uint& index);
    uint& getDOF_Index(const uint& dof);
    uint getNumOfDOF();

    //境界種類(Dirichlet | Neumann)
    void setBndType(const uint& bndType){ mnBndType= bndType;}
    uint& getBndType(){ return mnBndType;}

    //階層Level
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}

    //最大Level(BNode境界条件の領域確保に使用)
    void setMaxMGLevel(const uint& maxLevel){ mMaxMGLevel= maxLevel;}
    uint& getMaxMGLevel(){ return mMaxMGLevel;}



    //境界Node :=> BNode
    void resizeBNode(const uint& res_size){ mvBNode.resize(res_size);}
    void setBNode(const uint& index, CBoundaryNode *pBNode);
    void addBNode(CBoundaryNode* pBNode);
    uint getNumOfBNode(){ return mvBNode.size();}
    CBoundaryNode* getBNodeIX(const uint& index){ return mvBNode[index];}
    CBoundaryNode* getBNodeID(const uint& id){
        uint index= mmBNodeID2Index[id];
        return mvBNode[index];
    }

protected:
    virtual void distNeumannValue()=0;
    void distDirichletValue();
    virtual void distDirichletValue_at_CGrid()=0;//Level==0(coarse grid)
    virtual void distDirichletValue_at_FGrid()=0;//Level>=1(fine grid)

public:
    void distValueBNode();//distNeumannValue, distDirichletValue呼び出し
};
#endif	/* _BOUNDARYMESH_H */
}


