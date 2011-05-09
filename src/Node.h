//
// Node.h <- Vertex
//
//
//			2008.05.25
//			2008.12.18
//			k.Takeda
#include "TypeDef.h"

#include "NodeType.h"
#include "Logger.h"

#include "Vertex.h"
namespace pmw{
#ifndef NODE_HH_
#define NODE_HH_
class CNode:public CVertex{
public:
    // NodeID==VertexID
    CNode();
    virtual ~CNode();
    
protected:

    //// Arbitrary DOF (任意自由度) 変数コンテナー
    // vvdouble mvvArbitVar;
    //
    // Solid :tx,ty,tz
    // Shell :tx,ty,tz, rx,ry,rz
    //
    // Solid: x, y, z, xy, yz, zx
    // Shell: Mx,My,Mxy,Qx,Qy
    //
    // ---> VectorNode,ScalarNode,VectorScalarNode に移設.


    uint mMGLevel;//MultiGrid Level

    
    // prolongater プロロンゲーター用
    // 
    vector<CNode*> mvParentNode;//refine時にNodeの親となったNode

    // restriction リストリクター用
    vector<CNode*> mvChildNode;//refine時に生成した子のNode
    
    
public:
    // Node ID == Vertex ID

    //// DOF initialize
    //void InitializeNode3();// Solid Node 3
    //void InitializeNode6();// Shell Node 6
    //void InitializeNodeADOF(const vint& vParam, const uint& num_of_param);// Arbitrary DOF Node

    // MultiGrid Level
    void  setMGLevel(const uint& level){ mMGLevel=level;}
    uint& getMGLevel(){return mMGLevel;}

    // Nodeタイプ
    virtual uint& getType()=0;


    // parameter accessor
    //
    //virtual void reserveScalar(const uint& res_size)=0;
    virtual void resizeScalar(const uint& res_size)=0;
    //virtual void reserveVector(const uint& res_size)=0;
    virtual void resizeVector(const uint& res_size)=0;

    //virtual void setScalar(const double& val)=0;
    virtual void setScalar(const double& val, const uint& index)=0;
    //virtual void setVector(const vdouble& vVal)=0;
    //virtual void setVector(const double& val)=0;
    virtual void setVector(const double& val, const uint& index)=0;

    virtual double& getScalar(const uint& i)=0;
    virtual vdouble& getVector()=0;

    
    virtual uint numOfScalarParam()=0;
    virtual uint numOfVectorParam()=0;
    virtual uint numOfTotalParam()=0;

    
    // for prolongate
    vector<CNode*>& getParentNode(){ return mvParentNode;}
    CNode* getParentNode(const uint& index){ return mvParentNode[index];}
    void reserveParentNode(const uint& res_size){ mvParentNode.reserve(res_size);}
    void addParentNode(CNode* pNode){ mvParentNode.push_back(pNode);}
    uint getNumOfParentNode(){ return mvParentNode.size();}

    // for restriction
    vector<CNode*>& getChildNode(){ return mvChildNode;}
    CNode* getChildNode(const uint& index){ return mvChildNode[index];}
    void reserveChildNode(const uint& res_size){ mvChildNode.reserve(res_size);}
    void addChildNode(CNode* pNode){ mvChildNode.push_back(pNode);}
    uint getNumOfChildNode(){ return mvChildNode.size();}

};
#endif
}

