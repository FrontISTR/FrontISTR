//
// Node.h <- Vertex
//
//
//			2008.05.25
//			2008.12.18
//			k.Takeda

#ifndef NODE_HH_5DBF3A6F_8A65_4ddd_B389_2CB193B87807
#define NODE_HH_5DBF3A6F_8A65_4ddd_B389_2CB193B87807

#include "CommonStd.h"
#include "TypeDef.h"

#include "NodeType.h"
#include "Logger.h"

#include "Vertex.h"
namespace pmw{
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

    // Nodeの接続先; Node, Element_Index
    //
    //vector<vector<CNode*> >   mvvConnNode;// ConnectionNodes, 1st_index => MGLevel
    //vector<vector<uint> > mvvConnElemIndex;// ConnectionNodes, 1st_index => MGlevel

    
    // prolongater プロロンゲーター用
    // 
    vector<CNode*> mParentsNode;//refine時にNodeの親となったNode
    
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
    virtual void reserveScalar(const uint& res_size)=0;
    virtual void reserveVector(const uint& res_size)=0;

    virtual void setScalar(const double& val)=0;
    virtual void setScalar(const double& val, const uint& index)=0;
    virtual void setVector(const vdouble& vVal)=0;
    virtual void setVector(const double& val)=0;
    virtual void setVector(const double& val, const uint& index)=0;

    virtual double& getScalar(const uint& i)=0;
    virtual vdouble& getVector()=0;

    
    virtual uint numOfScalarParam()=0;
    virtual uint numOfVectorParam()=0;

    
    // prolongater
    vector<CNode*>& getParentsNode(){ return mParentsNode;}
    CNode* getParentsNode(const uint& index){ return mParentsNode[index];}

    void reserveParentsNode(const uint& res_size){ mParentsNode.reserve(res_size);}
    void addParentsNode(CNode* pNode){ mParentsNode.push_back(pNode);}
    uint getNumOfParentsNode(){ return mParentsNode.size();}
};
}
#endif
