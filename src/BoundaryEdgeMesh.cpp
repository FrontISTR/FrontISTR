//
//  BoundaryEdgeMesh.cpp
//
//
//
//              2010.04.13
//              k.Takeda
#include <vector>

#include "BoundaryEdgeMesh.h"
#include "BoundaryNode.h"
using namespace pmw;


CBoundaryEdgeMesh::CBoundaryEdgeMesh()
{
    ;
}

CBoundaryEdgeMesh::~CBoundaryEdgeMesh()
{
    for_each(mvBEdge.begin(), mvBEdge.end(), DeleteObject());
}


void CBoundaryEdgeMesh::resizeEdge(const uint& res_size)
{
    mvBEdge.resize(res_size);
}


void CBoundaryEdgeMesh::setBEdge(const uint& index, CBoundaryEdge* pBEdge)
{
    uint id;
    id= pBEdge->getID();

    mmBEdgeID2Index[id]= index;

    mvBEdge[index]= pBEdge;
}


void CBoundaryEdgeMesh::addBEdge(CBoundaryEdge* pBEdge)
{
    mvBEdge.push_back(pBEdge);
    //    mvbMarkingEdge.push_back(false);//マーキング初期値:false

    uint id;
    id= pBEdge->getID();
    mmBEdgeID2Index[id]= mvBEdge.size()-1;
}


// 辺-集合 領域確保
//
void CBoundaryEdgeMesh::resizeAggEdge()
{
    uint res_size = mvBNode.size();

    mvAggregateEdge.resize(res_size);
}

// 辺-集合 計算
//
void CBoundaryEdgeMesh::setupAggEdge()
{
    
    //下位グリッドのAggデータをクリア
    //----
    uint numOfBNode= mvBNode.size();
    CBoundaryNode *pBNode;
    uint ibnode;
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pBNode->clearAggElemID();
        pBNode->clearNeibElemVert();
    };
    
    //初期設定=> 頂点に辺IDをセット
    //----
    uint numOfEdge= mvBEdge.size();
    CBoundaryEdge *pBEdge;
    uint iedge;
    for(iedge=0; iedge < numOfEdge; iedge++){
        pBEdge= mvBEdge[iedge];
        
        pBEdge->setupVertexElemID();//頂点に 辺-ID をセット
    };

    //Aggデータに移す
    //----
    for(ibnode=0; ibnode < numOfBNode; ibnode++){

        pBNode= mvBNode[ibnode];

        uint numOfAgg= pBNode->getNumOfAggElem();
        uint iagg;
        for(iagg=0; iagg < numOfAgg; iagg++){
            mvAggregateEdge[ibnode].push_back(pBNode->getAggElemID(iagg));
        };
    };
}

// 辺ノード生成
//
void CBoundaryEdgeMesh::GeneEdgeBNode()
{
    uint countID= mvBNode.size();// ID生成のベース番号

    uint numOfEdge= mvBEdge.size();
    uint iedge;
    
    CBoundaryEdge *pBEdge;

    mvBEdgeBNode.reserve(numOfEdge);
    CBoundaryNode *pBNode;

    for(iedge=0; iedge < numOfEdge; iedge++){

        pBEdge= mvBEdge[iedge];

        //新BNode(辺のBNode)
        pBNode = new CBoundaryNode;//<<<<<<<<<<<<< new

        countID += iedge;
        pBNode->setID(countID);        //新BNodeのID
        pBNode->setMGLevel(mMGLevel+1);//新BNodeのLevel
        pBNode->resizeValue(mMaxMGLevel-mMGLevel);

        pBEdge->setEdgeBNode(pBNode);
        mvBEdgeBNode.push_back(pBNode);
        
        
        //辺"BNode"に,辺"Node"をセット
        pBEdge->setupNode();
    };
}


// Refine 辺の再分割  => progBMeshにセット
// ----
void CBoundaryEdgeMesh::refine(CBoundaryEdgeMesh* pProgEdgeMesh)
{
    CBoundaryEdge *pBEdge;
    vector<CBoundaryEdge*> vBEdge;
    
    uint countID(0);//新BEdgeのID
    uint iedge, numOfBEdge= mvBEdge.size();
    for(iedge=0; iedge < numOfBEdge; iedge++){
        pBEdge= mvBEdge[iedge];

        pBEdge->refine(countID, mvDOF);// Refine

        vBEdge.clear();
        vBEdge= pBEdge->getProgParts();

        uint iprog;
        for(iprog=0; iprog < vBEdge.size(); iprog++){
            pProgEdgeMesh->addBEdge(vBEdge[iprog]);
        };
    };

    uint numOfBNode= mvBNode.size();
    uint numOfEdgeBNode= mvBEdgeBNode.size();
    uint numOfProgBNode= numOfBNode + numOfEdgeBNode;

    pProgEdgeMesh->resizeBNode(numOfProgBNode);

    uint ibnode;
    CBoundaryNode *pBNode;
    // コースグリッドのBNodeをセット
    // --
    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pProgEdgeMesh->setBNode(ibnode, pBNode);
    };
    // 新BNode( 辺BNode )をセット
    // --
    for(ibnode=numOfBNode; ibnode < numOfProgBNode; ibnode++){
        pBNode= mvBEdgeBNode[ibnode-numOfBNode];

        pProgEdgeMesh->setBNode(ibnode, pBNode);
    };
}

// ノイマン型境界値のBNodeへの再配分
//
void CBoundaryEdgeMesh::distNeumannValue()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ノイマン条件が設定されていなければ,Errorを返す.
    // --
    if(mnBndType != BoundaryType::Neumann){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryEdgeMesh::distNeumannValue");
        return;
    }
    
    uint idof, dof, numOfDOF;
    uint inode,numOfBNode=mvBNode.size();
    CBoundaryNode *pBNode;
    
    //節点力-ゼロクリア
    for(inode=0; inode < numOfBNode; inode++){
        pBNode= mvBNode[inode];

        numOfDOF= getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++){
            dof= getDOF(idof);
            pBNode->initValue(dof, mMGLevel);//mgLevel別-自由度別の値をゼロクリア
        };
    };
    
    //等価節点力配分
    CShapeLine *pShLine = CShapeLine::Instance();//形状関数
    CBoundaryEdge *pBEdge;
    uint  iedge, ivert, numOfEdge=mvBEdge.size();
    double entVal,integVal,nodalVal;
    
    for(iedge=0; iedge < numOfEdge; iedge++){
        pBEdge = mvBEdge[iedge];
        
        uint numOfDOF = getNumOfDOF();
        for(idof=0; idof < numOfDOF; idof++){
            
            dof = getDOF(idof);
            entVal = pBEdge->getBndValue(dof);

            //debug
            cout << "CBoundaryEdgeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", entVal=" << entVal << endl;
            
            switch(pBEdge->getBEdgeShape()){
                case(ElementType::Beam):
                    for(ivert=0; ivert < 2; ivert++){
                        integVal= pShLine->getIntegValue2(ivert);
                        nodalVal= integVal * entVal;//等価節点力

                        //debug
                        cout << "CBoundaryEdgeMesh::distNeumannValue, mgLevel=" << mMGLevel << ", dof=" << dof << ", nodalVal=" << nodalVal << endl;

                        pBNode= pBEdge->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//加算
                    };
                    break;
                case(ElementType::Beam2):
                    for(ivert=0; ivert < 3; ivert++){
                        integVal= pShLine->getIntegValue3(ivert);
                        nodalVal= integVal * entVal;//等価節点力

                        pBNode= pBEdge->getBNode(ivert);
                        pBNode->addValue(dof, mMGLevel, nodalVal);//加算
                    };
                    break;
                default:
                    break;
            }//switch
        };//idof loop
    };//iedge loop
}

// ディレクレ型境界値のBNodeへの再配分(Level==0の場合の最初の処理)
//
void CBoundaryEdgeMesh::distDirichletValue_at_CGrid()
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    // ノイマン条件が設定されていなければ,Errorを返す.
    // --
    if(mnBndType != BoundaryType::Dirichlet){
        pLogger->Info(Utility::LoggerMode::Error, "BoundaryType Error,  CBoundaryEdgeMesh::distDirichletValue");
        return;
    }

    //節点集合
    uint inode, numOfBNode=mvBNode.size();
    for(inode=0; inode < numOfBNode; inode++){
        
        CBoundaryNode *pBNode = mvBNode[inode];
        vuint vAggID = mvAggregateEdge[inode];
        
        double dLength, dVal, dDirichletVal;
        CBoundaryEdge *pBEdge;
        uint numOfDOF = getNumOfDOF();
        uint idof, dof;
        for(idof=0; idof < numOfDOF; idof++){

            dof = getDOF(idof);

            dLength=0.0; dVal=0.0;

            uint iedge, numOfEdge=vAggID.size();
            for(iedge=0; iedge < numOfEdge; iedge++){
                uint nEdgeID = vAggID[iedge];
                uint nEdgeIndex = mmBEdgeID2Index[nEdgeID];

                pBEdge = mvBEdge[nEdgeIndex];

                dVal    += pBEdge->getLength() * pBEdge->getBndValue(dof);
                dLength += pBEdge->getLength();
            };
            dDirichletVal = dVal/dLength;// 節点値 = Σ(距離×値)/Σ(距離)

            pBNode->setValue(dof, mMGLevel, dDirichletVal);
        };
    };
}

// Fine Grid
// Level>=1:BNode間の平均
//
void CBoundaryEdgeMesh::distDirichletValue_at_FGrid()
{
    //上位のグリッドの値をきめる
    //
    uint iedge, numOfEdge=mvBEdge.size();
    CBoundaryEdge *pBEdge;
    for(iedge=0; iedge < numOfEdge; iedge++){
        pBEdge = mvBEdge[iedge];

        uint idof, dof;
        for(idof=0; idof < getNumOfDOF(); idof++){
            dof = getDOF(idof);
            pBEdge->distDirichletVal(dof, mMGLevel);
        };
    };
}


// Refine時のmvBEdgeBNodeの解放
//
void CBoundaryEdgeMesh::deleteProgData()
{
    vector<CBoundaryNode*>().swap(mvBEdgeBNode);
}






