//
//  FileWriterAVS.cpp
//
//          2011.04.09
//          k.Takeda
#include "FileWriterAVS.h"
using namespace FileIO;


CFileWriterAVS::CFileWriterAVS()
{
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(0);

    uiint nNumOfMesh = pAssy->getNumOfMesh();

    mvLabel.resize(nNumOfMesh); //ラベル名の集合
    mmUnit.resize(nNumOfMesh);  //ラベル毎の単位名
    mmDOF.resize(nNumOfMesh);   //ラベル毎のDOF
    mmVecVal.resize(nNumOfMesh);//ラベル毎の値の配列(節点×自由度)
}
CFileWriterAVS::~CFileWriterAVS()
{
    ;
}

void CFileWriterAVS::WriteDebug(ofstream& ofs, const uiint& mgLevel)
{
    Utility::CLogger* pLogger= Utility::CLogger::Instance();

    pLogger->Info(Utility::LoggerMode::Warn, "invalid method FileWriterAVS::WriteDebug");
}

//
// メッシュ部分出力
//
void CFileWriterAVS::WriteMesh(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    // ----
    // MicroAVS UCD format 出力 (Mesh)
    // ----
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CAssyVector *pSolAssyVector= pAssy->getSolutionAssyVector(iMesh);
    
    //出力幅
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    
    uiint nIWidth= nLength, nDWidth=15;
    //----print:Node数 要素数 節点パラメータ数 要素パラメータ数 モデルパラメータ数
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfNode()    << " ";
    ofs << setw(nIWidth) << right << dec << pAssy->getMesh(iMesh)->getNumOfElement() << " ";
    ofs << setw(nIWidth) << right << dec << pSolAssyVector->getDOF() << " " << 0 << " " << 0 << endl;

    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++)
    {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();
        double x = pAssy->getMesh(iMesh)->getNodeIX(i)->getX();
        double y = pAssy->getMesh(iMesh)->getNodeIX(i)->getY();
        double z = pAssy->getMesh(iMesh)->getNodeIX(i)->getZ();

        //----print:Node X Y Z
        ofs << setw(nIWidth) << right << dec << id << " ";
        ofs << setw(nDWidth) << right << scientific  << x << " " << y << " " << z << endl;

    }//Nodeループ

    vector<pmw::CNode*> vNode;
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfElement(); i++)
    {
        vNode.clear();
        pmw::CElement *pElem = pAssy->getMesh(iMesh)->getElementIX(i);

        uiint id = pElem->getID();
        vNode = pElem->getNode();

        if(pElem->getType() == pmw::ElementType::Hexa){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Hexa2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " hex2  ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Hexa2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Tetra2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tet2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Tetra2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Prism2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " prism2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Prism2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Quad2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " quad2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Quad2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Triangle2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " tri2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Triangle2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
        if(pElem->getType() == pmw::ElementType::Beam2){
            //----print:要素構成Node
            ofs << setw(nIWidth) << right << dec << id << " ";
            ofs << setw(1) << right << dec << 0 << " line2 ";
            for(uiint iENode=0; iENode < pmw::NumberOfNode::Beam2(); iENode++){
                ofs << setw(nIWidth) << right << dec << vNode[iENode]->getID() << " ";
            };
            ofs << endl;
        }
    };//Elementループ
}

//
// 基礎変数の出力
//
void CFileWriterAVS::WriteBasis(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    // ----
    // MicroAVS UCD format 出力
    // ----
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    pmw::CAssyVector *pSolAssyVector= pAssy->getSolutionAssyVector(iMesh);

    CFileWriterAVS::WriteMesh(ofs, iMesh, mgLevel);// メッシュ部分出力

    //出力幅
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    
    uiint nIWidth= nLength, nDWidth=15;
    uiint nDOF = pSolAssyVector->getDOF();// 解ベクトルのDOF

    //----print:成分数 成分構成数
    ofs << 1 << " " << nDOF << endl;
    //----print:ラベル,単位
    ofs << "basis, variable" << endl;

    // 計算結果 : 節点ごとの値
    for(uiint i=0; i < pAssy->getMesh(iMesh)->getNumOfNode(); i++)
    {
        pmw::CNode *pNode = pAssy->getMesh(iMesh)->getNodeIX(i);
        uiint id = pNode->getID();

        //----print: NodeID
        ofs << setw(nIWidth) << right << dec << id << " ";
        for(uiint idof=0; idof < nDOF; idof++){
            //----print: 解ベクトルの値
            ofs << setw(nDWidth) << right << scientific << pSolAssyVector->getValue(iMesh, i, idof) << " ";
        }
        ofs << endl;
    };//Nodeループ
    
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Info, "Output UCD file");
}


//
// 変数登録 pValue[nMaxNum]
//
void CFileWriterAVS::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    string sLabel=cLabel;
    mvLabel[iMesh].push_back(sLabel);

    string sUnit=cUnit;
    mmUnit[iMesh][sLabel]= sUnit;

    mmDOF[iMesh][sLabel]= nNumOfDOF;
    
}
void CFileWriterAVS::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    string sLabel=cLabel;

    uiint nDOF= mmDOF[iMesh][sLabel];
    vdouble vValue;
    vValue.resize(nNumOfNode*nDOF);
    for(uiint i=0; i < vValue.size(); i++) vValue[i]=pvValue[i];

    mmVecVal[iMesh][sLabel]= vValue;
}
//
// 登録変数の出力
//
void CFileWriterAVS::WriteFEM(ofstream& ofs, const uiint& iMesh, const uiint& mgLevel)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    
    
    CFileWriterAVS::WriteMesh(ofs, iMesh, mgLevel);//メッシュ部分出力
    
    pmw::CAssyModel *pAssy= mpGMGModel->getAssyModel(mgLevel);
    
    //出力幅
    pmw::CMesh *pMesh= pAssy->getMesh(iMesh);
    uiint nNumOfNode= pMesh->getNumOfNode();
    stringstream ss;
    ss << nNumOfNode;
    uiint nLength= ss.str().length();
    uiint nIWidth= nLength, nDWidth=15;

//    vector<vstring> mvLabel;//ラベル名の集合
//    vector<map<string, string> >  mmUnit;//ラベル毎の単位名
//    vector<map<string, uiint> >   mmDOF; //ラベル毎のDOF
//    vector<map<string, vdouble> > mmVecVal;//ラベル毎の値の配列(節点×自由度)
    
    //----print:成分数 
    uiint nNumOfCompo = mvLabel[iMesh].size();
    ofs << nNumOfCompo << " ";
    for(uiint i=0; i < nNumOfCompo; i++){
        string sLabel= mvLabel[iMesh][i];
        uiint nDOF= mmDOF[iMesh][sLabel];
        //-----print:成分構成数
        ofs << nDOF  << " ";
    }
    ofs << endl;
    
    //----print:ラベル,単位
    for(uiint i=0; i < nNumOfCompo; i++){
        string sLabel= mvLabel[iMesh][i];
        ofs << mvLabel[iMesh][i] << ", " << mmUnit[iMesh][sLabel] << endl;
    }
    
    
    // 計算結果 : 節点ごとの値
    for(uiint iNode=0; iNode < nNumOfNode; iNode++)
    {
        pmw::CNode *pNode = pMesh->getNodeIX(iNode);
        uiint id = pNode->getID();

        //----print:NodeID
        ofs << setw(nIWidth) << right << dec << id << " ";
        //----print: 値 値 値 値 ・・・・
        for(uiint icompo=0; icompo < nNumOfCompo; icompo++){
            string sLabel= mvLabel[iMesh][icompo];
            uiint nDOF= mmDOF[iMesh][sLabel];
            for(uiint idof=0; idof < nDOF; idof++)
            ofs << setw(nDWidth) << right << scientific << mmVecVal[iMesh][sLabel][iNode*nDOF + idof] << " ";
        }
        ofs << endl;
    };//Nodeループ

    pLogger->Info(Utility::LoggerMode::Info, "Output UCD file");
}











