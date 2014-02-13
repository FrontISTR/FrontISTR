/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/HEC_MW3.cpp
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
#include "Mesh.h"
#include "ShapePrism.h"
#include "ShapeHexa.h"
#include "Jacobian.h"
#include "AssyModel.h"
#include "GMGModel.h"
#include "AssyMatrix.h"
#include "AssyVector.h"
#include "SolverCG.h"
#include "SolverBiCGSTAB.h"
#include "SolverGPBiCG.h"
#include "SolverGMRES.h"
#include "SolverMG.h"

#include "MatrixBCRS.h"
#include "HEC_MPI.h"

#ifdef MSVC
#include "HEC_MW3.hxx"
#else
#include "HEC_MW3.h"
#endif


using namespace pmw;
#define MW_ERROR 0
#define MW_SUCCESS 1
CMW::CMW(void)
{
    mpMPI = CHecMPI::Instance();
}
CMW::~CMW(void)
{
    ;
}
uiint CMW::Initialize(int argc, char** argv)
{
    mpMPI->Initialize(argc, argv);
    uiint rank= mpMPI->getRank();


    mbRefine_Use = false;

    mpGMGModel= CGMGModel::Instance();

    mnSolutionType = SolutionType::FEM;//----

    mpFactory = CMeshFactory::Instance();
    mpFactory->setGMGModel(mpGMGModel);
    mpFactory->setSolutionType(mnSolutionType);//----

    mpFileIO = FileIO::CFileIO::Instance();
    mb_file = false;
    mpFileIO->setFactory(mpFactory);
    mpFileIO->setSolutionType(mnSolutionType);//----


    mpLogger = Utility::CLogger::Instance();

    mpShapeHexa= pmw::CShapeHexa::Instance();
    mpShapeHexaNic= pmw::CShapeHexaNic::Instance();
    mpShapeTetra= pmw::CShapeTetra::Instance();
    mpShapePrism= pmw::CShapePrism::Instance();
    mpShapeQuad= pmw::CShapeQuad::Instance();
    mpShapeTriangle= pmw::CShapeTriangle::Instance();
    mpShapeLine= pmw::CShapeLine::Instance();
    mpShapeCatalog= pmw::CShapeFunctionCatalog::Instance();
    
    mpISTR2Edge= pmw::CISTR2Edge::Instance();
    mpEdge2ISTR= pmw::CEdge2ISTR::Instance();


    mpGMGModel->initAssyModel();//AssyModel(Level=0)の生成
    mpFactory->setMGLevel(0);

    mpLogger->setMode(Utility::LoggerMode::MWDebug);
    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn,    Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info,    Utility::LoggerDevice::Display);
    
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");
    mpLogger->initializeLogFile(rank);

    // API_Fortran 定数
    mnMPI_UIINT= mpMPI->MPI_UIINT();
    mnMPI_IINT = mpMPI->MPI_IINT();
    mnMyRank   = rank;
    mnNumOfProcess= mpMPI->getNumOfProcess();
    
    return rank;
}
uiint CMW::BaseName_BCast(uiint& nLength, string& sName, uiint nType)
{
    if(mpMPI->getRank()==0){
        sName= mpFileIO->getFstr_FileName(nType);//--- fstr全体制御データ記載の各filename文字列長をPE=0から取得
        nLength= sName.length();
    }
    MPI_Datatype MPI_UIINT= mpMPI->MPI_UIINT();
    mpMPI->Barrier(MPI_COMM_WORLD);
    mpMPI->Bcast(&nLength, 1, MPI_UIINT, 0, MPI_COMM_WORLD);
    
    char *cBaseName;
    cBaseName = (char*)malloc(sizeof(char) * nLength);

    if(mpMPI->getRank()==0){
        for(uiint i=0; i < nLength; i++) cBaseName[i] = sName[i];
    }
    mpMPI->Barrier(MPI_COMM_WORLD);
    mpMPI->Bcast((void*)cBaseName, nLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    if(mpMPI->getRank()!=0){
        mpFileIO->setFstr_FileName(cBaseName, nLength, nType);//--- fstr全体制御ファイル内容をBCAST
    }
    free(cBaseName);

    return nLength;
}
uiint CMW::Initialize_fstr(int argc, char** argv, string& ctrlname)
{
    CMW::Initialize(argc, argv);
    
    mpFileIO->markingFstrStyle();
    //--
    // 全体制御ファイルに記載の各ファイル名をPE=0から取得 && PE > 0に配布
    //--
    if(mpMPI->getRank()==0){
        bool bCntSuccess = mpFileIO->Read_fstr_CntFile(ctrlname);// hecmw_ctrl.dat 読み込み
        if(!bCntSuccess){
            mpLogger->Info(Utility::LoggerMode::Error, "MW::Initialize_fstr interruption");
            return MW_ERROR;
        }
    }
    uiint nLength;
    string sBaseName;
    for(uiint iType=0; iType < FileIO::FileTypeMW2::Limit; iType++){
        CMW::BaseName_BCast(nLength, sBaseName, iType);//-------全体制御データのfilenameを全プロセスで共有
    };
    msOutputFileName  += mpFileIO->getFstr_MeshFileName()    + ".";
    msInputFileName   += mpFileIO->getFstr_MeshFileName()    + ".";
    msResFileName     += mpFileIO->getFstr_RestartFileName() + ".";
    msRltFileName     += mpFileIO->getFstr_ResultFileName()  + ".";
    msPartOutFileName += mpFileIO->getFstr_PartFileName_OUT()+ ".";

    uiint rank= mpMPI->getRank();
    stringstream ss;
    ss << rank;

    msOutputFileName  += ss.str();
    msInputFileName   += ss.str();
    msResFileName     += ss.str();
    msRltFileName     += ss.str();
    msPartOutFileName += ss.str();
    msOutputFileName += ".out";

    
    //--
    // 全体制御データ記載のRefine数、RefineTypeをPE=0から取得 && PE > 0 に配布
    //--
    uiint nRefineNum(0), nRefineType(1);
    MPI_Datatype MPI_UIINT= mpMPI->MPI_UIINT();
    
    if(mpMPI->getRank()==0){
        nRefineNum= mpFileIO->getFstr_RefineNum();  //-- fstr全体制御データ記載のRefineNumをPE=0から取得
        nRefineType= mpFileIO->getFstr_RefineType();//-- fstr全体制御データ記載のRefineTypeをPE=0から取得
    }
    // PE=0のRefine数、RefineTypeをBCastで配布
    mpMPI->Barrier(MPI_COMM_WORLD);
    mpMPI->Bcast((void*)&nRefineNum, 1, MPI_UIINT, 0, MPI_COMM_WORLD);
    mpMPI->Bcast((void*)&nRefineType,1, MPI_UIINT, 0, MPI_COMM_WORLD);
    if(mpMPI->getRank()!=0){
        mpFileIO->setFstr_RefineNum(nRefineNum);  //-- RefineNumを PE > 0 に配布
        mpFileIO->setFstr_RefineType(nRefineType);//-- RefineTypeを PE > 0 に配布
    }

    return MW_SUCCESS;
}
string CMW::getFstr_FileName_Mesh()
{
    return msInputFileName;
}
string CMW::getFstr_FileName_Control()
{
    return mpFileIO->getFstr_ControlFileName();
}
string CMW::getFstr_FileName_Result()
{
    return msRltFileName;
}
string CMW::getFstr_FileName_Restart()
{
    return msResFileName;
}
string CMW::getFstr_FileName_PartIn()
{
    return mpFileIO->getFstr_PartFileName_IN();
}
string CMW::getFstr_FileName_PartOut()
{
    return msPartOutFileName;
}
string CMW::getFstr_FileName_VisMesh()
{
    return mpFileIO->getFstr_VisFileName_Mesh();
}
string CMW::getFstr_FileName_VisIn()
{
    return mpFileIO->getFstr_VisFileName_IN();
}
string CMW::getFstr_FileName_VisOut()
{
    return mpFileIO->getFstr_VisFileName_OUT();
}
string CMW::getFstr_RefineCADFitName()
{
    return mpFileIO->getFstr_RefineCADFitName();
}
uiint  CMW::getFstr_RefineNum()
{
    return mpFileIO->getFstr_RefineNum();
}
uiint  CMW::getFstr_RefineType()
{
    return mpFileIO->getFstr_RefineType();
}


void CMW::Banner_Display()
{
    mpLogger->InfoDisplay();
}


#ifdef REVOCAP_REFINE //==========================================REVOCAP_REFINE
//----
// REVOCAP_Refiner
//----
uiint CMW::RevocapRefine(string& filename, const uiint& nRefine)
{
    if(mbRefine_Use){
        mpLogger->Info(Utility::LoggerMode::Info, "You have already used the Refine, RevocapRefine_Function can not be used.");
        return MW_ERROR;
    }else{
        mbRefine_Use = true;
    }

    rcapInitRefiner(0,0);//初期化

    //CADファイルセット
    if(filename.length() > 4){
        rcapSetCADFilename(filename.c_str());
        mpLogger->Info(Utility::LoggerMode::Info, " CAD file was set. ", filename);
    }else{
        mpLogger->Info(Utility::LoggerMode::Info, " CAD file has not been set.");
    }


    // 階層数分のAssyModelを用意
    mpFactory->GeneAssyModel(nRefine+1);//(Level > 0)のAssyModel生成
    mpFactory->setMGLevel(nRefine);
    mpFactory->setupNodeGridSize();//Level:0のNodeにMG階層数:全てのMeshに実行
    
    // コースグリッド && メッシュパーツ一個
    mpAssy = mpGMGModel->getAssyModel(0);//(Level=0)initで生成済み


    //============================== パーツを特定：ここから
    mpMesh = mpAssy->getMesh(0);
    CIndexBucket *pBucket = mpMesh->getBucket();// コースグリッド Bucket
    mpMesh->setupAggregate(0);// コースグリッドAggElement
    //============================== パーツを特定：ここまで
    
    
    
    // 節点登録
    size_t nodeCount = mpMesh->getNumOfNode();
    float64_t* coords = (float64_t*)calloc( nodeCount*3, sizeof(float64_t));
    int32_t* globalIDs = (int32_t*)calloc( nodeCount, sizeof(int32_t));
    int32_t* localIDs  = (int32_t*)calloc( nodeCount, sizeof(int32_t));

    for(uiint i=0; i < nodeCount; i++){
        CNode *pNode= mpMesh->getNodeIX(i);
        coords[i*3]    = pNode->getX();
        coords[i*3 + 1]= pNode->getY();
        coords[i*3 + 2]= pNode->getZ();

        uiint id = pNode->getID();
        globalIDs[i]= id;                       // Node ID
        localIDs[i] = pBucket->getIndexNode(id);// Node index

        ////cout << "MW::RevocapRefine ---- globalID:" << globalIDs[i] << "  localID:" << localIDs[i] << "  rank:" << mpMPI->getRank() << endl;// debug
    };
    rcapSetNode64(nodeCount, coords, globalIDs, localIDs);//Node番号のセット
    

    ////cout << "MW::RevocapRefine ----------------- A rank:" << mpMPI->getRank() << endl;

    // 要素別 節点数
    size_t    nNumNodeHexa = NumberOfNode::Hexa(); // 要素別のノード数
    size_t    nNumNodeHexa2 = NumberOfNode::Hexa2();
    size_t    nNumNodeTetra = NumberOfNode::Tetra();
    size_t    nNumNodeTetra2 = NumberOfNode::Tetra2();
    size_t    nNumNodePrism = NumberOfNode::Prism();
    size_t    nNumNodePrism2 = NumberOfNode::Prism2();
    size_t    nNumNodeQuad = NumberOfNode::Quad();
    size_t    nNumNodeQuad2 = NumberOfNode::Quad2();
    size_t    nNumNodeTri = NumberOfNode::Triangle();
    size_t    nNumNodeTri2= NumberOfNode::Triangle2();
    size_t    nNumNodeBeam = NumberOfNode::Beam();
    size_t    nNumNodeBeam2 = NumberOfNode::Beam2();
    size_t    nNumNodePoint = 1;

    ////cout << "MW::RevocapRefine ----------------- B rank:" << mpMPI->getRank() << endl;

    
    // 要素 (初期値)
    size_t nElemCount = mpMesh->getNumOfElement();

    size_t hexaCount(0), hexa2Count(0), tetCount(0), tet2Count(0);
    size_t prismCount(0), prism2Count(0), quadCount(0), quad2Count(0);
    size_t triCount(0), tri2Count(0), beamCount(0), beam2Count(0);
    vector<CElement*> vHexaElem,  vHexa2Elem;
    vector<CElement*> vTetraElem, vTetra2Elem;
    vector<CElement*> vPrismElem, vPrism2Elem;
    vector<CElement*> vQuadElem,  vQuad2Elem;
    vector<CElement*> vTriElem,   vTri2Elem;
    vector<CElement*> vBeamElem,  vBeam2Elem;
    // Elementタイプ別 集合
    for(size_t i=0; i < nElemCount; i++){
        CElement* pElem= mpMesh->getElementIX(i);
        uiint nType = pElem->getType();
        switch(nType){
            case(ElementType::Hexa):  hexaCount++;  vHexaElem.push_back(pElem); break;
            case(ElementType::Hexa2): hexa2Count++; vHexa2Elem.push_back(pElem); break;
            case(ElementType::Tetra):  tetCount++;   vTetraElem.push_back(pElem); break;
            case(ElementType::Tetra2): tet2Count++;  vTetra2Elem.push_back(pElem);  break;
            case(ElementType::Prism):  prismCount++;   vPrismElem.push_back(pElem); break;
            case(ElementType::Prism2): prism2Count++;  vPrism2Elem.push_back(pElem); break;
            case(ElementType::Quad):  quadCount++;   vQuadElem.push_back(pElem);  break;
            case(ElementType::Quad2): quad2Count++;  vQuad2Elem.push_back(pElem);  break;
            case(ElementType::Triangle):  triCount++;  vTriElem.push_back(pElem);  break;
            case(ElementType::Triangle2): tri2Count++; vTri2Elem.push_back(pElem); break;
            case(ElementType::Beam):   beamCount++;  vBeamElem.push_back(pElem); break;
            case(ElementType::Beam2):  beam2Count++; vBeam2Elem.push_back(pElem); break;
        }
    };
    // --
    // 体積Gr,面Gr,辺Gr (初期登録) --> Quad,Triangle要素を面Grデータとする.
    // --
    // !境界条件
    // !通信界面
    // ×接合面 : Mesh間をまたぐのでREVOCAPに登録先がない.(同じ要素番号が登録できない)
    // --
    uiint nNumOfVolMesh = mpMesh->getNumOfBoundaryVolumeMesh();
    uiint nNumOfFaceMesh = mpMesh->getNumOfBoundaryFaceMesh();
    uiint nNumOfEdgeMesh = mpMesh->getNumOfBoundaryEdgeMesh();
    uiint nNumOfCommMesh2 = mpMesh->getCommMesh2Size();

    ////cout << " nNumOfVolMesh  : " << nNumOfVolMesh   << " rank:" << mpMPI->getRank() << endl;
    ////cout << " nNumOfFaceMesh : " << nNumOfFaceMesh  << " rank:" << mpMPI->getRank() << endl;
    ////cout << " nNumOfEdgeMesh : " << nNumOfEdgeMesh  << " rank:" << mpMPI->getRank() << endl;
    ////cout << " nNumOfCommMesh2: " << nNumOfCommMesh2 << " rank:" << mpMPI->getRank() << endl;

    // 境界番号ごとの個数 : 其々形状別
    size_t *bndHexaCount, *bndHexa2Count, *bndTetCount, *bndTet2Count, *bndPrismCount, *bndPrism2Count;
    size_t *bndQuadCount, *bndQuad2Count, *bndTriCount, *bndTri2Count, *bndBeamCount, *bndBeam2Count;
    
    // 境界の節点配列
    int32_t **bndHexaNodes, **bndHexa2Nodes, **bndTetNodes, **bndTet2Nodes, **bndPrismNodes, **bndPrism2Nodes;
    int32_t **bndQuadNodes, **bndQuad2Nodes, **bndTriNodes, **bndTri2Nodes, **bndBeamNodes, **bndBeam2Nodes;
    


    // 境界番-形状、  通信番号-形状
    vuint vBndVolType;  //境界番号->形状(体積)
    vuint vBndFaceType; //境界番号->形状(面)
    vuint vBndEdgeType; //境界番号->形状(辺)
//////////////    vuint vCommFaceType;//通信番号->形状(面、辺)
    // 形状別-面数,Vol数,辺数
    vuint vBndHexaCount, vBndHexa2Count, vBndTetCount, vBndTet2Count, vBndPrismCount, vBndPrism2Count;
    vuint vBndQuadCount, vBndQuad2Count, vBndTriCount, vBndTri2Count, vBndBeamCount, vBndBeam2Count;
//////////////    
    // 境界番号-形状の何番目
    vuint vBndVolTypeSeq;
    vuint vBndFaceTypeSeq;
    vuint vBndEdgeTypeSeq;
////////////////////////    vuint vCommFaceTypeSeq;
    // 形状[]=境界番号,  形状[]=通信番号
    vuint vBndHexaBndID, vBndHexa2BndID, vBndTetBndID, vBndTet2BndID, vBndPrismBndID, vBndPrism2BndID;
    vuint vBndQuadBndID, vBndQuad2BndID, vBndTriBndID, vBndTri2BndID, vBndBeamBndID, vBndBeam2BndID;
////////    vuint vCommQuadCommID, vCommQuad2CommID, vCommTriCommID, vCommTri2CommID, vCommBeamCommID, vCommBeam2CommID, vCommPointCommID;
    
    // 境界名、境界タイプ、境界自由度数, 自由度数に対応した自由度番号
    vstring vBndHexaName, vBndHexa2Name, vBndTetName, vBndTet2Name, vBndPrismName, vBndPrism2Name, vBndQuadName, vBndQuad2Name, vBndTriName, vBndTri2Name, vBndBeamName, vBndBeam2Name;
    vuint   vBndHexaBType, vBndHexa2BType, vBndTetBType, vBndTet2BType, vBndPrismBType, vBndPrism2BType, vBndQuadBType, vBndQuad2BType, vBndTriBType, vBndTri2BType, vBndBeamBType, vBndBeam2BType;
    vuint   vBndHexaDOF, vBndHexa2DOF, vBndTetDOF, vBndTet2DOF, vBndPrismDOF, vBndPrism2DOF, vBndQuadDOF, vBndQuad2DOF, vBndTriDOF, vBndTri2DOF, vBndBeamDOF, vBndBeam2DOF;
    vvuint  vvBndHexaDOF, vvBndHexa2DOF, vvBndTetDOF, vvBndTet2DOF, vvBndPrismDOF, vvBndPrism2DOF, vvBndQuadDOF, vvBndQuad2DOF, vvBndTriDOF, vvBndTri2DOF, vvBndBeamDOF, vvBndBeam2DOF;


    //形状の何番目かをカウントする；シーケンス番号
    uiint hexaSeq, hexa2Seq, tetSeq, tet2Seq, prismSeq, prism2Seq, quadSeq, quad2Seq, triSeq, tri2Seq, beamSeq, beam2Seq, pointSeq;//pointSeqは通信用途(通信点は1個しかない)
    hexaSeq=0, hexa2Seq=0, tetSeq=0, tet2Seq=0, prismSeq=0, prism2Seq=0,
    quadSeq=0, quad2Seq=0, triSeq=0, tri2Seq=0, beamSeq=0, beam2Seq=0;

    // 境界体積:BoundaryVol
    for(uiint i=0; i < nNumOfVolMesh; i++){
        CBoundaryVolumeMesh *pVolMesh= mpMesh->getBndVolumeMeshIX(i);

        pVolMesh->setMaxMGLevel(nRefine);
        pVolMesh->resizeCGrid_BNodeValue(nRefine);

        uiint nNumOfVol= pVolMesh->getNumOfVolume();
        if(nNumOfVol > 0){
            //--
            // 1個ひろって全ての代表とする：全て同一形状の必要あり(境界条件別に同じ形状の要素の必要あり)
            //--
            CBoundaryVolume *pVol= pVolMesh->getBVolumeIX(0);

            if(pVol->getElemType()==ElementType::Hexa)
                Rcap_BndParamSetup(ElementType::Hexa, nNumOfVol, vBndHexaBndID, vBndHexaName, vBndHexaBType,
                                   vBndHexaDOF, vvBndHexaDOF, vBndVolType, vBndHexaCount, vBndVolTypeSeq, hexaSeq, pVolMesh);
            if(pVol->getElemType()==ElementType::Hexa2)
                Rcap_BndParamSetup(ElementType::Hexa2, nNumOfVol, vBndHexa2BndID, vBndHexa2Name, vBndHexa2BType,
                                   vBndHexa2DOF, vvBndHexa2DOF, vBndVolType, vBndHexa2Count, vBndVolTypeSeq, hexa2Seq, pVolMesh);
            if(pVol->getElemType()==ElementType::Tetra)
                Rcap_BndParamSetup(ElementType::Tetra, nNumOfVol, vBndTetBndID, vBndTetName, vBndTetBType,
                                   vBndTetDOF, vvBndTetDOF, vBndVolType, vBndTetCount, vBndVolTypeSeq, tetSeq, pVolMesh);
            if(pVol->getElemType()==ElementType::Tetra2)
                Rcap_BndParamSetup(ElementType::Tetra2, nNumOfVol, vBndTet2BndID, vBndTet2Name, vBndTet2BType,
                                   vBndTet2DOF, vvBndTet2DOF, vBndVolType, vBndTet2Count, vBndVolTypeSeq, tet2Seq, pVolMesh);
            if(pVol->getElemType()==ElementType::Prism)
                Rcap_BndParamSetup(ElementType::Prism, nNumOfVol, vBndPrismBndID, vBndPrismName, vBndPrismBType,
                                   vBndPrismDOF, vvBndPrismDOF, vBndVolType, vBndPrismCount, vBndVolTypeSeq, prismSeq, pVolMesh);
            if(pVol->getElemType()==ElementType::Prism2)
                Rcap_BndParamSetup(ElementType::Prism2, nNumOfVol, vBndPrism2BndID, vBndPrism2Name, vBndPrism2BType,
                                   vBndPrism2DOF, vvBndPrism2DOF, vBndVolType, vBndPrism2Count, vBndVolTypeSeq, prism2Seq, pVolMesh);
        }
    };

    // 境界面:BoundaryFace
    for(uiint i=0; i < nNumOfFaceMesh; i++){
        CBoundaryFaceMesh *pFaceMesh= mpMesh->getBndFaceMeshIX(i);

        pFaceMesh->setMaxMGLevel(nRefine);
        pFaceMesh->resizeCGrid_BNodeValue(nRefine);

        uiint nNumOfFace= pFaceMesh->getNumOfBFace();
        if(nNumOfFace > 0){
            //一個ひろって全ての代表とする：全て同一形状の必要あり
            CBoundaryFace *pFace= pFaceMesh->getBFaceIX(0);

            // 境界番号別の形状:先頭の形状を取得して代表とする:全て同一形状の必要あり
            if(pFace->getBFaceShape()==ElementType::Quad){
                vBndQuadBndID.push_back(pFaceMesh->getID());     //境界ID
                vBndQuadName.push_back(pFaceMesh->getName());    //境界名
                vBndQuadBType.push_back(pFaceMesh->getBndType());//境界種類
                vBndQuadDOF.push_back(pFaceMesh->getNumOfDOF()); //境界DOF数
                
                vuint vDOF; vDOF.resize(pFaceMesh->getNumOfDOF());
                for(uiint idof=0; idof < pFaceMesh->getNumOfDOF(); idof++){
                    vDOF[idof]= pFaceMesh->getDOF(idof);
                    ////uiint dof = pFaceMesh->getDOF(idof);
                    ////cout << "MW::RevocapRefine ---------- CG pFace Value:" << pFace->getBndValue(dof) << "  dof:" << dof << " rank:" << mpMPI->getRank() << endl;
                }
                vvBndQuadDOF.push_back(vDOF);

                vBndFaceType.push_back(ElementType::Quad);
                vBndQuadCount.push_back(nNumOfFace);
                vBndFaceTypeSeq.push_back(quadSeq);
                quadSeq++;
            }
            if(pFace->getBFaceShape()==ElementType::Quad2){
                vBndQuad2BndID.push_back(pFaceMesh->getID());
                vBndQuad2Name.push_back(pFaceMesh->getName());
                vBndQuad2BType.push_back(pFaceMesh->getBndType());
                vBndQuad2DOF.push_back(pFaceMesh->getNumOfDOF());

                vuint vDOF; vDOF.resize(pFaceMesh->getNumOfDOF());
                for(uiint idof=0; idof < pFaceMesh->getNumOfDOF(); idof++)
                    vDOF[idof]= pFaceMesh->getDOF(idof);
                vvBndQuad2DOF.push_back(vDOF);

                vBndFaceType.push_back(ElementType::Quad2);
                vBndQuad2Count.push_back(nNumOfFace);
                vBndFaceTypeSeq.push_back(quad2Seq);
                quad2Seq++;
            }
            if(pFace->getBFaceShape()==ElementType::Triangle){
                vBndTriBndID.push_back(pFaceMesh->getID());
                vBndTriName.push_back(pFaceMesh->getName());
                vBndTriBType.push_back(pFaceMesh->getBndType());
                vBndTriDOF.push_back(pFaceMesh->getNumOfDOF());

                vuint vDOF; vDOF.resize(pFaceMesh->getNumOfDOF());
                for(uiint idof=0; idof < pFaceMesh->getNumOfDOF(); idof++)
                    vDOF[idof]= pFaceMesh->getDOF(idof);
                vvBndTriDOF.push_back(vDOF);

                vBndFaceType.push_back(ElementType::Triangle);
                vBndTriCount.push_back(nNumOfFace);
                vBndFaceTypeSeq.push_back(triSeq);
                triSeq++;
            }
            if(pFace->getBFaceShape()==ElementType::Triangle2){
                vBndTri2BndID.push_back(pFaceMesh->getID());
                vBndTri2Name.push_back(pFaceMesh->getName());
                vBndTri2BType.push_back(pFaceMesh->getBndType());
                vBndTri2DOF.push_back(pFaceMesh->getNumOfDOF());

                vuint vDOF; vDOF.resize(pFaceMesh->getNumOfDOF());
                for(uiint idof=0; idof < pFaceMesh->getNumOfDOF(); idof++)
                    vDOF[idof]= pFaceMesh->getDOF(idof);
                vvBndTri2DOF.push_back(vDOF);

                vBndFaceType.push_back(ElementType::Triangle2);
                vBndTri2Count.push_back(nNumOfFace);
                vBndFaceTypeSeq.push_back(tri2Seq);
                tri2Seq++;
            }
        };
    };

    // 境界辺:BoundaryEdge
    for(uiint i=0; i < nNumOfEdgeMesh; i++){
        CBoundaryEdgeMesh *pEdgeMesh= mpMesh->getBndEdgeMeshIX(i);

        pEdgeMesh->setMaxMGLevel(nRefine);
        pEdgeMesh->resizeCGrid_BNodeValue(nRefine);

        uiint nNumOfBeam= pEdgeMesh->getNumOfEdge();
        if(nNumOfBeam > 0){
            CBoundaryEdge *pEdge= pEdgeMesh->getBEdgeIX(0);
            
            if(pEdge->getBEdgeShape()==ElementType::Beam){
                vBndBeamBndID.push_back(pEdgeMesh->getID());
                vBndBeamName.push_back(pEdgeMesh->getName());
                vBndBeamBType.push_back(pEdgeMesh->getBndType());
                vBndBeamDOF.push_back(pEdgeMesh->getNumOfDOF());

                vuint vDOF; vDOF.resize(pEdgeMesh->getNumOfDOF());
                for(uiint idof=0; idof < pEdgeMesh->getNumOfDOF(); idof++)
                    vDOF[idof]= pEdgeMesh->getDOF(idof);
                vvBndBeamDOF.push_back(vDOF);

                vBndEdgeType.push_back(ElementType::Beam);
                vBndBeamCount.push_back(nNumOfBeam);
                vBndEdgeTypeSeq.push_back(beamSeq);
                beamSeq++;
            }
            if(pEdge->getBEdgeShape()==ElementType::Beam2){
                vBndBeam2BndID.push_back(pEdgeMesh->getID());
                vBndBeam2Name.push_back(pEdgeMesh->getName());
                vBndBeam2BType.push_back(pEdgeMesh->getBndType());
                vBndBeam2DOF.push_back(pEdgeMesh->getNumOfDOF());

                vuint vDOF; vDOF.resize(pEdgeMesh->getNumOfDOF());
                for(uiint idof=0; idof < pEdgeMesh->getNumOfDOF(); idof++)
                    vDOF[idof]= pEdgeMesh->getDOF(idof);
                vvBndBeam2DOF.push_back(vDOF);

                vBndEdgeType.push_back(ElementType::Beam2);
                vBndBeam2Count.push_back(nNumOfBeam);
                vBndEdgeTypeSeq.push_back(beam2Seq);
                beam2Seq++;
            }
        };
    };

    ////cout << "MW::RevocapRefine ----------------- C rank:" << mpMPI->getRank() << endl;

    bndHexaCount = (size_t*)calloc(vBndHexaCount.size(), sizeof(size_t));
    bndHexa2Count= (size_t*)calloc(vBndHexa2Count.size(), sizeof(size_t));
    bndTetCount  = (size_t*)calloc(vBndTetCount.size(), sizeof(size_t));
    bndTet2Count = (size_t*)calloc(vBndTet2Count.size(), sizeof(size_t));
    bndPrismCount= (size_t*)calloc(vBndPrismCount.size(), sizeof(size_t));
    bndPrism2Count= (size_t*)calloc(vBndPrism2Count.size(), sizeof(size_t));
    
    bndQuadCount = (size_t*)calloc(vBndQuadCount.size(), sizeof(size_t));
    bndQuad2Count= (size_t*)calloc(vBndQuad2Count.size(), sizeof(size_t));
    bndTriCount  = (size_t*)calloc(vBndTriCount.size(), sizeof(size_t));
    bndTri2Count = (size_t*)calloc(vBndTri2Count.size(), sizeof(size_t));
    bndBeamCount = (size_t*)calloc(vBndBeamCount.size(), sizeof(size_t));
    bndBeam2Count= (size_t*)calloc(vBndBeam2Count.size(), sizeof(size_t));

    
    for(uiint i=0; i < vBndHexaCount.size(); i++) bndHexaCount[i]= vBndHexaCount[i];//境界番号、形状別の境界要素数
    for(uiint i=0; i < vBndHexa2Count.size(); i++) bndHexa2Count[i]= vBndHexa2Count[i];
    for(uiint i=0; i < vBndTetCount.size(); i++) bndTetCount[i]= vBndTetCount[i];
    for(uiint i=0; i < vBndTet2Count.size(); i++) bndTet2Count[i]= vBndTet2Count[i];
    for(uiint i=0; i < vBndPrismCount.size(); i++) bndPrismCount[i]= vBndPrismCount[i];
    for(uiint i=0; i < vBndPrism2Count.size(); i++) bndPrism2Count[i]= vBndPrism2Count[i];

    for(uiint i=0; i < vBndQuadCount.size(); i++) bndQuadCount[i]= vBndQuadCount[i];
    for(uiint i=0; i < vBndQuad2Count.size(); i++) bndQuad2Count[i]= vBndQuad2Count[i];
    for(uiint i=0; i < vBndTriCount.size(); i++) bndTriCount[i]= vBndTriCount[i];
    for(uiint i=0; i < vBndTri2Count.size(); i++) bndTri2Count[i]= vBndTri2Count[i];
    for(uiint i=0; i < vBndBeamCount.size(); i++) bndBeamCount[i]= vBndBeamCount[i];
    for(uiint i=0; i < vBndBeam2Count.size(); i++) bndBeam2Count[i]= vBndBeam2Count[i];


    bndHexaNodes = (int32_t**)calloc(vBndHexaCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndHexaCount.size(); i++) bndHexaNodes[i] = (int32_t*)calloc(vBndHexaCount[i]*8, sizeof(int32_t));
    bndHexa2Nodes = (int32_t**)calloc(vBndHexa2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndHexa2Count.size(); i++) bndHexa2Nodes[i] = (int32_t*)calloc(vBndHexa2Count[i]*20, sizeof(int32_t));
    bndTetNodes = (int32_t**)calloc(vBndTetCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndTetCount.size(); i++) bndTetNodes[i] = (int32_t*)calloc(vBndTetCount[i]*4, sizeof(int32_t));
    bndTet2Nodes = (int32_t**)calloc(vBndTet2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndTet2Count.size(); i++) bndTet2Nodes[i] = (int32_t*)calloc(vBndTet2Count[i]*10, sizeof(int32_t));
    bndPrismNodes = (int32_t**)calloc(vBndPrismCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndPrismCount.size(); i++) bndPrismNodes[i] = (int32_t*)calloc(vBndPrismCount[i]*6, sizeof(int32_t));
    bndPrism2Nodes = (int32_t**)calloc(vBndPrism2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndPrism2Count.size(); i++) bndPrism2Nodes[i] = (int32_t*)calloc(vBndPrism2Count[i]*15, sizeof(int32_t));


    bndQuadNodes = (int32_t**)calloc(vBndQuadCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndQuadCount.size(); i++) bndQuadNodes[i]= (int32_t*)calloc(vBndQuadCount[i]*4, sizeof(int32_t));
    bndQuad2Nodes = (int32_t**)calloc(vBndQuad2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndQuad2Count.size(); i++) bndQuad2Nodes[i]= (int32_t*)calloc(vBndQuad2Count[i]*8, sizeof(int32_t));
    bndTriNodes = (int32_t**)calloc(vBndTriCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndTriCount.size(); i++) bndTriNodes[i]= (int32_t*)calloc(vBndTriCount[i]*3, sizeof(int32_t));
    bndTri2Nodes = (int32_t**)calloc(vBndTri2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndTri2Count.size(); i++) bndTri2Nodes[i]= (int32_t*)calloc(vBndTri2Count[i]*6, sizeof(int32_t));
    bndBeamNodes = (int32_t**)calloc(vBndBeamCount.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndBeamCount.size(); i++) bndBeamNodes[i]= (int32_t*)calloc(vBndBeamCount[i]*2, sizeof(int32_t));
    bndBeam2Nodes = (int32_t**)calloc(vBndBeam2Count.size(), sizeof(int32_t*));
    for(uiint i=0; i < vBndBeam2Count.size(); i++) bndBeam2Nodes[i]= (int32_t*)calloc(vBndBeam2Count[i]*3, sizeof(int32_t));

    ////cout << "MW::RevocapRefine ----------------- C0 rank:" << mpMPI->getRank() << endl;

    //---- REVOCAP引数用の〜Nodes[]に節点Index番号を代入:Vol ----
    //--
    // 対応Node番号代入(tempNodes[]を介在して代入)
    //--
    for(uiint i=0; i < nNumOfVolMesh; i++){
        CBoundaryVolumeMesh *pVolMesh= mpMesh->getBndVolumeMeshIX(i);

        int32_t *tempNodes;// 介在ポインタ
        // 0個の場合に対応
        if( pVolMesh->getNumOfVolume() > 0 ){
            if(vBndVolType[i]==ElementType::Hexa){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndHexaNodes[seq];
            }
            if(vBndVolType[i]==ElementType::Hexa2){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndHexa2Nodes[seq];
            }
            if(vBndVolType[i]==ElementType::Tetra){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndTetNodes[seq];
            }
            if(vBndVolType[i]==ElementType::Tetra2){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndTet2Nodes[seq];
            }
            if(vBndVolType[i]==ElementType::Prism){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndPrismNodes[seq];
            }
            if(vBndVolType[i]==ElementType::Prism2){
                uiint seq = vBndVolTypeSeq[i];
                tempNodes = bndPrism2Nodes[seq];
            }
        }
        // Bnd : Node配列をセット
        uiint nNumOfVol= pVolMesh->getNumOfVolume();
        uiint nNumOfDOF= pVolMesh->getNumOfDOF();
        for(uiint ivol=0; ivol < nNumOfVol; ivol++){

            CBoundaryVolume *pVol = pVolMesh->getBVolumeIX(ivol);
            uiint nNumOfNode= pVol->getNumOfBNode();

            for(uiint ibnode=0; ibnode < nNumOfNode; ibnode++){
                CBoundaryNode *pBNode= pVol->getBNode(ibnode);
                CNode *pNode= pBNode->getNode();
                uiint id= pNode->getID();
                uiint index= pBucket->getIndexNode(id);

                tempNodes[ivol*nNumOfNode + ibnode]= index;//対応アドレスに節点Indexをセット
            };
            if( pVolMesh->getBndType()==BoundaryType::Neumann ){
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof= pVolMesh->getDOF(idof);
                    //----------------------------------------------------------------------------------
                    Rcap_EquivalentNodalForce(0, dof, pVolMesh, pVol);//---ついで:コースグリッドの等価節点力を分配
                    //----------------------------------------------------------------------------------
                };//idof
            }
        };
    };

    ////cout << "MW::RevocapRefine ----------------- C1 rank:" << mpMPI->getRank() << endl;

    //---- REVOCAP引数用の〜Nodes[]に節点Index番号を代入:Face ----
    //--
    // 対応Node番号代入(tempNodes[]を介在して代入)
    //--
    for(uiint i=0; i < nNumOfFaceMesh; i++){
        CBoundaryFaceMesh *pFaceMesh= mpMesh->getBndFaceMeshIX(i);
        
        int32_t *tempNodes;// 介在ポインタ
        //0個の場合に対応
        if( pFaceMesh->getNumOfBFace() > 0 ){
            // 境界番号-形状-形状の何番目?かを探して=> 対応するアドレスをtempにセット
            
            if(vBndFaceType[i]==ElementType::Quad){
                uiint seq= vBndFaceTypeSeq[i];
                tempNodes = bndQuadNodes[seq];
            }
            if(vBndFaceType[i]==ElementType::Quad2){
                uiint seq= vBndFaceTypeSeq[i];
                tempNodes = bndQuad2Nodes[seq];
            }
            if(vBndFaceType[i]==ElementType::Triangle){
                uiint seq= vBndFaceTypeSeq[i];
                tempNodes = bndTriNodes[seq];
            }
            if(vBndFaceType[i]==ElementType::Triangle2){
                uiint seq= vBndFaceTypeSeq[i];
                tempNodes = bndTri2Nodes[seq];
            }
        }
        // Bnd : Node配列をセット
        uiint nNumOfFace= pFaceMesh->getNumOfBFace();
        uiint nNumOfDOF= pFaceMesh->getNumOfDOF();
        for(uiint iface=0; iface < nNumOfFace; iface++){

            CBoundaryFace *pFace= pFaceMesh->getBFaceIX(iface);
            uiint nNumOfNode= pFace->getNumOfBNode();

            for(uiint ibnode=0; ibnode < nNumOfNode; ibnode++){
                CBoundaryNode *pBNode= pFace->getBNode(ibnode);
                CNode *pNode= pBNode->getNode();
                uiint id= pNode->getID();
                uiint index= pBucket->getIndexNode(id);
                
                tempNodes[iface*nNumOfNode + ibnode]= index;//対応アドレスに節点Indexをセット
            };
            if( pFaceMesh->getBndType()==BoundaryType::Neumann ){
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof= pFaceMesh->getDOF(idof);
                    //-----------------------------------------------------------------------------------
                    Rcap_EquivalentNodalForce(0, dof, pFaceMesh, pFace);// ついで:コースグリッドの等価節点力を分配
                    //-----------------------------------------------------------------------------------
                };//idof
            }
        };
    };

    ////cout << "MW::RevocapRefine ----------------- C2 rank:" << mpMPI->getRank() << endl;

    //--
    // 対応Node番号代入(tempNodes[]を介在して代入)
    //--
    for(uiint i=0; i < nNumOfEdgeMesh; i++){
        CBoundaryEdgeMesh *pEdgeMesh= mpMesh->getBndEdgeMeshIX(i);
        
        int32_t *tempNodes;// 介在ポインタ
        // 0個の場合に対応
        if( pEdgeMesh->getNumOfEdge() > 0 ){
            // 境界番号-形状-形状の何番目?かを探して=> 対応するアドレスをtempにセット
            //
            if(vBndEdgeType[i]==ElementType::Beam){
                uiint seq= vBndEdgeTypeSeq[i];
                tempNodes = bndBeamNodes[seq];
            }
            if(vBndEdgeType[i]==ElementType::Beam2){
                uiint seq= vBndEdgeTypeSeq[i];
                tempNodes = bndBeam2Nodes[seq];
            }
        }
        // Bnd : Node配列をセット
        uiint nNumOfEdge= pEdgeMesh->getNumOfEdge();
        uiint nNumOfDOF= pEdgeMesh->getNumOfDOF();
        for(uiint iedge=0; iedge < nNumOfEdge; iedge++){

            CBoundaryEdge *pEdge= pEdgeMesh->getBEdgeIX(iedge);
            uiint nNumOfNode= pEdge->getNumOfBNode();

            for(uiint ibnode=0; ibnode < nNumOfNode; ibnode++){
                CBoundaryNode *pBNode= pEdge->getBNode(ibnode);
                CNode *pNode= pBNode->getNode();
                uiint id= pNode->getID();
                uiint index= pBucket->getIndexNode(id);

                tempNodes[iedge*nNumOfNode + ibnode]= index;//対応アドレスに節点Indexをセット
            };
            if( pEdgeMesh->getBndType()==BoundaryType::Neumann ){
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof= pEdgeMesh->getDOF(idof);
                    //------------------------------------------------------------------------------------
                    Rcap_EquivalentNodalForce(0, dof, pEdgeMesh, pEdge);//---ついで:コースグリッドの等価節点力を分配
                    //------------------------------------------------------------------------------------
                };//dof
            }
        };
    };
    
    ////cout << "MW::RevocapRefine ----------------- C3 rank:" << mpMPI->getRank() << endl;


    // 通信面:CommFace の形状別カウント数
    map<uiint,uiint> mCommQuadCount, mCommQuad2Count, mCommTriCount, mCommTri2Count, mCommBeamCount, mCommBeam2Count, mCommPointCount;
    map<uiint,uiint> mCommID;
    vvuint vvCommPointNIDs; vvCommPointNIDs.resize(nNumOfCommMesh2);//commPointのノードID確保用

    for(uiint i=0; i < nNumOfCommMesh2; i++){
        CCommMesh2 *pCommMesh= mpMesh->getCommMesh2IX(i);
        uiint nCommID= pCommMesh->getID();
        mCommID[i]= nCommID;

        uiint nNumOfCommFace= pCommMesh->getCommFaceSize();
        uiint quadSeq(0),quad2Seq(0),triSeq(0),tri2Seq(0),beamSeq(0),beam2Seq(0),pointSeq(0);

        for(uiint iface=0; iface < nNumOfCommFace; iface++){
            CCommFace *pCommFace= pCommMesh->getCommFaceIX(iface);
            
            if(pCommFace->getType()==ElementType::Quad)   quadSeq++;
            if(pCommFace->getType()==ElementType::Quad2)  quad2Seq++;
            if(pCommFace->getType()==ElementType::Triangle)  triSeq++;
            if(pCommFace->getType()==ElementType::Triangle2) tri2Seq++;
            if(pCommFace->getType()==ElementType::Beam)  beamSeq++;
            if(pCommFace->getType()==ElementType::Beam2) beam2Seq++;
            // commPointはここでノードIDを取得しておく(リファインを行わない=>NodeIDを確保するのはここだけ)
            if(pCommFace->getType()==ElementType::Point){
                pointSeq++;
                CCommNode *pCommNode= pCommFace->getCommNode(0);
                vvCommPointNIDs[i].push_back( pCommNode->getNodeID() );//---------------commPoint NodeID取得
            }
        };
        mCommQuadCount[nCommID]=quadSeq; mCommQuad2Count[nCommID]=quad2Seq; mCommTriCount[nCommID]=triSeq;
        mCommTri2Count[nCommID]=tri2Seq; mCommBeamCount[nCommID]=beamSeq; mCommBeam2Count[nCommID]=beam2Seq;
        mCommPointCount[nCommID]=pointSeq;
    };
    //通信メッシュごとの個数 
    size_t *commQuadCount, *commQuad2Count, *commTriCount, *commTri2Count, *commBeamCount, *commBeam2Count, *commPointCount;

    //通信の節点配列
    int32_t **commQuadNodes, **commQuad2Nodes, **commTriNodes, **commTri2Nodes, **commBeamNodes, **commBeam2Nodes, **commPointNodes;
    
    ////cout << "MW::RevocapRefine ----------------- D rank:" << mpMPI->getRank() << endl;

    commQuadCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commQuad2Count = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commTriCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commTri2Count = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commBeamCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commBeam2Count = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    commPointCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));

    for(uiint i=0; i < nNumOfCommMesh2; i++){
        CCommMesh2 *pCommMesh= mpMesh->getCommMesh2IX(i);
        uiint nCommID= pCommMesh->getID();
        
        commQuadCount[i]= mCommQuadCount[nCommID];
        commQuad2Count[i]= mCommQuad2Count[nCommID];
        commTriCount[i]= mCommTriCount[nCommID];

        ////cout << "MW::RevocapRefine ----------------- D0  icom:" << i << " commTriCount:" << commTriCount[i] << " rank:" << mpMPI->getRank() << endl;

        commTri2Count[i]= mCommTri2Count[nCommID];
        commBeamCount[i]= mCommBeamCount[nCommID];
        commBeam2Count[i]= mCommBeam2Count[nCommID];
        commPointCount[i]= mCommPointCount[nCommID];
    }

    commQuadNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commQuad2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commTriNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commTri2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commBeamNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commBeam2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    commPointNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    
    for(uiint i=0; i < nNumOfCommMesh2; i++){
        commQuadNodes[i]= (int32_t*)calloc( commQuadCount[i]*4, sizeof(int32_t));
        commQuad2Nodes[i]= (int32_t*)calloc( commQuad2Count[i]*8, sizeof(int32_t));
        commTriNodes[i]= (int32_t*)calloc( commTriCount[i]*3, sizeof(int32_t));
        commTri2Nodes[i]= (int32_t*)calloc( commTri2Count[i]*6, sizeof(int32_t));
        commBeamNodes[i]= (int32_t*)calloc( commBeamCount[i]*2, sizeof(int32_t));
        commBeam2Nodes[i]= (int32_t*)calloc( commBeam2Count[i]*3, sizeof(int32_t));
        commPointNodes[i]= (int32_t*)calloc( commPointCount[i]*1, sizeof(int32_t));
    };
    // commPointのノードIDは、コースグリッドをそのままファイングリッドで利用するので、ここでノードIDを代入しておく.
    for(uiint i=0; i < nNumOfCommMesh2; i++){
        for( uiint ipoint=0 ; ipoint < commPointCount[i]; ipoint++){
            uiint id = vvCommPointNIDs[i][ipoint];//-----commPointのNodeID(MeshのNodeID番号)
            uiint index= pBucket->getIndexNode(id);
            commPointNodes[i][ipoint]= index;//----------index番号(REVOCAP_Refinerの節点番号)
        };
    };
    //----
    // commNode -> 対応Node番号代入(BndはtempNodes[]を介在していたが、Commは直接代入)
    //----
    for(uiint i=0; i < nNumOfCommMesh2; i++){
        CCommMesh2 *pCommMesh= mpMesh->getCommMesh2IX(i);

        quadSeq=0, quad2Seq=0, triSeq=0, tri2Seq=0, beamSeq=0, beam2Seq=0, pointSeq=0;//面別のカウント

        // Comm : Node配列をセット
        uiint nNumOfFace= pCommMesh->getCommFaceSize();
        for(uiint iface=0; iface < nNumOfFace; iface++){

            CCommFace *pCommFace= pCommMesh->getCommFaceIX(iface);
            uiint nNumOfNode= pCommFace->getCommNodeSize();
            
            // 型ごとにNode番号(REVOCAP_Refinerの節点番号 : MW3の節点Index番号)
            if(pCommFace->getType()==ElementType::Quad){
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commQuadNodes[i][quadSeq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                quadSeq++;
            }
            if(pCommFace->getType()==ElementType::Quad2){
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commQuad2Nodes[i][quad2Seq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                quad2Seq++;
            }
            if(pCommFace->getType()==ElementType::Triangle){
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commTriNodes[i][triSeq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                triSeq++;
            }
            if(pCommFace->getType()==ElementType::Triangle2){
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commTri2Nodes[i][tri2Seq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                tri2Seq++;
            }
            if(pCommFace->getType()==ElementType::Beam){  
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commBeamNodes[i][beamSeq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                beamSeq++;
            }
            if(pCommFace->getType()==ElementType::Beam2){
                for(uiint icnode=0; icnode < nNumOfNode; icnode++){
                    CCommNode *pCommNode = pCommFace->getCommNode(icnode);
                    CNode *pNode= pCommNode->getNode();
                    uiint id= pNode->getID();
                    uiint index= pBucket->getIndexNode(id);

                    commBeam2Nodes[i][beam2Seq*nNumOfNode + icnode]= index;//---対応アドレスに節点Indexをセット
                };
                beam2Seq++;
            }
        };
    };
    
    
    //初期 要素ノード配列 (コースグリッドのElementNodes配列 )
    int32_t*  hexaElemNodes = (int32_t*)calloc(nNumNodeHexa*hexaCount, sizeof(int32_t)); // 8節点/1要素
    int32_t*  hexa2ElemNodes= (int32_t*)calloc(nNumNodeHexa2*hexa2Count, sizeof(int32_t));
    int32_t*  tetraElemNodes = (int32_t*)calloc(nNumNodeTetra*tetCount, sizeof(int32_t));
    int32_t*  tetra2ElemNodes = (int32_t*)calloc(nNumNodeTetra2*tet2Count, sizeof(int32_t));
    int32_t*  prismElemNodes = (int32_t*)calloc(nNumNodePrism*prismCount, sizeof(int32_t));
    int32_t*  prism2ElemNodes = (int32_t*)calloc(nNumNodePrism2*prism2Count, sizeof(int32_t));
    int32_t*  quadElemNodes = (int32_t*)calloc(nNumNodeQuad*quadCount, sizeof(int32_t));
    int32_t*  quad2ElemNodes = (int32_t*)calloc(nNumNodeQuad2*quad2Count, sizeof(int32_t));
    int32_t*  triElemNodes = (int32_t*)calloc(nNumNodeTri*triCount, sizeof(int32_t));
    int32_t*  tri2ElemNodes = (int32_t*)calloc(nNumNodeTri2*tri2Count, sizeof(int32_t));
    int32_t*  beamElemNodes = (int32_t*)calloc(nNumNodeBeam*beamCount, sizeof(int32_t));
    int32_t*  beam2ElemNodes = (int32_t*)calloc(nNumNodeBeam2*beam2Count, sizeof(int32_t));
    
    Rcap_CrsElemNodes( nNumNodeHexa,  vHexaElem, pBucket,  hexaElemNodes); //Hexa
    Rcap_CrsElemNodes( nNumNodeHexa2, vHexa2Elem, pBucket, hexa2ElemNodes);//Hexa2
    Rcap_CrsElemNodes( nNumNodeTetra,  vTetraElem, pBucket,  tetraElemNodes); //Tetra
    Rcap_CrsElemNodes( nNumNodeTetra2, vTetra2Elem, pBucket, tetra2ElemNodes);//Tetra2
    Rcap_CrsElemNodes( nNumNodePrism,  vPrismElem, pBucket,  prismElemNodes); //Prism
    Rcap_CrsElemNodes( nNumNodePrism2, vPrism2Elem, pBucket, prism2ElemNodes);//Prism2
    Rcap_CrsElemNodes( nNumNodeQuad,  vQuadElem, pBucket,  quadElemNodes); //Quad
    Rcap_CrsElemNodes( nNumNodeQuad2, vQuad2Elem, pBucket, quad2ElemNodes);//Quad2
    Rcap_CrsElemNodes( nNumNodeTri,  vTriElem, pBucket, triElemNodes);  //Triangle
    Rcap_CrsElemNodes( nNumNodeTri2, vTri2Elem, pBucket, tri2ElemNodes);//Triangle2
    Rcap_CrsElemNodes( nNumNodeBeam,  vBeamElem, pBucket, beamElemNodes);  //Beam
    Rcap_CrsElemNodes( nNumNodeBeam2, vBeam2Elem, pBucket, beam2ElemNodes);//Beam2
    
    
    //////debug
    ////cout << " hexaCount :" << hexaCount  << " hexa2Count :" << hexa2Count  << endl;
    ////cout << " tetCount  :" << tetCount   << " tet2Count  :" << tet2Count   << endl;
    ////cout << " prismCount:" << prismCount << " prism2Count:" << prism2Count << endl;
    ////cout << " quadCount :" << quadCount  << " quad2Count :" << quad2Count  << endl;
    ////cout << " triCount  :" << triCount   << " tri2Count  :" << tri2Count   << endl;
    ////cout << " beamCount :" << beamCount  << " beam2Count :" << beam2Count  << endl;

    ////cout << "MW::RevocapRefine ----------------- E rank:" << mpMPI->getRank() << endl;

    // ----
    // リファイン:Mesh
    // ----
    size_t   refineHexaCount,  refineHexa2Count, refineTetCount, refineTet2Count;
    size_t   refinePrismCount, refinePrism2Count;
    size_t   refineQuadCount, refineQuad2Count, refineTriCount, refineTri2Count;
    size_t   refineBeamCount, refineBeam2Count;
    int32_t *refineHexaNodes,  *refineHexa2Nodes,  *refineTetNodes,  *refineTet2Nodes;
    int32_t *refinePrismNodes, *refinePrism2Nodes, *refineQuadNodes, *refineQuad2Nodes;
    int32_t *refineTriNodes,   *refineTri2Nodes,   *refineBeamNodes,   *refineBeam2Nodes;

    float64_t *resultCoord; int32_t refineNodeCount; int32_t  nodeOffset=0;
    // ----
    // リファイン:Bnd && Comm
    // ----
    // 境界番号ごとの個数, 通信メッシュごとの個数 : 其々形状別
    size_t *refineBndHexaCount, *refineBndHexa2Count, *refineBndTetCount, *refineBndTet2Count, *refineBndPrismCount, *refineBndPrism2Count;
    size_t *refineBndQuadCount, *refineBndQuad2Count, *refineBndTriCount, *refineBndTri2Count, *refineBndBeamCount, *refineBndBeam2Count;
    size_t *refineCommQuadCount, *refineCommQuad2Count, *refineCommTriCount, *refineCommTri2Count, *refineCommBeamCount, *refineCommBeam2Count;
    // 境界、通信の節点配列
    int32_t **refineBndHexaNodes, **refineBndHexa2Nodes, **refineBndTetNodes, **refineBndTet2Nodes, **refineBndPrismNodes, **refineBndPrism2Nodes;
    int32_t **refineBndQuadNodes, **refineBndQuad2Nodes, **refineBndTriNodes, **refineBndTri2Nodes, **refineBndBeamNodes, **refineBndBeam2Nodes;
    int32_t **refineCommQuadNodes, **refineCommQuad2Nodes, **refineCommTriNodes, **refineCommTri2Nodes, **refineCommBeamNodes, **refineCommBeam2Nodes;
    
    refineBndHexaCount = (size_t*)calloc(vBndHexaCount.size(), sizeof(size_t));   refineBndHexa2Count = (size_t*)calloc(vBndHexa2Count.size(), sizeof(size_t));
    refineBndTetCount = (size_t*)calloc(vBndTetCount.size(), sizeof(size_t));     refineBndTet2Count = (size_t*)calloc(vBndTet2Count.size(), sizeof(size_t));
    refineBndPrismCount = (size_t*)calloc(vBndPrismCount.size(), sizeof(size_t)); refineBndPrism2Count = (size_t*)calloc(vBndPrism2Count.size(), sizeof(size_t));
    refineBndQuadCount = (size_t*)calloc(vBndQuadCount.size(), sizeof(size_t)); refineBndQuad2Count= (size_t*)calloc(vBndQuad2Count.size(), sizeof(size_t));
    refineBndTriCount = (size_t*)calloc(vBndTriCount.size(), sizeof(size_t));   refineBndTri2Count = (size_t*)calloc(vBndTri2Count.size(), sizeof(size_t));
    refineBndBeamCount = (size_t*)calloc(vBndBeamCount.size(), sizeof(size_t)); refineBndBeam2Count = (size_t*)calloc(vBndBeam2Count.size(), sizeof(size_t));

    refineCommQuadCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t)); refineCommQuad2Count= (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    refineCommTriCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));   refineCommTri2Count = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));
    refineCommBeamCount = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t)); refineCommBeam2Count = (size_t*)calloc( nNumOfCommMesh2, sizeof(size_t));

    refineBndHexaNodes = (int32_t**)calloc(vBndHexaCount.size(), sizeof(int32_t*)); refineBndHexa2Nodes = (int32_t**)calloc(vBndHexa2Count.size(), sizeof(int32_t*));
    refineBndTetNodes = (int32_t**)calloc(vBndTetCount.size(), sizeof(int32_t*));   refineBndTet2Nodes = (int32_t**)calloc(vBndTet2Count.size(), sizeof(int32_t*));
    refineBndPrismNodes = (int32_t**)calloc(vBndPrismCount.size(), sizeof(int32_t*));refineBndPrism2Nodes = (int32_t**)calloc(vBndPrism2Count.size(), sizeof(int32_t*));
    refineBndQuadNodes = (int32_t**)calloc(vBndQuadCount.size(), sizeof(int32_t*)); refineBndQuad2Nodes = (int32_t**)calloc(vBndQuad2Count.size(), sizeof(int32_t*));
    refineBndTriNodes = (int32_t**)calloc(vBndTriCount.size(), sizeof(int32_t*));   refineBndTri2Nodes = (int32_t**)calloc(vBndTri2Count.size(), sizeof(int32_t*));
    refineBndBeamNodes = (int32_t**)calloc(vBndBeamCount.size(), sizeof(int32_t*)); refineBndBeam2Nodes = (int32_t**)calloc(vBndBeam2Count.size(), sizeof(int32_t*));

    refineCommQuadNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*)); refineCommQuad2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    refineCommTriNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));   refineCommTri2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    refineCommBeamNodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*)); refineCommBeam2Nodes = (int32_t**)calloc( nNumOfCommMesh2, sizeof(int32_t*));
    

    //--------------------------------------------------------------------------
    //-------------- リファイン ループ -------------------------------------------
    //--------------------------------------------------------------------------
    for(uiint i=0; i < nRefine; i++){
        if(i > 0){
            // Mesh : 次ステップのためコースグリッドElementNodesを入れ替え
            free( hexaElemNodes );   hexaElemNodes = refineHexaNodes;    hexaCount= refineHexaCount;
            free( hexa2ElemNodes );  hexa2ElemNodes= refineHexa2Nodes;   hexa2Count=refineHexa2Count;
            free( tetraElemNodes );  tetraElemNodes= refineTetNodes;     tetCount= refineTetCount;
            free( tetra2ElemNodes ); tetra2ElemNodes= refineTet2Nodes;   tet2Count= refineTet2Count;
            free( prismElemNodes );  prismElemNodes= refinePrismNodes;   prismCount= refinePrismCount;
            free( prism2ElemNodes ); prism2ElemNodes= refinePrism2Nodes; prism2Count= refinePrism2Count;
            free( quadElemNodes );   quadElemNodes= refineQuadNodes;     quadCount= refineQuadCount;
            free( quad2ElemNodes );  quad2ElemNodes= refineQuad2Nodes;   quad2Count= refineQuad2Count;
            free( triElemNodes );    triElemNodes= refineTriNodes;       triCount= refineTriCount;
            free( tri2ElemNodes );   tri2ElemNodes= refineTri2Nodes;     tri2Count= refineTri2Count;
            free( beamElemNodes );   beamElemNodes= refineBeamNodes;     beamCount= refineBeamCount;
            free( beam2ElemNodes );  beam2ElemNodes= refineBeam2Nodes;   beam2Count= refineBeam2Count;
            // Bnd
            for(uiint ibnd=0; ibnd < vBndHexaCount.size(); ibnd++){ free( bndHexaNodes[ibnd] ); bndHexaNodes[ibnd]= refineBndHexaNodes[ibnd]; bndHexaCount[ibnd]= refineBndHexaCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndHexa2Count.size(); ibnd++){ free( bndHexa2Nodes[ibnd] ); bndHexa2Nodes[ibnd]= refineBndHexa2Nodes[ibnd]; bndHexa2Count[ibnd]= refineBndHexa2Count[ibnd];}
            for(uiint ibnd=0; ibnd < vBndTetCount.size(); ibnd++){ free( bndTetNodes[ibnd] ); bndTetNodes[ibnd]= refineBndTetNodes[ibnd]; bndTetCount[ibnd]= refineBndTetCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndTet2Count.size(); ibnd++){ free( bndTet2Nodes[ibnd] ); bndTet2Nodes[ibnd]= refineBndTet2Nodes[ibnd]; bndTet2Count[ibnd]= refineBndTet2Count[ibnd];}
            for(uiint ibnd=0; ibnd < vBndPrismCount.size(); ibnd++){ free( bndPrismNodes[ibnd] ); bndPrismNodes[ibnd]= refineBndPrismNodes[ibnd]; bndPrismCount[ibnd]= refineBndPrismCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndPrism2Count.size(); ibnd++){ free( bndPrism2Nodes[ibnd] ); bndPrism2Nodes[ibnd]= refineBndPrism2Nodes[ibnd]; bndPrism2Count[ibnd]= refineBndPrism2Count[ibnd];}

            for(uiint ibnd=0; ibnd < vBndQuadCount.size(); ibnd++){  free( bndQuadNodes[ibnd] ); bndQuadNodes[ibnd]= refineBndQuadNodes[ibnd]; bndQuadCount[ibnd]= refineBndQuadCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndQuad2Count.size(); ibnd++){ free( bndQuad2Nodes[ibnd] ); bndQuad2Nodes[ibnd]= refineBndQuad2Nodes[ibnd]; bndQuad2Count[ibnd]= refineBndQuad2Count[ibnd];}
            for(uiint ibnd=0; ibnd < vBndTriCount.size(); ibnd++){ free( bndTriNodes[ibnd] ); bndTriNodes[ibnd]= refineBndTriNodes[ibnd]; bndTriCount[ibnd]= refineBndTriCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndTri2Count.size(); ibnd++){ free( bndTri2Nodes[ibnd] ); bndTri2Nodes[ibnd]= refineBndTri2Nodes[ibnd]; bndTri2Count[ibnd]= refineBndTri2Count[ibnd];}
            for(uiint ibnd=0; ibnd < vBndBeamCount.size(); ibnd++){ free( bndBeamNodes[ibnd] ); bndBeamNodes[ibnd]= refineBndBeamNodes[ibnd]; bndBeamCount[ibnd]= refineBndBeamCount[ibnd];}
            for(uiint ibnd=0; ibnd < vBndBeam2Count.size(); ibnd++){ free( bndBeam2Nodes[ibnd] ); bndBeam2Nodes[ibnd]= refineBndBeam2Nodes[ibnd]; bndBeam2Count[ibnd]= refineBndBeam2Count[ibnd];}
            // Comm : Pointはリファイン自体をしないので、入れ替えも無し.
            for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){
                free( commQuadNodes[icom] );  commQuadNodes[icom]= refineCommQuadNodes[icom];   commQuadCount[icom]= refineCommQuadCount[icom];
                free( commQuad2Nodes[icom] ); commQuad2Nodes[icom]= refineCommQuad2Nodes[icom]; commQuad2Count[icom]= refineCommQuad2Count[icom];
                free( commTriNodes[icom] );   commTriNodes[icom]= refineCommTriNodes[icom];     commTriCount[icom]= refineCommTriCount[icom];
                free( commTri2Nodes[icom] );  commTri2Nodes[icom]= refineCommTri2Nodes[icom];   commTri2Count[icom]= refineCommTri2Count[icom];
                free( commBeamNodes[icom] );  commBeamNodes[icom]= refineCommBeamNodes[icom];   commBeamCount[icom]= refineCommBeamCount[icom];
                free( commBeam2Nodes[icom] ); commBeam2Nodes[icom]= refineCommBeam2Nodes[icom]; commBeam2Count[icom]= refineCommBeam2Count[icom];
            };
        }
        //Mesh:refine時の要素数を取得 => refineNodes領域確保(calloc)
        refineHexaCount = rcapRefineElement( hexaCount, RCAP_HEXAHEDRON, hexaElemNodes, NULL);    refineHexaNodes = (int32_t*)calloc( nNumNodeHexa*refineHexaCount, sizeof(int32_t));//Hexa
        refineHexa2Count= rcapRefineElement( hexa2Count, RCAP_HEXAHEDRON2, hexa2ElemNodes,NULL);  refineHexa2Nodes = (int32_t*)calloc( nNumNodeHexa2*refineHexa2Count, sizeof(int32_t));//Hexa2
        refineTetCount = rcapRefineElement( tetCount, RCAP_TETRAHEDRON, tetraElemNodes, NULL);    refineTetNodes = (int32_t*)calloc( nNumNodeTetra*refineTetCount, sizeof(int32_t));//Tetra
        refineTet2Count = rcapRefineElement( tet2Count, RCAP_TETRAHEDRON2, tetra2ElemNodes, NULL);refineTet2Nodes = (int32_t*)calloc( nNumNodeTetra2*refineTet2Count, sizeof(int32_t));//Tetra2 
        refinePrismCount= rcapRefineElement( prismCount, RCAP_WEDGE, prismElemNodes, NULL);       refinePrismNodes= (int32_t*)calloc( nNumNodePrism*refinePrismCount, sizeof(int32_t));//Prism 
        refinePrism2Count= rcapRefineElement( prism2Count, RCAP_WEDGE2, prism2ElemNodes, NULL);   refinePrism2Nodes= (int32_t*)calloc( nNumNodePrism2*refinePrism2Count, sizeof(int32_t));//Prism2 
        refineQuadCount= rcapRefineElement( quadCount, RCAP_QUAD, quadElemNodes, NULL);           refineQuadNodes= (int32_t*)calloc( nNumNodeQuad*refineQuadCount, sizeof(int32_t));//Quad 
        refineQuad2Count= rcapRefineElement( quad2Count, RCAP_QUAD2, quad2ElemNodes, NULL);       refineQuad2Nodes= (int32_t*)calloc( nNumNodeQuad2*refineQuad2Count, sizeof(int32_t));//Quad2 
        refineTriCount= rcapRefineElement( triCount, RCAP_TRIANGLE, triElemNodes, NULL);          refineTriNodes= (int32_t*)calloc( nNumNodeTri*refineTriCount, sizeof(int32_t));//Tri 
        refineTri2Count= rcapRefineElement( tri2Count, RCAP_TRIANGLE2, tri2ElemNodes, NULL);      refineTri2Nodes= (int32_t*)calloc( nNumNodeTri2*refineTri2Count, sizeof(int32_t));//Tri2 
        refineBeamCount= rcapRefineElement( beamCount, RCAP_SEGMENT, beamElemNodes, NULL);        refineBeamNodes= (int32_t*)calloc( nNumNodeBeam*refineBeamCount, sizeof(int32_t));//Beam 
        refineBeam2Count= rcapRefineElement( beam2Count, RCAP_SEGMENT2, beam2ElemNodes, NULL);    refineBeam2Nodes= (int32_t*)calloc( nNumNodeBeam2*refineBeam2Count, sizeof(int32_t));//Beam2 
        //Bnd:refine時の要素数を取得 => refineNodesの領域確保(calloc)
        for(uiint ibnd=0; ibnd < vBndHexaCount.size(); ibnd++){ refineBndHexaCount[ibnd]= rcapRefineElement( bndHexaCount[ibnd], RCAP_HEXAHEDRON, bndHexaNodes[ibnd], NULL); refineBndHexaNodes[ibnd]=(int32_t*)calloc( nNumNodeHexa*refineBndHexaCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndHexa2Count.size(); ibnd++){ refineBndHexa2Count[ibnd]= rcapRefineElement( bndHexa2Count[ibnd], RCAP_HEXAHEDRON2, bndHexa2Nodes[ibnd], NULL); refineBndHexa2Nodes[ibnd]=(int32_t*)calloc( nNumNodeHexa2*refineBndHexa2Count[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndTetCount.size(); ibnd++){ refineBndTetCount[ibnd]= rcapRefineElement( bndTetCount[ibnd], RCAP_TETRAHEDRON, bndTetNodes[ibnd], NULL); refineBndTetNodes[ibnd]=(int32_t*)calloc( nNumNodeTetra*refineBndTetCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndTet2Count.size(); ibnd++){ refineBndTet2Count[ibnd]= rcapRefineElement( bndTet2Count[ibnd], RCAP_TETRAHEDRON2, bndTet2Nodes[ibnd], NULL); refineBndTet2Nodes[ibnd]=(int32_t*)calloc( nNumNodeTetra2*refineBndTet2Count[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndPrismCount.size(); ibnd++){ refineBndPrismCount[ibnd]= rcapRefineElement( bndPrismCount[ibnd], RCAP_WEDGE, bndPrismNodes[ibnd], NULL); refineBndPrismNodes[ibnd]=(int32_t*)calloc( nNumNodePrism*refineBndPrismCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndPrism2Count.size(); ibnd++){ refineBndPrism2Count[ibnd]= rcapRefineElement( bndPrism2Count[ibnd], RCAP_WEDGE2, bndPrism2Nodes[ibnd], NULL); refineBndPrism2Nodes[ibnd]=(int32_t*)calloc( nNumNodePrism2*refineBndPrism2Count[ibnd], sizeof(int32_t));}

        for(uiint ibnd=0; ibnd < vBndQuadCount.size(); ibnd++){ refineBndQuadCount[ibnd]= rcapRefineElement( bndQuadCount[ibnd], RCAP_QUAD, bndQuadNodes[ibnd], NULL);    refineBndQuadNodes[ibnd]=(int32_t*)calloc( nNumNodeQuad*refineBndQuadCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndQuad2Count.size(); ibnd++){  refineBndQuad2Count[ibnd]= rcapRefineElement( bndQuad2Count[ibnd], RCAP_QUAD2, bndQuad2Nodes[ibnd], NULL); refineBndQuad2Nodes[ibnd] = (int32_t*)calloc( nNumNodeQuad2*refineBndQuad2Count[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndTriCount.size(); ibnd++){ refineBndTriCount[ibnd]= rcapRefineElement( bndTriCount[ibnd], RCAP_TRIANGLE, bndTriNodes[ibnd], NULL); refineBndTriNodes[ibnd] = (int32_t*)calloc( nNumNodeTri*refineBndTriCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndTri2Count.size(); ibnd++){  refineBndTri2Count[ibnd]= rcapRefineElement( bndTri2Count[ibnd], RCAP_TRIANGLE2, bndTri2Nodes[ibnd], NULL); refineBndTri2Nodes[ibnd] = (int32_t*)calloc( nNumNodeTri2*refineBndTri2Count[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndBeamCount.size(); ibnd++){  refineBndBeamCount[ibnd]= rcapRefineElement( bndBeamCount[ibnd], RCAP_SEGMENT, bndBeamNodes[ibnd], NULL); refineBndBeamNodes[ibnd] = (int32_t*)calloc( nNumNodeBeam*refineBndBeamCount[ibnd], sizeof(int32_t));}
        for(uiint ibnd=0; ibnd < vBndBeam2Count.size(); ibnd++){  refineBndBeam2Count[ibnd]= rcapRefineElement( bndBeam2Count[ibnd], RCAP_SEGMENT2, bndBeam2Nodes[ibnd], NULL); refineBndBeam2Nodes[ibnd] = (int32_t*)calloc( nNumNodeBeam2*refineBndBeam2Count[ibnd], sizeof(int32_t));}
        //Comm:refine時の要素数を取得 => refineNodesの領域確保(calloc) : Pointはリファイン自体をしない
        for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
            refineCommQuadCount[icom]= rcapRefineElement( commQuadCount[icom], RCAP_QUAD, commQuadNodes[icom], NULL);   refineCommQuadNodes[icom] = (int32_t*)calloc( nNumNodeQuad*refineCommQuadCount[icom], sizeof(int32_t));
            refineCommQuad2Count[icom]= rcapRefineElement( commQuad2Count[icom], RCAP_QUAD2, commQuad2Nodes[icom], NULL); refineCommQuad2Nodes[icom] = (int32_t*)calloc( nNumNodeQuad2*refineCommQuad2Count[icom], sizeof(int32_t));
            refineCommTriCount[icom]= rcapRefineElement( commTriCount[icom], RCAP_TRIANGLE, commTriNodes[icom], NULL);  refineCommTriNodes[icom] = (int32_t*)calloc( nNumNodeTri*refineCommTriCount[icom], sizeof(int32_t));
            refineCommTri2Count[icom]= rcapRefineElement( commTri2Count[icom], RCAP_TRIANGLE2, commTri2Nodes[icom], NULL); refineCommTri2Nodes[icom] = (int32_t*)calloc( nNumNodeTri2*refineCommTri2Count[icom], sizeof(int32_t));
            refineCommBeamCount[icom]= rcapRefineElement( commBeamCount[icom], RCAP_SEGMENT, commBeamNodes[icom], NULL); refineCommBeamNodes[icom] = (int32_t*)calloc( nNumNodeBeam*refineCommBeamCount[icom], sizeof(int32_t));
            refineCommBeam2Count[icom]= rcapRefineElement( commBeam2Count[icom], RCAP_SEGMENT2, commBeam2Nodes[icom], NULL);refineCommBeam2Nodes[icom] = (int32_t*)calloc( nNumNodeBeam2*refineCommBeam2Count[icom], sizeof(int32_t));
        };
        

        // Element Refine
        refineHexaCount = rcapRefineElement( hexaCount, RCAP_HEXAHEDRON, hexaElemNodes, refineHexaNodes);
        refineHexa2Count = rcapRefineElement( hexa2Count, RCAP_HEXAHEDRON2, hexa2ElemNodes, refineHexa2Nodes);
        refineTetCount = rcapRefineElement( tetCount, RCAP_TETRAHEDRON, tetraElemNodes, refineTetNodes);
        refineTet2Count = rcapRefineElement( tet2Count, RCAP_TETRAHEDRON2, tetra2ElemNodes, refineTet2Nodes);
        refinePrismCount = rcapRefineElement( prismCount, RCAP_WEDGE, prismElemNodes, refinePrismNodes);
        refinePrism2Count = rcapRefineElement( prism2Count, RCAP_WEDGE2, prism2ElemNodes, refinePrism2Nodes);
        refineQuadCount = rcapRefineElement( quadCount, RCAP_QUAD, quadElemNodes, refineQuadNodes);
        refineQuad2Count = rcapRefineElement( quad2Count, RCAP_QUAD2, quad2ElemNodes, refineQuad2Nodes);
        refineTriCount = rcapRefineElement( triCount, RCAP_TRIANGLE, triElemNodes, refineTriNodes);
        refineTri2Count = rcapRefineElement( tri2Count, RCAP_TRIANGLE2, tri2ElemNodes, refineTri2Nodes);
        refineBeamCount = rcapRefineElement( beamCount, RCAP_SEGMENT, beamElemNodes, refineBeamNodes);
        refineBeam2Count= rcapRefineElement( beam2Count, RCAP_SEGMENT2, beam2ElemNodes, refineBeam2Nodes);
        
        ////cout << "i:" << i << " REVOCAP refineHexaCount : " << refineHexaCount
        ////                  << " refineNodeCount : " << rcapGetNodeCount() << " rank:" << mpMPI->getRank() << endl;

        //Bnd refine
        for(uiint ibnd=0; ibnd < vBndHexaCount.size(); ibnd++)  refineBndHexaCount[ibnd] = rcapRefineElement( bndHexaCount[ibnd], RCAP_HEXAHEDRON, bndHexaNodes[ibnd], refineBndHexaNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndHexa2Count.size(); ibnd++)  refineBndHexa2Count[ibnd] = rcapRefineElement( bndHexa2Count[ibnd], RCAP_HEXAHEDRON2, bndHexa2Nodes[ibnd], refineBndHexa2Nodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndTetCount.size(); ibnd++)  refineBndTetCount[ibnd] = rcapRefineElement( bndTetCount[ibnd], RCAP_TETRAHEDRON, bndTetNodes[ibnd], refineBndTetNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndTet2Count.size(); ibnd++)  refineBndTet2Count[ibnd] = rcapRefineElement( bndTet2Count[ibnd], RCAP_TETRAHEDRON2, bndTet2Nodes[ibnd], refineBndTet2Nodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndPrismCount.size(); ibnd++)  refineBndPrismCount[ibnd] = rcapRefineElement( bndPrismCount[ibnd], RCAP_WEDGE, bndPrismNodes[ibnd], refineBndPrismNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndPrism2Count.size(); ibnd++)  refineBndPrism2Count[ibnd] = rcapRefineElement( bndPrism2Count[ibnd], RCAP_WEDGE2, bndPrism2Nodes[ibnd], refineBndPrism2Nodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndQuadCount.size(); ibnd++)  refineBndQuadCount[ibnd] = rcapRefineElement( bndQuadCount[ibnd], RCAP_QUAD, bndQuadNodes[ibnd], refineBndQuadNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndQuad2Count.size(); ibnd++) refineBndQuad2Count[ibnd] = rcapRefineElement( bndQuad2Count[ibnd], RCAP_QUAD2, bndQuad2Nodes[ibnd], refineBndQuad2Nodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndTriCount.size(); ibnd++)   refineBndTriCount[ibnd] = rcapRefineElement( bndTriCount[ibnd], RCAP_TRIANGLE, bndTriNodes[ibnd], refineBndTriNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndTri2Count.size(); ibnd++)  refineBndTri2Count[ibnd] = rcapRefineElement( bndTri2Count[ibnd], RCAP_TRIANGLE2, bndTri2Nodes[ibnd], refineBndTri2Nodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndBeamCount.size(); ibnd++)  refineBndBeamCount[ibnd] = rcapRefineElement( bndBeamCount[ibnd], RCAP_SEGMENT, bndBeamNodes[ibnd], refineBndBeamNodes[ibnd]);
        for(uiint ibnd=0; ibnd < vBndBeam2Count.size(); ibnd++) refineBndBeam2Count[ibnd] = rcapRefineElement( bndBeam2Count[ibnd], RCAP_SEGMENT2, bndBeam2Nodes[ibnd], refineBndBeam2Nodes[ibnd]);
        
        //Comm refine : Pointはリファイン自体をしない.
        for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
            refineCommQuadCount[icom] = rcapRefineElement( commQuadCount[icom], RCAP_QUAD, commQuadNodes[icom], refineCommQuadNodes[icom]);
            refineCommQuad2Count[icom] = rcapRefineElement( commQuad2Count[icom], RCAP_QUAD2, commQuad2Nodes[icom], refineCommQuad2Nodes[icom]);
            refineCommTriCount[icom] = rcapRefineElement( commTriCount[icom], RCAP_TRIANGLE, commTriNodes[icom], refineCommTriNodes[icom]);
            refineCommTri2Count[icom] = rcapRefineElement( commTri2Count[icom], RCAP_TRIANGLE2, commTri2Nodes[icom], refineCommTri2Nodes[icom]);
            refineCommBeamCount[icom] = rcapRefineElement( commBeamCount[icom], RCAP_SEGMENT, commBeamNodes[icom], refineCommBeamNodes[icom]);
            refineCommBeam2Count[icom] = rcapRefineElement( commBeam2Count[icom], RCAP_SEGMENT2, commBeam2Nodes[icom], refineCommBeam2Nodes[icom]);
        
            ////cout << "refineCommCount icom:" << icom << " tri_count:" << refineCommTriCount[icom] << " rank:" << mpMPI->getRank() << endl;
        };

        rcapCommit();//----------------------------- 唯一のCommit

        //---
        // コースグリッドMesh
        //---
        uiint nMeshID= mpMesh->getMeshID();
        uiint crsNumOfNode= mpMesh->getNumOfNode();
        CBNodeMeshGrp* pBNodeMeshGrp= mpMesh->getBNodeMeshGrp();//コースグリッドのBNodeMesh(点Gr)
        CNode *pCrsLastNode = mpMesh->getNodeIX(crsNumOfNode-1);
        uiint nCrsMaxNodeID= pCrsLastNode->getID();
        //---
        // ファイングリッドMesh
        //---
        //　AssyModel
        CAssyModel *progAssy;
        progAssy= mpGMGModel->getAssyModel(i+1);//------- Level:i+1 --- FineGrid AssyModel

        ////cout << "MW::RevocapRefine ----------------- F0 rank:" << mpMPI->getRank() << endl;

        progAssy->resizeMesh(1);
        //　Mesh
        CMesh *progMesh;
        progMesh = new CMesh();//---------------------------------------new Mesh
        progMesh->setMeshID(nMeshID);
        progMesh->setMGLevel(i+1);//--------- Level:i+1
        progMesh->setMaxMGLevel(nRefine);
        progMesh->setSolutionType(SolutionType::FEM);
        progMesh->setProp(0);

        //--
        // BNodeMeshを上位グリッドにセット
        //--
        if(pBNodeMeshGrp) progMesh->setBNodeMeshGrp(pBNodeMeshGrp);//-------------コースグリッドと同じBNodeMeshをセット(点Gr)
        
        ////cout << "MW::RevocapRefine ----------------- F1 rank:" << mpMPI->getRank() << endl;
        mpFactory->setupBucketMesh(i+1, nMeshID, nMeshID);
        ////cout << "MW::RevocapRefine ----------------- F2 rank:" << mpMPI->getRank() << endl;
        progAssy->setMesh(progMesh, 0);//-------------------------AssyModelにセット
        
        ////cout << "MW::RevocapRefine ----------------- G0 rank:" << mpMPI->getRank() << endl;
        //----
        // 新Node : ファイングリッドMesh
        //----
        if(i > 0) free( resultCoord );
        refineNodeCount = rcapGetNodeCount();
        resultCoord = (float64_t*)calloc( refineNodeCount*3, sizeof(float64_t));
        rcapGetNodeSeq64( refineNodeCount, nodeOffset, resultCoord );
        
        // Node生成
        Rcap_NodeGene( i+1,  nRefine, crsNumOfNode, refineNodeCount, resultCoord, mpMesh, progMesh);

        progMesh->setupNumOfNode();
        
        // Bucket for Node
        uiint nProgMaxNodeID = nCrsMaxNodeID + (refineNodeCount-crsNumOfNode);
        progMesh->initBucketNode(nProgMaxNodeID, 0);//maxID = コースグリッドMaxID + (増加個数), minID=0
        progMesh->setupBucketNode(); 

        ////cout << "MW::RevocapRefine ----------------- G1 rank:" << mpMPI->getRank() << endl;

        //----
        // 新Element : ファイングリッドMesh
        //----
        uiint refineTotalCount = refineHexaCount + refineHexa2Count + refineTetCount + refineTet2Count
                               + refinePrismCount + refinePrism2Count + refineQuadCount + refineQuad2Count
                               + refineTriCount + refineTri2Count + refineBeamCount + refineBeam2Count;
        // Element生成
        uiint nIDCount(0);// 新Element ID : 全ての種類の新ID
        Rcap_ElemGene( RCAP_HEXAHEDRON, nIDCount, refineHexaCount, refineHexaNodes, progMesh );// Hexa
        Rcap_ElemGene( RCAP_HEXAHEDRON2, nIDCount, refineHexa2Count, refineHexa2Nodes, progMesh );// Hexa2
        Rcap_ElemGene( RCAP_TETRAHEDRON, nIDCount, refineTetCount, refineTetNodes, progMesh );// Tetra
        Rcap_ElemGene( RCAP_TETRAHEDRON, nIDCount, refineTet2Count, refineTet2Nodes, progMesh );// Tetra2
        Rcap_ElemGene( RCAP_WEDGE, nIDCount, refinePrismCount, refinePrismNodes, progMesh );// Prism
        Rcap_ElemGene( RCAP_WEDGE2, nIDCount, refinePrism2Count, refinePrism2Nodes, progMesh );// Prism2
        Rcap_ElemGene( RCAP_QUAD, nIDCount, refineQuadCount, refineQuadNodes, progMesh );// Quad
        Rcap_ElemGene( RCAP_QUAD2, nIDCount, refineQuad2Count, refineQuad2Nodes, progMesh );// Quad2
        Rcap_ElemGene( RCAP_TRIANGLE, nIDCount, refineTriCount, refineTriNodes, progMesh );// Triangle
        Rcap_ElemGene( RCAP_TRIANGLE2, nIDCount, refineTri2Count, refineTri2Nodes, progMesh );// Triangle2
        Rcap_ElemGene( RCAP_SEGMENT, nIDCount, refineBeamCount, refineBeamNodes, progMesh );// Beam
        Rcap_ElemGene( RCAP_SEGMENT2, nIDCount, refineBeam2Count, refineBeam2Nodes, progMesh );// Beam2
        
        ////cout << "MW::RevocapRefine ----------------- G2 rank:" << mpMPI->getRank() << endl;
        
        progMesh->setupNumOfElement();
        // Bucket for Element
        progMesh->initBucketElement(refineTotalCount-1, 0);//maxID=要素数-1, minID=0
        progMesh->setupBucketElement();
        
        ////cout << "MW::RevocapRefine ----------------- G3 rank:" << mpMPI->getRank() << endl;

        // ----
        // AggregeteElementの設定
        // ----
        progMesh->resizeAggregate( progMesh->getNumOfNode() );
        mpFactory->GeneAggregate(i+1, 0, progMesh->getNumOfNode() );//Level=i+1, Mesh=0
        progMesh->setupAggregate(i+1);//ファインLevel=i+1

        ////cout << "MW::RevocapRefine ----------------- G4 rank:" << mpMPI->getRank() << endl;

        //---
        // 新Bnd:ファイングリッド
        //---
        Rcap_SetBVolMesh( i+1, nRefine, vBndHexaBndID, vBndHexaBType, vBndHexaName, vBndHexaDOF, vvBndHexaDOF, progMesh, mpMesh);
        Rcap_SetBVolMesh( i+1, nRefine, vBndHexa2BndID, vBndHexa2BType, vBndHexa2Name, vBndHexa2DOF, vvBndHexa2DOF, progMesh, mpMesh);
        Rcap_SetBVolMesh( i+1, nRefine, vBndTetBndID, vBndTetBType, vBndTetName, vBndTetDOF, vvBndTetDOF, progMesh, mpMesh);
        Rcap_SetBVolMesh( i+1, nRefine, vBndTet2BndID, vBndTet2BType, vBndTet2Name, vBndTet2DOF, vvBndTet2DOF, progMesh, mpMesh);
        Rcap_SetBVolMesh( i+1, nRefine, vBndPrismBndID, vBndPrismBType, vBndPrismName, vBndPrismDOF, vvBndPrismDOF, progMesh, mpMesh);
        Rcap_SetBVolMesh( i+1, nRefine, vBndPrism2BndID, vBndPrism2BType, vBndPrism2Name, vBndPrism2DOF, vvBndPrism2DOF, progMesh, mpMesh);

        Rcap_SetBFaceMesh( i+1, nRefine, vBndQuadBndID, vBndQuadBType, vBndQuadName, vBndQuadDOF, vvBndQuadDOF, progMesh, mpMesh);
        Rcap_SetBFaceMesh( i+1, nRefine, vBndQuad2BndID, vBndQuad2BType, vBndQuad2Name, vBndQuad2DOF, vvBndQuad2DOF, progMesh, mpMesh);
        Rcap_SetBFaceMesh( i+1, nRefine, vBndTriBndID, vBndTriBType, vBndTriName, vBndTriDOF, vvBndTriDOF, progMesh, mpMesh);
        Rcap_SetBFaceMesh( i+1, nRefine, vBndTri2BndID, vBndTri2BType, vBndTri2Name, vBndTri2DOF, vvBndTri2DOF, progMesh, mpMesh);

        Rcap_SetBEdgeMesh( i+1, nRefine, vBndBeamBndID, vBndBeamBType, vBndBeamName, vBndBeamDOF, vvBndBeamDOF, progMesh, mpMesh);
        Rcap_SetBEdgeMesh( i+1, nRefine, vBndBeam2BndID, vBndBeam2BType, vBndBeam2Name, vBndBeam2DOF, vvBndBeam2DOF, progMesh, mpMesh);

        ////cout << "MW::RevocapRefine ----------------- G5 rank:" << mpMPI->getRank() << endl;

        //BNode生成 : refineBndNodesの節点番号をソート・マージしてから"BNode生成"
        // 1 節点番号をソート・マージ
        vvuint vvHexaNodeIndex = Rcap_NodeIX_SortMerge(vBndHexaBndID.size(), refineBndHexaCount, 8, refineBndHexaNodes);
        vvuint vvHexa2NodeIndex = Rcap_NodeIX_SortMerge(vBndHexa2BndID.size(), refineBndHexa2Count, 20, refineBndHexa2Nodes);
        vvuint vvTetNodeIndex = Rcap_NodeIX_SortMerge(vBndTetBndID.size(), refineBndTetCount, 4, refineBndTetNodes);
        vvuint vvTet2NodeIndex = Rcap_NodeIX_SortMerge(vBndTet2BndID.size(), refineBndTet2Count, 10, refineBndTet2Nodes);
        vvuint vvPrismNodeIndex = Rcap_NodeIX_SortMerge(vBndPrismBndID.size(), refineBndPrismCount, 6, refineBndPrismNodes);
        vvuint vvPrism2NodeIndex = Rcap_NodeIX_SortMerge(vBndPrism2BndID.size(), refineBndPrism2Count, 15, refineBndPrism2Nodes);
        
        vvuint vvQuadNodeIndex  = Rcap_NodeIX_SortMerge(vBndQuadBndID.size(), refineBndQuadCount, 4, refineBndQuadNodes);
        vvuint vvQuad2NodeIndex = Rcap_NodeIX_SortMerge(vBndQuad2BndID.size(), refineBndQuad2Count, 8, refineBndQuad2Nodes);
        vvuint vvTriNodeIndex = Rcap_NodeIX_SortMerge(vBndTriBndID.size(), refineBndTriCount, 3, refineBndTriNodes);
        vvuint vvTri2NodeIndex = Rcap_NodeIX_SortMerge(vBndTri2BndID.size(), refineBndTri2Count, 6, refineBndTri2Nodes);
        vvuint vvBeamNodeIndex = Rcap_NodeIX_SortMerge(vBndBeamBndID.size(), refineBndBeamCount, 2, refineBndBeamNodes);
        vvuint vvBeam2NodeIndex = Rcap_NodeIX_SortMerge(vBndBeam2BndID.size(), refineBndBeam2Count, 3, refineBndBeam2Nodes);

        // 2 ソートマージされた節点番号に沿って ----> "BNode"生成
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndHexaBndID, mpMesh, progMesh, vvHexaNodeIndex);
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndHexa2BndID, mpMesh, progMesh, vvHexa2NodeIndex);
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndTetBndID, mpMesh, progMesh, vvTetNodeIndex);
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndTet2BndID, mpMesh, progMesh, vvTet2NodeIndex);
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndPrismBndID, mpMesh, progMesh, vvPrismNodeIndex);
        Rcap_BNodeGene_VolMesh( i+1, nRefine, vBndPrism2BndID, mpMesh, progMesh, vvPrism2NodeIndex);

        Rcap_BNodeGene_FaceMesh( i+1, nRefine, vBndQuadBndID, mpMesh, progMesh, vvQuadNodeIndex);
        Rcap_BNodeGene_FaceMesh( i+1, nRefine, vBndQuad2BndID, mpMesh, progMesh, vvQuad2NodeIndex);
        Rcap_BNodeGene_FaceMesh( i+1, nRefine, vBndTriBndID, mpMesh, progMesh, vvTriNodeIndex);
        Rcap_BNodeGene_FaceMesh( i+1, nRefine, vBndTri2BndID, mpMesh, progMesh, vvTri2NodeIndex);

        Rcap_BNodeGene_EdgeMesh( i+1, nRefine, vBndBeamBndID, mpMesh, progMesh, vvBeamNodeIndex);
        Rcap_BNodeGene_EdgeMesh( i+1, nRefine, vBndBeam2BndID, mpMesh, progMesh, vvBeam2NodeIndex);

        ////cout << "MW::RevocapRefine ----------------- G6 rank:" << mpMPI->getRank() << endl;

        // 逆引き : Volume構成Node番号から、NodeIndexの配列番号(インデックス)を取得  #境界Grの節点配列のインデックス==BNodeのインデックス
        //         Face構成Node番号から、NodeIndexの配列番号(インデックス)を取得　  #境界Grの節点配列のインデックス==BNodeのインデックス
        //         BVolume生成, BFace生成, BEdge生成で利用
        vector<map<uiint, uiint> > vmID2IndexHexa, vmID2IndexHexa2, vmID2IndexTet, vmID2IndexTet2, vmID2IndexPrism, vmID2IndexPrism2;
        vector<map<uiint, uiint> > vmID2IndexQuad, vmID2IndexQuad2, vmID2IndexTri, vmID2IndexTri2, vmID2IndexBeam, vmID2IndexBeam2;
        Rcap_NodeNum2Index(vmID2IndexHexa, vvHexaNodeIndex);  Rcap_NodeNum2Index(vmID2IndexHexa2, vvHexa2NodeIndex);
        Rcap_NodeNum2Index(vmID2IndexTet, vvTetNodeIndex);  Rcap_NodeNum2Index(vmID2IndexTet2, vvTet2NodeIndex);
        Rcap_NodeNum2Index(vmID2IndexPrism, vvPrismNodeIndex);Rcap_NodeNum2Index(vmID2IndexPrism2, vvPrism2NodeIndex);
        Rcap_NodeNum2Index(vmID2IndexQuad, vvQuadNodeIndex);  Rcap_NodeNum2Index(vmID2IndexQuad2, vvQuad2NodeIndex);
        Rcap_NodeNum2Index(vmID2IndexTri,  vvTriNodeIndex);   Rcap_NodeNum2Index(vmID2IndexTri2,  vvTri2NodeIndex);
        Rcap_NodeNum2Index(vmID2IndexBeam, vvBeamNodeIndex);  Rcap_NodeNum2Index(vmID2IndexBeam2, vvBeam2NodeIndex);
        
        //BVol生成, BFace生成, BEdge生成
        Rcap_BVolGene( i+1, vBndHexaBndID, vmID2IndexHexa, refineBndHexaCount, ElementType::Hexa, nNumNodeHexa, refineBndHexaNodes, progMesh, mpMesh);
        Rcap_BVolGene( i+1, vBndHexa2BndID, vmID2IndexHexa2, refineBndHexa2Count, ElementType::Hexa2, nNumNodeHexa2, refineBndHexa2Nodes, progMesh, mpMesh);
        Rcap_BVolGene( i+1, vBndTetBndID, vmID2IndexTet, refineBndTetCount, ElementType::Tetra, nNumNodeTetra, refineBndTetNodes, progMesh, mpMesh);
        Rcap_BVolGene( i+1, vBndTet2BndID, vmID2IndexTet2, refineBndTet2Count, ElementType::Tetra2, nNumNodeTetra2, refineBndTet2Nodes, progMesh, mpMesh);
        Rcap_BVolGene( i+1, vBndPrismBndID, vmID2IndexPrism, refineBndPrismCount, ElementType::Prism, nNumNodePrism, refineBndPrismNodes, progMesh, mpMesh);
        Rcap_BVolGene( i+1, vBndPrism2BndID, vmID2IndexPrism2, refineBndPrism2Count, ElementType::Prism2, nNumNodePrism2, refineBndPrism2Nodes, progMesh, mpMesh);
        
        Rcap_BFaceGene( i+1, vBndQuadBndID, vmID2IndexQuad, refineBndQuadCount, ElementType::Quad, nNumNodeQuad, refineBndQuadNodes, progMesh, mpMesh);
        Rcap_BFaceGene( i+1, vBndQuad2BndID,vmID2IndexQuad2,refineBndQuad2Count,ElementType::Quad2, nNumNodeQuad2,refineBndQuad2Nodes, progMesh, mpMesh);
        Rcap_BFaceGene( i+1, vBndTriBndID,  vmID2IndexTri,  refineBndTriCount,  ElementType::Triangle, nNumNodeTri, refineBndTriNodes, progMesh, mpMesh);
        Rcap_BFaceGene( i+1, vBndTri2BndID, vmID2IndexTri2, refineBndTri2Count, ElementType::Triangle2, nNumNodeTri2, refineBndTri2Nodes, progMesh, mpMesh);

        Rcap_BEdgeGene( i+1, vBndBeamBndID,  vmID2IndexBeam,  refineBndBeamCount, ElementType::Beam, nNumNodeBeam, refineBndBeamNodes, progMesh, mpMesh);
        Rcap_BEdgeGene( i+1, vBndBeam2BndID, vmID2IndexBeam2, refineBndBeam2Count, ElementType::Beam2, nNumNodeBeam2, refineBndBeam2Nodes, progMesh, mpMesh);
        
        //debug:BoundaryFace確認
        Rcap_BFaceMeshDebug(vBndQuadBndID, progMesh);

        ////cout << "MW::RevocapRefine ----------------- G7 rank:" << mpMPI->getRank() << endl;

        //--------------
        //    新Com
        //--------------
        // CommNode生成の前処理
        // 1. refineNodes を 便宜上vvuintに移し替え
        vvuint vvCommQuadRefineNodes, vvCommQuad2RefineNodes, vvCommTriRefineNodes, vvCommTri2RefineNodes;
        vvuint vvCommBeamRefineNodes, vvCommBeam2RefineNodes, vvCommPointRefineNodes;
        //---
        Rcap_NodesPoint2vuint( mCommID, 4, refineCommQuadCount, refineCommQuadNodes, vvCommQuadRefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 8, refineCommQuad2Count, refineCommQuad2Nodes, vvCommQuad2RefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 3, refineCommTriCount, refineCommTriNodes, vvCommTriRefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 6, refineCommTri2Count, refineCommTri2Nodes, vvCommTri2RefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 2, refineCommBeamCount, refineCommBeamNodes, vvCommBeamRefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 3, refineCommBeam2Count, refineCommBeam2Nodes, vvCommBeam2RefineNodes);
        Rcap_NodesPoint2vuint( mCommID, 1, commPointCount, commPointNodes, vvCommPointRefineNodes);

        ////cout << "MW::RevocapRefine ----------------- G8 rank:" << mpMPI->getRank() << endl;

        // 2. CommNodeを各階層で独立して生成(コースグリッドCommNodeは不使用)
        //    # vmaElemNIndex, maNNum2CommNNum は、Face生成で利用するためCommNode生成時にデータ作成
        //---
        // vvNum4EachType:コースグリッドでの面、辺の個数
        //---
        vvuint vvNum4EachType;
        Rcap_CommFaceCount( vvNum4EachType, mCommID, commQuadCount, commQuad2Count, commTriCount, commTri2Count, commBeamCount, commBeam2Count, commPointCount);
        //--
        // ファイングリッドにCommMesh2を生成
        //--
        Rcap_SetCommMesh2( i+1, mCommID, vvNum4EachType, mpMesh, progMesh);

        ////cout << "MW::RevocapRefine ----------------- G9 rank:" << mpMPI->getRank() << endl;

        //---
        // 全て足し合わせたRefineNodes
        //---
        vvuint vvCommRefineNodes;
        Rcap_SumRefineNodes( vvCommRefineNodes, vvCommQuadRefineNodes, vvCommQuad2RefineNodes, vvCommTriRefineNodes, vvCommTri2RefineNodes, vvCommBeamRefineNodes, vvCommBeam2RefineNodes, vvCommPointRefineNodes);
        //--
        // CommNode生成
        //--
        vvvuint vvvElemNum;
        vector<map<uiint, vuint> > vmaElemNIndex;
        vector<map<uiint,uiint> > maNNum2CommNNum;
        Rcap_CommNodeGene( mCommID, vvNum4EachType, vvvElemNum, vmaElemNIndex, progMesh, vvCommRefineNodes, maNNum2CommNNum );

        ////cout << "MW::RevocapRefine ----------------- G10 rank:" << mpMPI->getRank() << endl;

        // 3. CommFace生成 : mapデータを利用　Node番号からCommNode番号が取得できる
        //    # CommFaceは計算では使用しない、デバッグ用途( Rcapの場合 )
        //
        Rcap_CommFaceGene( i+1, mCommID,  vvvElemNum, vmaElemNIndex, vvCommRefineNodes, progMesh, maNNum2CommNNum);


        ////cout << "MW::RevocapRefine ----------------- G11 rank:" << mpMPI->getRank() << endl;

        //---
        // カレントMeshの入れ替え
        //---
        mpMesh=progMesh;
        mpAssy=progAssy;

        ////cout << "progMesh NumOfNode:    " << progMesh->getNumOfNode() << endl;
        ////cout << "progMesh NumOfElement: " << progMesh->getNumOfElement() << endl;

        rcapClearRefiner();//--------------------- REVOCAP Clear

    };//=============================================( i < nRefine ) ループ　終端

    
    // ----
    // 境界条件(Dirichlet,Neumann)マーキング：MG境界処理
    // ----
    mpFactory->setupBNodeMarking();

    // ----
    // Rank大の通信Nodeマーキング:並列計算での荷重処理
    // ----
    mpFactory->setupLargeRankCommNode_Marking();

    ////cout << "MW::RevocapRefine ----------------- Refine_End rank:" << mpMPI->getRank() << endl;
    
    rcapTermRefiner();//終了
    // ----
    // REVOCAP関連メモリー解放
    // ----
    free( coords );
    ////cout << "MW::RevocapRefine ----------------- H0" << endl;
    if(globalIDs) free( globalIDs );
    ////cout << "MW::RevocapRefine ----------------- H1" << endl;
    if(localIDs) free( localIDs );
    ////cout << "MW::RevocapRefine ----------------- H2" << endl;
    
    // Mesh free
    free( hexaElemNodes );  free( hexa2ElemNodes );
    free( tetraElemNodes ); free( tetra2ElemNodes );
    free( prismElemNodes ); free( prism2ElemNodes );
    free( quadElemNodes );  free( quad2ElemNodes );
    free( triElemNodes );   free( tri2ElemNodes );
    free( beamElemNodes );  free( beam2ElemNodes);
    ////cout << "MW::RevocapRefine ----------------- H3" << endl;
    if(nRefine > 0){
        free( refineHexaNodes ); free( refineHexa2Nodes );
        free( refineTetNodes );  free( refineTet2Nodes );
        free( refinePrismNodes );free( refinePrism2Nodes );
        free( refineQuadNodes ); free( refineQuad2Nodes );
        free( refineTriNodes );  free( refineTri2Nodes );
        free( refineBeamNodes ); free( refineBeam2Nodes );
        ////cout << "MW::RevocapRefine ----------------- H4" << endl;
        free( resultCoord );
        ////cout << "MW::RevocapRefine ----------------- H5" << endl;
    }

    // Bnd free
    free( bndHexaCount );  if(nRefine > 0) free( refineBndHexaCount);
    free( bndHexa2Count ); if(nRefine > 0) free( refineBndHexa2Count);
    free( bndTetCount );   if(nRefine > 0) free( refineBndTetCount);
    free( bndTet2Count );  if(nRefine > 0) free( refineBndTet2Count);
    free( bndPrismCount ); if(nRefine > 0) free( refineBndPrismCount);
    free( bndPrism2Count ); if(nRefine > 0) free( refineBndPrism2Count);
    free( bndQuadCount );  if(nRefine > 0) free( refineBndQuadCount );
    free( bndQuad2Count ); if(nRefine > 0) free( refineBndQuad2Count );
    free( bndTriCount );   if(nRefine > 0) free( refineBndTriCount );
    free( bndTri2Count );  if(nRefine > 0) free( refineBndTri2Count );
    free( bndBeamCount );  if(nRefine > 0) free( refineBndBeamCount );
    free( bndBeam2Count ); if(nRefine > 0) free( refineBndBeam2Count );
    ////cout << "MW::RevocapRefine ----------------- H6" << endl;
    for(uiint ibnd=0; ibnd < vBndHexaCount.size(); ibnd++){  free( bndHexaNodes[ibnd] ); if(nRefine > 0) free( refineBndHexaNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndHexa2Count.size(); ibnd++){ free( bndHexa2Nodes[ibnd] ); if(nRefine > 0) free( refineBndHexa2Nodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndTetCount.size(); ibnd++){  free( bndTetNodes[ibnd] ); if(nRefine > 0) free( refineBndTetNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndTet2Count.size(); ibnd++){ free( bndTet2Nodes[ibnd] ); if(nRefine > 0) free( refineBndTet2Nodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndPrismCount.size(); ibnd++){  free( bndPrismNodes[ibnd] ); if(nRefine > 0) free( refineBndPrismNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndPrism2Count.size(); ibnd++){ free( bndPrism2Nodes[ibnd] ); if(nRefine > 0) free( refineBndPrism2Nodes[ibnd]); }

    for(uiint ibnd=0; ibnd < vBndQuadCount.size(); ibnd++){  free( bndQuadNodes[ibnd] ); if(nRefine > 0) free( refineBndQuadNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndQuad2Count.size(); ibnd++){ free( bndQuad2Nodes[ibnd] ); if(nRefine > 0) free( refineBndQuad2Nodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndTriCount.size(); ibnd++){ free( bndTriNodes[ibnd] ); if(nRefine > 0) free( refineBndTriNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndTri2Count.size(); ibnd++){ free( bndTri2Nodes[ibnd] ); if(nRefine > 0) free( refineBndTri2Nodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndBeamCount.size(); ibnd++){ free( bndBeamNodes[ibnd] ); if(nRefine > 0) free( refineBndBeamNodes[ibnd]); }
    for(uiint ibnd=0; ibnd < vBndBeam2Count.size(); ibnd++){ free( bndBeam2Nodes[ibnd] ); if(nRefine > 0) free( refineBndBeam2Nodes[ibnd]); }
    ////cout << "MW::RevocapRefine ----------------- H7" << endl;
    free(bndHexaNodes);  if(nRefine > 0) free(refineBndHexaNodes);
    free(bndHexa2Nodes); if(nRefine > 0) free(refineBndHexa2Nodes);
    free(bndTetNodes);  if(nRefine > 0) free(refineBndTetNodes);
    free(bndTet2Nodes); if(nRefine > 0) free(refineBndTet2Nodes);
    free(bndPrismNodes);  if(nRefine > 0) free(refineBndPrismNodes);
    free(bndPrism2Nodes); if(nRefine > 0) free(refineBndPrism2Nodes);

    free(bndQuadNodes); if(nRefine > 0) free(refineBndQuadNodes);
    free(bndQuad2Nodes);if(nRefine > 0) free(refineBndQuad2Nodes);
    free(bndTriNodes);  if(nRefine > 0) free(refineBndTriNodes);
    free(bndTri2Nodes); if(nRefine > 0) free(refineBndTri2Nodes);
    free(bndBeamNodes); if(nRefine > 0) free(refineBndBeamNodes);
    free(bndBeam2Nodes);if(nRefine > 0) free(refineBndBeam2Nodes);
    ////cout << "MW::RevocapRefine ----------------- H8" << endl;

    // Comm free
    free( commQuadCount);  if(nRefine > 0) free( refineCommQuadCount);
    free( commQuad2Count); if(nRefine > 0) free( refineCommQuad2Count);
    free( commTriCount);   if(nRefine > 0) free( refineCommTriCount);
    free( commTri2Count);  if(nRefine > 0) free( refineCommTri2Count);
    free( commBeamCount);  if(nRefine > 0) free( refineCommBeamCount);
    free( commBeam2Count); if(nRefine > 0) free( refineCommBeam2Count);
    free( commPointCount); 
    ////cout << "MW::RevocapRefine ----------------- H9" << endl;
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commQuadNodes[icom] ); if(nRefine > 0) free(refineCommQuadNodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commQuad2Nodes[icom] ); if(nRefine > 0) free(refineCommQuad2Nodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commTriNodes[icom] ); if(nRefine > 0) free(refineCommTriNodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commTri2Nodes[icom] ); if(nRefine > 0) free( refineCommTri2Nodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commBeamNodes[icom] ); if(nRefine > 0) free( refineCommBeamNodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commBeam2Nodes[icom] ); if(nRefine > 0) free( refineCommBeam2Nodes[icom]);}
    for(uiint icom=0; icom < nNumOfCommMesh2; icom++ ){ free( commPointNodes[icom] );}
    ////cout << "MW::RevocapRefine ----------------- H10" << endl;
    free(commQuadNodes);   if(nRefine > 0) free(refineCommQuadNodes);
    free(commQuad2Nodes);  if(nRefine > 0) free(refineCommQuad2Nodes);
    free(commTriNodes);    if(nRefine > 0) free(refineCommTriNodes);
    free(commTri2Nodes);   if(nRefine > 0) free(refineCommTri2Nodes);
    free(commBeamNodes);   if(nRefine > 0) free(refineCommBeamNodes);
    free(commBeam2Nodes);  if(nRefine > 0) free(refineCommBeam2Nodes);
    free(commPointNodes); 

    return MW_SUCCESS;
}
//境界条件管理変数セットアップ:Bnd Volume
void CMW::Rcap_BndParamSetup(uiint nElemType, uiint& nNumOfVol,vuint& vBndID, vstring& vBndName, vuint& vBndBType,
        vuint& vBndDOF, vvuint& vvBndDOF, vuint& vBndType, vuint& vBndCount, vuint& vBndTypeSeq, uiint& nSeq,
        CBoundaryVolumeMesh *pVolMesh)
{
    vBndID.push_back(pVolMesh->getID());     //境界ID
    vBndName.push_back(pVolMesh->getName());    //境界名
    vBndBType.push_back(pVolMesh->getBndType());//境界種類
    vBndDOF.push_back(pVolMesh->getNumOfDOF()); //境界DOF数

    vuint vDOF; vDOF.resize(pVolMesh->getNumOfDOF());
    for(uiint idof=0; idof < pVolMesh->getNumOfDOF(); idof++){
        vDOF[idof]= pVolMesh->getDOF(idof);
    }
    vvBndDOF.push_back(vDOF);

    vBndType.push_back(nElemType);//境界番号-形状タイプ
    vBndCount.push_back(nNumOfVol);
    vBndTypeSeq.push_back(nSeq);//境界番号-形状の何番目
    nSeq++;
}

//コースグリッドのElementNodes配列をセット
void CMW::Rcap_CrsElemNodes(size_t nNNode, vector<CElement*> vElem, CIndexBucket *pBucket, int32_t *elemNodes)
{
    for(size_t i; i < vElem.size(); i++){
        CElement *pElem= vElem[i];//親要素 Hexa
        for(size_t ii=0; ii < nNNode; ii++){
            CNode *pNode = pElem->getNode(ii);
            uiint id = pNode->getID();
            elemNodes[i*nNNode + ii] = pBucket->getIndexNode(id);
        };
    }
}

void CMW::Rcap_NodeGene(uiint iLevel, uiint nMaxLevel, size_t crsNumOfNode, size_t refineNodeCount, float64_t* resultCoord, CMesh *pCrsMesh, CMesh *pProgMesh)
{
    // コースグリッドNode 1点 取得
    CNode *pCrsNode= pCrsMesh->getNodeIX(0);
    uiint nNodeType= pCrsNode->getType();
    uiint nSDOF,nVDOF;
    uiint nMaxID, nCountID(0);
    // コースグリッドNodeを新Meshにセット
    for(size_t inode=0; inode < crsNumOfNode; inode++){
        CNode *pNode = pCrsMesh->getNodeIX(inode);
        pProgMesh->setNode(pNode);//-------------------ファインMeshへNodeをセット

        if(inode==crsNumOfNode-1) nMaxID=pNode->getID();//------- コースグリッド MaxIDの取得
    }
    // 新Nodeを新Meshにセット: #コースグリッド以降の番号が新Node
    for(size_t inode = crsNumOfNode; inode < refineNodeCount; inode++){
        CNode *pNode;
        if(nNodeType==NodeType::Scalar){
            pNode= new CScalarNode();//----------------------- new Node
            nSDOF= pCrsNode->getScalarDOF();
            pNode->setScalarDOF(nSDOF);
        }
        if(nNodeType==NodeType::Vector){
            pNode= new CVectorNode();//----------------------- new Node
            nVDOF= pCrsNode->getVectorDOF();
            pNode->setVectorDOF(nVDOF);
        }
        if(nNodeType==NodeType::ScalarVector){
            pNode= new CScalarVectorNode();//----------------- new Node
            nSDOF= pCrsNode->getScalarDOF();  nVDOF= pCrsNode->getVectorDOF();
            pNode->setScalarDOF(nSDOF);  pNode->setVectorDOF(nVDOF);
        }

        int8_t nType;//REVOCAPタイプ
        int32_t *crsNodes = (int32_t*)calloc(20, sizeof(int32_t));
        nType = rcapGetOriginal(inode, crsNodes);

        pNode->resizeGridLevel(nMaxLevel+1);//--------------- GridLevel 確保

        CNode *pParentNode;
        uiint nNNode;
        switch(nType){
            case(RCAP_HEXAHEDRON): nNNode = NumberOfNode::Hexa();   break;
            case(RCAP_HEXAHEDRON2): nNNode = NumberOfNode::Hexa2();   break;
            case(RCAP_TETRAHEDRON): nNNode = NumberOfNode::Tetra();   break;
            case(RCAP_TETRAHEDRON2): nNNode = NumberOfNode::Tetra2();   break;
            case(RCAP_WEDGE): nNNode = NumberOfNode::Prism();   break;
            case(RCAP_WEDGE2): nNNode = NumberOfNode::Prism2();   break;
            case(RCAP_QUAD): nNNode = NumberOfNode::Quad();   break;
            case(RCAP_QUAD2): nNNode = NumberOfNode::Quad2();    break;
            case(RCAP_TRIANGLE): nNNode = NumberOfNode::Triangle();    break;
            case(RCAP_TRIANGLE2): nNNode = NumberOfNode::Triangle2();   break;
            case(RCAP_SEGMENT): nNNode = NumberOfNode::Beam();     break;
            case(RCAP_SEGMENT2): nNNode = NumberOfNode::Beam2();    break;
        }
        for(uiint k=0; k < nNNode; k++){ 
            pParentNode= pCrsMesh->getNodeIX(crsNodes[k]);
            pNode->addParentNode(iLevel, pParentNode);
        };
        free( crsNodes );

        vdouble vCoord; vCoord.resize(3);
        vCoord[0]= resultCoord[inode*3]; vCoord[1]= resultCoord[inode*3 +1];  vCoord[2]= resultCoord[inode*3 +2];

        pNode->setCoord(vCoord);
        nCountID++;
        pNode->setID(nCountID + nMaxID);//ID:コースグリッドMaxID + カウント

        pProgMesh->setNode(pNode);//-------------------ファインMeshへNodeをセット
    };
}
void CMW::Rcap_ElemGene(int8_t nType, uiint& nIDCount, size_t refineCount, int32_t *refineNodes, CMesh *pProgMesh)
{
    size_t nNNode;
    CElement *pElem;
    // Element 節点数
    switch(nType){
        case(RCAP_HEXAHEDRON):  nNNode = NumberOfNode::Hexa();  break;
        case(RCAP_HEXAHEDRON2): nNNode = NumberOfNode::Hexa2(); break;
        case(RCAP_TETRAHEDRON):  nNNode = NumberOfNode::Tetra();  break;
        case(RCAP_TETRAHEDRON2): nNNode = NumberOfNode::Tetra2(); break;
        case(RCAP_WEDGE):      nNNode = NumberOfNode::Prism();  break;
        case(RCAP_WEDGE2):     nNNode = NumberOfNode::Prism2(); break;
        case(RCAP_QUAD):  nNNode = NumberOfNode::Quad();  break;
        case(RCAP_QUAD2): nNNode = NumberOfNode::Quad2(); break;
        case(RCAP_TRIANGLE):  nNNode = NumberOfNode::Triangle();  break;
        case(RCAP_TRIANGLE2): nNNode = NumberOfNode::Triangle2(); break;
        case(RCAP_SEGMENT):  nNNode = NumberOfNode::Beam();  break;
        case(RCAP_SEGMENT2): nNNode = NumberOfNode::Beam2(); break;
    }

    // Element 生成
    for(size_t ie=0; ie < refineCount; ie++){
        switch(nType){
            case(RCAP_HEXAHEDRON):  pElem = new CHexa;  break;
            case(RCAP_HEXAHEDRON2): pElem = new CHexa2; break;
            case(RCAP_TETRAHEDRON):  pElem = new CTetra;  break;
            case(RCAP_TETRAHEDRON2): pElem = new CTetra2; break;
            case(RCAP_WEDGE):   pElem = new CPrism;  break;
            case(RCAP_WEDGE2):  pElem = new CPrism2; break;
            case(RCAP_QUAD):  pElem = new CQuad;  break;
            case(RCAP_QUAD2): pElem = new CQuad2; break;
            case(RCAP_TRIANGLE):  pElem = new CTriangle;  break;
            case(RCAP_TRIANGLE2): pElem = new CTriangle2; break;
            case(RCAP_SEGMENT):  pElem = new CBeam;  break;
            case(RCAP_SEGMENT2): pElem = new CBeam2; break;
        }
        pProgMesh->setElement(pElem);
        pElem->initialize();
        pElem->setID(nIDCount);
        nIDCount++;// 新Element ID

        //Nodeセット
        for(size_t ii=0; ii < nNNode; ii++){
            uiint nNodeIX = refineNodes[nNNode*ie + ii];//通し番号をIDとして扱う
            CNode *pNode = pProgMesh->getNodeIX(nNodeIX);
            pElem->setNode(pNode, ii);

            //2次要素 辺ノード:REVOCAPからMWの辺番号へ変換
            uiint iedge;
            switch(nType){
                case(RCAP_HEXAHEDRON2):
                    if(ii >= 8){ iedge=ii-8; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
                case(RCAP_TETRAHEDRON2):
                    if(ii == 4){ iedge=1; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 5){ iedge=2; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 6){ iedge=0; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 7){ iedge=3; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 8){ iedge=4; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 9){ iedge=5; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
                case(RCAP_WEDGE2):
                    if(ii == 6){ iedge=2; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 7){ iedge=1; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 8){ iedge=0; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 9){ iedge=7; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii ==10){ iedge=8; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii ==11){ iedge=6; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii ==12){ iedge=3; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii ==13){ iedge=4; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii ==14){ iedge=5; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
                case(RCAP_QUAD2):
                    if(ii >= 4){ iedge=ii-4; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
                case(RCAP_TRIANGLE2):
                    if(ii == 3){ iedge=1; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 4){ iedge=2; pElem->setEdgeInterNode(pNode, iedge);}
                    if(ii == 5){ iedge=0; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
                case(RCAP_SEGMENT2):
                    if(ii == 2){ iedge=0; pElem->setEdgeInterNode(pNode, iedge);}
                    break;
            }
        };
        
    };//ie loop
}

void CMW::Rcap_SetBFaceMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        CBoundaryFaceMesh *pProgBFaceMesh= new CBoundaryFaceMesh;//----------------------- BoundaryFaceMesh 生成
        pProgBFaceMesh->setMaxMGLevel(nMaxLevel);
        pProgBFaceMesh->setMGLevel(iLevel);
        pProgBFaceMesh->setID(vBndID[ibnd]);
        pProgBFaceMesh->setBndType(vBndBType[ibnd]);//境界タイプ
        pProgBFaceMesh->setName(vBndName[ibnd]);    //境界名
        pProgBFaceMesh->resizeDOF(vNumDOF[ibnd]);   //境界自由度数
        uiint idof, dof;
        for(idof=0; idof < vNumDOF[ibnd]; idof++){
            dof = vvDOF[ibnd][idof];
            pProgBFaceMesh->setDOF(idof, dof);
        }

        //Poland
        CBoundaryFaceMesh *pBFaceMesh= pMesh->getBndFaceMeshID(vBndID[ibnd]);
        map<uiint,CPoland*> mPoland= pBFaceMesh->getPoland();
        pProgBFaceMesh->setPoland(mPoland);

        vuint vPolandDOF= pBFaceMesh->getPolandDOF();
        pProgBFaceMesh->setPolandDOF(vPolandDOF);

        pProgMesh->setBndFaceMesh(pProgBFaceMesh);
    }
}
void CMW::Rcap_SetBEdgeMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        CBoundaryEdgeMesh *pProgBEdgeMesh= new CBoundaryEdgeMesh;//----------------------- BoundaryEdgeMesh 生成
        pProgBEdgeMesh->setMaxMGLevel(nMaxLevel);
        pProgBEdgeMesh->setMGLevel(iLevel);
        pProgBEdgeMesh->setID(vBndID[ibnd]);
        pProgBEdgeMesh->setBndType(vBndBType[ibnd]);//境界タイプ
        pProgBEdgeMesh->setName(vBndName[ibnd]);    //境界名
        pProgBEdgeMesh->resizeDOF(vNumDOF[ibnd]);   //境界名
        uiint idof, dof;
        for(idof=0; idof < vNumDOF[ibnd]; idof++){
            dof = vvDOF[ibnd][idof];
            pProgBEdgeMesh->setDOF(idof, dof);
        }
        //Poland
        CBoundaryEdgeMesh *pBEdgeMesh= pMesh->getBndEdgeMeshID(vBndID[ibnd]);
        map<uiint,CPoland*> mPoland= pBEdgeMesh->getPoland();
        pProgBEdgeMesh->setPoland(mPoland);

        vuint vPolandDOF= pBEdgeMesh->getPolandDOF();
        pProgBEdgeMesh->setPolandDOF(vPolandDOF);

        pProgMesh->setBndEdgeMesh(pProgBEdgeMesh);
    }
}
void CMW::Rcap_SetBVolMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        CBoundaryVolumeMesh *pProgBVolMesh = new CBoundaryVolumeMesh;//------------- BoundaryVolumeMesh 生成
        pProgBVolMesh->setMaxMGLevel(nMaxLevel);
        pProgBVolMesh->setMGLevel(iLevel);
        pProgBVolMesh->setID(vBndID[ibnd]);
        pProgBVolMesh->setBndType(vBndBType[ibnd]);//境界タイプ
        pProgBVolMesh->setName(vBndName[ibnd]);    //境界名
        pProgBVolMesh->resizeDOF(vNumDOF[ibnd]);   //境界名
        uiint idof, dof;
        for(idof=0; idof < vNumDOF[ibnd]; idof++){
            dof = vvDOF[ibnd][idof];
            pProgBVolMesh->setDOF(idof, dof);
        }
        //Poland
        CBoundaryVolumeMesh *pBVolMesh= pMesh->getBndVolumeMeshID(vBndID[ibnd]);
        map<uiint,CPoland*> mPoland= pBVolMesh->getPoland();
        pProgBVolMesh->setPoland(mPoland);

        vuint vPolandDOF= pBVolMesh->getPolandDOF();
        pProgBVolMesh->setPolandDOF(vPolandDOF);

        pProgMesh->setBndVolumeMesh(pProgBVolMesh);
    }
}
vvuint CMW::Rcap_NodeIX_SortMerge(size_t nNumOfID, size_t *nEntityCount, size_t nNumOfEntityNode, int32_t **refineNodes)
{
    vvuint vvNodeIndex;
    for(size_t i=0; i < nNumOfID; i++){
        vuint vNodeIndex;
        for(size_t ient=0; ient < nEntityCount[i]; ient++){
            for(size_t inode=0; inode < nNumOfEntityNode; inode++){
                uiint nNodeIX = refineNodes[i][ient*nNumOfEntityNode + inode];
                vNodeIndex.push_back(nNodeIX);
            }
        }
        // ソート & マージ
        std::sort(vNodeIndex.begin(), vNodeIndex.end());
        std::vector<uiint>::iterator new_end = std::unique(vNodeIndex.begin(), vNodeIndex.end());
        vNodeIndex.erase(new_end, vNodeIndex.end());

        vvNodeIndex.push_back(vNodeIndex);
    }
    
    return vvNodeIndex;
}

void CMW::Rcap_BNodeGene_VolMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){

        uiint nBndID = vBndID[ibnd];
        //ファイングリッドのBoundaryMesh
        CBoundaryVolumeMesh *pBVolMesh= pProgMesh->getBndVolumeMeshID(nBndID);

        //コースグリッドのBoundaryMesh
        CBoundaryVolumeMesh *pCrsBVolMesh= pMesh->getBndVolumeMeshID(nBndID);
        uiint nCrsBNodeNum = pCrsBVolMesh->getNumOfBNode();

        if(nCrsBNodeNum >= vvNodeIndex[ibnd].size()) return;// Error

        map<uiint, uiint> NodeNum2BNodeNum;//-- コースグリッド節点番号(Index)から、コースグリッドBNode番号(Index)を取得

        //コースグリッドのBNodeは、そのまま利用する.
        for(uiint i=0; i < nCrsBNodeNum; i++){
            CBoundaryNode *pBNode= pCrsBVolMesh->getBNodeIX(i);//---- コースグリッドBNode

            pBVolMesh->addBNode(pBNode);

            uiint nNodeNum= vvNodeIndex[ibnd][i];
            NodeNum2BNodeNum[nNodeNum]= i;
        }
        //Poland
        map<uiint,CPoland*> mPoland= pBVolMesh->getPoland();
        CCalc* pCalc= CCalc::Instance();

        //--
        //ファイングリッドNode -> ファイングリッドBNode生成
        //--
        for(uiint i=nCrsBNodeNum; i < vvNodeIndex[ibnd].size(); i++){
            uiint nIndex= vvNodeIndex[ibnd][i];
            CNode *pNode= pProgMesh->getNodeIX(nIndex);

            //--------- Dirichlet値 分配 -----------//
            vdouble vVal;
            int32_t *crsNodes;
            if(pBVolMesh->getBndType()==BoundaryType::Dirichlet){
                vuint vBNodeIndex;//親BNodeインデックス番号
                // コースグリッドNode
                int8_t  nType;//REVOCAPタイプ
                int32_t progNodeNum = (int32_t)nIndex;//------------- Node通し番号: 型サイズ !問題
                crsNodes = (int32_t*)calloc(20, sizeof(int32_t));
                nType = rcapGetOriginal(progNodeNum, crsNodes);

                // コースグリッドBNode取得 => vVal計算
                // 要素中心のノード
                if(nType==RCAP_HEXAHEDRON)  Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 8, crsNodes);
                if(nType==RCAP_HEXAHEDRON2) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 20, crsNodes);
                if(nType==RCAP_TETRAHEDRON) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 4, crsNodes);
                if(nType==RCAP_TETRAHEDRON2) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 10, crsNodes);
                if(nType==RCAP_WEDGE)  Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 6, crsNodes);
                if(nType==RCAP_WEDGE2) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 15, crsNodes);
                // 面中心のノード
                if(nType==RCAP_QUAD)  Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 4, crsNodes);
                if(nType==RCAP_QUAD2)  Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 8, crsNodes);
                if(nType==RCAP_TRIANGLE) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 3, crsNodes);
                if(nType==RCAP_TRIANGLE2) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 6, crsNodes);
                // 辺中心のノード
                if(nType==RCAP_SEGMENT) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 2, crsNodes);
                if(nType==RCAP_SEGMENT2) Rcap_BNodeValueDist(iLevel, pBVolMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 3, crsNodes);

                free( crsNodes );//コースグリッドNode配列 解放
            }

            CBoundaryNode *pBNode = new CBoundaryNode;//------------- BNode 生成
            pBNode->setNode(pNode);
            pBNode->setID(i);
            pBNode->setMGLevel(iLevel);
            pBNode->resizeValue(nMaxLevel+1 - iLevel);//引数:Level差

            ////cout << "MW::Rcap_BNodeGene_VolMesh  i:" << i << " NodeID:" << pNode->getID() << endl;

            
            //--------- Dirichlet値セット -----------//
            if(pBVolMesh->getBndType()==BoundaryType::Dirichlet){
                for(uiint idof=0; idof < vVal.size(); idof++){
                    uiint dof = pBVolMesh->getDOF(idof);
                    
                    pBNode->initRcapBool(dof, nMaxLevel);//------------------ BOOL
                    pBNode->setEntValue(dof, iLevel, vVal[idof]);

                    ////cout << "MW::Rcap_BNodeGene_VolMesh  dof:" << dof
                    ////        << "  i:" << i << " NodeID:" << pNode->getID()
                    ////        << "  vVal:" << vVal[idof] << "  rank:" << mpMPI->getRank() << endl;

                    if( pBVolMesh->existPoland(dof) ){
                        double x=pBNode->getX(),y=pBNode->getY(),z=pBNode->getZ();
                        pCalc->setElementParam(vVal[idof], x,y,z);
                        double calcVal= pCalc->Exec(mPoland[dof]);
                        pBNode->setValue(dof, iLevel, calcVal);
                    }else{
                        pBNode->setValue(dof, iLevel, vVal[idof]);
                    }
                };
            }

            pBVolMesh->addBNode(pBNode);

        };//----------------------------------------ファイングリッドBNode loop
    };//ibnd loop
}

void CMW::Rcap_BNodeGene_FaceMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex)
{
    ////cout << "MW::Rcap_BNodeGene_FaceMesh  --- A  rank:" << mpMPI->getRank() << endl;

    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){

        uiint nBndID = vBndID[ibnd];
        //ファイングリッドのBoundaryMesh
        CBoundaryFaceMesh *pBFaceMesh= pProgMesh->getBndFaceMeshID(nBndID);

        //コースグリッドのBoundaryMesh
        CBoundaryFaceMesh *pCrsBFaceMesh= pMesh->getBndFaceMeshID(nBndID);
        uiint nCrsBNodeNum = pCrsBFaceMesh->getNumOfBNode();

        if(nCrsBNodeNum >= vvNodeIndex[ibnd].size()) return;// Error

        map<uiint, uiint> NodeNum2BNodeNum;//-- コースグリッド節点番号(Index)から、コースグリッドBNode番号(Index)を取得

        //コースグリッドのBNodeは、そのまま利用する.
        for(uiint i=0; i < nCrsBNodeNum; i++){
            CBoundaryNode *pBNode= pCrsBFaceMesh->getBNodeIX(i);//---- コースグリッドBNode

            pBFaceMesh->addBNode(pBNode);

            uiint nNodeNum= vvNodeIndex[ibnd][i];
            NodeNum2BNodeNum[nNodeNum]= i;
        }
        //Poland
        map<uiint,CPoland*> mPoland= pBFaceMesh->getPoland();
        CCalc *pCalc= CCalc::Instance();

        //--
        //ファイングリッドNode -> ファイングリッドBNode生成
        //--
        for(uiint i=nCrsBNodeNum; i < vvNodeIndex[ibnd].size(); i++){
            uiint nIndex= vvNodeIndex[ibnd][i];
            CNode *pNode= pProgMesh->getNodeIX(nIndex);

            //--------- Dirichlet値 分配 -----------//
            vdouble vVal;
            int32_t *crsNodes;
            if(pBFaceMesh->getBndType()==BoundaryType::Dirichlet){
                vuint vBNodeIndex;//親BNodeインデックス番号
                // コースグリッドNode
                int8_t  nType;//REVOCAPタイプ
                int32_t progNodeNum = (int32_t)nIndex;//------------- Node通し番号: 型サイズ !問題
                crsNodes = (int32_t*)calloc(8, sizeof(int32_t));
                nType = rcapGetOriginal(progNodeNum, crsNodes);
                //
                // コースグリッドBNode取得 => vVal計算(エンティティ値)
                // 
                // 面
                if(nType==RCAP_QUAD)  Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 4, crsNodes);
                if(nType==RCAP_QUAD2)  Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 8, crsNodes);
                if(nType==RCAP_TRIANGLE) Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 3, crsNodes);
                if(nType==RCAP_TRIANGLE2) Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 6, crsNodes);
                // 辺
                if(nType==RCAP_SEGMENT) Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 2, crsNodes);
                if(nType==RCAP_SEGMENT2) Rcap_BNodeValueDist(iLevel, pBFaceMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 3, crsNodes);
                
                free( crsNodes );//コースグリッドNode配列 解放
            }

            ////cout << "MW::Rcap_BNodeGene_FaceMesh  --- B  i:" << i << "  rank:" << mpMPI->getRank() << endl;
            
            CBoundaryNode *pBNode = new CBoundaryNode;//------------- BNode 生成
            pBNode->setNode(pNode);
            pBNode->setID(i);
            pBNode->setMGLevel(iLevel);
            pBNode->resizeValue(nMaxLevel+1 - iLevel);//引数:Level差
            
            //--------- Dirichlet値セット -----------//
            if(pBFaceMesh->getBndType()==BoundaryType::Dirichlet){

                ////cout << "MW::Rcap_BNodeGene_FaceMesh  --- C vVal.size::" << vVal.size() << endl;

                for(uiint idof=0; idof < vVal.size(); idof++){
                    uiint dof = pBFaceMesh->getDOF(idof);
                    
                    pBNode->initRcapBool(dof, nMaxLevel);//------------------ BOOL
                    pBNode->setEntValue(dof, iLevel, vVal[idof]);

                    ////cout << "MW::Rcap_BNodeGene_FaceMesh  dof:" << dof
                    ////        << "  i:" << i << " NodeID:" << pNode->getID()
                    ////        << "  vVal:" << vVal[idof] << "  rank:" << mpMPI->getRank() << endl;

                    if( pBFaceMesh->existPoland(dof) ){
                        double x=pBNode->getX(),y=pBNode->getY(),z=pBNode->getZ();
                        pCalc->setElementParam(vVal[idof], x,y,z);
                        double calcVal= pCalc->Exec(mPoland[dof]);

                        pBNode->setValue(dof, iLevel, calcVal);
                    }else{
                        pBNode->setValue(dof, iLevel, vVal[idof]);
                    }
                };
            }

            pBFaceMesh->addBNode(pBNode);
            
        };//----------------------- ファイングリッドBNode loop

    };//ibnd loop
}
// Dirichlet値 分配 : コースグリッドのBNodeから値を取得
// # VolMesh
void CMW::Rcap_BNodeValueDist(uiint iLevel, CBoundaryVolumeMesh *pBVolMesh, vdouble& vVal, vuint& vBNodeIndex,
                                map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfEntityNode, int32_t *crsNodes)
{
    //コースグリッドのBNode番号
    for(uiint icrs=0; icrs < nNumOfEntityNode; icrs++){
        uiint crsNodeNum = crsNodes[icrs]; // 親"Node"のインデックス番号
        uiint crsBNodeNum= NodeNum2BNodeNum[crsNodeNum];
        vBNodeIndex.push_back(crsBNodeNum);// 親"BNode"のインデックス番号
    };

    uiint nNumOfDOF = pBVolMesh->getNumOfDOF();

    ////cout << "MW::Rcap_BNodeValueDist  (pBVolMesh) nNumOfDOF:" << nNumOfDOF << endl;

    vVal.resize(nNumOfDOF);
    for(uiint i=0; i < nNumOfDOF; i++) vVal[i]=0.0;

    // 親"BNode"のDirichlet値 加算集積
    for(uiint ibnode=0; ibnode < nNumOfEntityNode; ibnode++){
        uiint index = vBNodeIndex[ibnode];// 親"BNode"のインデックス番号
        CBoundaryNode *pCrsBNode = pBVolMesh->getBNodeIX(index);

        for(uiint idof=0; idof < nNumOfDOF; idof++){
            uiint dof = pBVolMesh->getDOF(idof);//自由度Indexに対する自由度番号
            double entVal= pCrsBNode->getEntValue(dof, iLevel-1);//////////////////////////// iLevel-1 '12.10.04
            vVal[idof] += entVal;

            //
            // コースグリッド節点　'12.10.04 追加
            //
            if( !pCrsBNode->isSetupValue(dof, iLevel) ){
                pCrsBNode->setEntValue(dof, iLevel, entVal);//---------------- コースグリッド iLevel:entVal
                if( pBVolMesh->existPoland(dof) ){
                    double x=pCrsBNode->getX(),y=pCrsBNode->getY(),z=pCrsBNode->getZ();
                    CPoland* pPoland=pBVolMesh->getPoland(dof);
                    CCalc*   pCalc=CCalc::Instance();
                    pCalc->setElementParam(entVal, x,y,z);
                    double calcVal= pCalc->Exec(pPoland);

                    pCrsBNode->setValue(dof, iLevel, calcVal);
                }else{
                    pCrsBNode->setValue(dof, iLevel, entVal);
                }
            }
        };
    };
    // 新Dirichlet値:加算した値を平均化
    for(uiint idof=0; idof < nNumOfDOF; idof++){
        vVal[idof] /= nNumOfEntityNode;
    };
}
// # FaceMesh
void CMW::Rcap_BNodeValueDist(uiint iLevel, CBoundaryFaceMesh *pBFaceMesh, vdouble& vVal, vuint& vBNodeIndex,
                                map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfFaceNode, int32_t *crsNodes)
{
    //コースグリッドのBNode番号
    for(uiint icrs=0; icrs < nNumOfFaceNode; icrs++){
        uiint crsNodeNum = crsNodes[icrs]; // 親"Node"のインデックス番号
        uiint crsBNodeNum= NodeNum2BNodeNum[crsNodeNum];
        vBNodeIndex.push_back(crsBNodeNum);// 親"BNode"のインデックス番号
    };

    uiint nNumOfDOF = pBFaceMesh->getNumOfDOF();

    ////cout << "MW::Rcap_BNodeValueDist --- nNumOfDOF:" << nNumOfDOF << endl;

    vVal.resize(nNumOfDOF);
    for(uiint i=0; i < nNumOfDOF; i++) vVal[i]=0.0;

    // 親"BNode"のDirichlet値 加算集積
    for(uiint ibnode=0; ibnode < nNumOfFaceNode; ibnode++){
        uiint index = vBNodeIndex[ibnode];// 親"BNode"のインデックス番号
        CBoundaryNode *pCrsBNode = pBFaceMesh->getBNodeIX(index);

        for(uiint idof=0; idof < nNumOfDOF; idof++){
            uiint dof = pBFaceMesh->getDOF(idof);//自由度Indexに対する自由度番号
            double entVal= pCrsBNode->getEntValue(dof, iLevel-1);//////////////////////////// iLevel-1 '12.10.04
            vVal[idof] += entVal;

            ////cout << "MW::Rcap_BNodeValueDist(Face)  CrsBNode index:" << index
            ////        << "  NodeID:" << pCrsBNode->getNode()->getID()
            ////        << "  EntValue:" << pCrsBNode->getEntValue(dof, iLevel-1)
            ////        << "  iLevel:" << iLevel << "  rank:" << mpMPI->getRank() << endl;

            //
            // コースグリッド節点　'12.10.04 追加
            //
            if( !pCrsBNode->isSetupValue(dof, iLevel) ){
                pCrsBNode->setEntValue(dof, iLevel, entVal);//---------------- コースグリッド iLevel:entVal
                if( pBFaceMesh->existPoland(dof) ){
                    double x=pCrsBNode->getX(),y=pCrsBNode->getY(),z=pCrsBNode->getZ();
                    CPoland* pPoland=pBFaceMesh->getPoland(dof);
                    CCalc*   pCalc=CCalc::Instance();
                    pCalc->setElementParam(entVal, x,y,z);
                    double calcVal= pCalc->Exec(pPoland);

                    pCrsBNode->setValue(dof, iLevel, calcVal);
                }else{
                    pCrsBNode->setValue(dof, iLevel, entVal);
                }
            }
        };
    };
    // 新Dirichlet値(エンティティ値):'12.10.04
    for(uiint idof=0; idof < nNumOfDOF; idof++){
        
        vVal[idof] /= nNumOfFaceNode;

        ////uiint dof= pBFaceMesh->getDOF(idof);//自由度Indexに対する自由度番号
        ////cout << "MW::Rcap_BNodeValueDist(Face) dof:" << dof
        ////        << " vVal:" << vVal[idof] << " nFaceNode数:" << nNumOfFaceNode
        ////        << " rank:" << mpMPI->getRank() << endl;
    };
}

void CMW::Rcap_BNodeGene_EdgeMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex)
{

    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){

        uiint nBndID = vBndID[ibnd];
        //ファイングリッドのBoundaryMesh
        CBoundaryEdgeMesh *pBEdgeMesh= pProgMesh->getBndEdgeMeshID(nBndID);

        //コースグリッドのBoundaryMesh
        CBoundaryEdgeMesh *pCrsBEdgeMesh= pMesh->getBndEdgeMeshID(nBndID);
        uiint nCrsBNodeNum = pCrsBEdgeMesh->getNumOfBNode();

        if(nCrsBNodeNum >= vvNodeIndex[ibnd].size()) return;// Error

        map<uiint, uiint> NodeNum2BNodeNum;//-- コースグリッド節点番号(Index)から、コースグリッドBNode番号(Index)を取得

        //コースグリッドのBNodeは、そのまま利用する.
        for(uiint i=0; i < nCrsBNodeNum; i++){
            CBoundaryNode *pBNode= pCrsBEdgeMesh->getBNodeIX(i);//---- コースグリッドBNode

            pBEdgeMesh->addBNode(pBNode);

            uiint nNodeNum= vvNodeIndex[ibnd][i];
            NodeNum2BNodeNum[nNodeNum]= i;
        }
        //Poland
        map<uiint,CPoland*> mPoland= pBEdgeMesh->getPoland();
        CCalc *pCalc= CCalc::Instance();
        
        //新BNode
        for(uiint i=nCrsBNodeNum; i < vvNodeIndex[ibnd].size(); i++){
            uiint nIndex= vvNodeIndex[ibnd][i];
            CNode *pNode= pProgMesh->getNodeIX(nIndex);

            CBoundaryNode *pBNode = new CBoundaryNode;//------------- BNode 生成
            pBNode->setNode(pNode);
            pBNode->setID(i);
            pBNode->setMGLevel(iLevel);
            pBNode->resizeValue(nMaxLevel+1 - iLevel);//引数:Level差
            
            //--------- Dirichlet値 分配 -----------//
            vdouble vVal;
            int32_t *crsNodes;
            if(pBEdgeMesh->getBndType()==BoundaryType::Dirichlet){
                vuint vBNodeIndex;//親BNodeインデックス番号
                // コースグリッドNode
                int8_t  nType;//REVOCAPタイプ
                int32_t progNodeNum = (int32_t)nIndex;//------------- Node通し番号: 型サイズ !問題
                crsNodes = (int32_t*)calloc(2, sizeof(int32_t));
                nType = rcapGetOriginal(progNodeNum, crsNodes);

                // コースグリッドBNode => vVal計算
                if(nType==RCAP_SEGMENT){
                    Rcap_BNodeValueDist(iLevel, pBEdgeMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 2, crsNodes);
                }
                if(nType==RCAP_SEGMENT2){
                    Rcap_BNodeValueDist(iLevel, pBEdgeMesh, vVal, vBNodeIndex, NodeNum2BNodeNum, 3, crsNodes);
                }
                free( crsNodes );//コースグリッドNode配列 解放
            }
            //--------- Dirichlet値セット -----------//
            if(pBEdgeMesh->getBndType()==BoundaryType::Dirichlet){
                for(uiint idof=0; idof < vVal.size(); idof++){
                    uiint dof = pBEdgeMesh->getDOF(idof);
                    
                    pBNode->initRcapBool(dof, nMaxLevel);//------------------ BOOL
                    pBNode->setEntValue(dof, iLevel, vVal[idof]);

                    if( pBEdgeMesh->existPoland(dof) ){
                        double x=pBNode->getX(),y=pBNode->getY(),z=pBNode->getZ();
                        pCalc->setElementParam(vVal[idof], x,y,z);
                        double calcVal= pCalc->Exec(mPoland[dof]);
                        pBNode->setValue(dof, iLevel, calcVal);
                    }else{
                        pBNode->setValue(dof, iLevel, vVal[idof]);
                    }
                };
            }

            pBEdgeMesh->addBNode(pBNode);
        }
    }
}
void CMW::Rcap_BNodeValueDist(uiint iLevel, CBoundaryEdgeMesh *pBEdgeMesh, vdouble& vVal, vuint& vBNodeIndex,
                                map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfEdgeNode, int32_t *crsNodes)
{
    //コースグリッドのBNode番号
    for(uiint icrs=0; icrs < nNumOfEdgeNode; icrs++){
        uiint crsNodeNum = crsNodes[icrs]; // 親"Node"のインデックス番号
        uiint crsBNodeNum= NodeNum2BNodeNum[crsNodeNum];
        vBNodeIndex.push_back(crsBNodeNum);// 親"BNode"のインデックス番号
    };

    uiint nNumOfDOF = pBEdgeMesh->getNumOfDOF();

    vVal.resize(nNumOfDOF);
    for(uiint i=0; i < nNumOfDOF; i++) vVal[i]=0.0;
    
    // 親"BNode"のDirichlet値 加算集積
    for(uiint ibnode=0; ibnode < nNumOfEdgeNode; ibnode++){
        uiint index = vBNodeIndex[ibnode];// 親"BNode"のインデックス番号
        CBoundaryNode *pCrsBNode = pBEdgeMesh->getBNodeIX(index);

        for(uiint idof=0; idof < nNumOfDOF; idof++){
            uiint dof = pBEdgeMesh->getDOF(idof);//自由度Indexに対する自由度番号
            double entVal= pCrsBNode->getEntValue(dof, iLevel-1);//////////////////////////// iLevel-1 '12.10.04
            vVal[idof] += entVal;

            //
            // コースグリッド節点　'12.10.04 追加
            //
            if( !pCrsBNode->isSetupValue(dof, iLevel) ){
                pCrsBNode->setEntValue(dof, iLevel, entVal);//---------------- コースグリッド iLevel:entVal
                if( pBEdgeMesh->existPoland(dof) ){
                    double x=pCrsBNode->getX(),y=pCrsBNode->getY(),z=pCrsBNode->getZ();
                    CPoland* pPoland=pBEdgeMesh->getPoland(dof);
                    CCalc*   pCalc=CCalc::Instance();
                    pCalc->setElementParam(entVal, x,y,z);
                    double calcVal= pCalc->Exec(pPoland);

                    pCrsBNode->setValue(dof, iLevel, calcVal);
                }else{
                    pCrsBNode->setValue(dof, iLevel, entVal);
                }
            }
        };
    };
    // 新Dirichlet値(エンティティ値)
    for(uiint idof=0; idof < nNumOfDOF; idof++){
        vVal[idof] /= nNumOfEdgeNode;
    };
}
//
// 逆引き生成 : Meshの節点番号 => 境界条件GrのNode配列のインデックス番号 || Vol,Face生成で利用::"BNodeは、境界条件Node配列インデックス順に、生成されている"
//
//   # Meshの節点番号は、Mesh内での通し番号(Index番号)
//
void CMW::Rcap_NodeNum2Index(vector<map<uiint, uiint> >& vmNodeNum2Index, vvuint& NodeIndex)
{
    // [BndID]map[節点番号] == インデックス番号 を生成 //

    uiint nNumOfBnd= NodeIndex.size();
    vmNodeNum2Index.resize(nNumOfBnd);

    for(uiint ibnd=0; ibnd < nNumOfBnd; ibnd++){
        uiint nNumOfNode= NodeIndex[ibnd].size();

        for(uiint index=0; index < nNumOfNode; index++){
            uiint nNodeNum = NodeIndex[ibnd][index];
            vmNodeNum2Index[ibnd][nNodeNum]=index;
        };
    };
}
void CMW::Rcap_BVolGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t *nVolCount, uiint nShape, size_t& nNumOfVolNode,
                        int32_t **refineNodes, CMesh *pProgMesh, CMesh *pCrsMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        uiint nBndID= vBndID[ibnd];
        CBoundaryVolumeMesh *pBVolMesh = pProgMesh->getBndVolumeMeshID(nBndID);
        CBoundaryVolumeMesh *pCrsBVolMesh= pCrsMesh->getBndVolumeMeshID(nBndID);


        //Neumannの場合:Value初期化
        if(pBVolMesh->getBndType()==BoundaryType::Neumann){
            uiint nNumOfDOF = pBVolMesh->getNumOfDOF();
            uiint nNumOfBNode = pBVolMesh->getNumOfBNode();

            for(uiint inode=0; inode < nNumOfBNode; inode++){
                CBoundaryNode *pBNode= pBVolMesh->getBNodeIX(inode);
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof = pBVolMesh->getDOF(idof);

                    pBNode->initValue(dof, iLevel);//----------------Neumann初期化
                };
            };
        }

        //Volume構成BNodeの取得
        for(uiint ivol=0; ivol < nVolCount[ibnd]; ivol++){
            CBoundaryVolume *pBVol;
            //------------------------------------------------------- BoundaryVolume 生成
            switch(nShape){
                case(ElementType::Hexa):   pBVol = new CBoundaryHexa();  pBVol->setOrder(ElementOrder::First); break;
                case(ElementType::Hexa2):  pBVol = new CBoundaryHexa();  pBVol->setOrder(ElementOrder::Second); break;
                case(ElementType::Tetra):  pBVol = new CBoundaryTetra(); pBVol->setOrder(ElementOrder::First); break;
                case(ElementType::Tetra2): pBVol = new CBoundaryTetra(); pBVol->setOrder(ElementOrder::Second); break;
                case(ElementType::Prism):  pBVol = new CBoundaryPrism(); pBVol->setOrder(ElementOrder::First); break;
                case(ElementType::Prism2): pBVol = new CBoundaryPrism(); pBVol->setOrder(ElementOrder::Second); break;
                default: mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_BVolGene nShape tpye error");
            }
            pBVol->resizeBNode(nNumOfVolNode);

            vector<CNode*> vNode;//要素検索 パラメータ
            uiint nElementID;    //要素検索検索 解

            for(uiint inode=0; inode < nNumOfVolNode; inode++){
                uiint nNodeNum = refineNodes[ibnd][nNumOfVolNode*ivol + inode];
                uiint index= vmID2Index[ibnd][nNodeNum];

                CBoundaryNode *pBNode= pBVolMesh->getBNodeIX(index);

                if(nShape==ElementType::Hexa || nShape==ElementType::Hexa2){
                    if(inode < 8){
                        pBVol->setBNode(inode, pBNode);
                    }else{
                        uiint iedge = inode - 8;//---- MWとREVOCAPの2次Nodeの番号付け同じ
                        pBVol->setEdgeBNode(iedge, pBNode);
                    }
                }
                if(nShape==ElementType::Tetra || nShape==ElementType::Tetra2){
                    if(inode < 4){
                        pBVol->setBNode(inode, pBNode);
                    }else{
                        uiint iedge;
                        if(inode == 4){ iedge=1; }//----
                        if(inode == 5){ iedge=2; }//
                        if(inode == 6){ iedge=0; }// REVOCAP -> MWの辺番号へ変換
                        if(inode == 7){ iedge=3; }//
                        if(inode == 8){ iedge=4; }//
                        if(inode == 9){ iedge=5; }//----
                        pBVol->setEdgeBNode(iedge, pBNode);
                    }
                }
                if(nShape==ElementType::Prism || nShape==ElementType::Prism2){
                    if(inode < 6){
                        pBVol->setBNode(inode, pBNode);
                    }else{
                        uiint iedge;
                        if(inode == 6){ iedge=2; }//----
                        if(inode == 7){ iedge=1; }//
                        if(inode == 8){ iedge=0; }//
                        if(inode == 9){ iedge=7; }//
                        if(inode ==10){ iedge=8; }// REVOCAP -> MWの辺番号へ変換
                        if(inode ==11){ iedge=6; }//
                        if(inode ==12){ iedge=3; }//
                        if(inode ==13){ iedge=4; }//
                        if(inode ==14){ iedge=5; }//----
                        pBVol->setEdgeBNode(iedge, pBNode);
                    }
                }
                CNode *pNode = pBNode->getNode();
                vNode.push_back(pNode);
            };//inode loop
            
            Rcap_ElemSearch(nElementID, vNode, pProgMesh);

            pBVol->setID(ivol);
            pBVol->setElementID(nElementID);
            CElement *pElem = pProgMesh->getElement(nElementID);
            pBVol->setElement(pElem);

            pBVolMesh->addBVolume(pBVol);

            double fnCubic = pBVol->calcVolume();

            // コースグリッドのBVolのVal => ファイングリッドBVolへセット
            // -- Neumann --
            if(pBVolMesh->getBndType()==BoundaryType::Neumann){
                uiint nCrsIndex = ivol/8;//------- コースグリッド要素番号:リファインするとき8分割されるので,整数を8で除算 #Hexa,Tetra,Prism全て8分割

                CBoundaryVolume *pCrsBVol = pCrsBVolMesh->getBVolumeIX(nCrsIndex);
                double crsCubic = pCrsBVol->getCubicVolume();

                for(uiint idof=0; idof < pBVolMesh->getNumOfDOF(); idof++){
                    uiint dof = pBVolMesh->getDOF(idof);
                    double val = pCrsBVol->getBndValue(dof);

                    val = val*(fnCubic/crsCubic);
                    pBVol->setBndValue(dof, val);//------- ファイングリッドに値をセット

                    Rcap_EquivalentNodalForce(iLevel, dof, pBVolMesh, pBVol);//------ 等価節点力 配分
                };
            }

        };//ivol loop

    };//ibnd loop
}
void CMW::Rcap_BFaceGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t* nFaceCount, uiint nShape, size_t& nNumOfFaceNode,
                         int32_t** refineNodes, CMesh* pProgMesh, CMesh *pCrsMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        
        uiint nBndID= vBndID[ibnd];
        CBoundaryFaceMesh *pBFaceMesh = pProgMesh->getBndFaceMeshID(nBndID);
        CBoundaryFaceMesh *pCrsBFaceMesh= pCrsMesh->getBndFaceMeshID(nBndID);

        
        //Neumannの場合:Value初期化
        if(pBFaceMesh->getBndType()==BoundaryType::Neumann){
            uiint nNumOfDOF = pBFaceMesh->getNumOfDOF();
            uiint nNumOfBNode = pBFaceMesh->getNumOfBNode();

            for(uiint inode=0; inode < nNumOfBNode; inode++){
                CBoundaryNode *pBNode= pBFaceMesh->getBNodeIX(inode);
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof = pBFaceMesh->getDOF(idof);

                    pBNode->initValue(dof, iLevel);//-----------Neumann初期化
                };
            };
        }

        //FaceNodeの取得
        for(uiint iface=0; iface < nFaceCount[ibnd]; iface++){

            CBoundaryFace *pBFace = new CBoundaryFace();//------------------ BoundaryFace 生成
            pBFace->resizeBNode(nNumOfFaceNode);
            
            vector<CNode*> vNode;     //面検索 パラメータ
            uiint nElementID, nFaceID;//面検索 解

            for(uiint inode=0; inode < nNumOfFaceNode; inode++){
                uiint nNodeNum = refineNodes[ibnd][nNumOfFaceNode*iface + inode];
                uiint index= vmID2Index[ibnd][nNodeNum];

                CBoundaryNode *pBNode= pBFaceMesh->getBNodeIX(index);

                if(nShape==ElementType::Quad || nShape==ElementType::Quad2){
                    if(inode < 4){
                        pBFace->setBNode(inode, pBNode);
                    }else{
                        uiint iedge= inode-4;// MWとREVOCAPのQuad2次の番号の付け方は同一
                        pBFace->setEdgeBNode(iedge, pBNode);
                    }
                }
                if(nShape==ElementType::Triangle || nShape==ElementType::Triangle2){
                    if(inode < 3){
                        pBFace->setBNode(inode, pBNode);
                    }else{
                        uiint iedge;
                        if(inode==3) iedge=1;// --
                        if(inode==4) iedge=2;// REVOCAPのTri2次の番号 => MWの付け方に変更
                        if(inode==5) iedge=0;// --
                        pBFace->setEdgeBNode(iedge, pBNode);
                    }
                }
                CNode *pNode = pBNode->getNode();
                vNode.push_back(pNode);
            };//inode loop

            Rcap_ElemFaceSearch(nElementID, nFaceID, vNode, pProgMesh);//----- 要素番号、面番号の検索

            pBFace->setID(iface);
            pBFace->setBFaceShape(nShape);
            pBFace->setElementID(nElementID);
            pBFace->setElementFaceID(nFaceID);

            CElement *pElem = pProgMesh->getElement(nElementID);
            pBFace->setElement(pElem);
            
            double fnArea = pBFace->calcArea();

            pBFaceMesh->addBFace(pBFace);

            // コースグリッドのBFaceのVal => ファイングリッドBFaceへセット
            // -- Neumann --
            if(pBFaceMesh->getBndType()==BoundaryType::Neumann){

                uiint nCrsIndex = iface/4;//------- コースグリッド面番号:三角形、四辺形ともリファイン時に4分割されるので、整数を4で除算

                CBoundaryFace *pCrsBFace = pCrsBFaceMesh->getBFaceIX(nCrsIndex);
                double crsArea = pCrsBFace->getArea();

                for(uiint idof=0; idof < pBFaceMesh->getNumOfDOF(); idof++){
                    uiint dof = pBFaceMesh->getDOF(idof);
                    double val = pCrsBFace->getBndValue(dof);

                    val = val*(fnArea/crsArea);
                    pBFace->setBndValue(dof, val);//------- ファイングリッドに値をセット

                    Rcap_EquivalentNodalForce(iLevel, dof, pBFaceMesh, pBFace);//------ 等価節点力 配分
                };
            }

        };//iface loop
    };//ibnd loop
}
// 要素番号、面番号検索 : 外表面であるから要素は一個に限定される.
void CMW::Rcap_ElemFaceSearch(uiint& nElementID, uiint& nFaceID, vector<CNode*>& vNode, CMesh *pProgMesh)
{
    if(vNode.size() < 3){
        mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_ElemFaceSearch, vNode size error");
        return;
    }

    CAggregateElement *pAggElem0, *pAggElem1, *pAggElem2;
    CNode *pNode;
    
    pNode = vNode[0];
    pAggElem0 = pProgMesh->getAggElem(pNode->getID());
    
    pNode = vNode[1];
    pAggElem1 = pProgMesh->getAggElem(pNode->getID());
    
    pNode = vNode[2];
    pAggElem2 = pProgMesh->getAggElem(pNode->getID());


    vuint vShareID;
    uiint nTargetID;
    uiint nNNElem0 = pAggElem0->getNumOfElement();

    //1. Agg0 && Agg1 = vShareID
    for(uiint i=0; i < nNNElem0; i++){
        nTargetID = pAggElem0->get(i)->getID();

        uiint nNNElem1 = pAggElem1->getNumOfElement();
        for(uiint ii=0; ii < nNNElem1; ii++){
            if(pAggElem1->get(ii)->getID() == nTargetID){ vShareID.push_back(nTargetID); break;}
        };
    };

    //2. vShareID && Agg2 = vAnsID
    vuint vAnsID;
    for(uiint i=0; i < vShareID.size(); i++){
        nTargetID = vShareID[i];

        uiint nNNElem2 = pAggElem2->getNumOfElement();
        for(uiint ii=0; ii < nNNElem2; ii++){
            if(pAggElem2->get(ii)->getID() == nTargetID){ vAnsID.push_back(nTargetID); break;}
        };
    };

    // !要素番号:外表面の面集合であるから要素は一個に限定される.
    if(vAnsID.size() == 1){
        nElementID = vAnsID[0];//----------------- 面の属する要素ID
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_ElemFaceSearch ----- Unable to determine the Element_ID");
    }
    
    ////cout << "MW::Rcap_ElemFaceSearch ----- nElementID :" << nElementID << endl;// debug

    CElement *pElem = pProgMesh->getElement(nElementID);
    nFaceID = pElem->getFaceIndex(vNode[0], vNode[1], vNode[2]);//-------------------- 面番号
    
    ////cout << "MW::Rcap_ElemFaceSearch ----- nFaceID :" << nFaceID << endl;// debug
}
// 要素番号、辺番号検索 : 外Edgeであるから要素は一個に限定される.
void CMW::Rcap_ElemEdgeSearch(uiint& nElementID, uiint& nEdgeID, vector<CNode*>& vNode, CMesh *pProgMesh)
{
    if(vNode.size() < 2){
        mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_ElemEdgeSearch, vNode size error");
        return;
    }

    CAggregateElement *pAggElem0, *pAggElem1;
    CNode *pNode;
    
    pNode = vNode[0];
    pAggElem0 = pProgMesh->getAggElem(pNode->getID());

    pNode = vNode[1];
    pAggElem1 = pProgMesh->getAggElem(pNode->getID());

    vuint vAnsID;
    uiint nTargetID;
    uiint nNNElem0 = pAggElem0->getNumOfElement();

    // Agg0 && Agg1 = vAnsID
    for(uiint i=0; i < nNNElem0; i++){
        nTargetID = pAggElem0->get(i)->getID();

        uiint nNNElem1 = pAggElem1->getNumOfElement();
        for(uiint ii=0; ii < nNNElem1; ii++){
            if(pAggElem1->get(ii)->getID() == nTargetID){ vAnsID.push_back(nTargetID); break;}
        };
    };

    // !要素番号:辺に集合する要素は複数だが、一つを選べば済む.
    if(vAnsID.size() > 0){
        nElementID = vAnsID[0];//----------------- 辺の属する要素ID
    }else{
        uiint rank = mpMPI->getRank();
        mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_ElemEdgeSearch ----- Unable to determine the Element_ID,  rank:", rank);
    }
    
    ////cout << "MW::Rcap_ElemEdgeSearch ----- nElementID :" << nElementID << endl;// debug

    CElement *pElem = pProgMesh->getElement(nElementID);
    nEdgeID = pElem->getEdgeIndex(vNode[0], vNode[1]);

    ////cout << "MW::Rcap_ElemEdgeSearch ----- nEdgeID :" << nEdgeID << endl;// debug
}
// 要素番号、局所節点番号検索 : 表面の角の一点であるから一個の要素に限定される(1点通信テーブル)
void CMW::Rcap_ElemPointSearch(uiint& nElementID, uiint& nLocalID, CNode* pNode, CMesh* pProgMesh)
{
    CAggregateElement *pAggElem= pProgMesh->getAggElem(pNode->getID());

    if(pAggElem->getNumOfElement() == 0) cout << "MW::Rcap_ElemPointSearch ----- AggElement size Error" << endl;

    CElement *pElem= pAggElem->get(0);//----- 複数所有している中の一つの要素を取得

    nElementID= pElem->getID();//------------------------- 解:要素ID
    nLocalID= pElem->getLocalVertNum(pNode->getID());//--- 解:局所節点番号
}

// 要素番号検索
void CMW::Rcap_ElemSearch(uiint& nElementID, vector<CNode*>& vNode, CMesh *pProgMesh)
{
    vector<CAggregateElement*> vAggElem;
    CNode *pNode;
    CAggregateElement *pAggElem;

    //全Nodeの要素集合を集める
    for(uiint i=0; i < vNode.size(); i++){
        pNode = vNode[i];
        pAggElem = pProgMesh->getAggElem(pNode->getID());

        vAggElem.push_back(pAggElem);
    };

    vuint vPrevID;//節点が共有するElementIDを入れていく.Prev
    vuint vCurrID;//節点が共有するElementIDを入れていく.Current

    // vPrevID 初期化
    pAggElem = vAggElem[0];
    uiint nNumOfElem = pAggElem->getNumOfElement();
    for(uiint i=0; i < nNumOfElem; i++){
        vPrevID.push_back(pAggElem->get(i)->getID());//最初の節点のElementID集合
    };

    // vPrevID && AggElem = vCurrID 
    for(uiint i=1; i < vAggElem.size(); i++){//-------- i=1で始めてる
        pAggElem = vAggElem[i];//Current AggElement
        nNumOfElem = pAggElem->getNumOfElement();
        
        for(uiint ipre=0; ipre < vPrevID.size(); ipre++)
        for(uiint ii=0; ii < nNumOfElem; ii++)
            if(vPrevID[ipre]==pAggElem->get(ii)->getID()) vCurrID.push_back(vPrevID[ipre]);
        
        //Prev を Curr に入れ替え
        if(i != vAggElem.size()-1){
            vPrevID.clear();
            for(uiint ii=0; ii < vCurrID.size(); ii++) vPrevID.push_back(vCurrID[ii]);
            vCurrID.clear();
        }
    };

    // 解
    if(vCurrID.size()==1){
        nElementID = vCurrID[0];
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "CMW::Rcap_ElemSearch ----- Unable to determine the Element_ID");
    }
    
    ////cout << "MW::Rcap_ElemSearch ----- nElementID :" << nElementID << endl;// debug
}

//Neumann境界の等価節点力 配分 : Face
void CMW::Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryFaceMesh *pBFaceMesh, CBoundaryFace *pBFace)
{
    uiint ivert;
    CBoundaryNode *pBNode;
    double entVal, nodalVal, integVal, calcVal;

    // 等価節点力の初期化: Rcap_BFaceGene で処理済み

    CCalc *pCalc= CCalc::Instance();
    CPoland *pNumForm=NULL;
    //
    //等価節点力 分配加算
    //
    entVal = pBFace->getBndValue(dof);//------計算式の基礎データ
    double x=pBFace->getCenterX(), y=pBFace->getCenterY(), z=pBFace->getCenterZ();

    if( pBFaceMesh->existPoland(dof) ){
        pNumForm= pBFaceMesh->getPoland(dof);//Poland[dof]

        pCalc->setElementParam(entVal, x,y,z);
        calcVal = pCalc->Exec(pNumForm);

        cout << "MW::Rcap_EquivalentNodalForce  dof:" << dof << "  calcVal:" << calcVal << " rank:" << mpMPI->getRank() << endl;
    }

    switch(pBFace->getBFaceShape()){
        case(ElementType::Quad):
            for(ivert=0; ivert < 4; ivert++){
                integVal= mpShapeQuad->getIntegValue4(ivert);
                //nodalVal= entVal * integVal;
                if( pBFaceMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBFace->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Quad2):
            for(ivert=0; ivert < 8; ivert++){
                integVal= mpShapeQuad->getIntegValue8(ivert);
                //nodalVal= entVal * integVal;
                if( pBFaceMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBFace->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Triangle):
            for(ivert=0; ivert < 3; ivert++){
                integVal= mpShapeTriangle->getIntegValue3(ivert);
                //nodalVal= entVal * integVal;
                if( pBFaceMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBFace->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
        case(ElementType::Triangle2):
            for(ivert=0; ivert < 6; ivert++){
                integVal= mpShapeTriangle->getIntegValue6(ivert);
                //nodalVal= entVal * integVal;
                if( pBFaceMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBFace->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
    }//switch end
        
}
//Neumann境界の等価節点力 配分 : Edge
void CMW::Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryEdgeMesh *pBEdgeMesh, CBoundaryEdge *pBEdge)
{
    uiint ivert;
    CBoundaryNode *pBNode;
    double entVal, nodalVal, integVal, calcVal;


    // 等価節点力の初期化: Rcap_BEdgeGene で処理済み

    CCalc *pCalc= CCalc::Instance();
    CPoland *pNumForm=NULL;
    //
    //等価節点力 分配加算
    //
    entVal = pBEdge->getBndValue(dof);//境界値:基礎データ
    double x=pBEdge->getCenterX(),y=pBEdge->getCenterY(),z=pBEdge->getCenterZ();

    if( pBEdgeMesh->existPoland(dof) ){
        pNumForm= pBEdgeMesh->getPoland(dof);//Poland[dof]
        pCalc->setElementParam(entVal, x,y,z);
        calcVal = pCalc->Exec(pNumForm);
    }

    switch(pBEdge->getBEdgeShape()){
        case(ElementType::Beam):case(ElementType::Line):
            for(ivert=0; ivert < 2; ivert++){
                integVal= mpShapeLine->getIntegValue2(ivert);
                //nodalVal= entVal * integVal;
                if( pBEdgeMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBEdge->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Beam2):case(ElementType::Line2):
            for(ivert=0; ivert < 3; ivert++){
                integVal= mpShapeLine->getIntegValue3(ivert);
                //nodalVal= entVal * integVal;
                if( pBEdgeMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBEdge->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
    }//switch end

}
//Neumann境界の等価節点力 配分 : Volume
void CMW::Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryVolumeMesh *pBVolMesh, CBoundaryVolume *pBVol)
{
    uiint ivert;
    CBoundaryNode *pBNode;
    double entVal, nodalVal, integVal, calcVal;

    // 等価節点力の初期化: Rcap_BEdgeGene で処理済み

    CCalc *pCalc= CCalc::Instance();
    CPoland *pNumForm=NULL;
    //
    //等価節点力 分配加算
    //
    entVal = pBVol->getBndValue(dof);
    double x=pBVol->getCenterX(),y=pBVol->getCenterY(),z=pBVol->getCenterZ();

    if( pBVolMesh->existPoland(dof) ){
        pNumForm= pBVolMesh->getPoland(dof);//Poland[dof]
        pCalc->setElementParam(entVal, x,y,z);
        calcVal = pCalc->Exec(pNumForm);
    }

    switch(pBVol->getElemType()){
        case(ElementType::Hexa):
            for(ivert=0; ivert < 8; ivert++){
                integVal= mpShapeHexa->getIntegralValue8(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Hexa2):
            for(ivert=0; ivert < 20; ivert++){
                integVal= mpShapeHexa->getIntegralValue20(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
        case(ElementType::Tetra):
            for(ivert=0; ivert < 4; ivert++){
                integVal= mpShapeTetra->getIntegValue4(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Tetra2):
            for(ivert=0; ivert < 10; ivert++){
                integVal= mpShapeTetra->getIntegValue10(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
        case(ElementType::Prism):
            for(ivert=0; ivert < 6; ivert++){
                integVal= mpShapePrism->getIntegValue6(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            };
            break;
        case(ElementType::Prism2):
            for(ivert=0; ivert < 15; ivert++){
                integVal= mpShapePrism->getIntegValue15(ivert);
                //nodalVal= entVal * integVal;
                if( pBVolMesh->existPoland(dof) ){
                    nodalVal= calcVal * integVal;
                }else{
                    nodalVal= entVal * integVal;
                }
                pBNode= pBVol->getBNode(ivert);
                pBNode->addValue(dof, iLevel, nodalVal);
            }
            break;
    }//switch end
    
}

void CMW::Rcap_BFaceMeshDebug(vuint& vBndID, CMesh *pProgMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){
        
        uiint nBndID= vBndID[ibnd];
        CBoundaryFaceMesh *pBFaceMesh= pProgMesh->getBndFaceMeshID(nBndID);
        uiint nNumOfFace= pBFaceMesh->getNumOfBFace();
        
        ////cout << " BFaceMesh ID: " << nBndID  << endl;
        
        for(uiint iface=0; iface < nNumOfFace; iface++){
            CBoundaryFace *pBFace= pBFaceMesh->getBFaceIX(iface);
            uiint nNumOfBNode= pBFace->getNumOfBNode();

            string sStr0, sStr1;
            stringstream ss0, ss1;
            for(uiint ibnode=0; ibnode < nNumOfBNode; ibnode++){
                CBoundaryNode *pBNode= pBFace->getBNode(ibnode);
                CNode *pNode= pBNode->getNode();

                ss0 << pBNode->getID();           ss1 << pNode->getID();
                sStr0 +=  " " + ss0.str();  sStr1 +=  " " + ss1.str();
                ss0.str("");                      ss1.str("");
            }
            ////cout << ":" << iface << " BNode" << sStr0 << " | Node" << sStr1 << endl;
            sStr0.clear();  sStr1.clear();
        }
    }
}
void CMW::Rcap_BEdgeGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t *nEdgeCount, uiint nShape, size_t& nNumOfEdgeNode,
                         int32_t **refineNodes, CMesh *pProgMesh, CMesh *pCrsMesh)
{
    for(uiint ibnd=0; ibnd < vBndID.size(); ibnd++){

        uiint nBndID= vBndID[ibnd];
        CBoundaryEdgeMesh *pBEdgeMesh = pProgMesh->getBndEdgeMeshID(nBndID);
        CBoundaryEdgeMesh *pCrsBEdgeMesh= pCrsMesh->getBndEdgeMeshID(nBndID);

        //Neumannの場:Value初期化
        if(pBEdgeMesh->getBndType()==BoundaryType::Neumann){
            uiint nNumOfDOF = pBEdgeMesh->getNumOfDOF();
            uiint nNumOfBNode = pBEdgeMesh->getNumOfBNode();

            for(uiint inode=0; inode < nNumOfBNode; inode++){
                CBoundaryNode *pBNode= pBEdgeMesh->getBNodeIX(inode);
                for(uiint idof=0; idof < nNumOfDOF; idof++){
                    uiint dof = pBEdgeMesh->getDOF(idof);

                    pBNode->initValue(dof, iLevel);//-------------Neumann初期化
                };
            };
        }

        //EdgeNodeの取得
        for(uiint iedge=0; iedge < nEdgeCount[ibnd]; iedge++){

            CBoundaryEdge *pBEdge = new CBoundaryEdge();//-------------- 生成
            pBEdge->resizeBNode(nNumOfEdgeNode);

            vector<CNode*> vNode;     // 要素,辺 : 検索パラメータ
            uiint nElementID, nEdgeID;// 要素,辺 : 解

            for(uiint inode=0; inode < nNumOfEdgeNode; inode++){
                uiint nNodeNum = refineNodes[ibnd][nNumOfEdgeNode*iedge + inode];
                uiint index= vmID2Index[ibnd][nNodeNum];

                CBoundaryNode *pBNode= pBEdgeMesh->getBNodeIX(index);

                if(nShape==ElementType::Beam || nShape==ElementType::Beam2){
                    if(inode < 2){
                        pBEdge->setBNode(inode, pBNode);
                    }else{
                        uiint iedge= inode-2;
                        pBEdge->setEdgeBNode(pBNode);
                    }
                }
                CNode *pNode = pBNode->getNode();
                vNode.push_back(pNode);

            };//inode loop

            Rcap_ElemEdgeSearch(nElementID, nEdgeID, vNode, pProgMesh);//---------- 要素、辺 検索
            
            pBEdge->setID(iedge);
            pBEdge->setBEdgeShape(nShape);
            pBEdge->setElementID(nElementID);
            pBEdge->setElementEdgeID(nEdgeID);

            CElement *pElem = pProgMesh->getElement(nElementID);
            pBEdge->setElement(pElem);

            double fnLength = pBEdge->calcLength();

            pBEdgeMesh->addBEdge(pBEdge);

            // コースグリッドのBFaceのVal => ファイングリッドBFaceへセット
            // -- Neumann --
            if(pBEdgeMesh->getBndType()==BoundaryType::Neumann){

                uiint nCrsIndex = iedge/2;//------- コースグリッド辺番号:2分割なので、整数を2で除算

                CBoundaryEdge *pCrsBEdge = pCrsBEdgeMesh->getBEdgeIX(nCrsIndex);
                double crsLength = pCrsBEdge->getLength();

                for(uiint idof=0; idof < pBEdgeMesh->getNumOfDOF(); idof++){
                    uiint dof = pBEdgeMesh->getDOF(idof);
                    double val = pCrsBEdge->getBndValue(dof);

                    val = val*(fnLength/crsLength);
                    pBEdge->setBndValue(dof, val);//------- ファイングリッドに値をセット

                    Rcap_EquivalentNodalForce(iLevel, dof, pBEdgeMesh, pBEdge);//------ 等価節点力配分
                };
            }

        };//iedge loop
    };//ibnd loop
}

void CMW::Rcap_SetCommMesh2(uiint iLevel, map<uiint,uiint> mComID, vvuint& vvNum4EachType, CMesh *pCrsMesh, CMesh *pProgMesh)
{
    uiint nNumOfType= 7;
    
    for(uiint icom=0; icom < mComID.size(); icom++){

        CCommMesh2 *pCommMesh2 = new CCommMesh2;//--------------CommMesh2の生成
        
        uiint nComID = mComID[icom];
        uiint nFaceCount=0;
        uiint nNumOfCommNode=0;
        
        for(uiint itype=0; itype < nNumOfType; itype++){
            nFaceCount += vvNum4EachType[icom][itype];

            // 節点数が重複しているが、そのままリザーブに利用
            if(itype==0) nNumOfCommNode += vvNum4EachType[icom][itype]*4;//Quad
            if(itype==1) nNumOfCommNode += vvNum4EachType[icom][itype]*8;//Quad2
            if(itype==2) nNumOfCommNode += vvNum4EachType[icom][itype]*3;//Tri
            if(itype==3) nNumOfCommNode += vvNum4EachType[icom][itype]*6;//Tri2
            if(itype==4) nNumOfCommNode += vvNum4EachType[icom][itype]*2;//Beam
            if(itype==5) nNumOfCommNode += vvNum4EachType[icom][itype]*3;//Beam2
            if(itype==6) nNumOfCommNode += vvNum4EachType[icom][itype]*1;//Point
        };
        
        pCommMesh2->setLevel(iLevel);
        pCommMesh2->setID(nComID);
        pCommMesh2->reserveCommFace(nFaceCount);
        pCommMesh2->reserveCommNode(nNumOfCommNode);

        CCommMesh2 *pCrsCommMesh2 = pCrsMesh->getCommMesh2(nComID);
        uiint nRank = pCrsCommMesh2->getRank();
        uiint nTransRank= pCrsCommMesh2->getTrasmitRank();

        pCommMesh2->setRank(nRank);
        pCommMesh2->setTransmitRank(nTransRank);

        pProgMesh->setCommMesh2(pCommMesh2);
    }
}
//--
// int32_t** => vvuint 移し替え
//--
void CMW::Rcap_NodesPoint2vuint(map<uiint,uiint>& mComID, uiint nNumOfFaceNode, size_t *nFaceCount, int32_t **refineNodes, vvuint& vvRefineNodes)
{
    vvRefineNodes.resize(mComID.size());

    for(uiint icom=0; icom < mComID.size(); icom++){
        uiint nMax = nFaceCount[icom] * nNumOfFaceNode;
        vvRefineNodes[icom].resize(nMax);

        for(uiint index=0; index < nMax; index++)
            vvRefineNodes[icom][index] = refineNodes[icom][index];
        
        ////cout << "CMW::Rcap_NodesPoint2vuint icom:" << icom  << " face_count:" << nFaceCount[icom]
        ////        << " face_node数:" << nNumOfFaceNode << " size:" << vvRefineNodes[icom].size() << endl;
    };
}
//--
// 形状別リファイン・ノード配列を一つにまとめる:順序=Quad,Quad2,Tri,Tri2,Beam,Beam2,Point
//--
void CMW::Rcap_SumRefineNodes(vvuint& vvCommRefineNodes, vvuint& vvQuad, vvuint& vvQuad2, vvuint& vvTri, vvuint& vvTri2,
                            vvuint& vvBeam, vvuint& vvBeam2, vvuint& vvPoint)
{
    uiint nNumOfType=7;

    vvCommRefineNodes.resize(vvQuad.size());//CommMesh2数==通信テーブル数

    for(uiint icom=0; icom < vvQuad.size(); icom++){
        for(uiint itype=0; itype < nNumOfType; itype++){
            if(itype==0){
                for(uiint index=0; index < vvQuad[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvQuad[icom][index]);
            }else if(itype==1){
                for(uiint index=0; index < vvQuad2[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvQuad2[icom][index]);
            }else if(itype==2){
                for(uiint index=0; index < vvTri[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvTri[icom][index]);
            }else if(itype==3){
                for(uiint index=0; index < vvTri2[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvTri2[icom][index]);
            }else if(itype==4){
                for(uiint index=0; index < vvBeam[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvBeam[icom][index]);
            }else if(itype==5){
                for(uiint index=0; index < vvBeam2[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvBeam2[icom][index]);
            }else if(itype==6){
                for(uiint index=0; index < vvPoint[icom].size(); index++)
                    vvCommRefineNodes[icom].push_back(vvPoint[icom][index]);
            }
        };//for(itype)end
        
        ////cout << "CMW::Rcap_SumRefineNodes  icom:" << icom << " size:" << vvCommRefineNodes[icom].size() << endl;
        
    };//for(icom)end
    
    
}
// Type別の面数 :コースグリッド
void CMW::Rcap_CommFaceCount( vvuint& vvNum4EachType, map<uiint,uiint> mComID, size_t* commQuadCount, size_t* commQuad2Count,
                              size_t* commTriCount, size_t* commTri2Count, size_t* commBeamCount, size_t* commBeam2Count,
                              size_t* commPointCount)
{
    uiint nNumOfType=7;

    vvNum4EachType.resize(mComID.size());
    for(uiint icom=0; icom < mComID.size(); icom++){
        vvNum4EachType[icom].resize(nNumOfType);
        for(uiint itype=0; itype < nNumOfType; itype++){
            if(itype==0) vvNum4EachType[icom][itype]= commQuadCount[icom];
            if(itype==1) vvNum4EachType[icom][itype]= commQuad2Count[icom];
            if(itype==2) vvNum4EachType[icom][itype]= commTriCount[icom];
            if(itype==3) vvNum4EachType[icom][itype]= commTri2Count[icom];
            if(itype==4) vvNum4EachType[icom][itype]= commBeamCount[icom];
            if(itype==5) vvNum4EachType[icom][itype]= commBeam2Count[icom];
            if(itype==6) vvNum4EachType[icom][itype]= commPointCount[icom];
        };
    };
}
//-----------------------
// 全タイプのCommFaceを処理
//-----------------------
// Quad,Quad2,Triangle,Triangle2,Beam,Beam2,Point を一度に処理する.
// # vNum4EachType:Type別の要素数 => 最後に次ステップのためにリファイン数に変更
//--
void CMW::Rcap_CommNodeGene( map<uiint,uiint> mComID, vvuint& vvNum4EachType, vvvuint& vvvElemNum, vector<map<uiint, vuint> >& vmaElemNIndex,
                             CMesh *pProgMesh, vvuint& vvRefineNodes,
                             vector<map<uiint, uiint> >& maNNum2CommNNum)
{
    //面の倍数、辺の倍数
    uiint nFaceMultiple=4, nEdgeMultiple=2;

    ////if( mpMPI->getRank()==0)
    ////    cout << "CMW::Rcap_CommNodeGene ------ A " << endl;//<< " rank:" << mpMPI->getRank() << endl;

    uiint nNumOfType=7;

    //--
    // コースグリッドの面数を、リファイン後の数に変更
    //--
    // # vvNum4EachType:コースグリッドでの面数
    for(uiint icom=0; icom < mComID.size(); icom++){
        for(uiint itype=0; itype < nNumOfType; itype++){

            uiint nNumOfFace= vvNum4EachType[icom][itype];//タイプ別の面数(Quad数、…Triangle数、… Beam数、…Point数)
            uiint nRefineNumOfFace;

            if(itype==0) nRefineNumOfFace=nNumOfFace*nFaceMultiple;//QUAD  (RCAP_QUAD)
            if(itype==1) nRefineNumOfFace=nNumOfFace*nFaceMultiple;//QUAD2 (RCAP_QUAD2)
            if(itype==2) nRefineNumOfFace=nNumOfFace*nFaceMultiple;//Triangle  (RCAP_TRIANGLE)
            if(itype==3) nRefineNumOfFace=nNumOfFace*nFaceMultiple;//Triangle2 (RCAP_TRIANGLE2)
            if(itype==4) nRefineNumOfFace=nNumOfFace*nEdgeMultiple;//Beam  (RCAP_SEGMENT)
            if(itype==5) nRefineNumOfFace=nNumOfFace*nEdgeMultiple;//Beam2 (RCAP_SEGMENT2)
            if(itype==6) nRefineNumOfFace=nNumOfFace*1;//Point (RCAP_VERTEX) <--------------- 点はリファインしない

            vvNum4EachType[icom][itype]=nRefineNumOfFace;//リファイン後のFace数
        };
    };

    ////if( mpMPI->getRank()==0)
    ////    cout << "CMW::Rcap_CommNodeGene ------ B " << endl; //<< " rank:" << mpMPI->getRank() << endl;
    
    //--
    // 面ごとに(仮)要素番号を割り振る: [icom][itype][面数]=(仮)要素番号
    //--
    vvvElemNum;//[icom][itype][面数]=(仮)要素番号  : CommNode & CommFace生成の両方で利用.
    vvvElemNum.resize(mComID.size());
    
    for(uiint icom=0; icom < mComID.size(); icom++){
        uiint nElemCount=0;
        vvvElemNum[icom].resize(nNumOfType);
        for(uiint itype=0; itype < nNumOfType; itype++){
            uiint nNumOfFace= vvNum4EachType[icom][itype];
            vvvElemNum[icom][itype].resize(nNumOfFace);
            
            for(uiint iface=0; iface < nNumOfFace; iface++){
                vvvElemNum[icom][itype][iface]=nElemCount++;//-----全体の面数をカウントして(仮)要素番号を割り振る.
            };
        };
    };

    ////if( mpMPI->getRank()==0)
    ////    cout << "CMW::Rcap_CommNodeGene ------ C " << endl;//<< " rank:" << mpMPI->getRank() << endl;

    //--
    // Nodeの属性として形状別の構成ノード数を与える
    //--
    // # どの形状タイプに属しているか => 形状別の構成ノード数
    // #注意) 節点は全て要素別に並んでいる => 重複して並んでいるので、必ずどこかの面に属している
    //
    vvuint vvNodesBaseNum;//------------- Node index ごとのBaseNum:形状の構成ノード数
    vvNodesBaseNum.resize(mComID.size());

    for(uiint icom=0; icom < mComID.size(); icom++) vvNodesBaseNum[icom].resize(vvRefineNodes[icom].size());
    
    vuint vQuadEnd,vQuad2End,vTriEnd,vTri2End,vBeamEnd,vBeam2End,vPointEnd;//icom別 index終端

    for(uiint icom=0; icom < mComID.size(); icom++){
        uiint nQuadEnd,nQuad2End,nTriEnd,nTri2End,nBeamEnd,nBeam2End,nPointEnd;//各形状のindex終端
        uiint nNodeStart, nNodeEnd;

        ////if( mpMPI->getRank()==0)
        ////    cout << "CMW::Rcap_CommNodeGene ------ icom:" << icom << " comID:" << mComID[icom] << " size:" << vvRefineNodes[icom].size() << endl;//<< " rank:" << mpMPI->getRank() << endl;
        
        for(uiint itype=0; itype < nNumOfType; itype++){

            uiint nNumOfFace= vvNum4EachType[icom][itype];//タイプ別の面数(Quad数、…Triangle数、… Beam数、…Point数)
            
            if(itype==0){
                nNodeStart=0; nNodeEnd= nNumOfFace*4;
                nQuadEnd=nNodeEnd;
                vQuadEnd.push_back(nQuadEnd);
                ////cout << "CMW::Rcap_CommNodeGene ------ Quad  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 4;//Quad:4
            }else if(itype==1){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*8;
                nQuad2End=nNodeEnd;
                vQuad2End.push_back(nQuad2End);
                ////cout << "CMW::Rcap_CommNodeGene ------ Quad2  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 8;//Quad2:8
            }else if(itype==2){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*3;
                nTriEnd=nNodeEnd;
                vTriEnd.push_back(nTriEnd);
                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene ------ Tri  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " comID:" << mComID[icom] << endl;//<< " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 3;//Triangle:3
            }else if(itype==3){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*6;
                nTri2End=nNodeEnd;
                vTri2End.push_back(nTri2End);
                ////cout << "CMW::Rcap_CommNodeGene ------ Tri2  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 6;//Triangle2:6
            }else if(itype==4){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*2;
                nBeamEnd=nNodeEnd;
                vBeamEnd.push_back(nBeamEnd);
                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene ------ Beam  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " comID:" << mComID[icom] << endl;//<< " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 2;//Beam:2
            }else if(itype==5){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*3;
                nBeam2End=nNodeEnd;
                vBeam2End.push_back(nBeam2End);
                ////cout << "CMW::Rcap_CommNodeGene ------ Beam2  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 3;//Beam2:3
            }else if(itype==6){
                nNodeStart=nNodeEnd; nNodeEnd= nNodeStart + nNumOfFace*1;
                nPointEnd=nNodeEnd;
                vPointEnd.push_back(nPointEnd);
                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene ------ Point  NodeStart:" << nNodeStart << " NodeEnd:" << nNodeEnd << " comID:" << mComID[icom] << endl;//<< " rank:" << mpMPI->getRank() << endl;
                for(uiint index=nNodeStart; index < nNodeEnd; index++) vvNodesBaseNum[icom][index]= 1;//Point:1
            }else{
                //ne-yo
                mpLogger->Info(Utility::LoggerMode::Error, "MW::Rcap_CommNodeGene, vvNodesType setting error.");
            }
        };//for(itype)end
    };//for(icom)end

    ////if( mpMPI->getRank()==0)
    ////    cout << "CMW::Rcap_CommNodeGene ------ D" << endl;//<< " rank:" << mpMPI->getRank() << endl;

    //--
    // 参照引数で渡されたvector<map>のリサイズ
    //--
    maNNum2CommNNum.resize(mComID.size());

    //--
    //
    //--
    for(uiint icom=0; icom < mComID.size(); icom++){
        uiint nComID = mComID[icom];
        CCommMesh2 *pCommMesh2 = pProgMesh->getCommMesh2(nComID);
        
        //ファイングリッドのNode番号の重複番号を削除
        // ---
        // # vvRefineNodes[icom]は(仮)要素並びでNode番号を並べている => Hashデータ生成 :本物要素ではなくて、面を要素と見立てた(仮)要素.
        // ---
        uiint nBaseNum;                // CommFaceの構成ノード数
        map< uiint, vuint > maNodeElem;  // [REVOCAP_RefineのNode番号][Node所属要素数] = 要素番号
        map< uiint, vuint > maElemNIndex;// (仮)要素番号ごとのindex配列 :: [(仮)要素番号][構成ノード数] = index : indexは,RefineNodes[index]のindex
        uiint nNumOfIndex= vvRefineNodes[icom].size();

        for(uiint index=0; index < nNumOfIndex; index++){
            //
            // [icom][itype][iface]=(仮)要素番号 :=> Node番号別のデータ(maNodeElem)に変換
            //
            nBaseNum= vvNodesBaseNum[icom][index];

            uiint nElemNumT;// (仮仮)要素番号 : 面,辺を要素とみなした場合の要素番号
            uiint itype, iface;

            if(index < vQuadEnd[icom]){//------------------------Quad
                itype= 0;
                iface= index/nBaseNum;// Quad 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];

                maElemNIndex[nElemNumT].push_back(index);// (仮)要素ごとのindexが出る変数
            }
            if(index >= vQuadEnd[icom] && index < vQuad2End[icom]){//--Quad2
                itype = 1;
                iface = (index-vQuadEnd[icom])/nBaseNum;//Quad2 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];
                
                maElemNIndex[nElemNumT].push_back(index);
            }
            if(index >= vQuad2End[icom] && index < vTriEnd[icom]){//--Triangle
                itype = 2;
                iface = (index-vQuad2End[icom])/nBaseNum;//Tri 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];

                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene -------D0 Tri elem_num:" << nElemNumT
                ////            << " comID:" << nComID << " index:" << index << " iface:" << iface << endl;//<< " rank:" << mpMPI->getRank() << endl;
                
                maElemNIndex[nElemNumT].push_back(index);
            }
            if(index >= vTriEnd[icom] && index < vTri2End[icom]){//--Triangle2
                itype = 3;
                iface = (index-vTriEnd[icom])/nBaseNum;//Tri2 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];
                
                maElemNIndex[nElemNumT].push_back(index);
            }
            if(index >= vTri2End[icom] && index < vBeamEnd[icom]){//--Beam
                itype = 4;
                iface = (index-vTri2End[icom])/nBaseNum;//Beam 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];

                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene -------D0 Beam elem_num:" << nElemNumT
                ////            << " comID:" << nComID << " index:" << index << " iface:" << iface << endl;//<< " rank:" << mpMPI->getRank() << endl;

                maElemNIndex[nElemNumT].push_back(index);
            }
            if(index >= vBeamEnd[icom] && index < vBeam2End[icom]){//--Beam2
                itype = 5;
                iface = (index-vBeamEnd[icom])/nBaseNum;//Beam2 面番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];
                
                maElemNIndex[nElemNumT].push_back(index);
            }
            if(index >= vBeam2End[icom]){//--Point
                itype = 6;
                iface = index-vBeam2End[icom];//Point 面番号 番号(0番から)
                nElemNumT = vvvElemNum[icom][itype][iface];

                ////if( mpMPI->getRank()==0)
                ////    cout << "CMW::Rcap_CommNodeGene -------D0 Point elem_num:" << nElemNumT
                ////            << " comID:" << nComID << " index:" << index << " iface:" << iface << endl;//<< " rank:" << mpMPI->getRank() << endl;

                maElemNIndex[nElemNumT].push_back(index);
            }

            uiint nNodeNum = vvRefineNodes[icom][index];//REVOCAP_RefinerのMesh節点番号(通し番号)
            maNodeElem[nNodeNum].push_back(nElemNumT);  //REVOCAP_Refineの節点番号:=> (仮)要素番号

        };//for(index)end
        

        vmaElemNIndex.push_back(maElemNIndex);//-------- CommFace生成だけで利用(CommNode生成は無関係)

        
        bool* vMarking = new bool[nNumOfIndex];
        for(uiint index=0; index < nNumOfIndex; index++)
            vMarking[index]=true;//--------------------------------- マーキング初期化
        //--
        // 重複節点のマーキング
        //--
        for(uiint index=0; index < nNumOfIndex; index++){

            nBaseNum= vvNodesBaseNum[icom][index];//------(仮)要素の構成ノード数
            uiint nNodeNum = vvRefineNodes[icom][index];//節点番号(REVOCAP_RefineのMesh節点番号)

            if(vMarking[index]){
                for(uiint ielem=0; ielem < maNodeElem[nNodeNum].size(); ielem++){
                    uiint nElemNum = maNodeElem[nNodeNum][ielem];//Nodeが所属する(仮)要素番号
                    //--
                    // 複数の要素に所属している場合に無効マーキング
                    //--
                    // ! (仮)要素ごとのindexが出る変数: maElemNIndex : [(仮)要素番号]=vuint : 構成ノードindex配列 : indexはRefineNodes[index]のindex
                    if(ielem > 0){
                        for(uiint k=0; k < maElemNIndex[nElemNum].size(); k++ ){
                            uiint ix = maElemNIndex[nElemNum][k];
                            if(vvRefineNodes[icom][ix] == nNodeNum ){
                                vMarking[ix] = false;//--- 同一Node番号を無効マーキング
                            }
                        };
                    }//if(ielem > 0)
                };
            }//if(Marking==true)
        };

        ////if( mpMPI->getRank()==0)
        ////    cout << "CMW::Rcap_CommNodeGene ------ F comID:" << nComID << endl;//<< " rank:" << mpMPI->getRank() << endl;

        //----
        // 各階層で個別にCommNode配列を生成
        // # コースグリッドのCommNodeは無視!!!!!!!
        //----
        uiint nComNodeNum=0;//CommNode IDカウンター
        for(uiint index=0; index < nNumOfIndex; index++){

            uiint nNodeNum = vvRefineNodes[icom][index];//REVOCAP_RefineのMesh節点番号 == MW3の節点Index番号

            if(vMarking[index]){
                CCommNode *pCommNode = new CCommNode;//----------- CCommNode生成
                pCommNode->setID(nComNodeNum);
                maNNum2CommNNum[icom][nNodeNum]=nComNodeNum;//-------- NodeNum => CommNodeNum : CommFace生成で利用する.
                nComNodeNum++;

                CNode *pNode= pProgMesh->getNodeIX(nNodeNum);
                pCommNode->setNode(pNode);
                pCommNode->setCoord(pNode->getCoord());

                pCommMesh2->addCommNode(pCommNode);
            }//if Marking
        };
        delete []vMarking;
        
        ////if( mpMPI->getRank()==0)
        ////    cout << "CMW::Rcap_CommNodeGene ------ G comID:" << nComID << endl;//<< " rank:" << mpMPI->getRank() << endl;
        
    };// for(icom)end
    
}
//----
// CommFace生成 # Rcapの場合は通信に無関係 : データチェックのために生成
//----
// # 変更が生じて面倒な時は、CommFace生成を丸ごと削除しても計算には差し支えない.
//----
void CMW::Rcap_CommFaceGene( uiint iLevel, map<uiint,uiint>& mComID, 
                            vvvuint& vvvElemNum, vector<map<uiint, vuint> >& vmaElemNIndex, 
                            vvuint& vvRefineNodes, CMesh *pProgMesh,
                            vector<map<uiint, uiint> >& maNNum2CommNNum )
{
    uiint nNumOfType=7;

    for(uiint icom=0; icom < mComID.size(); icom++){

        uiint nComID = mComID[icom];
        CCommMesh2 *pCommMesh2 = pProgMesh->getCommMesh2(nComID);

        ////cout << "CMW::Rcap_CommFaceGene ------ A" << " rank:" << mpMPI->getRank() << endl;
        
        for(uiint itype=0; itype < nNumOfType; itype++){
            uiint nNumOfVert, nNumOfEdge, nOrder;
            switch(itype){
                case(0):  nNumOfVert=NumberOfVertex::Quad();     nNumOfEdge=NumberOfEdge::Quad();     nOrder=ElementOrder::First;  break;
                case(1):  nNumOfVert=NumberOfVertex::Quad();     nNumOfEdge=NumberOfEdge::Quad();     nOrder=ElementOrder::Second;  break;
                case(2):  nNumOfVert=NumberOfVertex::Triangle(); nNumOfEdge=NumberOfEdge::Triangle(); nOrder=ElementOrder::First;  break;
                case(3):  nNumOfVert=NumberOfVertex::Triangle(); nNumOfEdge=NumberOfEdge::Triangle(); nOrder=ElementOrder::Second;  break;
                case(4):  nNumOfVert=NumberOfVertex::Beam();     nNumOfEdge=NumberOfEdge::Beam();     nOrder=ElementOrder::First;  break;
                case(5):  nNumOfVert=NumberOfVertex::Beam();     nNumOfEdge=NumberOfEdge::Beam();     nOrder=ElementOrder::Second;  break;
                case(6):  nNumOfVert=NumberOfVertex::Point();    nNumOfEdge=NumberOfEdge::Point();    nOrder=ElementOrder::Zero;  break;
            }//switch(itype)end

            for(uiint iface=0; iface < vvvElemNum[icom][itype].size(); iface++){
                uiint nElemNumT;

                CCommFace *pCommFace = new CCommFace;//---------------------- CommFace生成
                pCommFace->initialize( nNumOfVert, nNumOfEdge, nOrder);
                
                nElemNumT = vvvElemNum[icom][itype][iface];//(仮)要素番号 => CommFace ID番号

                pCommFace->setID(nElemNumT);
                pCommFace->setMGLevel(iLevel);//------------------ 引数のiLevelは、ここだけ使用.

                vector<CNode*> vNode;     // 要素、面番号: 検索パラメータ
                uiint nElementID, nFaceID;// 要素、面番号(実際は形状による：面番号、辺番号、局所節点番号): 解

                uiint nNumOfFaceNode= vmaElemNIndex[icom][nElemNumT].size();
                //局所節点ループ
                for(uiint k=0; k < nNumOfFaceNode; k++){
                    uiint index = vmaElemNIndex[icom][nElemNumT][k];    //REVOCAP_Refine: RefineNodesのインデックス
                    uiint nNodeNum = vvRefineNodes[icom][index];        //REVOCAP_Refine: Meshの節点番号==MW3の節点INDEX
                    uiint nComNodeNum = maNNum2CommNNum[icom][nNodeNum];//CommNodeのID番号取得
                    
                    CCommNode *pCommNode = pCommMesh2->getCommNode(nComNodeNum);
                    pCommFace->setCommNode( k, pCommNode);
                    
                    //辺ノードの設定
                    if(k >= nNumOfVert){
                        uiint iedge;
                        if(itype==1)//--------------QUAD2
                            iedge = k - nNumOfVert;
                        if(itype==3){//-------------Triangle2
                            if(k==3) iedge=1;
                            if(k==4) iedge=2;
                            if(k==5) iedge=0;
                        }
                        if(itype==5) iedge=0;//-----Beam2
                        
                        pCommFace->setEdgeCommNode(pCommNode, iedge);//(2次の場合):辺上のCommNode
                    }

                    vNode.push_back(pCommNode->getNode());//--- vNode:要素検索に用いる.

                };//for(k)end
                
                pCommMesh2->addCommFace(pCommFace);//--------------- CommMesh2にCommFaceを保存(追加)

                //要素ID, 要素面番号(面番号、辺番号、局所節点番号)取得
                switch(itype){
                    case(0): case(1): case(2): case(3)://--Quad,Quad2,Triangle,Triangle2
                        Rcap_ElemFaceSearch(nElementID, nFaceID, vNode, pProgMesh);
                        break;
                    case(4): case(5)://--------------------Beam,Beam2
                        Rcap_ElemEdgeSearch(nElementID, nFaceID, vNode, pProgMesh);
                        break;
                    case(6)://-----------------------------Point
                        Rcap_ElemPointSearch(nElementID, nFaceID, vNode[0], pProgMesh);//-- pNode 1点
                        break;
                    default:
                        break;
                }
                pCommFace->setElementID(nElementID); //要素ID
                ////  pCommFace->setElementFaceID(nFaceID);// × 使わなくなった: 面:要素面番号、辺:要素辺番号、点:局所節点番号
                
            };//for(iface)end
        };//for(itype)end
        //----------------------------CommFace生成 (Rcapの場合は通信には無関係:データチェック用) END

        ////cout << "CMW::Rcap_CommFaceGene ------ C" << " rank:" << mpMPI->getRank() << endl;

    };// for(icom) end
}

////void CMW::Rcap_CommFaceGene(uiint iLevel, map<uiint,uiint>& mComID, int8_t nType, size_t *nFaceCount, int32_t **refineNodes, CMesh *pProgMesh, vector<map<uiint, uiint> >& maNNum2CommNNum )
////{
////    uiint nNumOfFaceNode, nNumOfVert, nNumOfEdge, nOrder;
////
////    switch(nType){
////        case(RCAP_QUAD):  nNumOfFaceNode=4; nNumOfVert=4; nNumOfEdge=4; nOrder=ElementOrder::First;  break;
////        case(RCAP_QUAD2): nNumOfFaceNode=8; nNumOfVert=4; nNumOfEdge=4; nOrder=ElementOrder::Second;  break;
////        case(RCAP_TRIANGLE):  nNumOfFaceNode=3; nNumOfVert=3; nNumOfEdge=3; nOrder=ElementOrder::First;  break;
////        case(RCAP_TRIANGLE2): nNumOfFaceNode=6; nNumOfVert=3; nNumOfEdge=3; nOrder=ElementOrder::Second;  break;
////        case(RCAP_SEGMENT):  nNumOfFaceNode=2; nNumOfVert=2; nNumOfEdge=1; nOrder=ElementOrder::First;  break;
////        case(RCAP_SEGMENT2): nNumOfFaceNode=3; nNumOfVert=2; nNumOfEdge=1; nOrder=ElementOrder::Second;  break;
////        case(RCAP_VERTEX): nNumOfFaceNode=1; nNumOfVert=1; nNumOfEdge=0; nOrder=ElementOrder::Zero;  break;
////    }
////
////    for(uiint icom=0; icom < mComID.size(); icom++){
////        uiint nComID = mComID[icom];
////
////        CCommMesh2 *pCommMesh2 = pProgMesh->getCommMesh2(nComID);
////
////        for(uiint iface=0; iface < nFaceCount[icom]; iface++){
////
////            CCommFace *pCommFace= new CCommFace;//----------- CommFace (通信エンティティ) 生成
////            pCommFace->setID(iface);
////            pCommFace->setMGLevel(iLevel);
////            pCommFace->initialize(nNumOfVert, nNumOfEdge, nOrder);
////
////            vector<CNode*> vNode;     // 要素、面番号: 検索パラメータ
////            uiint nElementID, nFaceID;// 要素、面番号(実際は形状による：面番号、辺番号、局所節点番号): 解
////
////            for(uiint inode=0; inode < nNumOfFaceNode; inode++){
////                uiint nNodeNum = refineNodes[icom][ iface*nNumOfFaceNode + inode];
////                uiint nComNodeNum = maNNum2CommNNum[icom][nNodeNum];
////
////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(nComNodeNum);
////
////                pCommFace->setCommNode(inode, pCommNode);//-------------- 面全てのCommNode
////
////                //辺ノードの設定
////                if(inode >= nNumOfVert){
////                    uiint iedge;
////                    if(nType==RCAP_QUAD2 || nType==RCAP_SEGMENT2)
////                        iedge = inode - nNumOfVert;
////                    if(nType==RCAP_TRIANGLE2){
////                        if(inode==3) iedge=1;
////                        if(inode==4) iedge=2;
////                        if(inode==5) iedge=0;
////                    }
////                    if(nType==RCAP_SEGMENT2) iedge=0;
////
////                    pCommFace->setEdgeCommNode(pCommNode, iedge);//(2次の場合):辺上のCommNode
////                }
////
////                vNode.push_back(pCommNode->getNode());//--- vNode
////
////            };//inode loop
////
////            //要素ID, 要素面番号(面番号、辺番号、局所節点番号)取得
////            switch(nType){
////                case(RCAP_QUAD): case(RCAP_QUAD2): case(RCAP_TRIANGLE): case(RCAP_TRIANGLE2):
////                    Rcap_ElemFaceSearch(nElementID, nFaceID, vNode, pProgMesh);
////                    break;
////                case(RCAP_SEGMENT): case(RCAP_SEGMENT2):
////                    Rcap_ElemEdgeSearch(nElementID, nFaceID, vNode, pProgMesh);
////                    break;
////                case(RCAP_VERTEX):
////                    Rcap_ElemPointSearch(nElementID, nFaceID, vNode[0], pProgMesh);//-- pNode 1点
////                    break;
////                default:
////                    break;
////            }
////
////            pCommFace->setElementID(nElementID); //要素ID
////            ////  pCommFace->setElementFaceID(nFaceID);// × 使わなくなった: 面:要素面番号、辺:要素辺番号、点:局所節点番号
////
////            pCommMesh2->addCommFace(pCommFace);
////
////        };//iface loop
////    };//icom loop
////}

       //=======================================================================
#endif //======================================================== REVOCAP_REFINE
       //=======================================================================

void CMW::SortMerge(vuint& vec)
{
    std::sort(vec.begin(), vec.end());
    std::vector<uiint>::iterator new_end = std::unique(vec.begin(), vec.end());
    vec.erase(new_end, vec.end());
}

void CMW::clearMW3(const uiint& mgLevel)
{
    mpAssy = mpGMGModel->getAssyModel(mgLevel);

    uiint nNumOfMesh = mpAssy->getNumOfMesh();
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
        mpMesh = mpAssy->getMesh(imesh);
        mpMesh->clear();
    };

    uiint nNumOfConMesh = mpAssy->getNumOfContactMesh();
    for(uiint icmesh=0; icmesh < nNumOfConMesh; icmesh++){
        mpConMesh = mpAssy->getContactMesh(icmesh);
        mpConMesh->clear();
    }
}

uiint CMW::Finalize()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();
    
    mpMPI->Finalize();

    return 1;
}
uiint CMW::FileRead(string& basename, bool bBinary)
{
    mpFileIO->setBaseName(basename);
    msInputFileName += basename + ".";
    msOutputFileName+= basename + ".";
    msResFileName   += basename + ".";
    msRltFileName   += basename + ".";
    uiint rank= mpMPI->getRank();
    stringstream ss;
    ss << rank;
    msInputFileName  += ss.str();
    msOutputFileName += ss.str();
    msResFileName    += ss.str();
    msRltFileName    += ss.str();
    msInputFileName  += ".msh";
    msOutputFileName += ".out";
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");
    mpFileIO->ReadFile(msInputFileName, bBinary);
    mb_file = true;
    return 1;
}
uiint CMW::FileRead_fstr(bool bBinary)
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead_fstr");
    mpFileIO->ReadFile(msInputFileName, bBinary);
    mb_file = true;
    return 1;
}
uiint CMW::FileDebugWrite()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 DataCheck Write");
    uiint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile_Debug(msOutputFileName, nLevel+1);
    return 1;
}
void CMW::PrintRlt_Start(const uiint& nStep, bool bBinary)
{
    mpFileIO->PrintResult_Start(nStep, msRltFileName, bBinary);
}
//--
// width=出力行の幅{一行の文字数}:  format= %d:int32, %f:double(fixed), %e:double(scientific), %s:const char*
//--
void CMW::PrintRlt(const uiint& width, const char* format, ...)
{
    uiint nLength = strlen(format);
    
    vector<void*> param;
    vint    vnVal;
    vdouble vdVal, veVal;
    vstring vsVal;
    
    va_list list;
    va_start( list, format);
    for(uiint i=0; i < nLength; i++){
        if(format[i] == '%'){
            ++i;
            switch( format[i] ){
            case('d'):
                { iint nVal = va_arg(list, iint);  vnVal.push_back(nVal);}
                break;
            case('f'):
                { double dVal= va_arg(list, double); vdVal.push_back(dVal);}
                break;
            case('e'):
                { double dVal= va_arg(list, double); veVal.push_back(dVal);}
                break;
            case('s'):
                { string sVal= va_arg(list, const char*); vsVal.push_back(sVal);}
                break;
            default:
                break;
            }
        }
    };
    va_end( list );

    uiint icase_d(0), icase_f(0), icase_e(0), icase_s(0);
    for(uiint i=0; i < nLength; i++){
        if(format[i] == '%'){
            ++i;
            switch( format[i] ){
            case('d'):
                vnVal[icase_d];
                param.push_back(&vnVal[icase_d]);
                ++icase_d;
                break;
            case('f'):
                vdVal[icase_f];
                param.push_back(&vdVal[icase_f]);
                ++icase_f;
                break;
            case('e'):
                veVal[icase_e];
                param.push_back(&veVal[icase_e]);
                ++icase_e;
                break;
            case('s'):
                vsVal[icase_s];
                param.push_back(&vsVal[icase_s][0]);
                ++icase_s;
                break;
            }
        }
    };
    
    mpFileIO->PrintResult(width, format, param);
}
void CMW::PrintRlt_End()
{
    mpFileIO->PrintResult_End();
}
void CMW::PrintMicroAVS_Basis(const uiint& ieq)
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    CMesh *pMesh;
    stringstream ss;
    string sFileName;
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        pMesh= mpAssy->getMesh(iMesh);
        ss.clear(); ss.str("");
        ss << pMesh->getMeshID();
        sFileName = msRltFileName + "." + ss.str() + ".inp";
        mpFileIO->WriteAVS_Basis(sFileName, ieq, iMesh, nMaxLevel);
    };
}
void CMW::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    mpFileIO->recAVS_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
void CMW::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    mpFileIO->recAVS_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
void CMW::PrintMicroAVS_FEM()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    CMesh *pMesh;
    stringstream ss;
    string sFileName;
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        pMesh= mpAssy->getMesh(iMesh);
        ss.clear(); ss.str("");
        ss << pMesh->getMeshID();
        sFileName = msRltFileName + "." + ss.str() + ".inp";
        mpFileIO->WriteAVS_FEM(sFileName, iMesh, nMaxLevel);
    };
}

//--
// vtk
//--
void CMW::recVTK_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    mpFileIO->recVTK_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
void CMW::recVTK_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    mpFileIO->recVTK_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
void CMW::PrintVTK_FEM()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    stringstream ss;
    string sFileName;
    string sFileNameRank0;
    string sBaseName = mpFileIO->getMeshFileBaseName();
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        ss.clear(); ss.str("");
        ss << iMesh;
        sFileName = msRltFileName + "." + ss.str() + ".vtu";
        mpFileIO->WriteVTK_FEM(sFileName, iMesh, nMaxLevel);
        uiint rank= mpMPI->getRank();
        if (rank == 0) {
            sFileNameRank0 = sBaseName + "." + ss.str() + ".pvtu";
            uiint nNumOfProcs = mpMPI->getNumOfProcess();
            mpFileIO->WriteVTK_FEM_Rank0(sFileNameRank0, iMesh, nMaxLevel, sBaseName, nNumOfProcs);
        }
    };
}
//--
// field view
//--
void CMW::recUNS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    mpFileIO->recUNS_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
void CMW::recUNS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    mpFileIO->recUNS_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
void CMW::PrintUNS_FEM()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    stringstream ss;
    string sFileName;
    string sBaseName = mpFileIO->getMeshFileBaseName();

    sFileName = msRltFileName + ".fvuns";
    mpFileIO->WriteUNS_FEM(sFileName, nNumOfMesh, nMaxLevel);
}

uiint CMW::FileWriteRes(const uiint& nStep, bool bBinary)
{
    mpLogger->Info(Utility::LoggerMode::Info, " HEC_MW3 RestartFile Write");
    mpFileIO->WriteResFile(nStep, msResFileName, bBinary);
    return 1;
}
//
// 引数:nAppNumOfLevel, nAppNumEquation はデータチェック用途
//
uiint CMW::SetRestart(const uiint& nStep, uiint& nAppNumOfLevel, uiint& nAppNumEquation, bool bBinary)
{
    uiint nNumOfLevel;  //--
    uiint nNumOfAlgebra;// Restart-Data から取得
    uiint nNumOfMesh;   //--

    if(mpFileIO->ReadAlgebraBlock(nStep, msResFileName, bBinary) ){

        nNumOfLevel= mpFileIO->getNumOfLevel();     // MGレベル数
        nNumOfAlgebra= mpFileIO->getNumOfEquation();// Ax=bの数
        nNumOfMesh= mpFileIO->getNumOfParts();      // ローカルMesh数(PEごとのMesh数):
                                                    //  ## 制御ファイルはグローバルMesh数なので比較チェックできない.
        //
        // compare "main define param" & "restart data"
        //
        if(nNumOfLevel != nAppNumOfLevel){
            mpLogger->Info(Utility::LoggerMode::Error, " MW::SetRestart, Num of Level mismatch");
            return MW_ERROR;
        }
        if(nNumOfAlgebra != nAppNumEquation){
            mpLogger->Info(Utility::LoggerMode::Error, " MW::SetRestart, Num of AlgebraEquation mismatch");
            return MW_ERROR;
        }
    }

    //リスタート・データのセット
    mpFileIO->ReadResBlock(nStep, msResFileName, bBinary);
    
    return MW_SUCCESS;
}

//--
// 線形方程式の生成 : 方程式の行列、ベクトルの生成, Film:接合面の方程式別の伝達率
//--
void CMW::GeneLinearAlgebra(const uiint& nNumOfAlgebra, const uiint& nGlobalNumOfMesh, uiint** vvDOF, double** vvTransCoeff)
{
    //--
    // 熱伝達率:[AssyMatrix]が生成される前にFilm値をセット
    //--
    CAssyModel *pAssyModel= mpGMGModel->getAssyModel(0);//--- コースグリッド・AssyModel

    uiint nNumOfContact= pAssyModel->getNumOfContactMesh();
    for(uiint iCont=0; iCont < nNumOfContact; iCont++){
        CContactMesh *pContactMesh= pAssyModel->getContactMesh(iCont);
        CFilm *pFilm= pContactMesh->getFilm();

        pFilm->initTransCoeff(nNumOfAlgebra);

        for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++){
            pFilm->setTransCoeff(vvTransCoeff[iCont][ieq], ieq);//-- 接合面-方程式番号の順
        };
    };

    CMW::GeneLinearAlgebra(nNumOfAlgebra, nGlobalNumOfMesh, vvDOF);
}
//--
// 線形方程式生成 Base
//--
void CMW::GeneLinearAlgebra(const uiint& nNumOfAlgebra, const uiint& nGlobalNumOfMesh, uiint** vvDOF)
{
    vvuint vecDOF;
    vecDOF.resize(nNumOfAlgebra);

    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++){
        vecDOF[ieq].resize(nGlobalNumOfMesh);

        for(uiint imesh=0; imesh < nGlobalNumOfMesh; imesh++){
            vecDOF[ieq][imesh]= vvDOF[ieq][imesh];
        };
    };

    //--
    // 方程式:Ax=b 生成
    //--
    CAssyModel *pAssyModel, *pCGridAssyModel;
    uiint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
        pAssyModel = mpGMGModel->getAssyModel(iLevel);
        if(iLevel > 0){
            pCGridAssyModel = mpGMGModel->getAssyModel(iLevel-1);
            pAssyModel->GeneLinearAlgebra(vecDOF, pCGridAssyModel);//----Ax=b生成
        }else{
            pAssyModel->GeneLinearAlgebra(vecDOF, NULL);//---------------Ax=b生成
        }
    };
}

//--
// 線形方程式の選択
//--
void CMW::SelectAlgebra(const uiint& iequ)
{
    mpAssyMatrix = mpAssy->getAssyMatrix(iequ);
    mpRHSAssyVector = mpAssy->getRHSAssyVector(iequ);
    mpSolAssyVector = mpAssy->getSolutionAssyVector(iequ);
}

//--
// 行列・ベクトル(線形方程式)
//--
uiint CMW::Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CMW::Matrix_Add_Elem %d %e \n", iElem, ElemMatrix[0]);
#endif

    mpAssyMatrix->Matrix_Add_Elem(mpAssy, iMesh, iElem, ElemMatrix);

#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMW::Matrix_Add_Elem \n");
#endif
    return 1;
}
uiint CMW::Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double* NodalMatrix)
{
    return mpAssyMatrix->Matrix_Add_Nodal(iMesh, inode, jnode, NodalMatrix);
}
void CMW::Matrix_Clear(const uiint& iMesh)
{
    mpAssyMatrix->Matrix_Clear(iMesh);
}
void CMW::Vector_Clear(const uiint& iMesh)
{
    mpRHSAssyVector->Vector_Clear(iMesh);
    mpSolAssyVector->Vector_Clear(iMesh);
}
void CMW::AssyMatrix_Clear()
{
    mpAssyMatrix->Matrix_Clear();
}
void CMW::AssyVector_Clear()
{
    mpRHSAssyVector->Vector_Clear();
    mpSolAssyVector->Vector_Clear();
}
//--
// Solid:ノイマン境界,  Flow:
//--
uiint CMW::Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);

    if(mpRHSAssyVector){
        // Solid
        if(pMesh->getProp()==CodeType::Solid){

            ////cout << " isDirichlet:" << pMesh->isDirichletBNode(iNode) << " NID:" << iNodeID << " iNode:" << iNode << " Val:" << value << endl;

            // Dirichlet境界と重ならない
            if(!pMesh->isDirichletBNode(iNode)) 
                // 一旦、右辺ベクトルにセットした後で、Rank大の値をRank小に加算する=> Solveの直前で、mpRHSAssyVectorをSumupする.
                mpRHSAssyVector->setValue(iMesh, iNode, nDOF, value);
        }
        // Flow
        if(pMesh->getProp()==CodeType::Flow){
            // 通常Nodeと、Rank小の通信Nodeのみに荷重を与える => Rank大の通信Nodeには荷重を与えない
            if(!pMesh->isLargeRankCommNode(iNode))
                mpRHSAssyVector->setValue(iMesh, iNode, nDOF, value);
        }
        
        return 1;
    }else{
        return 0;
    }
}
//--
// Solid:ノイマン境界,  Flow:
//--
uiint CMW::Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);
    
    if(mpRHSAssyVector){
        // Solid
        if(pMesh->getProp()==CodeType::Solid){
            // Dirichlet境界と重ならない
            if(!pMesh->isDirichletBNode(iNode))   
                // 一旦、右辺ベクトルにセットした後で、Rank大の値をRank小に加算する=> Solveの直前で、mpRHSAssyVectorをSumupする.
                mpRHSAssyVector->addValue(iMesh, iNode, nDOF, value);
        }
        // Flow
        if(pMesh->getProp()==CodeType::Flow){
            // 通常Nodeと、Rank小の通信Nodeのみに荷重を与える => Rank大の通信Nodeには荷重を与えない
            if(!pMesh->isLargeRankCommNode(iNode))
                mpRHSAssyVector->addValue(iMesh, iNode, nDOF, value);
        }
            
        return 1;
    }else{
        return 0;
    }
}

//--
// 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.
//--
uiint CMW::Set_BC_NL_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);

    if(mpRHSAssyVector){
        // Rank小 && Dirichletでない
        if( !pMesh->isLargeRankCommNode(iNode) ){// && !pMesh->isDirichletBNode(iNode) ){
            mpRHSAssyVector->setValue(iMesh, iNode, nDOF, value);
        }
        return 1;
    }else{
        return 0;
    }
}
//--
// 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.
//--
uiint CMW::Add_BC_NL_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);

    if(mpRHSAssyVector){
        // Rank小 && Dirichletでない
        if( !pMesh->isLargeRankCommNode(iNode) ){//&& !pMesh->isDirichletBNode(iNode) ){
            mpRHSAssyVector->addValue(iMesh, iNode, nDOF, value);
        }
        return 1;
    }else{
        return 0;
    }
}
//
// ディレクレ処理のために、行列を"1"、"0"で払う
//
uiint CMW::Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diagValue, double& solValue)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);

    mpRHSAssyVector->setValue(iMesh, iNode, nDOF, solValue);

    mpAssyMatrix->setZero_NonDiag(iMesh, iNode, nDOF, mpRHSAssyVector, solValue);    
    mpAssyMatrix->setValue(iMesh, iNode, nDOF, diagValue, mpRHSAssyVector, solValue);

    return 1;
}
uiint CMW::Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diagValue, double& rhsValue)
{
#ifdef ADVANCESOFT_DEBUG
    printf("enter CMW::Set_BC (Mat_RHS) %d %d %d %e %e \n", iMesh, iNode, nDOF, diagValue, rhsValue);
#endif

    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    uiint iNode = pBucket->getIndexNode(iNodeID);
    
    mpAssyMatrix->setValue_D(iMesh, iNode, nDOF, diagValue);
    mpRHSAssyVector->setValue(iMesh, iNode, nDOF, rhsValue);

#ifdef ADVANCESOFT_DEBUG
    printf("exit CMW::Set_BC (Mat_RHS) \n");
#endif
	return 1;
}
//--
// 線形ソルバー
//--
uiint CMW::Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition)
{
#ifdef ADVANCESOFT_DEBUG
    printf(" enter CMW::Solve %d %e \n", iter_max, tolerance);
#endif

    //--
    // 通常Nodeと、Rank小の通信Nodeのみに荷重を与える => Rank大の通信Nodeには荷重を与えない
    // # 一旦、右辺ベクトルにセットした後、Rank大の値をRank小に加算する
    // # 選択されているmpRHSAssyVectorのみをsumup.(MGLevelと, AlgebraEquation は Select****で選択済み)
    //--
    mpRHSAssyVector->sumupCommBoundary();

  
    bool flag_iter_log = false;
    bool flag_time_log = false;
    uiint  success;
    //////--
    ////// MG + solver type
    //////--
    ////if( precondition==2 ){
    ////    switch( method ){
    ////        case(1): method= 5; break;
    ////        case(2): method= 6; break;
    ////        case(3): method= 7; break;
    ////        case(4):
    ////            mpLogger->Info(Utility::LoggerMode::Warn,"this combination not supported. change to precondition:1");
    ////            precondition=1;
    ////            break;
    ////        case(5): case(6): case(7):
    ////            mpLogger->Info(Utility::LoggerMode::Warn,"solver&pre combination mismatch. change to precondition:1");
    ////            precondition=1;
    ////            break;
    ////        default:
    ////            break;
    ////    }
    ////}
    //--
    // solver type switch
    //--
    switch( method ){
        case( 1 ):{
            CSolverCG *solver = new CSolverCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
            success = solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);
            }
            break;
        case( 2 ):{
            CSolverBiCGSTAB *solver = new CSolverBiCGSTAB(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
            success = solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector); 
            }
            break;
        case( 3 ):{ 
            CSolverGPBiCG *solver = new CSolverGPBiCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
            success = solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector); 
            }
            break;
        case( 4 ):{ 
            CSolverGMRES *solver = new CSolverGMRES(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
            success = solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);
            }
            break;
        case( 5 ):case( 6 ):case( 7 ):{ 
            // 5:CG-MG, 6:BiCGSTAB-MG, 7:GPBiCG-MG
            CSolverMG *solver = new CSolverMG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
            success = solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);
            }
            break;
        default:
           break;
    }
    if( success==0 ) {
      std::cout << " Fails in solver! " << std::endl;
	  return 0;
    }

#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMW::Solve \n");
#endif
    return 1;
}

void CMW::GetSolution_Vector(double* buf, const uiint& imesh)
{
    CVector *pVec = mpSolAssyVector->getVector(imesh);
    CMesh *pMesh = mpAssy->getMesh(imesh);
    uiint nNumOfNode = pMesh->getNumOfNode();
    CVector *pSolVec= mpSolAssyVector->getVector(imesh);
    uiint nDOF = pSolVec->getDOF();
    
    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            buf[ inode*nDOF + idof ] = pVec->getValue(inode, idof);
        };
    };
}
void CMW::GetSolution_AssyVector(double* buf)
{
    uiint nNumOfMesh = mpAssy->getNumOfMesh();
    
    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){
        CMesh *pMesh= mpAssy->getMesh(iMesh);
        CVector *pSolVec= mpSolAssyVector->getVector(iMesh);
        
        uiint nNumOfNode = pMesh->getNodeSize();
        uiint nDOF= pSolVec->getDOF();
        uiint nMeshSize= iMesh * nNumOfNode;
        
        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                buf[nMeshSize + iNode*nDOF + idof]= mpSolAssyVector->getValue(iMesh, iNode, idof);
            };
        };
    };
}
void CMW::GetRHS_Vector(double* buf, const uiint& imesh)
{
    CVector *pRHSVec= mpRHSAssyVector->getVector(imesh);
    CMesh *pMesh= mpAssy->getMesh(imesh);
    
    uiint nNumOfNode= pMesh->getNumOfNode();
    uiint nDOF= pRHSVec->getDOF();
    
    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            buf[ inode*nDOF + idof ]= pRHSVec->getValue(inode, idof);
        };
    };
}
void CMW::GetRHS_AssyVector(double* buf)
{
    uiint nNumOfMesh = mpAssy->getNumOfMesh();
    
    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){
        CMesh *pMesh= mpAssy->getMesh(iMesh);
        CVector* pRHSVec= mpRHSAssyVector->getVector(iMesh);
        
        uiint nNumOfNode= pMesh->getNodeSize();
        uiint nDOF= pRHSVec->getDOF();
        uiint nMeshSize= iMesh * nNumOfNode;
        
        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                buf[nMeshSize + iNode*nDOF + idof]= mpRHSAssyVector->getValue(iMesh, iNode, idof);
            };
        };
    };
}
//--
// 非線形構造解析の残差力:右辺ベクトルをsumupしたベクトル
//--
void CMW::GetRHS_Load(double* buf, const uiint& imesh)
{
    mpRHSAssyVector->sumupCommBoundary();

    GetRHS_Vector(buf, imesh);
}
//--
// 非線形構造解析の残差力:右辺ベクトルをsumupしたベクトル
//--
void CMW::GetRHS_AssyLoad(double* buf)
{
    mpRHSAssyVector->sumupCommBoundary();

    GetRHS_AssyVector(buf);
}
double& CMW::GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof)
{
    return mpSolAssyVector->getValue(imesh, inode, idof);
}
double& CMW::GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof)
{
    return mpRHSAssyVector->getValue(imesh, inode, idof);
}
uiint& CMW::GetSolutionVector_DOF(const uiint& imesh)
{
    CVector *pSolVec= mpSolAssyVector->getVector(imesh);
    return pSolVec->getDOF();
}
uiint& CMW::GetRHSVector_DOF(const uiint& imesh)
{
    CVector *pRHSVec= mpRHSAssyVector->getVector(imesh);
    return pRHSVec->getDOF();
}
uiint& CMW::GetSolutionVector_DOF(const uiint& iLevel, const uiint& imesh, const uiint& ieq)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CAssyVector *pAssyVec = pAssy->getSolutionAssyVector(ieq);
    CVector *pVec = pAssyVec->getVector(imesh);

    return pVec->getDOF();
}
uiint& CMW::GetRHSVector_DOF(const uiint& iLevel, const uiint& imesh, const uiint& ieq)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CAssyVector *pAssyVec = pAssy->getRHSAssyVector(ieq);
    CVector *pVec = pAssyVec->getVector(imesh);

    return pVec->getDOF();
}
//--
// debug用 行列ダンプ(CRT表示)
//--
void CMW::dump_AssyMatrix()
{
    mpAssyMatrix->dump();
}

void CMW::dump_RHSAssyVector()
{
    mpRHSAssyVector->dump();
}


//////--
////// 通信界面のベクトル値 加算
////// -------------------
////// # ユーザー定義ベクトルの通信界面値の加算.(MW3ソルバーには関与しない)
////// -------------------
//////--
////void CMW::Update_Add(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////    CMesh *pMesh = pAssy->getMesh(imesh);
////    CHecMPI *pMPI= CHecMPI::Instance();
////
////    CIndexBucket *pBucket = pMesh->getBucket();
////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////
////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////
////        uiint nNumOfCNode, transRank;
////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////        transRank = pCommMesh2->getTrasmitRank();
////
////        // データ交換
////        double *buff;
////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////            CNode *pNode = pCommNode->getNode();
////            uiint id = pNode->getID();
////            uiint index = pBucket->getIndexNode(id);
////
////            for(uiint idof=0; idof < nDOF; idof++)
////                buff[icnode*nDOF + idof] = vec[index*nDOF + idof];//----------- 交換用 buffに代入
////
////        };//icnode loop: 通信節点ループ
////
////        Send_Recv_R(buff, nNumOfCNode, nDOF, transRank);//-----------データ交換: send=buff, buff=recv
////
////        // データ加算
////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////            CNode *pNode = pCommNode->getNode();
////            uiint id = pNode->getID();
////            uiint index = pBucket->getIndexNode(id);
////
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[index*nDOF + idof] += buff[icnode*nDOF + idof];//--------------- vecに加算
////
////        };//icnode loop: 通信節点ループ
////
////        free( buff );
////
////    };//icom loop: 通信メッシュループ
////
////}

//--
// 通信界面のベクトル値 コピー : MatVecで利用
//--
// # 解ベクトルを想定している(ユーザー定義ベクトルの更新)
//
void CMW::Update(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CMesh *pMesh = pAssy->getMesh(imesh);

    CIndexBucket *pBucket = pMesh->getBucket();
    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();

    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);

        uiint nNumOfCNode, myRank, transRank;
        nNumOfCNode = pCommMesh2->getCommNodeSize();
        myRank = pCommMesh2->getRank();
        transRank = pCommMesh2->getTrasmitRank();

        if(mpMPI->getRank()!= myRank){ return;}//---------------------ERROR

        // バッファー確保: send=buff, buff=recv
        double *buff;
        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));

        // データ送信
        if(myRank < transRank){
            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
                CNode *pNode = pCommNode->getNode();
                uiint id = pNode->getID();
                uiint index = pBucket->getIndexNode(id);

                for(uiint idof=0; idof < nDOF; idof++)
                    buff[icnode*nDOF + idof] = vec[index*nDOF + idof];//----------- 送信用 buffに代入

            };//icnode loop: 通信節点ループ

            mpMPI->Send(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+myRank, MPI_COMM_WORLD);//-----------データ送信
        }
        // データ受信して代入
        if(myRank > transRank){

            MPI_Status stat;
            mpMPI->Recv(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+transRank, MPI_COMM_WORLD, &stat);//----データ受信

            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
                CNode *pNode = pCommNode->getNode();
                uiint id = pNode->getID();
                uiint index = pBucket->getIndexNode(id);

                for(uiint idof=0; idof < nDOF; idof++)
                    vec[index*nDOF + idof] = buff[icnode*nDOF + idof];//------------- 代入
            };//icnode loop: 通信節点ループ
        }
        //バッファー解放
        free( buff );

    };//icom loop: 通信メッシュループ
}
//--
// 通信界面のベクトル値の加算(rank大から小に値を加算): MatVecで利用
//--
// # 右辺ベクトルを想定している.(ユーザー定義ベクトルのSumup)
//
void CMW::Sumup(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CMesh *pMesh = pAssy->getMesh(imesh);

    CIndexBucket *pBucket = pMesh->getBucket();
    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();

    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);

        uiint nNumOfCNode, myRank, transRank;
        nNumOfCNode = pCommMesh2->getCommNodeSize();
        myRank = pCommMesh2->getRank();
        transRank = pCommMesh2->getTrasmitRank();

        if(mpMPI->getRank()!= myRank){ return;}//---------------------ERROR

        // バッファー確保: send=buff, buff=recv
        double *buff;
        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));

        // データ送信
        if(myRank > transRank){
            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
                CNode *pNode = pCommNode->getNode();
                uiint id = pNode->getID();
                uiint index = pBucket->getIndexNode(id);

                for(uiint idof=0; idof < nDOF; idof++){
                    buff[icnode*nDOF + idof] = vec[index*nDOF + idof];//----------- 送信用 buffに代入
                    vec[index*nDOF + idof]= 0.0;//--------------------------------- rank大の通信節点は:"0.0"に設定
                };

            };//icnode loop: 通信節点ループ

            mpMPI->Send(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+myRank, MPI_COMM_WORLD);//-----------データ送信
        }
        // データ受信して加算
        if(myRank < transRank){
            MPI_Status stat;
            mpMPI->Recv(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+transRank, MPI_COMM_WORLD, &stat);//----データ受信

            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
                CNode *pNode = pCommNode->getNode();
                uiint id = pNode->getID();
                uiint index = pBucket->getIndexNode(id);

                for(uiint idof=0; idof < nDOF; idof++)
                    vec[index*nDOF + idof] += buff[icnode*nDOF + idof];//---------------加算
            };//icnode loop: 通信節点ループ
        }
        //バッファー解放
        free( buff );

    };//icom loop: 通信メッシュループ
}
////
////// 通信界面のベクトル値 平均 : アセンブルモデル間の通信 => 未実装
//////
//////  # ユーザー定義ベクトルの通信界面値の平均
//////
////void CMW::Update_Average(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////    CMesh *pMesh = pAssy->getMesh(imesh);
////    CHecMPI *pMPI= CHecMPI::Instance();
////
////    CIndexBucket *pBucket = pMesh->getBucket();
////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////
////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////
////        uiint nNumOfCNode, transRank;
////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////        transRank = pCommMesh2->getTrasmitRank();
////
////        // データ交換
////        double *buff;
////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////            CNode *pNode = pCommNode->getNode();
////            uiint id = pNode->getID();
////            uiint index = pBucket->getIndexNode(id);
////
////            for(uiint idof=0; idof < nDOF; idof++)
////                buff[icnode*nDOF + idof] = vec[index*nDOF + idof];//----------- 交換用 buffに代入
////
////        };//icnode loop: 通信節点ループ
////
////        Send_Recv_R(buff, nNumOfCNode, nDOF, transRank);//-----------データ交換: send=buff, buff=recv
////
////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////            CNode *pNode = pCommNode->getNode();
////            uiint id = pNode->getID();
////            uiint index = pBucket->getIndexNode(id);
////
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[index*nDOF + idof] = 0.5*( vec[index*nDOF + idof] + buff[icnode*nDOF + idof] );//--- 平均:0.5*(buff+vec)に代入
////
////        };//icnode loop: 通信節点ループ
////
////        free( buff );
////
////    };//icom loop: 通信メッシュループ
////}
//////////
////////// 通信界面の
////////// 　右辺ベクトルに使用することを前提としている => 自身Rank小のベクトル値が更新される. '12.06.07
//////////
////////void CMW::Update_Add(const uiint& iLevel, const uiint& imesh, const uiint& ieq, const uiint& nDOF, uiint nVecType)
////////{
////////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////////    CMesh *pMesh = pAssy->getMesh(imesh);
////////    CHecMPI *pMPI= CHecMPI::Instance();
////////
////////    CAssyVector *pAssyVec;
////////    switch(nVecType){
////////        case(SOL_VEC): pAssyVec = pAssy->getSolutionAssyVector(ieq); break;
////////        case(RHS_VEC): pAssyVec = pAssy->getRHSAssyVector(ieq); break;
////////        default: mpLogger->Info(Utility::LoggerMode::Error, "CMW::Update_Average, fifth arg value is incorrect "); return;
////////    }
////////    CVector *pVec = pAssyVec->getVector(imesh);
////////
////////    CIndexBucket *pBucket = pMesh->getBucket();
////////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////////
////////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////////
////////        uiint nNumOfCNode, transRank;
////////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////////        transRank = pCommMesh2->getTrasmitRank();
////////
////////        // データ交換
////////        double *buff;
////////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////            CNode *pNode = pCommNode->getNode();
////////            uiint id = pNode->getID();
////////            uiint index = pBucket->getIndexNode(id);
////////
////////            for(uiint idof=0; idof < nDOF; idof++)
////////                buff[icnode*nDOF + idof] = pVec->getValue(index, idof);//----------- 交換用 buffに代入
////////
////////        };//icnode loop: 通信節点ループ
////////
////////        Send_Recv_R(buff, nNumOfCNode, nDOF, transRank);//-----------データ交換: send=buff, buff=recv
////////
////////        // データ加算: 自身Rankが小の場合に加算  # 右辺ベクトルに使用することを前提としている。
////////        //                                   # 通信節点のRank大の右辺ベクトル値は:0.0
////////        if(pMPI->getRank() < transRank ){
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++){
////////                    double value = pVec->getValue(index, idof) + buff[icnode*nDOF + idof];//--------------- vecに加算
////////                    pVec->setValue(index, idof, value);
////////                };
////////
////////            };//icnode loop: 通信節点ループ
////////        }//if(myRank小)
////////
////////        free( buff );
////////
////////    };//icom loop: 通信メッシュループ
////////}
//////////
////////// 通信節点の値を同値にする。
//////////
////////// 解ベクトルに使用することを想定している.
//////////
////////void CMW::Update(const uiint& iLevel, const uiint& imesh, const uiint& ieq, const uiint& nDOF, uiint nVecType)
////////{
////////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////////    CMesh *pMesh = pAssy->getMesh(imesh);
////////
////////    CAssyVector *pAssyVec;
////////    switch(nVecType){
////////        case(SOL_VEC): pAssyVec = pAssy->getSolutionAssyVector(ieq); break;
////////        case(RHS_VEC): pAssyVec = pAssy->getRHSAssyVector(ieq); break;
////////        default: mpLogger->Info(Utility::LoggerMode::Error, "CMW::Update_Average, fifth arg value is incorrect "); return;
////////    }
////////    CVector *pVec = pAssyVec->getVector(imesh);
////////
////////    CIndexBucket *pBucket = pMesh->getBucket();
////////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////////
////////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////////
////////        uiint nNumOfCNode, myRank, transRank;
////////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////////        myRank = pCommMesh2->getRank();
////////        transRank = pCommMesh2->getTrasmitRank();
////////
////////        if(mpMPI->getRank()!= myRank){ return;}//---------------------ERROR
////////
////////        // バッファー確保: send=buff, buff=recv
////////        double *buff;
////////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////////
////////        // データ送信
////////        if(myRank < transRank){
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++)
////////                    buff[icnode*nDOF + idof] = pVec->getValue(index, idof);//----------- buffにvecを代入
////////
////////            };//icnode loop: 通信節点ループ
////////
////////            mpMPI->Send(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+myRank, MPI_COMM_WORLD);//-----------データ送信
////////        }
////////        // データ受信して代入
////////        if(myRank > transRank){
////////
////////            MPI_Status stat;
////////            mpMPI->Recv(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+transRank, MPI_COMM_WORLD, &stat);//----データ受信
////////
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++){
////////                    double value = buff[icnode*nDOF + idof];//--------------- vecにbuffを代入
////////                    pVec->setValue(index, idof, value);
////////                }
////////            };//icnode loop: 通信節点ループ
////////        }
////////        // バッファー解放
////////        free( buff );
////////
////////    };//icom loop: 通信メッシュループ
////////}
//////////
////////// ranK大の値をrank小に加算
//////////  # 右辺ベクトルを想定している.
//////////
////////void CMW::Sumup(const uiint& iLevel, const uiint& imesh, const uiint& ieq, const uiint& nDOF, uiint nVecType)
////////{
////////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////////    CMesh *pMesh = pAssy->getMesh(imesh);
////////
////////    CAssyVector *pAssyVec;
////////    switch(nVecType){
////////        case(SOL_VEC): pAssyVec = pAssy->getSolutionAssyVector(ieq); break;
////////        case(RHS_VEC): pAssyVec = pAssy->getRHSAssyVector(ieq); break;
////////        default: mpLogger->Info(Utility::LoggerMode::Error, "CMW::Update_Average, fifth arg value is incorrect "); return;
////////    }
////////    CVector *pVec = pAssyVec->getVector(imesh);
////////
////////    CIndexBucket *pBucket = pMesh->getBucket();
////////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////////
////////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////////
////////        uiint nNumOfCNode, myRank, transRank;
////////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////////        myRank = pCommMesh2->getRank();
////////        transRank = pCommMesh2->getTrasmitRank();
////////
////////        if(mpMPI->getRank()!= myRank){ return;}//---------------------ERROR
////////
////////        // バッファー確保: send=buff, buff=recv
////////        double *buff;
////////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////////
////////        // データ送信
////////        if(myRank > transRank){
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++){
////////                    buff[icnode*nDOF + idof] = pVec->getValue(index, idof);//----------- buffにvecを代入
////////                    pVec->setValue(index, idof, 0.0);//--------------------------------- rank大:値をクリア
////////                }
////////
////////            };//icnode loop: 通信節点ループ
////////
////////            mpMPI->Send(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+myRank, MPI_COMM_WORLD);//-----------データ送信
////////        }
////////        // データ受信して加算
////////        if(myRank < transRank){
////////
////////            MPI_Status stat;
////////            mpMPI->Recv(buff, nNumOfCNode*nDOF, MPI_DOUBLE, transRank, 100+transRank, MPI_COMM_WORLD, &stat);//----データ受信
////////
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++){
////////                    double value= pVec->getValue(index, idof);
////////                    value += buff[icnode*nDOF + idof];//--------------- vecにbuffを加算
////////                    pVec->setValue(index, idof, value);
////////                }
////////            };//icnode loop: 通信節点ループ
////////        }
////////        // バッファー解放
////////        free( buff );
////////
////////    };//icom loop: 通信メッシュループ
////////}
//////////
////////// #右辺ベクトルを想定している.
//////////
////////void CMW::Update_Average(const uiint& iLevel, const uiint& imesh, const uiint& ieq, const uiint& nDOF, uiint nVecType)
////////{
////////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////////    CMesh *pMesh = pAssy->getMesh(imesh);
////////    CHecMPI *pMPI= CHecMPI::Instance();
////////
////////    CAssyVector *pAssyVec;
////////    switch(nVecType){
////////        case(SOL_VEC): pAssyVec = pAssy->getSolutionAssyVector(ieq); break;
////////        case(RHS_VEC): pAssyVec = pAssy->getRHSAssyVector(ieq); break;
////////        default: mpLogger->Info(Utility::LoggerMode::Error, "CMW::Update_Average, fifth arg value is incorrect "); return;
////////    }
////////    CVector *pVec = pAssyVec->getVector(imesh);
////////
////////    CIndexBucket *pBucket = pMesh->getBucket();
////////    uiint nNumOfCommMesh2 = pMesh->getCommMesh2Size();
////////
////////    for(uiint icom=0; icom < nNumOfCommMesh2; icom++){
////////        CCommMesh2 *pCommMesh2 = pMesh->getCommMesh2IX(icom);
////////
////////        uiint nNumOfCNode, transRank;
////////        nNumOfCNode = pCommMesh2->getCommNodeSize();
////////        transRank = pCommMesh2->getTrasmitRank();
////////
////////        // データ交換
////////        double *buff;
////////        buff = (double*)calloc(nNumOfCNode*nDOF, sizeof(double));
////////        for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////            CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////            CNode *pNode = pCommNode->getNode();
////////            uiint id = pNode->getID();
////////            uiint index = pBucket->getIndexNode(id);
////////
////////            for(uiint idof=0; idof < nDOF; idof++)
////////                buff[icnode*nDOF + idof] = pVec->getValue(index, idof);//----------- buffに代入
////////
////////        };//icnode loop: 通信節点ループ
////////
////////        Send_Recv_R(buff, nNumOfCNode, nDOF, transRank);//-----------データ交換: send=buff, buff=recv
////////
////////        // データ代入(相手との平均)
////////        if(pMPI->getRank() < transRank ){
////////            for(uiint icnode=0; icnode < nNumOfCNode; icnode++){
////////                CCommNode *pCommNode = pCommMesh2->getCommNodeIX(icnode);
////////                CNode *pNode = pCommNode->getNode();
////////                uiint id = pNode->getID();
////////                uiint index = pBucket->getIndexNode(id);
////////
////////                for(uiint idof=0; idof < nDOF; idof++){
////////                    double value = 0.5*( pVec->getValue(index, idof) + buff[icnode*nDOF + idof] );//--- 平均:0.5*(buff+vec)に代入
////////                    pVec->setValue(index, idof, value);
////////                }
////////
////////            };//icnode loop: 通信節点ループ
////////        }//if(myRank小)
////////
////////        free( buff );
////////
////////    };//icom loop: 通信メッシュループ
////////}
////
//////--
////// ユーザー定義アセンブルベクトルのassy_vec のupdate_add(互いに加算)
//////--
////void CMW::Update_Assy_Add(const uiint& iLevel, double* assy_vec, uiint* vDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////
////    for(uiint imesh=0; imesh < pAssy->getNumOfMesh(); imesh++){
////        uiint nDOF = vDOF[imesh];
////        CMesh *pMesh = pAssy->getMesh(imesh);
////        uiint nNumOfNode = pMesh->getNumOfNode();
////
////        //パーツ別ベクトルのコピー
////        double *vec = (double*)malloc(nNumOfNode*nDOF*sizeof(double));
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[inode*nDOF + idof] = assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof];
////
////        Update_Add(iLevel, imesh, vec, nDOF);//パーツ別ベクトルのupdate_add
////
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof] = vec[inode*nDOF + idof];
////
////        free( vec );
////    };//imesh loop
////}
//////--
////// ユーザー定義アセンブルベクトルのassy_vec のupdate(rank小から大へ値を更新)
//////--
////void CMW::Update_Assy(const uiint& iLevel, double* assy_vec, uiint* vDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////
////    for(uiint imesh=0; imesh < pAssy->getNumOfMesh(); imesh++){
////        uiint nDOF = vDOF[imesh];
////        CMesh *pMesh = pAssy->getMesh(imesh);
////        uiint nNumOfNode = pMesh->getNumOfNode();
////
////        //パーツ別ベクトルのコピー
////        double *vec = (double*)malloc(nNumOfNode*nDOF*sizeof(double));
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[inode*nDOF + idof] = assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof];
////
////        Update(iLevel, imesh, vec, nDOF);//パーツ別ベクトルのupdate
////
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof] = vec[inode*nDOF + idof];
////
////        free( vec );
////    };//imesh loop
////}
//////--
////// ユーザー定義アセンブルベクトルのassy_vec のsumup(rank大から小へ加算)
//////--
////void CMW::Sumup_Assy(const uiint& iLevel, double* assy_vec, uiint* vDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////
////    for(uiint imesh=0; imesh < pAssy->getNumOfMesh(); imesh++){
////        uiint nDOF = vDOF[imesh];
////        CMesh *pMesh = pAssy->getMesh(imesh);
////        uiint nNumOfNode = pMesh->getNumOfNode();
////
////        //パーツ別ベクトルのコピー
////        double *vec = (double*)malloc(nNumOfNode*nDOF*sizeof(double));
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[inode*nDOF + idof] = assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof];
////
////        Sumup(iLevel, imesh, vec, nDOF);//パーツ別ベクトルのsumup
////
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof] = vec[inode*nDOF + idof];
////
////        free( vec );
////    };//imesh loop
////}
//////--
////// ユーザー定義アセンブルベクトルのassy_vec の平均
//////--
////void CMW::Update_Assy_Average(const uiint& iLevel, double* assy_vec, uiint* vDOF)
////{
////    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
////
////    for(uiint imesh=0; imesh < pAssy->getNumOfMesh(); imesh++){
////        uiint nDOF = vDOF[imesh];
////        CMesh *pMesh = pAssy->getMesh(imesh);
////        uiint nNumOfNode = pMesh->getNumOfNode();
////
////        //パーツ別ベクトルのコピー
////        double *vec = (double*)malloc(nNumOfNode*nDOF*sizeof(double));
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                vec[inode*nDOF + idof] = assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof];
////
////        Update_Average(iLevel, imesh, vec, nDOF);//パーツ別ベクトルのupdate_average
////
////        for(uiint inode=0; inode < nNumOfNode; inode++)
////            for(uiint idof=0; idof < nDOF; idof++)
////                assy_vec[imesh*nNumOfNode*nDOF + inode*nDOF + idof] = vec[inode*nDOF + idof];
////
////        free( vec );
////    };//imesh loop
////}


//--
// AssyMatrix版 Ax=y
//--
void CMW::MatVec_Assy(const uiint& iLevel, const uiint& ieq, double* assy_x, double* assy_y)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CAssyMatrix *pAssyMat = pAssy->getAssyMatrix(ieq);//方程式番号:ieq

    uiint nNumOfMesh = pAssy->getNumOfMesh();

    uiint *vDOF = (uiint*)malloc(nNumOfMesh*sizeof(uiint));
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++) vDOF[imesh]= pAssyMat->getMatrix(imesh)->getDOF();


    // ベクトル初期化
    CAssyVector *pX = pAssy->getSolutionAssyVector(ieq);//初期化に利用
    CAssyVector X(pX), Y(pX);

    for(uiint imesh=0; imesh < nNumOfMesh; imesh++){ X.Vector_Clear(imesh);  Y.Vector_Clear(imesh);}

    // Xベクトルにassy_xベクトルをコピー
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
        CMesh *pMesh = pAssy->getMesh(imesh);
        uiint nNumOfNode = pMesh->getNumOfNode();
        uiint nDOF = vDOF[imesh];
        for(uiint inode=0; inode < nNumOfNode; inode++){
            for(uiint idof=0; idof < nDOF; idof++){
                double val = assy_x[imesh*nNumOfNode*nDOF + inode*nDOF + idof];
                X.addValue(imesh, inode, idof, val);//assy_x から Xベクトルへ値をセット
            };//idof
        };//inode
    };//imesh

    pAssyMat->multVector(&X, &Y, NULL);// Y = AssyMat * X  : multVector内で通信(update, sumup)

    //Y => assy_y にコピー
    for(uiint imesh=0; imesh < nNumOfMesh; imesh++){
        CMesh *pMesh = pAssy->getMesh(imesh);
        uiint nNumOfNode = pMesh->getNumOfNode();
        uiint nDOF = vDOF[imesh];
        for(uiint inode=0; inode < nNumOfNode; inode++){
            for(uiint idof=0; idof < nDOF; idof++){
                double val = Y.getValue(imesh, inode, idof);
                assy_y[imesh*nNumOfNode*nDOF + inode*nDOF + idof]= val;
            };//idof
        };//inode
    };//imesh

    //// Update_Assy(iLevel, assy_y, vDOF);// yベクトルの更新(並列処理には適さない)

    free( vDOF );
}
//--
// MatrixBCRS版 Ax=y : xベクトルのupdate後、行列ベクトル演算、演算後にyベクトルをsumup
//--
void CMW::MatVec(const uiint& iLevel, const uiint& imesh, const uiint& ieq, double* x, double* y)
{
    CAssyModel *pAssy = mpGMGModel->getAssyModel(iLevel);
    CAssyMatrix *pAssyMat = pAssy->getAssyMatrix(ieq);//方程式番号:ieq
    CMatrixBCRS *pMat = pAssyMat->getMatrix(imesh);

    uiint nDOF = pMat->getDOF();

    Update(iLevel, imesh, x, nDOF);// xベクトル_update(通信1/2)

    //xベクトル から Vector Xに値をコピー
    CMesh *pMesh = pAssy->getMesh(imesh);
    CVector X(pMesh, nDOF), Y(pMesh, nDOF);
    uiint nNumOfNode = pMesh->getNumOfNode();
    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            double val = x[inode*nDOF + idof];
            X.setValue(inode, idof, val);
        };//idof
    };//inode

    pMat->multVector(&X, &Y);// Y = Mat * X  : パーツ別のmultVector

    //Vector Y から yベクトルに値をコピー
    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            double val = Y.getValue(inode, idof);
            y[inode*nDOF + idof] = val;
        };//idof
    };//inode

    Sumup(iLevel, imesh, y, nDOF); // yベクトル_sumup(通信2/2)

    //// Update(iLevel, imesh, y, nDOF);// yベクトルの更新(並列処理には適さない)
}

//--
// MW3 Refine
//--
uiint CMW::Refine(const uiint& nNumOfRefine)
{
    if(mbRefine_Use){
        mpLogger->Info(Utility::LoggerMode::Info, "You have already used the Refine, Refine_Function can not be used.");
        return MW_ERROR;
    }else{
        mbRefine_Use = true;
    }

    CHecMPI* pMPI=CHecMPI::Instance();

    uiint ilevel,mgLevel;
    CAssyModel *pAssy;
    uiint icon,numOfConMesh;
    CContactMesh *pConMesh;

    mpFactory->GeneAssyModel(nNumOfRefine+1);//Level=0 以外のAssyModelを生成
    mpFactory->setMGLevel(nNumOfRefine);
    mpFactory->setupNodeGridSize();//Level:0のNodeにMG階層数

    pAssy = mpGMGModel->getAssyModel(0);
    uiint iMesh, nNumOfMesh= pAssy->getNumOfMesh();

    // BoundaryMesh前処理
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        CMesh *pMesh= pAssy->getMesh(iMesh);
        uiint iBMesh, nNumOfBMesh= pMesh->getNumOfBoundaryFaceMesh();
        for(iBMesh=0; iBMesh < nNumOfBMesh; iBMesh++){
            CBoundaryFaceMesh *pBFMesh= pMesh->getBndFaceMeshIX(iBMesh);
            pBFMesh->setMaxMGLevel(nNumOfRefine);
            pBFMesh->resizeCGrid_BNodeValue(nNumOfRefine);
        };
        nNumOfBMesh= pMesh->getNumOfBoundaryEdgeMesh();
        for(iBMesh=0; iBMesh < nNumOfBMesh; iBMesh++){
            CBoundaryEdgeMesh *pBEMesh= pMesh->getBndEdgeMeshIX(iBMesh);
            pBEMesh->setMaxMGLevel(nNumOfRefine);
            pBEMesh->resizeCGrid_BNodeValue(nNumOfRefine);
        };
        nNumOfBMesh= pMesh->getNumOfBoundaryVolumeMesh();
        for(iBMesh=0; iBMesh < nNumOfBMesh; iBMesh++){
            CBoundaryVolumeMesh *pBVMesh= pMesh->getBndVolumeMeshIX(iBMesh);
            pBVMesh->setMaxMGLevel(nNumOfRefine);
            pBVMesh->resizeCGrid_BNodeValue(nNumOfRefine);
        }
    };
    //--
    // ContactMesh 準備(生成)
    //--
    uiint contactID, myRank, nNumTransRank, nProp, nNumOfAlgebra;
    vdouble vTransCoeff;//------------------------ 方程式別の伝達率
    CAssyModel *pPrevAssy= mpGMGModel->getAssyModel(0);
    uiint nNumOfCont= pPrevAssy->getNumOfContactMesh();
    //--
    // Level:1〜Max ,カレントレベル==ファイングリッド
    //--
    for(ilevel=1; ilevel< nNumOfRefine+1; ilevel++){
        pAssy= mpGMGModel->getAssyModel(ilevel);
        for(uiint icont=0; icont < nNumOfCont; icont++){
            CContactMesh *pPrevContactMesh= pPrevAssy->getContactMesh(icont);
            CContactMesh *pContactMesh = new CContactMesh;//--------------------ContactMesh生成

            contactID = pPrevContactMesh->getID();
            pContactMesh->setID(contactID);
            pContactMesh->setLevel(ilevel);
            myRank= pPrevContactMesh->getRank();
            nProp= pPrevContactMesh->getProp();

            CFilm *pFilm= pPrevContactMesh->getFilm();//--- コースグリッドのFilm

            pContactMesh->setRank(myRank);
            pContactMesh->setProp(nProp);
            pContactMesh->setFilm(pFilm);//--- 熱伝達率 : 各階層全て同じFilm
            
            pAssy->addContactMesh(pContactMesh, contactID);
        };
    };
    //--
    // Refine : Mesh, ContactMesh, CommMesh2, BoundaryMesh
    //--
    if(mb_file){

        ////cout << "MW::Refine ---------- A rank:" << pMPI->getRank() << endl;

        mpFactory->refineMesh();

        ////cout << "MW::Refine ---------- B rank:" << pMPI->getRank() << endl;

        mpFactory->refineContactMesh();
        
        ////cout << "MW::Refine ---------- C rank:" << pMPI->getRank() << endl;

        mgLevel= mpFactory->getMGLevel();

        for(ilevel=0; ilevel < mgLevel+1; ilevel++){
            pAssy= mpGMGModel->getAssyModel(ilevel);
            numOfConMesh= pAssy->getNumOfContactMesh();
            for(icon=0; icon < numOfConMesh; icon++){
                pConMesh= pAssy->getContactMesh(icon);
                pConMesh->setupSPointOnMFace();
                pConMesh->setupMPC_Coef();     
            };
        };

        ////cout << "MW::Refine ---------- D rank:" << pMPI->getRank() << endl;

        mpFactory->refineCommMesh2();

        ////cout << "MW::Refine ---------- E rank:" << pMPI->getRank() << endl;

        mpFactory->refineBoundary(); 

        ////cout << "MW::Refine ---------- F rank:" << pMPI->getRank() << endl;

        mpFactory->setupBNodeMarking();//--------------- 境界条件マーキング：MG境界処理
        mpFactory->setupLargeRankCommNode_Marking();//-- Rank大の通信ノードマーキング:並列計算での荷重処理
        return 1;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error," Not read file, MW::Refine()");
        return 0;
    }
    
}
void CMW::FinalizeRefine()
{
    CAssyModel *pAssy;
    uiint ilevel, mgLevel = mpFactory->getMGLevel();
    for(ilevel=0; ilevel < mgLevel+1; ilevel++){
        pAssy = mpGMGModel->getAssyModel(ilevel);
        CMesh *pMesh;
        uiint imesh, nNumOfMesh = pAssy->getNumOfMesh();
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            pMesh = pAssy->getMesh(imesh);
            pMesh->deleteProgData();
            pMesh->deleteAggregate_on_Node();
            CCommMesh2 *pCommMesh;
            uiint icom, nNumOfCommMesh = pMesh->getCommMesh2Size();
            for(icom=0; icom < nNumOfCommMesh; icom++){
                pCommMesh = pMesh->getCommMesh2IX(icom);
                pCommMesh->deleteProgData();
            };
            CBoundaryVolumeMesh *pBVolMesh;
            uiint ibvol, nNumOfBVolMesh = pMesh->getNumOfBoundaryVolumeMesh();
            for(ibvol=0; ibvol < nNumOfBVolMesh; ibvol++){
                pBVolMesh = pMesh->getBndVolumeMeshIX(ibvol);
                pBVolMesh->deleteProgData();
            };
            CBoundaryFaceMesh *pBFaceMesh;
            uiint ibface, nNumOfBFaceMesh = pMesh->getNumOfBoundaryFaceMesh();
            for(ibface=0; ibface < nNumOfBFaceMesh; ibface++){
                pBFaceMesh = pMesh->getBndFaceMeshIX(ibface);
                pBFaceMesh->deleteProgData();
            };
            CBoundaryEdgeMesh *pBEdgeMesh;
            uiint ibedge, nNumOfBEdgeMesh = pMesh->getNumOfBoundaryEdgeMesh();
            for(ibedge=0; ibedge < nNumOfBEdgeMesh; ibedge++){
                pBEdgeMesh = pMesh->getBndEdgeMeshIX(ibedge);
                pBEdgeMesh->deleteProgData();
            };
        };
        CContactMesh *pConMesh;
        uiint icmesh, nNumOfConMesh = pAssy->getNumOfContactMesh();
        for(icmesh=0; icmesh < nNumOfConMesh; icmesh++){
            pConMesh = pAssy->getContactMesh(icmesh);
            pConMesh->deleteProgData();
        };
    };
}
uiint CMW::GetNumOfAssembleModel()
{
    return mpGMGModel->getNumOfLevel();
}
void CMW::SelectAssembleModel(const uiint& mgLevel)
{
    mpAssy= mpGMGModel->getAssyModel(mgLevel);
}
uiint CMW::GetNumOfMeshPart()
{
    return mpAssy->getNumOfMesh();
}
void CMW::SelectMeshPart_ID(const uiint& mesh_id)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh_ID(mesh_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
void CMW::SelectMeshPart_IX(const uiint& index)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
uiint CMW::GetMeshID_Num(const uiint& index)//---------- Mesh 通し番号(Index)から、ID
{
    return mpAssy->getID_of_Mesh(index);
}
uiint CMW::GetMeshIndex_Num(const uiint& id)//---------- Mesh ID番号から、Index
{
    return mpAssy->getIndex_of_Mesh(id);
}

void CMW::SelectElement_ID(const uiint& elem_id)
{
    if(mpMesh){
        mpElement= mpMesh->getElement(elem_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
void CMW::SelectElement_IX(const uiint& index)
{
    if(mpMesh){
        mpElement= mpMesh->getElementIX(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
uiint CMW::GetElementType()
{
    return mpElement->getType();
}
uiint CMW::GetNumOfElementVert()
{
    return mpElement->getNumOfNode();
}
void CMW::GetElementVertNodeID(iint* vNodeID)
{
    uiint nNumOfNode;
    nNumOfNode= mpElement->getNumOfNode();
    CNode* pNode;
    uiint ivert;
    for(ivert=0; ivert< nNumOfNode; ivert++){
        pNode= mpElement->getNode(ivert);
        vNodeID[ivert]= pNode->getID();
    };
}
void CMW::GetElementVertNodeIndex(iint* vNodeIndex)
{
    CIndexBucket *pBucket = mpMesh->getBucket();

    uiint nNumOfNode;
    nNumOfNode= mpElement->getNumOfNode();
    CNode* pNode;
    uiint ivert;
    for(ivert=0; ivert< nNumOfNode; ivert++){
        pNode= mpElement->getNode(ivert);
        uiint id = pNode->getID();

        vNodeIndex[ivert]= pBucket->getIndexNode(id);
    };
}
uiint CMW::GetNumOfElementEdge()
{
    return mpElement->getNumOfEdge();
}
uiint CMW::GetNumOfElementFace()
{
    return mpElement->getNumOfFace();
}
void CMW::GetElementEdgeNodeID(iint* vNodeID)
{
    uiint numOfEdge;
    numOfEdge= mpElement->getNumOfEdge();
    CNode *pNode;
    uiint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++){
        pNode= mpElement->getEdgeInterNode(iedge);
        vNodeID[iedge]= pNode->getID();
    };
}
uiint CMW::GetElementFaceElementID(uiint faceIndex)
{
    CElement* elem = mpElement->getFaceElement(faceIndex);
    if(!elem){
        return NULL;
    }else{
        uiint id = elem->getID();
        return id;
    }
}
uiint CMW::GetNumOfElementEdgeElement(uiint edgeIndex)
{
    vector<CElement*> elements = mpElement->getEdgeElement(edgeIndex);
    return elements.size();

}
uiint CMW::GetElementEdgeElementID(uiint edgeIndex, uiint i)
{
    vector<CElement*> elements = mpElement->getEdgeElement(edgeIndex);
    if(i >= elements.size()){
        return NULL;
    }
    return elements[i]->getID();

}
////////////////////////////////////////
uiint CMW::getNodeSize()
{
    return mpMesh->getNodeSize();
}
uiint CMW::getElementSize()
{
    return mpMesh->getElementSize();
}
uiint CMW::getNodeSize(uiint iMesh)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    return pMesh->getNodeSize();
}
uiint CMW::getElementSize(uiint iMesh)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    return pMesh->getElementSize();
}
/////////////////////////////////////////
uiint& CMW::getNodeID(const uiint& index)
{
    CNode* pNode = mpMesh->getNodeIX(index);
    return pNode->getID();
}
uiint& CMW::getElementID(const uiint& index)
{
    CElement* pElem = mpMesh->getElementIX(index);
    return pElem->getID();
}
uiint& CMW::getNodeIndex(const uiint& id)
{
    CIndexBucket *pBucket = mpMesh->getBucket();
    return pBucket->getIndexNode(id);
}
uiint& CMW::getElementIndex(const uiint& id)
{
    CIndexBucket *pBucket = mpMesh->getBucket();
    return pBucket->getIndexElement(id);
}

uiint CMW::getNumOfParentNode(const uiint& id, const uiint& nLevel)
{
    CNode *pNode= mpMesh->getNode(id);

    uiint nNumP = pNode->getNumOfParentNode(nLevel);

    return nNumP;
}
uiint CMW::getParentNodeID(const uiint& id, const uiint& nLevel, const uiint& index)
{
    CNode *pNode= mpMesh->getNode(id);
    CNode *pPareNode= pNode->getParentNode(nLevel, index);

    return pPareNode->getID();
}

void CMW::constructNodeConnectFEM(const uiint& node_id)
{
    mv_ItemL.clear();
    mv_ItemU.clear();
    CIndexBucket *pBucket= mpMesh->getBucket();
    CElement *pElement;
    CAggregateElement *pAggElement= mpMesh->getAggElem(node_id);
    uiint nNumOfElement= pAggElement->getNumOfElement();
    uiint i_elem;
    for(i_elem=0; i_elem < nNumOfElement; i_elem++){
        pElement= pAggElement->get(i_elem);
        uiint k_node, nNumOfElemNode= pElement->getNumOfNode();
        for (k_node = 0; k_node < nNumOfElemNode; k_node++) {
            uiint k_node_id = pElement->getNode(k_node)->getID();
            uiint k_index = pBucket->getIndexNode(k_node_id);
            if(k_node_id < node_id){
                mv_ItemL.push_back(k_index);
            }else if(node_id < k_node_id){
                mv_ItemU.push_back(k_index);
            }
        };
    };
    sort(mv_ItemL.begin(), mv_ItemL.end());
    vector<uiint>::iterator new_end = unique(mv_ItemL.begin(), mv_ItemL.end());
    mv_ItemL.erase(new_end, mv_ItemL.end());
    sort(mv_ItemU.begin(), mv_ItemU.end());
    new_end = unique(mv_ItemU.begin(), mv_ItemU.end());
    mv_ItemU.erase(new_end, mv_ItemU.end());
}
void CMW::getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL)
{
    nNumOfItemU = mv_ItemU.size();
    nNumOfItemL = mv_ItemL.size();
}
void CMW::getNodeConnectFEM_Item(uiint itemU[], uiint itemL[])
{
    for(uiint i=0; i < mv_ItemL.size(); i++){
        itemL[i] = mv_ItemL[i];
    };
    for(uiint i=0; i < mv_ItemU.size(); i++){
        itemU[i] = mv_ItemU[i];
    };
}
void CMW::getNodeConnectFEM_Item_F(iint itemU[], iint itemL[])// Fortran
{
    uiint nMax;
    if(mv_ItemL.size() > IINT_MAX){
        nMax=IINT_MAX;
        mpLogger->Info(Utility::LoggerMode::Error, " getNodeConnectFEM_Item_F, over INT_MAX");
    }else{
        nMax=mv_ItemL.size();
    }
    for(uiint i=0; i < nMax; i++){
        itemL[i] = mv_ItemL[i];
    };
    if(mv_ItemL.size() > IINT_MAX){
        nMax=IINT_MAX;
        mpLogger->Info(Utility::LoggerMode::Error, " getNodeConnectFEM_Item_F, over INT_MAX");
    }else{
        nMax=mv_ItemU.size();
    }
    for(uiint i=0; i < nMax; i++){
        itemU[i] = mv_ItemU[i];
    };
}
uiint CMW::getNumOfAggregateElement(const uiint& node_id)
{
    CAggregateElement *pAggElement= mpMesh->getAggElem(node_id);
    return pAggElement->getNumOfElement();
}
uiint& CMW::getAggregateElementID(const uiint& node_id, const uiint& ielem)
{
    CAggregateElement *pAggElement= mpMesh->getAggElem(node_id);
    CElement *pElem = pAggElement->get(ielem);
    return pElem->getID();
}
void CMW::setupNeighbors()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    CMesh *pMesh;
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        pMesh= mpAssy->getMesh(iMesh);
        pMesh->setupAggregate(0);
        pMesh->setupEdgeElement(NULL, 0);
        pMesh->setupFaceElement2(NULL);
    };
}
void CMW::GetNodeCoord(const uiint& node_id, double& x, double& y, double& z)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    x= pNode->getX();
    y= pNode->getY();
    z= pNode->getZ();
}
uiint CMW::GetNumOfDOF(const uiint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    return pNode->getTotalDOF();
}
uiint CMW::GetNumOfScalar(const uiint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    return pNode->getScalarDOF();
}
uiint CMW::GetNumOfVector(const uiint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    return pNode->getVectorDOF();
}
uiint& CMW::GetNodeType(const uiint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    return pNode->getType();
}
uiint& CMW::NumOfIntegPoint(const uiint& shapeType)
{
    return mpShapeCatalog->NumOfIntegPoint(shapeType);
}
void CMW::ShapeFunc_on_pt(const uiint& shapeType, const uiint& igauss, vdouble& N)
{
    switch(shapeType){
        case(ShapeType::Hexa81):
            N= mpShapeHexa->N81(igauss);
            break;
        case(ShapeType::Hexa82):
            N= mpShapeHexa->N82(igauss);
            break;
        case(ShapeType::Hexa201):
            N= mpShapeHexa->N201(igauss);
            break;
        case(ShapeType::Hexa202):
            N= mpShapeHexa->N202(igauss);
            break;
        case(ShapeType::Hexa203):
            N= mpShapeHexa->N203(igauss);
            break;
        case(ShapeType::HexaNic111):
            N= mpShapeHexaNic->N111(igauss);
            break;
        case(ShapeType::HexaNic118):
            N= mpShapeHexaNic->N118(igauss);
            break;
        case(ShapeType::HexaNic1127):
            N= mpShapeHexaNic->N1127(igauss);
            break;
        case(ShapeType::Tetra41):
            N= mpShapeTetra->N41(igauss);
            break;
        case(ShapeType::Tetra101):
            N= mpShapeTetra->N101(igauss);
            break;
        case(ShapeType::Tetra104):
            N= mpShapeTetra->N104(igauss);
            break;
        case(ShapeType::Tetra1015):
            N= mpShapeTetra->N1015(igauss);
            break;
        case(ShapeType::Prism62):
            N= mpShapePrism->N62(igauss);
            break;
        case(ShapeType::Prism156):
            N= mpShapePrism->N156(igauss);
            break;
        case(ShapeType::Prism159):
            N= mpShapePrism->N159(igauss);
            break;
        case(ShapeType::Prism1518):
            N= mpShapePrism->N1518(igauss);
            break;
        case(ShapeType::Quad41):
            N= mpShapeQuad->N41(igauss);
            break;
        case(ShapeType::Quad84):
            N= mpShapeQuad->N84(igauss);
            break;
        case(ShapeType::Quad89):
            N= mpShapeQuad->N89(igauss);
            break;
        case(ShapeType::Triangle31):
            N= mpShapeTriangle->N31(igauss);
            break;
        case(ShapeType::Triangle63):
            N= mpShapeTriangle->N63(igauss);
            break;
        case(ShapeType::Line21):
            N= mpShapeLine->N21(igauss);
            break;
        case(ShapeType::Line32):
            N= mpShapeLine->N32(igauss);
            break;
        default:
            break;
    }
}
void CMW::ShapeFunc_on_pt(uiint shapeType, uiint igauss, double N[])
{
    vdouble vN;
    switch(shapeType){
        case(ShapeType::Hexa81):
            vN= mpShapeHexa->N81(igauss);
            break;
        case(ShapeType::Hexa82):
            vN= mpShapeHexa->N82(igauss);
            break;
        case(ShapeType::Hexa201):
            vN= mpShapeHexa->N201(igauss);
            break;
        case(ShapeType::Hexa202):
            vN= mpShapeHexa->N202(igauss);
            break;
        case(ShapeType::Hexa203):
            vN= mpShapeHexa->N203(igauss);
            break;
        case(ShapeType::HexaNic111):
            vN= mpShapeHexaNic->N111(igauss);
            break;
        case(ShapeType::HexaNic118):
            vN= mpShapeHexaNic->N118(igauss);
            break;
        case(ShapeType::HexaNic1127):
            vN= mpShapeHexaNic->N1127(igauss);
            break;
        case(ShapeType::Tetra41):
            vN= mpShapeTetra->N41(igauss);
            break;
        case(ShapeType::Tetra101):
            vN= mpShapeTetra->N101(igauss);
            break;
        case(ShapeType::Tetra104):
            vN= mpShapeTetra->N104(igauss);
            break;
        case(ShapeType::Tetra1015):
            vN= mpShapeTetra->N1015(igauss);
            break;
        case(ShapeType::Prism62):
            vN= mpShapePrism->N62(igauss);
            break;
        case(ShapeType::Prism156):
            vN= mpShapePrism->N156(igauss);
            break;
        case(ShapeType::Prism159):
            vN= mpShapePrism->N159(igauss);
            break;
        case(ShapeType::Prism1518):
            vN= mpShapePrism->N1518(igauss);
            break;
        case(ShapeType::Quad41):
            vN= mpShapeQuad->N41(igauss);
            break;
        case(ShapeType::Quad84):
            vN= mpShapeQuad->N84(igauss);
            break;
        case(ShapeType::Quad89):
            vN= mpShapeQuad->N89(igauss);
            break;
        case(ShapeType::Triangle31):
            vN= mpShapeTriangle->N31(igauss);
            break;
        case(ShapeType::Triangle63):
            vN= mpShapeTriangle->N63(igauss);
            break;
        case(ShapeType::Line21):
            vN= mpShapeLine->N21(igauss);
            break;
        case(ShapeType::Line32):
            vN= mpShapeLine->N32(igauss);
            break;
        default:
            break;
    }
    uiint i;
    for(i=0; i < vN.size(); i++) N[i]=vN[i];
}
void CMW::ShapeFunc(const uiint& shapeType, vvdouble& N)
{
    switch(shapeType){
        case(ShapeType::Hexa81):
            N= mpShapeHexa->N81();
            break;
        case(ShapeType::Hexa82):
            N= mpShapeHexa->N82();
            break;
        case(ShapeType::Hexa201):
            N= mpShapeHexa->N201();
            break;
        case(ShapeType::Hexa202):
            N= mpShapeHexa->N202();
            break;
        case(ShapeType::Hexa203):
            N= mpShapeHexa->N203();
            break;
        case(ShapeType::HexaNic111):
            N= mpShapeHexaNic->N111();
            break;
        case(ShapeType::HexaNic118):
            N= mpShapeHexaNic->N118();
            break;
        case(ShapeType::HexaNic1127):
            N= mpShapeHexaNic->N1127();
            break;
        case(ShapeType::Tetra41):
            N= mpShapeTetra->N41();
            break;
        case(ShapeType::Tetra101):
            N= mpShapeTetra->N101();
            break;
        case(ShapeType::Tetra104):
            N= mpShapeTetra->N104();
            break;
        case(ShapeType::Tetra1015):
            N= mpShapeTetra->N1015();
            break;
        case(ShapeType::Prism62):
            N= mpShapePrism->N62();
            break;
        case(ShapeType::Prism156):
            N= mpShapePrism->N156();
            break;
        case(ShapeType::Prism159):
            N= mpShapePrism->N159();
            break;
        case(ShapeType::Prism1518):
            N= mpShapePrism->N1518();
            break;
        case(ShapeType::Quad41):
            N= mpShapeQuad->N41();
            break;
        case(ShapeType::Quad84):
            N= mpShapeQuad->N84();
            break;
        case(ShapeType::Quad89):
            N= mpShapeQuad->N89();
            break;
        case(ShapeType::Triangle31):
            N= mpShapeTriangle->N31();
            break;
        case(ShapeType::Triangle63):
            N= mpShapeTriangle->N63();
            break;
        case(ShapeType::Line21):
            N= mpShapeLine->N21();
            break;
        case(ShapeType::Line32):
            N= mpShapeLine->N32();
            break;
        default:
            break;
    }
}
double& CMW::ShapeFunc_Hexa81(uiint igauss, uiint ishape)
{
    return mpShapeHexa->N81(igauss,ishape);
}
double& CMW::ShapeFunc_Hexa82(uiint igauss, uiint ishape)
{
    return mpShapeHexa->N82(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa201(uiint igauss, uiint ishape)
{
    return mpShapeHexa->N201(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa202(uiint igauss, uiint ishape)
{
    return mpShapeHexa->N202(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa203(uiint igauss, uiint ishape)
{
    return mpShapeHexa->N203(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra41(uiint igauss, uiint ishape)
{
    return mpShapeTetra->N41(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra101(uiint igauss, uiint ishape)
{
    return mpShapeTetra->N101(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra104(uiint igauss, uiint ishape)
{
    return mpShapeTetra->N104(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra1015(uiint igauss, uiint ishape)
{
    return mpShapeTetra->N1015(igauss, ishape);
}
double& CMW::ShapeFunc_Prism62(uiint igauss, uiint ishape)
{
    return mpShapePrism->N62(igauss, ishape);
}
double& CMW::ShapeFunc_Prism156(uiint igauss, uiint ishape)
{
    return mpShapePrism->N156(igauss, ishape);
}
double& CMW::ShapeFunc_Prism159(uiint igauss, uiint ishape)
{
    return mpShapePrism->N159(igauss, ishape);
}
double& CMW::ShapeFunc_Prism1518(uiint igauss, uiint ishape)
{
    return mpShapePrism->N1518(igauss, ishape);
}
double& CMW::ShapeFunc_Quad41(uiint igauss, uiint ishape)
{
    return mpShapeQuad->N41(igauss, ishape);
}
double& CMW::ShapeFunc_Quad84(uiint igauss, uiint ishape)
{
    return mpShapeQuad->N84(igauss, ishape);
}
double& CMW::ShapeFunc_Quad89(uiint igauss, uiint ishape)
{
    return mpShapeQuad->N89(igauss, ishape);
}
double& CMW::ShapeFunc_Triangle31(uiint igauss, uiint ishape)
{
    return mpShapeTriangle->N31(igauss, ishape);
}
double& CMW::ShapeFunc_Triangle63(uiint igauss, uiint ishape)
{
    return mpShapeTriangle->N63(igauss, ishape);
}
double& CMW::ShapeFunc_Line21(uiint igauss, uiint ishape)
{
    return mpShapeLine->N21(igauss, ishape);
}
double& CMW::ShapeFunc_Line32(uiint igauss, uiint ishape)
{
    return mpShapeLine->N32(igauss, ishape);
}
void CMW::dNdr_on_pt(const uiint& shapeType, const uiint& igauss, vvdouble& dNdr)
{
    switch(shapeType){
        case(ShapeType::Hexa81):
            dNdr= mpShapeHexa->dNdr81(igauss);
            break;
        case(ShapeType::Hexa82):
            dNdr= mpShapeHexa->dNdr82(igauss);
            break;
        case(ShapeType::Hexa201):
            dNdr= mpShapeHexa->dNdr201(igauss);
            break;
        case(ShapeType::Hexa202):
            dNdr= mpShapeHexa->dNdr202(igauss);
            break;
        case(ShapeType::Hexa203):
            dNdr= mpShapeHexa->dNdr203(igauss);
            break;
        case(ShapeType::HexaNic111):
            dNdr= mpShapeHexaNic->dNdr111(igauss);
            break;
        case(ShapeType::HexaNic118):
            dNdr= mpShapeHexaNic->dNdr118(igauss);
            break;
        case(ShapeType::HexaNic1127):
            dNdr= mpShapeHexaNic->dNdr1127(igauss);
            break;
        case(ShapeType::Tetra41):
            dNdr= mpShapeTetra->dNdr41(igauss);
            break;
        case(ShapeType::Tetra101):
            dNdr= mpShapeTetra->dNdr101(igauss);
            break;
        case(ShapeType::Tetra104):
            dNdr= mpShapeTetra->dNdr104(igauss);
            break;
        case(ShapeType::Tetra1015):
            dNdr= mpShapeTetra->dNdr1015(igauss);
            break;
        case(ShapeType::Prism62):
            dNdr= mpShapePrism->dNdr62(igauss);
            break;
        case(ShapeType::Prism156):
            dNdr= mpShapePrism->dNdr156(igauss);
            break;
        case(ShapeType::Prism159):
            dNdr= mpShapePrism->dNdr159(igauss);
            break;
        case(ShapeType::Prism1518):
            dNdr= mpShapePrism->dNdr1518(igauss);
            break;
        case(ShapeType::Quad41):
            dNdr= mpShapeQuad->dNdr41(igauss);
            break;
        case(ShapeType::Quad84):
            dNdr= mpShapeQuad->dNdr84(igauss);
            break;
        case(ShapeType::Quad89):
            dNdr= mpShapeQuad->dNdr89(igauss);
            break;
        case(ShapeType::Triangle31):
            dNdr= mpShapeTriangle->dNdr31(igauss);
            break;
        case(ShapeType::Triangle63):
            dNdr= mpShapeTriangle->dNdr63(igauss);
            break;
        case(ShapeType::Line21):
            dNdr= mpShapeLine->dNdr21(igauss);
            break;
        case(ShapeType::Line32):
            dNdr= mpShapeLine->dNdr32(igauss);
            break;
        default:
            break;
    }
}
void CMW::dNdr(const uiint& shapeType, vvvdouble& dNdr)
{
    switch(shapeType){
        case(ShapeType::Hexa81):
            dNdr= mpShapeHexa->dNdr81();
            break;
        case(ShapeType::Hexa82):
            dNdr= mpShapeHexa->dNdr82();
            break;
        case(ShapeType::Hexa201):
            dNdr= mpShapeHexa->dNdr201();
            break;
        case(ShapeType::Hexa202):
            dNdr= mpShapeHexa->dNdr202();
            break;
        case(ShapeType::Hexa203):
            dNdr= mpShapeHexa->dNdr203();
            break;
        case(ShapeType::Tetra41):
            dNdr= mpShapeTetra->dNdr41();
            break;
        case(ShapeType::Tetra101):
            dNdr= mpShapeTetra->dNdr101();
            break;
        case(ShapeType::Tetra104):
            dNdr= mpShapeTetra->dNdr104();
            break;
        case(ShapeType::Tetra1015):
            dNdr= mpShapeTetra->dNdr1015();
            break;
        case(ShapeType::Prism62):
            dNdr= mpShapePrism->dNdr62();
            break;
        case(ShapeType::Prism156):
            dNdr= mpShapePrism->dNdr156();
            break;
        case(ShapeType::Prism159):
            dNdr= mpShapePrism->dNdr159();
            break;
        case(ShapeType::Prism1518):
            dNdr= mpShapePrism->dNdr1518();
            break;
        case(ShapeType::Quad41):
            dNdr= mpShapeQuad->dNdr41();
            break;
        case(ShapeType::Quad84):
            dNdr= mpShapeQuad->dNdr84();
            break;
        case(ShapeType::Quad89):
            dNdr= mpShapeQuad->dNdr89();
            break;
        case(ShapeType::Triangle31):
            dNdr= mpShapeTriangle->dNdr31();
            break;
        case(ShapeType::Triangle63):
            dNdr= mpShapeTriangle->dNdr63();
            break;
        case(ShapeType::Line21):
            dNdr= mpShapeLine->dNdr21();
            break;
        case(ShapeType::Line32):
            dNdr= mpShapeLine->dNdr32();
            break;
        default:
            break;
    }
}
void CMW::dNdr(const uiint& shapeType, double dNdr[])
{
    vvvdouble vdNdr;
    uiint numOfInteg, numOfShape, numOfAxis;
    switch(shapeType){
        case(ShapeType::Hexa81):
            vdNdr= mpShapeHexa->dNdr81();
            numOfInteg = 1; numOfShape = 8; numOfAxis = 3;
            break;
        case(ShapeType::Hexa82):
            vdNdr= mpShapeHexa->dNdr82();
            numOfInteg = 8; numOfShape = 8; numOfAxis = 3;
            break;
        case(ShapeType::Hexa201):
            vdNdr= mpShapeHexa->dNdr201();
            numOfInteg = 1; numOfShape = 20; numOfAxis = 3;
            break;
        case(ShapeType::Hexa202):
            vdNdr= mpShapeHexa->dNdr202();
            numOfInteg = 8; numOfShape = 20; numOfAxis = 3;
            break;
        case(ShapeType::Hexa203):
            vdNdr= mpShapeHexa->dNdr203();
            numOfInteg = 27; numOfShape = 20; numOfAxis = 3;
            break;
        case(ShapeType::Tetra41):
            vdNdr= mpShapeTetra->dNdr41();
            numOfInteg = 1; numOfShape = 4; numOfAxis = 3;
            break;
        case(ShapeType::Tetra101):
            vdNdr= mpShapeTetra->dNdr101();
            numOfInteg = 1; numOfShape = 10; numOfAxis = 3;
            break;
        case(ShapeType::Tetra104):
            vdNdr= mpShapeTetra->dNdr104();
            numOfInteg = 4; numOfShape = 10; numOfAxis = 3;
            break;
        case(ShapeType::Tetra1015):
            vdNdr= mpShapeTetra->dNdr1015();
            numOfInteg = 15; numOfShape = 10; numOfAxis = 3;
            break;
        case(ShapeType::Prism62):
            vdNdr= mpShapePrism->dNdr62();
            numOfInteg = 2; numOfShape = 6; numOfAxis = 3;
            break;
        case(ShapeType::Prism156):
            vdNdr= mpShapePrism->dNdr156();
            numOfInteg = 6; numOfShape = 15; numOfAxis = 3;
            break;
        case(ShapeType::Prism159):
            vdNdr= mpShapePrism->dNdr159();
            numOfInteg = 9; numOfShape = 15; numOfAxis = 3;
            break;
        case(ShapeType::Prism1518):
            vdNdr= mpShapePrism->dNdr1518();
            numOfInteg = 18; numOfShape = 15; numOfAxis = 3;
            break;
        case(ShapeType::Quad41):
            vdNdr= mpShapeQuad->dNdr41();
            numOfInteg = 1; numOfShape = 4; numOfAxis = 2;
            break;
        case(ShapeType::Quad84):
            vdNdr= mpShapeQuad->dNdr84();
            numOfInteg = 4; numOfShape = 8; numOfAxis = 2;
            break;
        case(ShapeType::Quad89):
            vdNdr= mpShapeQuad->dNdr89();
            numOfInteg = 9; numOfShape = 8; numOfAxis = 2;
            break;
        case(ShapeType::Triangle31):
            vdNdr= mpShapeTriangle->dNdr31();
            numOfInteg = 1; numOfShape = 3; numOfAxis = 2;
            break;
        case(ShapeType::Triangle63):
            vdNdr= mpShapeTriangle->dNdr63();
            numOfInteg = 3; numOfShape = 6; numOfAxis = 2;
            break;
        case(ShapeType::Line21):
            vdNdr= mpShapeLine->dNdr21();
            numOfInteg = 1; numOfShape = 2; numOfAxis = 1;
            break;
        case(ShapeType::Line32):
            vdNdr= mpShapeLine->dNdr32();
            numOfInteg = 2; numOfShape = 3; numOfAxis = 1;
            break;
        default:
            break;
    }
    uiint igauss,ishape,iaxis;
    uiint pos=0;
    for(igauss=0; igauss < numOfInteg; igauss++){
        for(ishape=0; ishape < numOfShape; ishape++){
            for(iaxis=0; iaxis < numOfAxis; iaxis++){
                dNdr[pos] = vdNdr[igauss][ishape][iaxis];
                pos++;
            }
        }
    }
}
double& CMW::dNdr_Hexa81_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeHexa->dNdr81(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa82_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeHexa->dNdr82(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa201_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeHexa->dNdr201(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa202_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeHexa->dNdr202(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa203_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeHexa->dNdr203(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTetra->dNdr41(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra101_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTetra->dNdr101(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra104_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTetra->dNdr104(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra1015_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTetra->dNdr1015(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism62_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapePrism->dNdr62(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism156_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapePrism->dNdr156(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism159_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapePrism->dNdr159(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism1518_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapePrism->dNdr1518(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeQuad->dNdr41(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad84_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeQuad->dNdr84(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad89_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeQuad->dNdr89(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tri31_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTriangle->dNdr31(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tri63_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis)
{
    return mpShapeTriangle->dNdr63(igauss, ishape, iaxis);
}
double& CMW::dNdr_Line21_on_pt_on_shape(uiint igauss, uiint ishape)
{
    return mpShapeLine->dNdr21(igauss, ishape);
}
double& CMW::dNdr_Line32_on_pt_on_shape(uiint igauss, uiint ishape)
{
    return mpShapeLine->dNdr32(igauss, ishape);
}
void CMW::Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index)
{
    mvdNdx.clear();
    CElement *pElement;
    pElement= mpMesh->getElementIX(elem_index);
    switch(elemType){
        case(ElementType::Hexa):
            if(numOfInteg==1){
                mpShapeHexa->Calc_dNdx8(numOfInteg,pElement);
                mvdNdx= mpShapeHexa->dNdx81();
            }
            if(numOfInteg==8){
                mpShapeHexa->Calc_dNdx8(numOfInteg,pElement);
                mvdNdx= mpShapeHexa->dNdx82();
            }
            break;
        case(ElementType::Hexa2):
            if(numOfInteg==1){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                mvdNdx= mpShapeHexa->dNdx201();
            }
            if(numOfInteg==8){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                mvdNdx= mpShapeHexa->dNdx202();
            }
            if(numOfInteg==27){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                mvdNdx= mpShapeHexa->dNdx203();
            }
            break;
        case(ElementType::Tetra):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx4(numOfInteg,pElement);
                mvdNdx= mpShapeTetra->dNdx41();
            }
            break;
        case(ElementType::Tetra2):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                mvdNdx= mpShapeTetra->dNdx101();
            }
            if(numOfInteg==4){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                mvdNdx= mpShapeTetra->dNdx104();
            }
            if(numOfInteg==15){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                mvdNdx= mpShapeTetra->dNdx1015();
            }
            break;
        case(ElementType::Prism):
            if(numOfInteg==2){
                mpShapePrism->Calc_dNdx6(numOfInteg,pElement);
                mvdNdx= mpShapePrism->dNdx62();
            }
            break;
        case(ElementType::Prism2):
            if(numOfInteg==6){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                mvdNdx= mpShapePrism->dNdx156();
            }
            if(numOfInteg==9){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                mvdNdx= mpShapePrism->dNdx159();
            }
            if(numOfInteg==18){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                mvdNdx= mpShapePrism->dNdx1518();
            }
            break;
        default:
            break;
    }
}
void CMW::dNdx_on_pt(const uiint& igauss, vvdouble& dNdX)
{
    dNdX = mvdNdx[igauss];
}
void CMW::dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdX)
{
    CElement *pElement;
    pElement= mpMesh->getElementIX(elem_index);
    switch(elemType){
        case(ElementType::Hexa):
            if(numOfInteg==1){
                mpShapeHexa->Calc_dNdx8(numOfInteg,pElement);
                dNdX= mpShapeHexa->dNdx81();
            }
            if(numOfInteg==8){
                mpShapeHexa->Calc_dNdx8(numOfInteg,pElement);
                dNdX= mpShapeHexa->dNdx82();
            }
            break;
        case(ElementType::Hexa2):
            if(numOfInteg==1){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                dNdX= mpShapeHexa->dNdx201();
            }
            if(numOfInteg==8){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                dNdX= mpShapeHexa->dNdx202();
            }
            if(numOfInteg==27){
                mpShapeHexa->Calc_dNdx20(numOfInteg,pElement);
                dNdX= mpShapeHexa->dNdx203();
            }
            break;
        case(ElementType::Tetra):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx4(numOfInteg,pElement);
                dNdX= mpShapeTetra->dNdx41();
            }
            break;
        case(ElementType::Tetra2):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                dNdX= mpShapeTetra->dNdx101();
            }
            if(numOfInteg==4){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                dNdX= mpShapeTetra->dNdx104();
            }
            if(numOfInteg==15){
                mpShapeTetra->Calc_dNdx10(numOfInteg,pElement);
                dNdX= mpShapeTetra->dNdx1015();
            }
            break;
        case(ElementType::Prism):
            if(numOfInteg==2){
                mpShapePrism->Calc_dNdx6(numOfInteg,pElement);
                dNdX= mpShapePrism->dNdx62();
            }
            break;
        case(ElementType::Prism2):
            if(numOfInteg==6){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                dNdX= mpShapePrism->dNdx156();
            }
            if(numOfInteg==9){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                dNdX= mpShapePrism->dNdx159();
            }
            if(numOfInteg==18){
                mpShapePrism->Calc_dNdx15(numOfInteg,pElement);
                dNdX= mpShapePrism->dNdx1518();
            }
            break;
        default:
            break;
    }
}
void CMW::dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& ielem, double dNdx[])
{
    Calculate_dNdx(elemType, numOfInteg, ielem);
    uiint numOfShape;
    switch(elemType){
        case(ElementType::Hexa):
            numOfShape = NumberOfNode::Hexa();
            break;
        case(ElementType::Hexa2):
            numOfShape = NumberOfNode::Hexa2();
            break;
        case(ElementType::Tetra):
            numOfShape = NumberOfNode::Tetra();
            break;
        case(ElementType::Tetra2):
            numOfShape = NumberOfNode::Tetra2();
            break;
        case(ElementType::Prism):
            numOfShape = NumberOfNode::Prism();
            break;
        case(ElementType::Prism2):
            numOfShape = NumberOfNode::Prism2();
            break;
        default:
            break;
    }
    uiint igauss, ishape, iaxis;
    uiint pos=0;
    for(igauss=0; igauss < numOfInteg; igauss++){
        for(ishape=0; ishape < numOfShape; ishape++){
            for(iaxis=0; iaxis < 3; iaxis++){
                dNdx[pos] = mvdNdx[igauss][ishape][iaxis];
                pos++;
            };
        };
    };
}
void CMW::detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ)
{
    switch(elemType){
        case(ElementType::Hexa):
            if(numOfInteg==1){
                detJ= mpShapeHexa->detJ81(igauss);
            }
            if(numOfInteg==8){
                detJ= mpShapeHexa->detJ82(igauss);
            }
            break;
        case(ElementType::Hexa2):
            if(numOfInteg==1){
                detJ= mpShapeHexa->detJ201(igauss);
            }
            if(numOfInteg==8){
                detJ= mpShapeHexa->detJ202(igauss);
            }
            if(numOfInteg==27){
                detJ= mpShapeHexa->detJ203(igauss);
            }
            break;
        case(ElementType::Tetra):
            if(numOfInteg==1){
                detJ= mpShapeTetra->detJ41(igauss);
            }
            break;
        case(ElementType::Tetra2):
            if(numOfInteg==1){
                detJ= mpShapeTetra->detJ101(igauss);
            }
            if(numOfInteg==4){
                detJ= mpShapeTetra->detJ104(igauss);
            }
            if(numOfInteg==15){
                detJ= mpShapeTetra->detJ1015(igauss);
            }
            break;
        case(ElementType::Prism):
            if(numOfInteg==2){
                detJ= mpShapePrism->detJ62(igauss);
            }
            break;
        case(ElementType::Prism2):
            if(numOfInteg==6){
                detJ= mpShapePrism->detJ156(igauss);
            }
            if(numOfInteg==9){
                detJ= mpShapePrism->detJ159(igauss);
            }
            if(numOfInteg==18){
                detJ= mpShapePrism->detJ1518(igauss);
            }
            break;
        default:
            break;
    }
}
void CMW::Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w)
{
    switch(elemType){
        case(ElementType::Hexa):case(ElementType::Hexa2):
            if(numOfInteg==1){
                w= mpShapeHexa->Weight3dpt1();
            }
            if(numOfInteg==8){
                w= mpShapeHexa->Weight3dpt2(igauss);
            }
            if(numOfInteg==27){
                w= mpShapeHexa->Weight3dpt3(igauss);
            }
            break;
        case(ElementType::Tetra):case(ElementType::Tetra2):
            if(numOfInteg==1){
                w= mpShapeTetra->Weight_pt1();
            }
            if(numOfInteg==4){
                w= mpShapeTetra->Weight_pt4(igauss);
            }
            if(numOfInteg==15){
                w= mpShapeTetra->Weight_pt15(igauss);
            }
            break;
        case(ElementType::Prism):case(ElementType::Prism2):
            if(numOfInteg==2){
                w= mpShapePrism->Weight_pt2(igauss);
            }
            if(numOfInteg==6){
                w= mpShapePrism->Weight_pt6(igauss);
            }
            if(numOfInteg==9){
                w= mpShapePrism->Weight_pt9(igauss);
            }
            if(numOfInteg==18){
                w= mpShapePrism->Weight_pt18(igauss);
            }
            break;
        default:
            break;
    }
}
uiint CMW::nodetype_s(){ return NodeType::Scalar;}
uiint CMW::nodetype_v(){ return NodeType::Vector;}
uiint CMW::nodetype_sv(){ return NodeType::ScalarVector;}

uiint CMW::elemtype_hexa(){ return ElementType::Hexa;}
uiint CMW::elemtype_hexa2(){ return ElementType::Hexa2;}
uiint CMW::elemtype_tetra(){ return ElementType::Tetra;}
uiint CMW::elemtype_tetra2(){ return ElementType::Tetra2;}
uiint CMW::elemtype_prism(){ return ElementType::Prism;}
uiint CMW::elemtype_prism2(){ return ElementType::Prism2;}
uiint CMW::elemtype_quad(){ return ElementType::Quad;}
uiint CMW::elemtype_quad2(){ return ElementType::Quad2;}
uiint CMW::elemtype_triangle(){ return ElementType::Triangle;}
uiint CMW::elemtype_triangle2(){ return ElementType::Triangle2;}
uiint CMW::elemtype_line(){ return ElementType::Beam;}
uiint CMW::elemtype_line2(){ return ElementType::Beam2;}

uiint CMW::fistr_elemtype_hexa(){ return FistrElementType::Hexa;}
uiint CMW::fistr_elemtype_hexa2(){ return FistrElementType::Hexa2;}
uiint CMW::fistr_elemtype_tetra(){ return FistrElementType::Tetra;}
uiint CMW::fistr_elemtype_tetra2(){ return FistrElementType::Tetra2;}
uiint CMW::fistr_elemtype_prism(){ return FistrElementType::Prism;}
uiint CMW::fistr_elemtype_prism2(){ return FistrElementType::Prism2;}
uiint CMW::fistr_elemtype_quad(){ return FistrElementType::Quad;}
uiint CMW::fistr_elemtype_quad2(){ return FistrElementType::Quad2;}
uiint CMW::fistr_elemtype_triangle(){ return FistrElementType::Triangle;}
uiint CMW::fistr_elemtype_triangle2(){ return FistrElementType::Triangle2;}
uiint CMW::fistr_elemtype_line(){ return FistrElementType::Beam;}
uiint CMW::fistr_elemtype_line2(){ return FistrElementType::Beam2;}
uiint CMW::fistr_elemtype_to_mw3_elemtype(const uiint& fistr_elemtype)
{
    switch(fistr_elemtype){
        case(FistrElementType::Hexa):
            return ElementType::Hexa;
        case(FistrElementType::Hexa2):
            return ElementType::Hexa2;
        case(FistrElementType::Tetra):
            return ElementType::Tetra;
        case(FistrElementType::Tetra2):
            return ElementType::Tetra2;
        case(FistrElementType::Prism):
            return ElementType::Prism;
        case(FistrElementType::Prism2):
            return ElementType::Prism2;
        case(FistrElementType::Quad):
            return ElementType::Quad;
        case(FistrElementType::Quad2):
            return ElementType::Quad2;
        case(FistrElementType::Triangle):
            return ElementType::Triangle;
        case(FistrElementType::Triangle2):
            return ElementType::Triangle2;
        case(FistrElementType::Beam):
            return ElementType::Beam;
        case(FistrElementType::Beam2):
            return ElementType::Beam2;
        case(FistrElementType::Point):
            return ElementType::Point;
        case(FistrElementType::Limit):
            return ElementType::Limit;
    }
}
uiint CMW::mw3_elemtype_to_fistr_elemtype(const uiint& mw3_elemtype)
{
    switch(mw3_elemtype){
        case(ElementType::Hexa):
            return FistrElementType::Hexa;
        case(ElementType::Hexa2):
            return FistrElementType::Hexa2;
        case(ElementType::Tetra):
            return FistrElementType::Tetra;
        case(ElementType::Tetra2):
            return FistrElementType::Tetra2;
        case(ElementType::Prism):
            return FistrElementType::Prism;
        case(ElementType::Prism2):
            return FistrElementType::Prism2;
        case(ElementType::Quad):
            return FistrElementType::Quad;
        case(ElementType::Quad2):
            return FistrElementType::Quad2;
        case(ElementType::Triangle):
            return FistrElementType::Triangle;
        case(ElementType::Triangle2):
            return FistrElementType::Triangle2;
        case(ElementType::Beam):
            return FistrElementType::Beam;
        case(ElementType::Beam2):
            return FistrElementType::Beam2;
        case(ElementType::Point):
            return FistrElementType::Point;
        case(ElementType::Limit):
            return FistrElementType::Limit;
    }
}
uiint CMW::shapetype_hexa81(){ return ShapeType::Hexa81;}
uiint CMW::shapetype_hexa82(){ return ShapeType::Hexa82;}
uiint CMW::shapetype_hexa201(){ return ShapeType::Hexa201;}
uiint CMW::shapetype_hexa202(){ return ShapeType::Hexa202;}
uiint CMW::shapetype_hexa203(){ return ShapeType::Hexa203;}
uiint CMW::shapetype_tetra41(){ return ShapeType::Tetra41;}
uiint CMW::shapetype_tetra101(){ return ShapeType::Tetra101;}
uiint CMW::shapetype_tetra104(){ return ShapeType::Tetra104;}
uiint CMW::shapetype_tetra1015(){ return ShapeType::Tetra1015;}
uiint CMW::shapetype_prism62(){ return ShapeType::Prism62;}
uiint CMW::shapetype_prism156(){ return ShapeType::Prism156;}
uiint CMW::shapetype_prism159(){ return ShapeType::Prism159;}
uiint CMW::shapetype_prism1518(){ return ShapeType::Prism1518;}
uiint CMW::shapetype_quad41(){ return ShapeType::Quad41;}
uiint CMW::shapetype_quad84(){ return ShapeType::Quad84;}
uiint CMW::shapetype_quad89(){ return ShapeType::Quad89;}
uiint CMW::shapetype_tri31(){ return ShapeType::Triangle31;}
uiint CMW::shapetype_tri63(){ return ShapeType::Triangle63;}
uiint CMW::shapetype_line21(){ return ShapeType::Line21;}
uiint CMW::shapetype_line32(){ return ShapeType::Line32;}

uiint CMW::GetNumOfBoundaryNodeMesh(){ return mpMesh->getNumOfBoundaryNodeMesh();}
uiint CMW::GetNumOfBoundaryFaceMesh(){ return mpMesh->getNumOfBoundaryFaceMesh();}
uiint CMW::GetNumOfBoundaryEdgeMesh(){ return mpMesh->getNumOfBoundaryEdgeMesh();}
uiint CMW::GetNumOfBoundaryVolumeMesh(){ return mpMesh->getNumOfBoundaryVolumeMesh();}
uiint CMW::GetBNDType_BNodeMesh(const uiint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    return pBNodeMesh->getBndType();
}
uiint CMW::GetBNDType_BFaceMesh(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getBndType();
}
uiint CMW::GetBNDType_BEdgeMesh(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getBndType();
}
uiint CMW::GetBNDType_BVolumeMesh(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getBndType();
}
uiint CMW::getNeumannType(){ return BoundaryType::Neumann;}
uiint CMW::getDirichletType(){ return BoundaryType::Dirichlet;}
uiint CMW::GetNumOfBNode_BNodeMesh(const uiint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    return pBNodeMesh->getNumOfBNode();
}
uiint CMW::GetNumOfBNode_BFaceMesh(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getNumOfBNode();
}
uiint CMW::GetNumOfBNode_BEdgeMesh(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getNumOfBNode();
}
uiint CMW::GetNumOfBNode_BVolumeMesh(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getNumOfBNode();
}
uiint CMW::GetNumOfDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    return pBNode->getNumOfDOF();
}
uiint CMW::GetNumOfDOF_BFaceMesh(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getNumOfDOF();
}
uiint CMW::GetNumOfDOF_BEdgeMesh(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getNumOfDOF();
}
uiint CMW::GetNumOfDOF_BVolumeMesh(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getNumOfDOF();
}
uiint CMW::GetDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& idof)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    return pBNode->getDOF(idof);
}
uiint CMW::GetDOF_BFaceMesh(const uiint& ibmesh, const uiint& idof)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getDOF(idof);
}
uiint CMW::GetDOF_BEdgeMesh(const uiint& ibmesh, const uiint& idof)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getDOF(idof);
}
uiint CMW::GetDOF_BVolumeMesh(const uiint& ibmesh, const uiint& idof)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getDOF(idof);
}
double& CMW::GetBNodeValue_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    
    return pBNode->getValue(dof);
}
double& CMW::GetBNodeValue_BFaceMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);
    return pBNode->getValue(dof, mgLevel);
}
double& CMW::GetBNodeValue_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBEdgeMesh->getBNodeIX(ibnode);
    return pBNode->getValue(dof, mgLevel);
}
double& CMW::GetBNodeValue_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBVolMesh->getBNodeIX(ibnode);
    return pBNode->getValue(dof, mgLevel);
}
uiint& CMW::GetNodeID_BNode_BNodeMesh(const uiint& ibmesh, const uiint& ibnode)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    CNode *pNode;
    pNode = pBNode->getNode();
    return pNode->getID();
}
uiint& CMW::GetNodeID_BNode_BFaceMesh(const uiint& ibmesh, const uiint& ibnode)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);
    CNode *pNode;
    pNode = pBNode->getNode();
    return pNode->getID();
}
uiint& CMW::GetNodeID_BNode_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBEdgeMesh->getBNodeIX(ibnode);
    CNode *pNode;
    pNode = pBNode->getNode();
    return pNode->getID();
}
uiint& CMW::GetNodeID_BNode_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBVolMesh->getBNodeIX(ibnode);
    CNode *pNode;
    pNode = pBNode->getNode();
    return pNode->getID();
}

#ifndef MSVC
//--
// Test
//--
double& CMW::GetBNode_X_BFaceMesh(const uiint& ibmesh, const uiint& ibnode)//Test
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);
    
    return pBNode->getX();
}
double& CMW::GetBNode_Y_BFaceMesh(const uiint& ibmesh, const uiint& ibnode)//Test
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);

    return pBNode->getY();
}
double& CMW::GetBNode_Z_BFaceMesh(const uiint& ibmesh, const uiint& ibnode)//Test
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);

    return pBNode->getZ();
}
#endif

uiint CMW::GetNumOfBFace(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getNumOfBFace();
}
double& CMW::GetBFaceValue(const uiint& ibmesh, const uiint& ibface, const uiint& dof)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryFace *pBFace;
    pBFace = pBFaceMesh->getBFaceIX(ibface);
    return pBFace->getBndValue(dof);
}
uiint CMW::GetNumOfBEdge(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getNumOfEdge();
}
double& CMW::GetBEdgeValue(const uiint& ibmesh, const uiint& ibedge, const uiint& dof)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryEdge *pBEdge;
    pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);
    return pBEdge->getBndValue(dof);
}
uiint CMW::GetNumOfBVolume(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getNumOfVolume();
}
double& CMW::GetBVolumeValue(const uiint& ibmesh, const uiint& ibvol, const uiint& dof)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryVolume *pBVol;
    pBVol = pBVolMesh->getBVolumeIX(ibvol);
    return pBVol->getBndValue(dof);
}
uiint CMW::GetNumOfNode_BFace(const uiint& ibmesh, const uiint& ibface)
{
    CBoundaryFaceMesh *pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryFace *pBFace = pBFaceMesh->getBFaceIX(ibface);
    return pBFace->getNumOfBNode();
}
uiint& CMW::GetNodeID_BFace(const uiint& ibmesh, const uiint& ibface, const uiint& ibnode)
{
    CBoundaryFaceMesh *pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryFace *pBFace = pBFaceMesh->getBFaceIX(ibface);
    CBoundaryNode *pBNode = pBFace->getBNode(ibnode);
    CNode *pNode = pBNode->getNode();
    return pNode->getID();
}
uiint CMW::GetNumOfNode_BEdge(const uiint& ibmesh, const uiint& ibedge)
{
    CBoundaryEdgeMesh *pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryEdge *pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);
    return pBEdge->getNumOfBNode();
}
uiint& CMW::GetNodeID_BEdge(const uiint& ibmesh, const uiint& ibedge, const uiint& ibnode)
{
    CBoundaryEdgeMesh *pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryEdge *pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);
    CBoundaryNode *pBNode = pBEdge->getBNode(ibnode);
    CNode *pNode = pBNode->getNode();
    return pNode->getID();
}
uiint CMW::GetNumOfNode_BVolume(const uiint& ibmesh, const uiint& ibvol)
{
    CBoundaryVolumeMesh *pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryVolume *pBVol = pBVolMesh->getBVolumeIX(ibvol);
    return pBVol->getNumOfBNode();
}
uiint& CMW::GetNodeID_BVolume(const uiint& ibmesh, const uiint& ibvol, const uiint& ibnode)
{
    CBoundaryVolumeMesh *pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryVolume *pBVol = pBVolMesh->getBVolumeIX(ibvol);
    CBoundaryNode *pBNode = pBVol->getBNode(ibnode);
    CNode *pNode = pBNode->getNode();
    return pNode->getID();
}
uiint CMW::GetBNodeMesh_NameLength(const uiint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    string sBndName = pBNodeMesh->getName();
    return sBndName.length();
}
string& CMW::GetBNodeMesh_Name(const uiint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    return pBNodeMesh->getName();
}
uiint CMW::GetBFaceMesh_NameLength(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    string sBndName = pBFaceMesh->getName();
    return sBndName.length();
}
string& CMW::GetBFaceMesh_Name(const uiint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    return pBFaceMesh->getName();
}
uiint CMW::GetBVolumeMesh_NameLength(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    string sBndName = pBVolMesh->getName();
    return sBndName.length();
}
string& CMW::GetBVolumeMesh_Name(const uiint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    return pBVolMesh->getName();
}
uiint CMW::GetBEdgeMesh_NameLength(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    string sBndName = pBEdgeMesh->getName();
    return sBndName.length();
}
string& CMW::GetBEdgeMesh_Name(const uiint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    return pBEdgeMesh->getName();
}
uiint CMW::GetEdgeID_BEdge(const uiint& ibmesh, const uiint& ibedge)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryEdge *pBEdge;
    pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);
    return pBEdge->getElementEdgeID();
}
uiint CMW::GetElemID_BEdge(const uiint& ibmesh, const uiint& ibedge)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    CBoundaryEdge *pBEdge;
    pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);
    return pBEdge->getElementID();
}
uiint CMW::GetFaceID_BFace(const uiint& ibmesh, const uiint& ibface)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryFace *pBFace;
    pBFace = pBFaceMesh->getBFaceIX(ibface);
    return pBFace->getElementFaceID();
}
uiint CMW::GetElemID_BFace(const uiint& ibmesh, const uiint& ibface)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    CBoundaryFace *pBFace;
    pBFace = pBFaceMesh->getBFaceIX(ibface);
    return pBFace->getElementID();
}
uiint CMW::GetElemID_BVolume(const uiint& ibmesh, const uiint& ibvol)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    CBoundaryVolume *pBVol;
    pBVol = pBVolMesh->getBVolumeIX(ibvol);
    return pBVol->getElementID();
}

//--
// MPI
//--
int CMW::GetRank()
{
    return mnMyRank;
}
int CMW::GetNumOfProcess()
{
    return mnNumOfProcess;
}
iint CMW::AllReduce(void* sendbuf, void* recvbuf, iint buf_size, MPI_Datatype datatype, MPI_Op op, MPI_Comm commworld)
{
    return mpMPI->Allreduce(sendbuf, recvbuf, buf_size, datatype, op, commworld);
}
iint CMW::Barrier(MPI_Comm commworld)
{
    return mpMPI->Barrier(commworld);
}
iint CMW::Abort(MPI_Comm commworld, int error)
{
    return mpMPI->Abort(commworld, error);
}
iint CMW::AllGather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm)
{
    return mpMPI->Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
}
iint CMW::Gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    return mpMPI->Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm);
}
iint CMW::Scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    return mpMPI->Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}
iint CMW::Send(void* buf, iint count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
    return mpMPI->Send(buf, count, datatype, dest, tag, comm);
}
iint CMW::Recv(void* buf, iint count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status)
{
    return mpMPI->Recv(buf, count, datatype, source, tag, comm, status);
}
iint CMW::Bcast(void* buf, iint cnt, MPI_Datatype type, int root, MPI_Comm comm)
{
    return mpMPI->Bcast(buf, cnt, type, root, comm);
}
iint CMW::GetNumOfNeibPE(const uiint& imesh)
{
    CMesh *pMesh = mpAssy->getMesh(imesh);
    return pMesh->getCommMesh2Size();
}
iint CMW::GetTransRank(const uiint& imesh, const uiint& ipe)
{
    CMesh *pMesh = mpAssy->getMesh(imesh);
    CCommMesh2 *pCommMesh = pMesh->getCommMesh2IX(ipe);
    
    uiint nTransRank= pCommMesh->getTrasmitRank();
    
    if(nTransRank < IINT_MAX){
        return nTransRank;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, " MW::GetTransRank ");
        return MW_ERROR;
    }
}
void CMW::Send_Recv_R(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
{
    iint nCount = num_of_node * dof_size;
    MPI_Request Req[2];
    MPI_Status  Sta[2];

    double* sbuf = (double*)malloc(sizeof(double) * nCount);
    for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
    mpMPI->Isend(sbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD, &Req[0]);

    double* rbuf = (double*)malloc(sizeof(double) * nCount);
    mpMPI->Irecv(rbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD, &Req[1]);

    mpMPI->Waitall(1, &Req[1], &Sta[1]);
    mpMPI->Waitall(1, &Req[0], &Sta[0]);
    
    for(iint i=0; i < nCount; i++) buf[i]=rbuf[i];
    
    free( sbuf );
    free( rbuf );
}
void CMW::Send_Recv_I(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
{
    iint nCount = num_of_node * dof_size;
    MPI_Request Req[2];
    MPI_Status  Sta[2];

    MPI_Datatype MPI_IINT= mpMPI->MPI_IINT();

    iint* sbuf = (iint*)malloc(sizeof(iint) * nCount);
    for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
    mpMPI->Isend(sbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD, &Req[0]);
    
    iint* rbuf = (iint*)malloc(sizeof(iint) * nCount);
    mpMPI->Irecv(rbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD, &Req[1]);

    mpMPI->Waitall(1, &Req[1], &Sta[1]);
    mpMPI->Waitall(1, &Req[0], &Sta[0]);
    
    for(iint i=0; i < nCount; i++) buf[i]=rbuf[i];

    free( sbuf );
    free( rbuf );
}
void CMW::Sendrecv_R(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
{
    iint nCount = num_of_node * dof_size;
    MPI_Request Req[2];
    MPI_Status  Sta[2];

    double *x = (double*)malloc(sizeof(double) * nCount);
    double* sbuf = (double*)malloc(sizeof(double) * nCount);
    double* rbuf = (double*)malloc(sizeof(double) * nCount);

    if(mpMPI->getRank() < trans_rank ){
        for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
        mpMPI->Send(sbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD);

        mpMPI->Recv(rbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD, &Sta[0]);
        for(iint i=0; i < nCount; i++) x[i]=rbuf[i];
    }
    if(mpMPI->getRank() > trans_rank ){
        mpMPI->Recv(rbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD, &Sta[0]);
        for(iint i=0; i < nCount; i++) x[i]=rbuf[i];

        for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
        mpMPI->Send(sbuf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD);
    }

    for(iint i=0; i < nCount; i++) buf[i]=x[i];


    free( x );
    free( sbuf );
    free( rbuf );
}
void CMW::Sendrecv_I(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
{
    iint nCount = num_of_node * dof_size;
    MPI_Request Req[2];
    MPI_Status  Sta[2];

    MPI_Datatype MPI_IINT= mpMPI->MPI_IINT();

    iint *x = (iint*)malloc(sizeof(iint) * nCount);
    iint* sbuf = (iint*)malloc(sizeof(iint) * nCount);
    iint* rbuf = (iint*)malloc(sizeof(iint) * nCount);

    if(mpMPI->getRank() < trans_rank ){
        for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
        mpMPI->Send(sbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD);

        mpMPI->Recv(rbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD, &Sta[0]);
        for(iint i=0; i < nCount; i++) x[i]=rbuf[i];
    }
    if(mpMPI->getRank() > trans_rank ){
        mpMPI->Recv(rbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD, &Sta[0]);
        for(iint i=0; i < nCount; i++) x[i]=rbuf[i];

        for(iint i=0; i < nCount; i++) sbuf[i]=buf[i];
        mpMPI->Send(sbuf, nCount, MPI_IINT, trans_rank, 0, MPI_COMM_WORLD);
    }

    for(iint i=0; i < nCount; i++) buf[i]=x[i];


    free( x );
    free( sbuf );
    free( rbuf );
}
MPI_Datatype CMW::MPI_UIINT()
{
    return mnMPI_UIINT;
}
MPI_Datatype CMW::MPI_IINT()
{
    return mnMPI_IINT;
}


//--
// 通信テーブル:CommMesh2
//--
uiint CMW::GetNumOfCommMesh()
{
    return mpMesh->getCommMesh2Size();
}
uiint CMW::GetNumOfCommNode(const uiint& icmesh)
{
    mpComMesh= mpMesh->getCommMesh2IX(icmesh);
    return mpComMesh->getCommNodeSize();
}
uiint& CMW::GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode)
{
    mpComMesh= mpMesh->getCommMesh2IX(icmesh);
    CCommNode *pComNode;
    pComNode = mpComMesh->getCommNodeIX(icnode);
    return pComNode->getNodeID();
}

//--
// 接合面:ContactMesh
//--
uiint CMW::GetNumOfContactMesh()
{
    return mpAssy->getNumOfContactMesh();
}
uiint CMW::GetContactMeshID(const uiint& icont)
{
    return mpAssy->getContactID(icont);
}




uiint CMW::GetNumOfElementGroup()
{
    return mpMesh->getNumOfElemGrp();
}
uiint CMW::GetNumOfElementID(const uiint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);
    return pElemGrp->getNumOfElementID();
}
uiint& CMW::GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);
    return pElemGrp->getElementID(index);
}
uiint CMW::GetElementGroupName_Length(const uiint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);
    return pElemGrp->getNameLength();
}
string& CMW::GetElementGroupName(const uiint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);
    return pElemGrp->getName();
}



void CMW::LoggerMode(const uiint& mode)
{
    switch(mode){
        case(Utility::LoggerMode::Debug):case(Utility::LoggerMode::Error):
        case(Utility::LoggerMode::Warn): case(Utility::LoggerMode::Info):
            mpLogger->setMode(mode);
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid logger mode, CMW::LoggerMode");
            break;
    }
}
void CMW::LoggerDevice(const uiint& mode, const uiint& device)
{
    switch(mode){
        case(Utility::LoggerMode::Debug):case(Utility::LoggerMode::Error):
        case(Utility::LoggerMode::Warn): case(Utility::LoggerMode::Info):
            mpLogger->setProperty(mode, device);
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid logger mode, CMW::LoggerDevice");
            break;
    }
}

void CMW::LoggerInfoMssg(const uiint& mode, const char* message)
{
    string sMessage(message);
    switch(mode){
        case(Utility::LoggerMode::Debug):case(Utility::LoggerMode::Error):
        case(Utility::LoggerMode::Warn): case(Utility::LoggerMode::Info):
            mpLogger->Info(mode, sMessage);
            break;
        default:
            mpLogger->Info(Utility::LoggerMode::Error, "invalid logger mode, CMW::LoggerInfo");
            break;
    }
}
//--
// format= %d:iint, %f:double(fixed), %e:double(scientific), %s:const char*
//--
void CMW::LoggerInfo(const uiint& mode, const char* format, ... )
{
    vector<void*> param;
    vint    vnVal;
    vuint   vuVal;
    vdouble vdVal, veVal;
    vstring vsVal;
    
    if(mode==Utility::LoggerMode::Debug || mode==Utility::LoggerMode::Error || 
           mode==Utility::LoggerMode::Warn  || mode==Utility::LoggerMode::Info){
        
        uiint nLength = strlen(format);

        va_list list;
        va_start( list, format);
        for(uiint i=0; i < nLength; i++){
            if(format[i] == '%'){
                ++i;
                switch( format[i] ){
                case('d'):
                    { iint nVal = va_arg(list, iint);  vnVal.push_back(nVal);}
                    break;
                case('u'):
                    { uiint nVal = va_arg(list, uiint);  vuVal.push_back(nVal);}
                    break;
                case('f'):
                    { double dVal= va_arg( list, double ); vdVal.push_back(dVal);}
                    break;
                case('e'):
                    { double dVal= va_arg( list, double ); veVal.push_back(dVal);}
                    break;
                case('s'):
                    { string sVal= va_arg( list, const char* ); vsVal.push_back(sVal);}
                    break;
                default:
                    break;
                }
            }
        };
        va_end( list );

        uiint icase_d(0), icase_u(0),icase_f(0), icase_e(0), icase_s(0);
        for(uiint i=0; i < nLength; i++){
            if(format[i] == '%'){
                ++i;
                switch( format[i] ){
                case('d'):
                    vnVal[icase_d];
                    param.push_back(&vnVal[icase_d]);
                    ++icase_d;
                    break;
                case('u'):
                    vuVal[icase_u];
                    param.push_back(&vuVal[icase_u]);
                    ++icase_u;
                    break;
                case('f'):
                    vdVal[icase_f];
                    param.push_back(&vdVal[icase_f]);
                    ++icase_f;
                    break;
                case('e'):
                    veVal[icase_e];
                    param.push_back(&veVal[icase_e]);
                    ++icase_e;
                    break;
                case('s'):
                    vsVal[icase_s];
                    param.push_back(&vsVal[icase_s][0]);
                    //cout << "CMW::LoggerInfo  vsVal=" << vsVal[icase_s] << ", param=" << (char*)param[param.size()-1] << endl;
                    ++icase_s;
                    break;
                }
            }
        };

        mpLogger->Info(mode, format, param);//--------------- Logger出力

    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "invalid logger mode, CMW::LoggerInfo");
    }
}
uiint CMW::getErrorMode()
{
    return Utility::LoggerMode::Error;
}
uiint CMW::getWarnMode()
{
    return Utility::LoggerMode::Warn;
}
uiint CMW::getInfoMode()
{
    return Utility::LoggerMode::Info;
}
uiint CMW::getDebugMode()
{
    return Utility::LoggerMode::Debug;
}
uiint CMW::getDiskDevice()
{
    return Utility::LoggerDevice::Disk;
}
uiint CMW::getDisplayDevice()
{
    return Utility::LoggerDevice::Display;
}

//////--
////// extern "C"
//////--
////MPI_Datatype mpi_uiint()
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->MPI_UIINT();
////}
////MPI_Datatype mpi_iint()
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->MPI_IINT();
////}
////int get_rank()
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->GetRank();
////}
////int get_num_of_process()
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->GetNumOfProcess();
////}
////iint allreduce(void* sendbuf, void* recvbuf, iint buf_size, MPI_Datatype datatype, MPI_Op op, MPI_Comm commworld)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->AllReduce(sendbuf, recvbuf, buf_size, datatype, op, commworld);
////}
////iint barrier(MPI_Comm commworld)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Barrier(commworld);
////}
////iint mpi_abort(MPI_Comm commworld, int error)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Abort(commworld, error);
////}
////iint allgather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->AllGather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
////}
////iint gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm);
////}
////iint scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
////}
////iint recv(void* buf, iint count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Recv(buf, count, datatype, source, tag, comm, status);
////}
////iint send(void* buf, iint count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Send(buf, count, datatype, dest, tag, comm);
////}
////iint bcast(void* buf, iint cnt, MPI_Datatype type, int root, MPI_Comm comm)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->Bcast(buf, cnt, type, root, comm);
////}
////iint get_num_of_neib_pe(const uiint& imesh)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->GetNumOfNeibPE(imesh);
////}
////iint get_trans_rank(const uiint& imesh, const uiint& ipe)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    return pMW->GetTransRank(imesh, ipe);
////}
////
////void send_recv_r(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    pMW->Send_Recv_R(buf, num_of_node, dof_size, trans_rank);
////}
////void send_recv_i(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank )
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    pMW->Send_Recv_I(buf, num_of_node, dof_size, trans_rank);
////}
////void sendrecv_r(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    pMW->Sendrecv_R(buf, num_of_node, dof_size, trans_rank);
////}
////void sendrecv_i(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank)
////{
////    pmw::CMW* pMW= pmw::CMW::Instance();
////
////    pMW->Sendrecv_I(buf, num_of_node, dof_size, trans_rank);
////}



