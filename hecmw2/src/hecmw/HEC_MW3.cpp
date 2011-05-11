//
//  HEC_MW3.cpp
//
//			2009.04.20
//			2008.11.19
//			k.Takeda
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

#ifdef MSVC
#include "HEC_MW3.hxx"
#else
#include "HEC_MW3.h"
#endif

using namespace pmw;

#define ERROR 0
#define SUCCESS 1

// Constructor
//
CMW::CMW(void)
{
    // GMGModel
    mpGMGModel= CGMGModel::Instance();

    // SolutionType
    mnSolutionType = SolutionType::FEM;

    // Factory
    mpFactory = CMeshFactory::Instance();
    mpFactory->setGMGModel(mpGMGModel);
    mpFactory->setSolutionType(mnSolutionType);

    // FileIO
    mpFileIO = FileIO::CFileIO::Instance();
    mb_file = false;// ファイル読み込み後？

    // Factory => FileIO
    mpFileIO->setFactory(mpFactory);
    mpFileIO->setSolutionType(mnSolutionType);

    // MPI
    mpMPI = CHecMPI::Instance();

    mpLogger = Utility::CLogger::Instance();

    // N:形状関数,dN/dr:導関数
    mpShapeHexa= pmw::CShapeHexa::Instance();
    mpShapeHexaNic= pmw::CShapeHexaNic::Instance();
    mpShapeTetra= pmw::CShapeTetra::Instance();
    mpShapePrism= pmw::CShapePrism::Instance();
    mpShapeQuad= pmw::CShapeQuad::Instance();
    mpShapeTriangle= pmw::CShapeTriangle::Instance();
    mpShapeLine= pmw::CShapeLine::Instance();

    // 形状間数カタログ
    mpShapeCatalog= pmw::CShapeFunctionCatalog::Instance();

    // 形状関数 番号 <=> 辺番号 変換
    mpISTR2Edge= pmw::CISTR2Edge::Instance();
    mpEdge2ISTR= pmw::CEdge2ISTR::Instance();
    
}

// Destructor
//
CMW::~CMW(void)
{
    ;
}

//----
// Initialize (標準)：MPI, Logger の初期化
//----
uiint CMW::Initialize(int argc, char** argv)
{
    // GMG に一つのAssyModelを生成 : コースグリッド
    mpGMGModel->initAssyModel();
    mpFactory->setMGLevel(0);


    // Logger設定
    //mpLogger->setMode(Utility::LoggerMode::Info);
    mpLogger->setMode(Utility::LoggerMode::MWDebug);
    
    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn,    Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info,    Utility::LoggerDevice::Display);
    
    // MPI 初期化
    mpMPI->Initialize(argc, argv);
    
    // 表示
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");
    
    // ログファイル・イニシャライズ
    uiint rank= mpMPI->getRank();
    mpLogger->initializeLogFile(rank);

    // REVOCAP_Refiner
    rcapInitRefiner(0,0);
}

//
// nLength : ファイル名 長さ
// sName : ctrlファイルから取得したファイル名
// 
//
uiint CMW::BaseName_BCast(int& nLength, string& sName, int nType)
{
    // "Rank 0"
    if(mpMPI->getRank()==0){
        sName= mpFileIO->getFstr_FileName(nType);
        nLength= sName.length();
    }
    mpMPI->Barrier(MPI_COMM_WORLD);// 同期
    mpMPI->Bcast(&nLength, 1, MPI_INT, 0, MPI_COMM_WORLD);// 文字列長さ情報を共有

    char cBaseName[nLength];
    
    // "Rank 0"
    if(mpMPI->getRank()==0){
        for(uiint i=0; i < nLength; i++) cBaseName[i] = sName[i];
    }
    mpMPI->Barrier(MPI_COMM_WORLD);// 同期
    mpMPI->Bcast((void*)cBaseName, nLength, MPI_CHAR, 0, MPI_COMM_WORLD);// ファイル名情報を共有

    // "Rank > 0"
    if(mpMPI->getRank()!=0){
        mpFileIO->setFstr_FileName(cBaseName, nLength, nType);
    }
}

//----
// Initialize(標準) + ctrlファイルのパス設定 + ファイルベース名をMPI経由で配布(rank!=0)
//----
uiint CMW::Initialize_fstr(int argc, char** argv, string& ctrlname)
{
    CMW::Initialize(argc, argv);//Initialize(標準)

    mpFileIO->markingFstrStyle();//FrontISTRスタイルでリスタートファイル拡張子を付与
    
    // "Rank 0"
    if(mpMPI->getRank()==0){
        //// ctrlファイルのパス指定
        //mpFileIO->setCntPathName(cntpath.c_str());
        
        // hecmw_ctrl.dat読み込み
        bool bCntSuccess = mpFileIO->Read_fstr_CntFile(ctrlname);
        
        if(!bCntSuccess){
            mpLogger->Info(Utility::LoggerMode::Error, "MW::Initialize_fstr interruption");
            return ERROR;
        }
    }

    //// データ・ファイルのパス指定
    //mpFileIO->setDatPathName(datpath.c_str());
    
    // ----
    // プロセス０が取得したファイル名を各Rankで共有
    // ----
    int nLength;     //ファイル名の長さ, MPIなのでint
    string sBaseName;//ベース名(テンポラリー)

    for(int iType=0; iType < FileIO::FileTypeMW2::Limit; iType++){
        CMW::BaseName_BCast(nLength, sBaseName, iType);
    };

    // --
    // rank番号を付加するファイル : HEC_MW3の管理
    // --
    // ファイルパス
    //    msOutputFileName = mpFileIO->getDatPathName();//Data Check File(MW3)
    //    msInputFileName  = mpFileIO->getDatPathName();
    //    msResFileName    = mpFileIO->getDatPathName();
    //    msRltFileName    = mpFileIO->getDatPathName();
    //    msPartOutFileName= mpFileIO->getDatPathName();

    // ファイルパス + ベース名 + "."
    msOutputFileName  += mpFileIO->getFstr_MeshFileName()    + ".";//Data Check File(MW3)
    msInputFileName   += mpFileIO->getFstr_MeshFileName()    + ".";//fstr
    msResFileName     += mpFileIO->getFstr_RestartFileName() + ".";//fstr
    msRltFileName     += mpFileIO->getFstr_ResultFileName()  + ".";//fstr
    msPartOutFileName += mpFileIO->getFstr_PartFileName_OUT()+ ".";//fstr

    // rank 取得
    uiint rank= mpMPI->getRank();

    stringstream ss;
    ss << rank;
    // rankを追加
    msOutputFileName  += ss.str();//Data Check File(MW3)
    msInputFileName   += ss.str();
    msResFileName     += ss.str();
    msRltFileName     += ss.str();
    msPartOutFileName += ss.str();

    // 拡張子を追加
    msOutputFileName += ".out";//Data Check File(MW3)
    
    return SUCCESS;
}
// --
// ファイル名(fstr) : rank番号を付加するファイルは、HEC_MW3の管理
// --
string CMW::getFstr_FileName_Mesh()
{
    return msInputFileName;
}
string CMW::getFstr_FileName_Control()
{
    //string sName= mpFileIO->getDatPathName() + mpFileIO->getFstr_ControlFileName();

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
    //string sName= mpFileIO->getDatPathName() + mpFileIO->getFstr_PartFileName_IN();

    return mpFileIO->getFstr_PartFileName_IN();
}
string CMW::getFstr_FileName_PartOut()
{
    return msPartOutFileName;
}
string CMW::getFstr_FileName_VisMesh()
{
    //string sName= mpFileIO->getDatPathName() + mpFileIO->getFstr_VisFileName_Mesh();

    return mpFileIO->getFstr_VisFileName_Mesh();
}
string CMW::getFstr_FileName_VisIn()
{
    //string sName= mpFileIO->getDatPathName() + mpFileIO->getFstr_VisFileName_IN();

    return mpFileIO->getFstr_VisFileName_IN();
}
string CMW::getFstr_FileName_VisOut()
{
    //string sName= mpFileIO->getDatPathName() + mpFileIO->getFstr_VisFileName_OUT();
    
    return mpFileIO->getFstr_VisFileName_OUT();
}

// --
// 横断幕表示
// --
void CMW::Banner_Display()
{
    mpLogger->InfoDisplay();
}


//----
// REVOCAP_Refiner
//----
void CMW::set_RevocapMesh(const uiint& nRefine)//REVOCAP_Refinerへのメッシュの登録
{
    size_t nodeCount;
    float64_t *coords;
    int32_t *globalIDs;
    int32_t *localIDs;

    nodeCount = mpMesh->getNumOfNode();
    coords = new float64_t[nodeCount*3];
    globalIDs = new int32_t [nodeCount];
    localIDs  = new int32_t [nodeCount];

    for(uiint i=0; i < nodeCount; i++){
        CNode *pNode= mpMesh->getNodeIX(i);
        coords[i*3]    = pNode->getX();
        coords[i*3 + 1]= pNode->getY();
        coords[i*3 + 2]= pNode->getZ();

        globalIDs[i]= pNode->getID();// ID
        localIDs[i] = i;             // index
    };

    rcapSetNode64(nodeCount, coords, globalIDs, localIDs);//Nodeのセット

    size_t elementCount;
    size_t refineElemCount;
    int32_t* elemNodes;
    int32_t* refineNodes;
    
    elementCount = mpMesh->getNumOfElement();
    elemNodes = new int32_t[elementCount*4];

    for(uiint i=0; i < elementCount; i++){
        CElement *pElem = mpMesh->getElementIX(i);
        
        for(uiint ii=0; ii < 4; ii++){
            CNode *pNode = pElem->getNode(ii);
            elemNodes[i*4 + ii] = pNode->getID();
        };
    };

    for(uiint i=0; i < nRefine; i++) 
        refineElemCount = rcapRefineElement( elementCount, RCAP_TETRAHEDRON, elemNodes, refineNodes);
}


// --
// MPI, Logger 後処理
// --
uiint CMW::Finalize()
{
    // MPI
    mpMPI->Finalize();

    // Logger
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();//ログファイル-クローズ

    // REVOCAP_Refiner
    rcapTermRefiner();

    return 1;
}

// ---------------------
// FileIO内からFactoryコール.
// 1.Factory内でAssyModel(マルチグリッド階層数分の生成)
// 2.Mesh(Level=0)を生成
// ---------------------
// --
// ファイル読み込み : MW3標準
// --
uiint CMW::FileRead(string& basename, bool bBinary)//(一般)
{
    // データファイルパス && ベースネームを保存
    //    mpFileIO->setDatPathName(datpath.c_str());
    mpFileIO->setBaseName(basename);

    // データファイルパス
    //    msInputFileName = datpath;
    //    msOutputFileName= datpath;//Data Check File
    //    msResFileName   = datpath;
    
    // ファイルパス + ファイルベース名 + ピリオド
    msInputFileName += basename + ".";
    msOutputFileName+= basename + ".";//Data Check File
    msResFileName   += basename + ".";
    msRltFileName   += basename + ".";
    
    // rank 取得
    uiint rank= mpMPI->getRank();

    stringstream ss;
    ss << rank;
    // ファイル名に、rankを追加
    msInputFileName  += ss.str();
    msOutputFileName += ss.str();//Data Check File
    msResFileName    += ss.str();
    msRltFileName    += ss.str();

    // ファイル名に、拡張子を追加
    msInputFileName  += ".msh";
    msOutputFileName += ".out";//Data Check File
    ////msResFileName+= ".res";//File..Chunkに管理を移行
    ////msRltFileName+= ".rlt";//出力時に指定 .inp とか..

    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");

    mpFileIO->ReadFile(msInputFileName, bBinary);
    mb_file = true;
    
    return 1;
}
// --
// ファイル読み込み : FrontISTR仕様
// --
uiint CMW::FileRead_fstr(bool bBinary)
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead_fstr");

    mpFileIO->ReadFile(msInputFileName, bBinary);
    mb_file = true;

    return 1;
}

uiint CMW::FileDebugWrite()//Data Check File output
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 DataCheck Write");

    uiint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile_Debug(msOutputFileName, nLevel+1);//階層数+1 <= 出力数

    return 1;
}

// --
// リザルト・ファイル
// --
void CMW::PrintRlt_Start(const uiint& nStep, bool bBinary)
{
    mpFileIO->PrintResult_Start(nStep, msRltFileName, bBinary);
}
// リザルト出力(可変長引数:ポインター)
// 
void CMW::PrintRlt_P(const uiint& width, const char* format, ...)
{
    vector<void*> param;
    iint *npVal; double *dpVal; string sVal; char *cVal;
    
    va_list list;

    va_start( list, format);
    uiint nLength = strlen(format);

    for(uiint i=0; i < nLength; i++){
        if(format[i] == '%'){
            ++i;
            switch( format[i] ){
            case('d'):
                npVal= va_arg( list, iint* );
                param.push_back(npVal);
                break;
            case('f'):
                dpVal= va_arg( list, double* );
                param.push_back(dpVal);
                break;
            case('e'):
                dpVal= va_arg( list, double* );
                param.push_back(dpVal);
                break;
            case('s'):
                sVal= va_arg( list, const char* );
                cVal = new char[sVal.length()+1];

                strcpy(cVal, sVal.c_str());
                param.push_back(cVal);
                break;
            default:
                break;
            }
        }
    };
    va_end( list );

    char fmt[nLength+1];
    strcpy(fmt, format);

    mpFileIO->PrintResult(width, fmt, param);

}
// リザルト出力(可変長引数:値)
//
void CMW::PrintRlt_R(const uiint& width, const char* format, ...)
{
    vector<void*> param;
    iint nVal; double dVal;  string sVal; char *cVal;
    
    va_list list;

    va_start( list, format);
    uiint nLength = strlen(format);
    
    for(uiint i=0; i < nLength; i++){
        if(format[i] == '%'){
            ++i;
            switch( format[i] ){
            case('d'):
                nVal= va_arg( list, iint );
                param.push_back(&nVal);
                break;
            case('f'):
                dVal= va_arg( list, double );
                param.push_back(&dVal);
                break;
            case('e'):
                dVal= va_arg( list, double );
                param.push_back(&dVal);
                break;
            case('s'):
                sVal= va_arg( list, const char* );
                cVal = new char[sVal.length()+1];

                strcpy(cVal, sVal.c_str());
                param.push_back(cVal);
                break;
            default:
                break;
            }
        }
    };
    va_end( list );

    char fmt[nLength+1];
    strcpy(fmt, format);

    mpFileIO->PrintResult(width, fmt, param);

}
void CMW::PrintRlt_End()
{
    mpFileIO->PrintResult_End();
}
//
// *.inp 出力 : 基礎変数
//
void CMW::PrintMicroAVS_Basis()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    stringstream ss;
    string sFileName;
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        ss.clear(); ss.str("");
        ss << iMesh;
        sFileName = msRltFileName + "." + ss.str() + ".inp";//BaseName + rank + part + 拡張子
        
        mpFileIO->WriteAVS_Basis(sFileName, iMesh, nMaxLevel);//nLevel:出力する階層
    };
}
//出力ラベルの登録
void CMW::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    mpFileIO->recAVS_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
//出力変数の登録
void CMW::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    mpFileIO->recAVS_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
//登録変数の出力
void CMW::PrintMicroAVS_FEM()
{
    uiint nMaxLevel= mpFactory->getMGLevel();
    uiint iMesh, nNumOfMesh=mpAssy->getNumOfMesh();
    stringstream ss;
    string sFileName;
    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        ss.clear(); ss.str("");
        ss << iMesh;
        sFileName = msRltFileName + "." + ss.str() + ".inp";//BaseName + rank + part + 拡張子
        
        mpFileIO->WriteAVS_FEM(sFileName, iMesh, nMaxLevel);
    };
}

//
// Res(リスタート)出力
//
uiint CMW::FileWriteRes(const uiint& nStep, bool bBinary)
{
    mpLogger->Info(Utility::LoggerMode::Info, " HEC_MW3 RestartFile Write");

    mpFileIO->WriteResFile(nStep, msResFileName, bBinary);

    return 1;
}
//--
// Restart data のセット
//--
// Resファイル入力
// 1. 方程式生成 => FileReaderAlgebra
// 2. 値のセット => FileReaderRes
//--
uiint CMW::SetRestart(const uiint& nStep, bool bBinary)
{
    mvAlgebraEquation.clear();
    // ---
    // Algebraブロック
    // ---
    if(mpFileIO->ReadAlgebraBlock(nStep, msResFileName, bBinary) ){

        // 線形代数方程式の個数、各方程式のDOF
        uiint nNumOfAlgebra = mpFileIO->getNumOfEquation();
        mvAlgebraEquation.resize(nNumOfAlgebra);

        for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++){
            mvAlgebraEquation[ieq]= mpFileIO->getEquationDOF(ieq);
        };

        if(nNumOfAlgebra > 0)
            // 線形代数方程式の確保(CMW::GeneLinearAlgebra)
            GeneLinearAlgebra(nNumOfAlgebra, &mvAlgebraEquation[0]);
    }

    // ---
    // Resブロック
    // ---
    mpFileIO->ReadResBlock(nStep, msResFileName, bBinary);
    
    
    return mvAlgebraEquation.size();//線形方程式確保数
}


//
// 全Levelに方程式を生成(AssyMatrixのコースグリッドも、ここで設定される)
//
void CMW::GeneLinearAlgebra(const uiint& nNumOfAlgebra, uiint* vDOF)
{
    vuint vecDOF;
    vecDOF.reserve(nNumOfAlgebra);
    for(uiint idof=0; idof < nNumOfAlgebra; idof++){
        vecDOF.push_back(vDOF[idof]);
    };

    CAssyModel *pAssyModel, *pCGridAssyModel;
    uiint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
        pAssyModel = mpGMGModel->getAssyModel(iLevel);

        if(iLevel > 0){
            pCGridAssyModel = mpGMGModel->getAssyModel(iLevel-1);  //コースグリッド AssyModel
            pAssyModel->GeneLinearAlgebra(vecDOF, pCGridAssyModel);//線形方程式の生成とコースグリッドのセット.
        }else{
            pAssyModel->GeneLinearAlgebra(vecDOF, NULL);
        }
    };
}
//
// iequ番目の方程式をロード{ Levelは選択済み == AssyModelは選択済み }
//
void CMW::SelectAlgebra(const uiint& iequ)
{
    // AssyMatrixのコースグリッドはGeneLinearAlgebraで既に設定済み
    //
    mpAssyMatrix = mpAssy->getAssyMatrix(iequ);
    mpRHSAssyVector = mpAssy->getRHSAssyVector(iequ);
    mpSolAssyVector = mpAssy->getSolutionAssyVector(iequ);
}


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
// inode,jnode:節点インデックス
uiint CMW::Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double* NodalMatrix)
{
    return mpAssyMatrix->Matrix_Add_Nodal(iMesh, inode, jnode, NodalMatrix);
}

// Matrix 0 clear
//
void CMW::Matrix_Clear(const uiint& iMesh)
{
    mpAssyMatrix->Matrix_Clear(iMesh);
}
void CMW::Vector_Clear(const uiint& iMesh)
{
    mpRHSAssyVector->Vector_Clear(iMesh);
    mpSolAssyVector->Vector_Clear(iMesh);
}

void CMW::Sample_Set_BC(uiint iMesh)
{
#ifdef ADVANCESOFT_DEBUG
	printf(" enter CMW::Sample_Set_BC \n");
#endif

	double X, Y, Z, value0, value1, value2;
	uiint iDof0, iDof1, iDof2;
	uiint iNodeMax = getNodeSize( iMesh );

        //cout << " -- CMW::Sample_Set_BC -- " << endl;
        
	for( uiint iNode = 0; iNode < iNodeMax; iNode++){
		X = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getX();
		Y = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getY();
		Z = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getZ();

                ////debug
                //cout << "Z = " << Z << endl;

                // x=1.0 y=1.0 z=4.0 の一点荷重に変更'11.01.14 {以前、x=1.0 z=4.0 の辺の節点全て}
		if(abs( Z - 4.0 ) < 1.0e-5 && abs( X - 1.0 ) < 1.0e-5 && abs( Y - 1.0) < 1.0e-5 ) {

                    //cout << " Z=4.0  iNode:" << iNode << endl;

                    value0 = 1.0e6;
                    iDof0 = 0;
                    Set_BC_RHS( iMesh, iNode, iDof0, value0);
		};
		if( (abs( Z ) < 1.0e-5) || (abs( Z - 8.0 ) < 1.0e-5) ) {

                    //cout << " Z=0.0  iNode:" << iNode << endl;
                    
                    value1 = 1.0e15;
                    value2 = 0.0;
                    iDof0 = 0; iDof1 = 1; iDof2 = 2;
                    Set_BC_Mat_RHS( iMesh, iNode, iDof0, value1, value2);
                    Set_BC_Mat_RHS( iMesh, iNode, iDof1, value1, value2);
                    Set_BC_Mat_RHS( iMesh, iNode, iDof2, value1, value2);
		}
	};
#ifdef ADVANCESOFT_DEBUG
 	printf(" enter CMW::Sample_Set_BC \n");
#endif
}

//
// 右辺ベクトルへ境界値をセット
//
uiint CMW::Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
#ifdef ADVANCESOFT_DEBUG
    printf("enter CMW::Set_BC (RHS) %d %d %d %e \n", iMesh, iNode, nDOF, value);
#endif
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();

    uiint iNode = pBucket->getIndexNode(iNodeID);//NodeIDのインデックス番号

    if(mpRHSAssyVector){
        mpRHSAssyVector->setValue(iMesh, iNode, nDOF, value);
        return 1;
    }else{
        return 0;
    }

#ifdef ADVANCESOFT_DEBUG
    printf("exit CMW::Set_BC (RHS) \n");
#endif
	
}
//
// 右辺ベクトルへの境界値の加算
//
uiint CMW::Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();

    uiint iNode = pBucket->getIndexNode(iNodeID);//NodeIDのインデックス番号

    if(mpRHSAssyVector){
        mpRHSAssyVector->addValue(iMesh, iNode, nDOF, value);
        return 1;
    }else{
        return 0;
    }
}


//
// 行列を"1","0"で払って右辺を振替 ： 対角項=diagVal, 非対角項=0、右辺ベクトル
//
uiint CMW::Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diagValue, double& solValue)
{
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();
    
    uiint iNode = pBucket->getIndexNode(iNodeID);//NodeIDのインデックス番号

    mpRHSAssyVector->setValue(iMesh, iNode, nDOF, solValue);
    
    mpAssyMatrix->setZero_NonDiag(iMesh, iNode, nDOF, mpRHSAssyVector, solValue);    // 非対角項=0 処理
    mpAssyMatrix->setValue(iMesh, iNode, nDOF, diagValue, mpRHSAssyVector, solValue);// 対角項
    
    return 1;
}

//
// ペナルティ法：行列の対角項、右辺ベクトルへ境界値をセット
//
uiint CMW::Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diagValue, double& rhsValue)
{
#ifdef ADVANCESOFT_DEBUG
    printf("enter CMW::Set_BC (Mat_RHS) %d %d %d %e %e \n", iMesh, iNode, nDOF, diagValue, rhsValue);
#endif
    CMesh *pMesh = mpAssy->getMesh(iMesh);
    CIndexBucket *pBucket = pMesh->getBucket();

    uiint iNode = pBucket->getIndexNode(iNodeID);//NodeIDのインデックス番号

    mpAssyMatrix->setValue_D(iMesh, iNode, nDOF, diagValue);
    mpRHSAssyVector->setValue(iMesh, iNode, nDOF, rhsValue);

#ifdef ADVANCESOFT_DEBUG
    printf("exit CMW::Set_BC (Mat_RHS) \n");
#endif

	return 1;
}

uiint CMW::Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition)
{
#ifdef ADVANCESOFT_DEBUG
  printf(" enter CMW::Solve %d %e \n", iter_max, tolerance);
#endif

////    // 行列をダンプ
////    cout << "-- Assy マトリックス --" << endl;
////    mpAssyMatrix->dump();
////    // ----
////    // 列ベクトルをダンプ
////    cout << "--   Sol ベクトル   --" << endl;
////    mpSolAssyVector->dump();
////    // ----
////    // 列ベクトルをダンプ
////    cout << "--   RHS ベクトル   --" << endl;
////    mpRHSAssyVector->dump();

  bool flag_iter_log = false;
  bool flag_time_log = false;
  uiint  success;
  
    switch( method ){
        case( 1 ):
	   {CSolverCG *solver = new CSolverCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   success=solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 2 ):
	   {CSolverBiCGSTAB *solver = new CSolverBiCGSTAB(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   success=solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 3 ):
           {CSolverGPBiCG *solver = new CSolverGPBiCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   success=solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 4 ):
	   {CSolverGMRES *solver = new CSolverGMRES(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   success=solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        default:
           break;
    }
	
    if( success == 0 ) {
      std::cout << " Fails in solver! " << std::endl;
	  return 0;
    }
    
#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMW::Solve \n");
#endif
    return 1;
}


//--
// 解ベクトルをbufへコピー  (SelectされているLevel && Selectされている方程式番号)
//--
void CMW::GetSolution_Vector(double* buf, const uiint& imesh)
{
    CVector *pVec = mpSolAssyVector->getVector(imesh);

    CMesh *pMesh = mpAssy->getMesh(imesh);
    uiint nNumOfNode = pMesh->getNumOfNode();
    uiint nDOF = mpSolAssyVector->getDOF();

    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            buf[ inode*nDOF + idof ] = pVec->getValue(inode, idof);
        };
    };
}
void CMW::GetSolution_AssyVector(double* buf)
{
    uiint nDOF = mpSolAssyVector->getDOF();
    uiint nNumOfMesh = mpAssy->getNumOfMesh();

    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){
        CMesh *pMesh = mpAssy->getMesh(iMesh);

        uiint nNumOfNode = pMesh->getNodeSize();
        uiint nMeshSize = iMesh * nNumOfNode;

        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                buf[nMeshSize + iNode*nDOF + idof] = mpSolAssyVector->getValue(iMesh, iNode, idof);
            };
        };
    };
}
//--
// 右辺ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
//--
void CMW::GetRHS_Vector(double* buf, const uiint& imesh)
{
    CVector *pVec = mpRHSAssyVector->getVector(imesh);

    CMesh *pMesh = mpAssy->getMesh(imesh);
    uiint nNumOfNode = pMesh->getNumOfNode();
    uiint nDOF = mpRHSAssyVector->getDOF();

    for(uiint inode=0; inode < nNumOfNode; inode++){
        for(uiint idof=0; idof < nDOF; idof++){
            buf[ inode*nDOF + idof ] = pVec->getValue(inode, idof);
        };
    };
}
void CMW::GetRHS_AssyVector(double* buf)
{
    uiint nDOF = mpRHSAssyVector->getDOF();
    uiint nNumOfMesh = mpAssy->getNumOfMesh();

    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){
        CMesh *pMesh = mpAssy->getMesh(iMesh);

        uiint nNumOfNode = pMesh->getNodeSize();
        uiint nMeshSize = iMesh * nNumOfNode;

        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                buf[nMeshSize + iNode*nDOF + idof] = mpRHSAssyVector->getValue(iMesh, iNode, idof);
            };
        };
    };
}
//--
// ベクトル値を取得
//--
double& CMW::GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof)
{
    return mpSolAssyVector->getValue(imesh, inode, idof);
}
double& CMW::GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof)
{
    return mpRHSAssyVector->getValue(imesh, inode, idof);
}
//--
// ベクトルDOFを取得
//--
uiint& CMW::GetSolutionVector_DOF()
{
    return mpSolAssyVector->getDOF();
}
uiint& CMW::GetRHSVector_DOF()
{
    return mpRHSAssyVector->getDOF();
}

//----
// AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
//----
void CMW::multVector(double* vX, double* vB)
{
    CAssyVector vecX(mpSolAssyVector);
    CAssyVector vecB(mpRHSAssyVector);

    vecX.setZero();

    // 入力値：Xベクトル値をAssyVectorに代入
    //
    uiint nDOF = mpSolAssyVector->getDOF();
    uiint nNumOfMesh = mpAssy->getNumOfMesh();

    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){

        CMesh *pMesh = mpAssy->getMesh(iMesh);
        uiint nNumOfNode = pMesh->getNodeSize();
        uiint nMeshSize = iMesh * nNumOfNode;

        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                double val = vX[nMeshSize + iNode*nDOF + idof];
                vecX.setValue(iMesh, iNode, idof, val);
            };
        };
    };

    vecB.setZero();
    //
    // A * vecX = vecB
    //
    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){
        mpAssyMatrix->multVector(iMesh, &vecX, &vecB);
    };

    // 出力値：BベクトルにAssyVectorの値を代入
    //
    for(uiint iMesh=0; iMesh < nNumOfMesh; iMesh++){

        CMesh *pMesh = mpAssy->getMesh(iMesh);
        uiint nNumOfNode = pMesh->getNodeSize();
        uiint nMeshSize = iMesh * nNumOfNode;

        for(uiint iNode=0; iNode < nNumOfNode; iNode++){
            for(uiint idof=0; idof < nDOF; idof++){
                double val = vecB.getValue(iMesh, iNode, idof);
                vB[nMeshSize + iNode*nDOF + idof] = val;
            };
        };
    };
}



//
// 1.マルチグリッドのための,Refineメッシュ生成
// 2.MPCのため,接合Meshのマスター面にスレーブ点をセット
//
uiint CMW::Refine(const uiint& nNumOfRefine)
{
    uiint ilevel,mgLevel;
    CAssyModel *pAssy;
    uiint icon,numOfConMesh;
    CContactMesh *pConMesh;

    // pre 1. AssyModel  generate
    mpFactory->GeneAssyModel(nNumOfRefine+1);//"Refine段数+1" => マルチグリッド階層数に合わせてAssyModelを生成
    mpFactory->setMGLevel(nNumOfRefine);     // Refine段数のセット

    // pre 2. コースグリッドのBNode_Valueの領域確保
    pAssy = mpGMGModel->getAssyModel(0);
    uiint iMesh, nNumOfMesh= pAssy->getNumOfMesh();
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

    // pre 3. ContactMeshのドンガラを全Levelに生成しておく
    //        # Level=0だけファイル読み込み時に生成ずみ
    uiint contactID, myRank, transRank, nProp;
    CAssyModel *pPrevAssy= mpGMGModel->getAssyModel(0);//Level=0から情報取得
    uiint nNumOfCont= pPrevAssy->getNumOfContactMesh();

    for(ilevel=1; ilevel< nNumOfRefine+1; ilevel++){
        
        pAssy= mpGMGModel->getAssyModel(ilevel);
        
        for(uiint icont=0; icont < nNumOfCont; icont++){
            CContactMesh *pPrevContactMesh= pPrevAssy->getContactMesh(icont);
            CContactMesh *pContactMesh = new CContactMesh;//-------------- 接合メッシュの生成

            contactID = pPrevContactMesh->getID();
            pContactMesh->setID(contactID);

            pContactMesh->setLevel(ilevel);

            myRank= pPrevContactMesh->getRank();
            transRank= pPrevContactMesh->getTransmitRank();
            nProp= pPrevContactMesh->getProp();

            pContactMesh->setRank(myRank);
            pContactMesh->setTransmitRank(transRank);
            pContactMesh->setProp(nProp);

            pAssy->addContactMesh(pContactMesh, contactID);//------------ 各LevelのAssyModelに接合メッシュをセット
        };
    };

    
    if(mb_file){// ファイル読み込み後
        
        //cout << "MW::Refine, A" << endl;

        mpFactory->refineMesh();       //MultiGridデータの生成(通信Meshも含む)
        mpFactory->refineContactMesh();//接合Mesh(MPCメッシュ)のMultiGridデータ生成

        //cout << "MW::Refine, B" << endl;
        
        mgLevel= mpFactory->getMGLevel();
        // MPCのマスタースレーブ構成のセット
        // 
        for(ilevel=0; ilevel < mgLevel+1; ilevel++){
            pAssy= mpGMGModel->getAssyModel(ilevel);
            
            numOfConMesh= pAssy->getNumOfContactMesh();

            for(icon=0; icon < numOfConMesh; icon++){
                pConMesh= pAssy->getContactMesh(icon);

                pConMesh->setupSPointOnMFace();//マスター面にスレーブ点をセット
                pConMesh->setupMPC_Coef();     //スレーブ点のマスター面の頂点ごとのCoefを計算
            };
        };

        //cout << "MW::Refine, C" << endl;

        mpFactory->refineCommMesh2();//要素分割型(節点共有型)の通信Mesh(CommMesh2)のRefine
        mpFactory->refineBoundary(); //境界条件の階層化

        //cout << "MW::Refine, D" << endl;
        
        return 1;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error," Not read file @MWMain::Refine()");
        return 0;
    }
    return 1;
}

//
// メモリーの解放(Node::AggElemIDの解放):Mesh操作によるデータ構築の最後に実行
//
void CMW::FinalizeRefine()
{
    // Refineに用いたvectorの解放
    //
    CAssyModel *pAssy;
    uiint ilevel, mgLevel = mpFactory->getMGLevel();
    
    for(ilevel=0; ilevel < mgLevel+1; ilevel++){
        pAssy = mpGMGModel->getAssyModel(ilevel);

        CMesh *pMesh;
        uiint imesh, nNumOfMesh = pAssy->getNumOfMesh();

        // Mesh      メッシュ本体
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            pMesh = pAssy->getMesh(imesh);
            
            pMesh->deleteProgData();// Elementの各Edge,Face の Node*,Element* 配列 を解放
            pMesh->deleteAggregate_on_Node();// Node::AggElemIDの解放

            // CommMesh2 通信メッシュ
            CCommMesh2 *pCommMesh;
            uiint icom, nNumOfCommMesh = pMesh->getCommMesh2Size();
            for(icom=0; icom < nNumOfCommMesh; icom++){
                pCommMesh = pMesh->getCommMesh2IX(icom);
                
                pCommMesh->deleteProgData();
            };

            // BoundaryVolumeMesh 境界メッシュ(体積)
            CBoundaryVolumeMesh *pBVolMesh;
            uiint ibvol, nNumOfBVolMesh = pMesh->getNumOfBoundaryVolumeMesh();
            for(ibvol=0; ibvol < nNumOfBVolMesh; ibvol++){
                pBVolMesh = pMesh->getBndVolumeMeshIX(ibvol);

                pBVolMesh->deleteProgData();
            };

            // BoundaryFaceMesh 境界メッシュ(面)
            CBoundaryFaceMesh *pBFaceMesh;
            uiint ibface, nNumOfBFaceMesh = pMesh->getNumOfBoundaryFaceMesh();
            for(ibface=0; ibface < nNumOfBFaceMesh; ibface++){
                pBFaceMesh = pMesh->getBndFaceMeshIX(ibface);

                pBFaceMesh->deleteProgData();
            };

            // BoundaryEdgeMesh 境界メッシュ(辺)
            CBoundaryEdgeMesh *pBEdgeMesh;
            uiint ibedge, nNumOfBEdgeMesh = pMesh->getNumOfBoundaryEdgeMesh();
            for(ibedge=0; ibedge < nNumOfBEdgeMesh; ibedge++){
                pBEdgeMesh = pMesh->getBndEdgeMeshIX(ibedge);

                pBEdgeMesh->deleteProgData();
            };
        };

        // ContactMesh 接合面メッシュ
        CContactMesh *pConMesh;
        uiint icmesh, nNumOfConMesh = pAssy->getNumOfContactMesh();

        for(icmesh=0; icmesh < nNumOfConMesh; icmesh++){
            pConMesh = pAssy->getContactMesh(icmesh);

            pConMesh->deleteProgData();// SkinFaceの辺のConNode*配列を解放
        };

    };// iLevel
    
}




//-------------------------------
// Fortran用の作業用オブジェクト選択
//-------------------------------

// Assemble Modelの個数==階層数(mMGLevel+1)
uiint CMW::GetNumOfAssembleModel()
{
    return mpGMGModel->getNumOfLevel();
}

// Assemble Modelの選択
void CMW::SelectAssembleModel(const uiint& mgLevel)
{
    mpAssy= mpGMGModel->getAssyModel(mgLevel);
}

// Meshパーツの個数
uiint CMW::GetNumOfMeshPart()
{
    return mpAssy->getNumOfMesh();
}

// Meshの選択 : IDによる選択
void CMW::SelectMeshPart_ID(const uiint& mesh_id)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh_ID(mesh_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
// Meshの選択 : 配列Index番号による選択
void CMW::SelectMeshPart_IX(const uiint& index)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
// Elementの選択 : IDによる選択
void CMW::SelectElement_ID(const uiint& elem_id)
{
    if(mpMesh){
        mpElement= mpMesh->getElement(elem_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
// Elementの選択 : 配列Index番号による選択
void CMW::SelectElement_IX(const uiint& index)
{
    if(mpMesh){
        mpElement= mpMesh->getElementIX(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
// Elementのタイプ取得
uiint CMW::GetElementType()
{
    return mpElement->getType();
}

// Elementの頂点数
uiint CMW::GetNumOfElementVert()
{
    return mpElement->getNumOfNode();
}

// Elementの頂点のノードID
void CMW::GetElementVertNodeID(iint* vNodeID)
{
    uiint numOfNode;
    numOfNode= mpElement->getNumOfNode();

    CNode* pNode;
    uiint ivert;
    for(ivert=0; ivert< numOfNode; ivert++){
        pNode= mpElement->getNode(ivert);

        vNodeID[ivert]= pNode->getID();
    };
}

// Elementの辺の数
uiint CMW::GetNumOfElementEdge()
{
    return mpElement->getNumOfEdge();
}

// Elementの辺のノードID
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

// --
// ID  &&  index
// --
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

// --
// 節点コネクティビティ
// --
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
                mv_ItemL.push_back(k_index);//Itemは節点配列インデックス 2011.04.21

            }else if(node_id < k_node_id){
                mv_ItemU.push_back(k_index);//Itemは節点配列インデックス 2011.04.21
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
//node_id : 上三角のNodeID数, 下三角のNodeID数
void CMW::getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL)
{
    nNumOfItemU = mv_ItemU.size();
    nNumOfItemL = mv_ItemL.size();
}
//node_id : 上三角のNodeID配列, 下三角のNodeID配列
void CMW::getNodeConnectFEM_Item(uiint itemU[], uiint itemL[])
{
    for(uiint i=0; i < mv_ItemL.size(); i++){
        itemL[i] = mv_ItemL[i];
    };
    for(uiint i=0; i < mv_ItemU.size(); i++){
        itemU[i] = mv_ItemU[i];
    };
}
void CMW::getNodeConnectFEM_Item_F(iint itemU[], iint itemL[])
{
    //FortranAPI用 (引数タイプ : signed int )
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
//--
// 節点周囲 要素
//--
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

// --
// ノード
// --
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
////void CMW::SetNodeValue(const uiint& node_id, double value[])
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    uiint numOfDOF = pNode->numOfTotalParam();
////
////    uiint idof;
////    if(NodeType::Scalar==pNode->getType()){
////        for(idof=0; idof < numOfDOF; idof++){
////            pNode->setScalar(value[idof], idof);
////        };
////    }
////    if(NodeType::Vector==pNode->getType()){
////        for(idof=0; idof < numOfDOF; idof++){
////            pNode->setVector(value[idof], idof);
////        };
////    }
////    if(NodeType::ScalarVector==pNode->getType()){
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetNodeValue(const uint&, double[])");
////    }
////}
////void CMW::SetNodeValue(const uiint& node_id, const uiint& idof, const double& value)
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::Scalar==pNode->getType()){
////        pNode->setScalar(value, idof);
////    }
////    if(NodeType::Vector==pNode->getType()){
////        pNode->setVector(value, idof);
////    }
////    if(NodeType::ScalarVector==pNode->getType()){
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetNodeValue(const uint&, const uint&, double&)");
////    }
////}
////
////void CMW::GetNodeValue(const uiint& node_id, double value[])
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::Scalar==pNode->getType()){
////        uiint idof, numOfScalar=pNode->numOfScalarParam();
////        for(idof=0; idof < numOfScalar; idof++)
////            value[idof] = pNode->getScalar(idof);
////    }
////    if(NodeType::Vector==pNode->getType()){
////        uiint idof, numOfVector=pNode->numOfVectorParam();
////        for(idof=0; idof < numOfVector; idof++)
////            value[idof] = pNode->getVector(idof);
////    }
////    if(NodeType::ScalarVector==pNode->getType()){
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetNodeValue(uint&, double[])");
////    }
////}
////double& CMW::GetNodeValue(const uiint& node_id, const uiint& idof)
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::Scalar==pNode->getType()){
////        return pNode->getScalar(idof);
////    }
////    if(NodeType::Vector==pNode->getType()){
////        return pNode->getVector(idof);
////    }
////    if(NodeType::ScalarVector==pNode->getType()){
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetNodeValue(uint&, uint&)");
////        return pLogger->getDDummyValue();
////    }
////}
//
// ScalarVector Node
//
////void CMW::SetSVNodeValue(const uiint& node_id, double v_value[], double s_value[])
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::ScalarVector==pNode->getType()){
////        uiint idof;
////
////        uiint numOfScalar=pNode->numOfScalarParam();
////        for(idof=0; idof < numOfScalar; idof++) pNode->setScalar(s_value[idof], idof);
////
////        uiint numOfVector=pNode->numOfVectorParam();
////        for(idof=0; idof < numOfVector; idof++) pNode->setVector(v_value[idof], idof);
////
////    }else{
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetSVNodeValue(uint&, uint&, double&, uint&, double&)");
////    }
////}
////void CMW::SetSVNodeValue(const uiint& node_id,
////        const uiint& v_dof, const double& v_value, const uiint& s_dof, const double& s_value)
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::ScalarVector==pNode->getType()){
////        pNode->setScalar(s_value, s_dof);
////        pNode->setVector(v_value, v_dof);
////    }else{
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetSVNodeValue(uint&, uint&, double&, uint&, double&)");
////    }
////}
////void CMW::GetSVNodeValue(const uiint& node_id, double v_value[], double s_value[])
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::ScalarVector==pNode->getType()){
////        uiint idof;
////        uiint numOfScalar = pNode->numOfScalarParam();
////        uiint numOfVector = pNode->numOfVectorParam();
////
////        for(idof=0; idof < numOfScalar; idof++)
////            s_value[idof] = pNode->getScalar(idof);
////        for(idof=0; idof < numOfVector; idof++)
////            v_value[idof] = pNode->getVector(idof);
////    }else{
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetSVNodeValue(uint&, double [], double [])");
////    }
////}
////void CMW::GetSVNodeValue(const uiint& node_id,
////        const uiint& v_dof, double& v_value, const uiint& s_dof, double& s_value)
////{
////    CNode *pNode;
////    pNode= mpMesh->getNode(node_id);
////
////    if(NodeType::ScalarVector==pNode->getType()){
////        s_value = pNode->getScalar(s_dof);
////        v_value = pNode->getVector(v_dof);
////    }else{
////        Utility::CLogger *pLogger = Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetSVNodeValue(uint&, uint&, double&, uint&, double&)");
////    }
////}




//--
// 積分点 情報: ShapeType別の積分点数を返す.
//--
// numOfInteg := 積分点数
uiint& CMW::NumOfIntegPoint(const uiint& shapeType)
{
    return mpShapeCatalog->NumOfIntegPoint(shapeType);
}
//    switch(shapeType){
//        case(ShapeType::Hexa81):
//            break;
//        case(ShapeType::Hexa82):
//            break;
//        case(ShapeType::Hexa201):
//            break;
//        case(ShapeType::Hexa202):
//            break;
//        case(ShapeType::Hexa203):
//            break;
//        case(ShapeType::HexaNic111):
//            break;
//        case(ShapeType::HexaNic118):
//            break;
//        case(ShapeType::HexaNic1127):
//            break;
//        case(ShapeType::Tetra41):
//            break;
//        case(ShapeType::Tetra101):
//            break;
//        case(ShapeType::Tetra104):
//            break;
//        case(ShapeType::Tetra1015):
//            break;
//        case(ShapeType::Prism62):
//            break;
//        case(ShapeType::Prism156):
//            break;
//        case(ShapeType::Prism159):
//            break;
//        case(ShapeType::Prism1518):
//            break;
//        case(ShapeType::Quad41):
//            break;
//        case(ShapeType::Quad84):
//            break;
//        case(ShapeType::Quad89):
//            break;
//        case(ShapeType::Triangle31):
//            break;
//        case(ShapeType::Triangle63):
//            break;
//        case(ShapeType::Line21):
//            break; 
//        case(ShapeType::Line32):
//            break;
//        default:
//            break;
//    }

//--
// N :形状関数
//--
// N(積分点ごと)
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


// N(まるごと)
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
// ----
// Fortran API用 形状関数
// ----
// N[igauss][ishape]
//
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



// dN/dr(積分点ごと)
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
// dN/dr(まるごと)
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
//        case(ShapeType::HexaNic111):
//            dNdr= mpShapeHexaNic->dNdr111();
//            break;
//        case(ShapeType::HexaNic118):
//            dNdr= mpShapeHexaNic->dNdr118();
//            break;
//        case(ShapeType::HexaNic1127):
//            dNdr= mpShapeHexaNic->dNdr1127();
//            break;
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


// ----
// Fortran API用 導関数
// ----
// dNdr[igauss][ishape][iaxis]
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


// dN/dx 計算
//
void CMW::Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index)
{
    // clear...
    mvdNdx.clear();

    //AssyModel,Mesh(Part)は選択してあるとする.
    CElement *pElement;
    pElement= mpMesh->getElementIX(elem_index);

    // switch文は,要素クラスにdNdx関数を移設すれば,解消するが,しばらくこのまま.
    switch(elemType){
        // 1次要素 Hexa
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
        // 2次要素 Hexa
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
        // 1次要素 Tetra
        case(ElementType::Tetra):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx4(numOfInteg,pElement);
                mvdNdx= mpShapeTetra->dNdx41();
            }
            break;
        // 2次要素 Tetra
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

// dN/dx(積分点ごと) : Calculate_dNdx(....)コール後に使用
//
void CMW::dNdx_on_pt(const uiint& igauss, vvdouble& dNdX)
{
    dNdX = mvdNdx[igauss];
}
// dN/dx(まるごと)
//
void CMW::dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdX)
{
    //AssyModel,Mesh(Part)は選択してあるとする.
    CElement *pElement;
    pElement= mpMesh->getElementIX(elem_index);

    //switch文は,しばらくこのまま.
    switch(elemType){
        // 1次要素 Hexa
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
        // 2次要素 Hexa
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
        // 1次要素 Tetra
        case(ElementType::Tetra):
            if(numOfInteg==1){
                mpShapeTetra->Calc_dNdx4(numOfInteg,pElement);
                dNdX= mpShapeTetra->dNdx41();
            }
            break;
        // 2次要素 Tetra
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
//---
// dNdx Fortran用途
//---
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


//--
// ヤコビアン行列式
//--
void CMW::detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ)
{
    // #これを呼び出す直前に使用した,dNdx計算の｜J|の値
    // -------------------------------------------
    switch(elemType){
        // 1次要素 Hexa
        case(ElementType::Hexa):
            if(numOfInteg==1){
                detJ= mpShapeHexa->detJ81(igauss);
            }
            if(numOfInteg==8){
                detJ= mpShapeHexa->detJ82(igauss);
            }
            break;
        // 2次要素 Hexa
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
        // 1次要素 Tetra
        case(ElementType::Tetra):
            if(numOfInteg==1){
                detJ= mpShapeTetra->detJ41(igauss);
            }
            break;
        // 2次要素 Tetra
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

// ガウス積分点の重み : Weight
//
void CMW::Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w)
{
    switch(elemType){
        //  Hexa
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
        
        // Tetra
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

        // Prism
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
        //-------
        default:
            break;
    }
}
/*
    switch(elemType){
        // 1次要素 Hexa
        case(ElementType::Hexa):
            if(numOfInteg==1){
                ;
            }
            if(numOfInteg==8){
                ;
            }
            break;
        // 2次要素 Hexa
        case(ElementType::Hexa2):
            if(numOfInteg==1){
                ;
            }
            if(numOfInteg==8){
                ;
            }
            if(numOfInteg==27){
                ;
            }
            break;
        // 1次要素 Tetra
        case(ElementType::Tetra):
            if(numOfInteg==1){
                ;
            }
            break;
        // 2次要素 Tetra
        case(ElementType::Tetra2):
            if(numOfInteg==1){
                ;
            }
            if(numOfInteg==4){
                ;
            }
            if(numOfInteg==15){
                ;
            }
            break;
        case(ElementType::Prism):
            if(numOfInteg==2){
                ;
            }
            break;
        case(ElementType::Prism2):
            if(numOfInteg==6){
                ;
            }
            if(numOfInteg==9){
                ;
            }
            if(numOfInteg==18){
                ;
            }
            break;
        default:
            break;
    }
*/
//----
// 節点タイプ for Fortran
//----
uiint CMW::nodetype_s(){ return NodeType::Scalar;}
uiint CMW::nodetype_v(){ return NodeType::Vector;}
uiint CMW::nodetype_sv(){ return NodeType::ScalarVector;}
//----
// 要素タイプ for Fortran
//----
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
//----
// FrontISTR 要素タイプ
//----
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
//----
// FrontISTR 要素タイプ　=> MW3 要素タイプ 変換
//----
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
//----
// MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
//----
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


//----
// 形状関数タイプ for Fortran
//----
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


//--
// Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
//--
uiint CMW::GetNumOfBoundaryNodeMesh(){ return mpMesh->getNumOfBoundaryNodeMesh();}
uiint CMW::GetNumOfBoundaryFaceMesh(){ return mpMesh->getNumOfBoundaryFaceMesh();}
uiint CMW::GetNumOfBoundaryEdgeMesh(){ return mpMesh->getNumOfBoundaryEdgeMesh();}
uiint CMW::GetNumOfBoundaryVolumeMesh(){ return mpMesh->getNumOfBoundaryVolumeMesh();}

// BND Type (Neumann || Dirichlet)
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
// BoundaryTypeを表す型(定数)
uiint CMW::getNeumannType(){ return BoundaryType::Neumann;}
uiint CMW::getDirichletType(){ return BoundaryType::Dirichlet;}

//--
// BNode数
//--
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
//--
// DOF数
//--
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
//--
// Boundary DOF番号("DOFインデックス"に対応するDOF番号)
//--
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
//--
// BoundaryNode 境界値
//--
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
//--
// Boundary Node の NodeID
//--
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
//--
// Face, Edge, Volume の境界値
//--
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
//--
// Face, Edge, Volume のNodeID
//--
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
//--
// Boundaryの名称
//--
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
//--
// Entity_ID of BoundaryMesh (for FrontISTR)
//--
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
// MPI (直接使用)
//--
int& CMW::GetRank()//自分のプロセス-ランクを取得
{
    return mpMPI->getRank();
}
int& CMW::GetNumOfProcess()
{
    return mpMPI->getNumOfProcess();
}
int CMW::AllReduce(void* sendbuf, void* recvbuf, int buf_size, int datatype, int op, int commworld)
{
    return mpMPI->Allreduce(sendbuf, recvbuf, buf_size, datatype, op, commworld);
}
int CMW::Barrier(int commworld)
{
    return mpMPI->Barrier(commworld);
}
int CMW::Abort(int commworld, int error)
{
    return mpMPI->Abort(commworld, error);
}
int CMW::AllGather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm)
{
    return mpMPI->Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
}
int CMW::Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    return mpMPI->Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm);
}
int CMW::Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    return mpMPI->Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}
int CMW::Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
    return mpMPI->Send(buf, count, datatype, dest, tag, comm);
}
int CMW::Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status)
{
    return mpMPI->Recv(buf, count, datatype, source, tag, comm, status);
}
int CMW::Bcast(void* buf, int cnt, MPI_Datatype type, int root, MPI_Comm comm)
{
    return mpMPI->Bcast(buf, cnt, type, root, comm);
}

// 以下の３メソッドは、ペア
// -----
// 1.Meshパーツが通信する相手の数 (SelectされているLevel)
uiint CMW::GetNumOfNeibPE(const uiint& imesh)
{
    CMesh *pMesh = mpAssy->getMesh(imesh);
    
    return pMesh->getCommMesh2Size();//"CommMesh2"
}
// 2.通信Mesh毎のランク番号
uiint& CMW::GetTransRank(const uiint& imesh, const uiint& ipe)
{
    CMesh *pMesh = mpAssy->getMesh(imesh);
    CCommMesh2 *pCommMesh = pMesh->getCommMesh2IX(ipe);

    return pCommMesh->getTrasmitRank();
}
// 3.bufの値を送信、受信をbufに代入
//
void CMW::Send_Recv_R(double* buf, const int& num_of_node, const int& dof_size, const int& trans_rank)
{
    int nCount = num_of_node * dof_size;
    
    // 送信
    mpMPI->Send(buf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD);

    // 受信
    MPI_Status sta;
    mpMPI->Recv(buf, nCount, MPI_DOUBLE, trans_rank, 0, MPI_COMM_WORLD, &sta);

}
void CMW::Send_Recv_I(int* buf, const int& num_of_node, const int& dof_size, const int& trans_rank)
{
    int nCount = num_of_node * dof_size;

    // 送信
    mpMPI->Send(buf, nCount, MPI_INT, trans_rank, 0, MPI_COMM_WORLD);

    // 受信
    MPI_Status sta;
    mpMPI->Recv(buf, nCount, MPI_INT, trans_rank, 0, MPI_COMM_WORLD, &sta);
}

////// bufの値を送信 => 受信した値をbufに代入、Nodeにセット
////// # bufのサイズ == NumOfCommNode * dof_size
//////
////void CMW::Send_Recv_R(double* buf, const uiint& nDOF)
////{
////    uiint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
////    CAssyModel *pAssy;
////    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
////        pAssy = mpGMGModel->getAssyModel(iLevel);
////
////        uiint imesh, nNumOfMesh = pAssy->getNumOfMesh();
////        CMesh *pMesh;
////        for(imesh=0; imesh < nNumOfMesh; imesh++){
////            pMesh = pAssy->getMesh(imesh);
////
////            uiint icmesh, nNumOfCommMesh = pMesh->getNumOfCommMesh();
////            CCommMesh2 *pCommMesh;
////            for(icmesh=0; icmesh < nNumOfCommMesh; icmesh++){
////                pCommMesh = pMesh->getCommMesh2IX(icmesh);
////
////                int transRank = pCommMesh->getTrasmitRank();
////
////                uiint icnode, nNumOfCommNode = pCommMesh->getCommNodeSize();
////                CCommNode *pCommNode;
////                for(icnode=0; icnode < nNumOfCommNode; icnode++){
////                    pCommNode = pCommMesh->getCommNodeIX(icnode);
////
////                    double sendbuf[nDOF];
////                    double recvbuf[nDOF];
////
////                    uiint idof;
////                    // 送信buf へ代入
////                    for(idof=0; idof < nDOF; idof++) sendbuf[idof] = buf[icnode*nDOF + idof];
////
////                    // 送信
////                    mpMPI->Send(sendbuf, nDOF, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD);
////
////
////                    // 受信
////                    MPI_Status sta;
////                    mpMPI->Recv(recvbuf, nDOF, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD, &sta);
////
////                    CNode* pNode = pCommNode->getNode();
////                    //
////                    // ・受信bufから引数bufへ代入
////                    // ・受信値をNodeの変数に代入
////                    //
////                    for(idof=0; idof < nDOF; idof++){
////
////                        buf[icnode*nDOF + idof] = recvbuf[idof];
////
////                        if(pNode->getType()==NodeType::Scalar) pNode->setScalar(recvbuf[idof], idof);
////                        if(pNode->getType()==NodeType::Vector) pNode->setVector(recvbuf[idof], idof);
////                        if(pNode->getType()==NodeType::ScalarVector){
////                            uiint nNumScalar = pNode->numOfScalarParam();
////                            if(idof < nNumScalar){
////                                pNode->setScalar(recvbuf[idof], idof);
////                            }else{
////                                pNode->setVector(recvbuf[idof-nNumScalar], idof);
////                            }
////                        }
////                    };// idof (Total DOF)
////                };// icnode (CommNode)
////            };// icmesh (CommMesh)
////        };// imesh (Mesh)
////    };// iLevel
////}
////
//////
////// 通信Nodeの値を入れ替えて更新
//////
////void CMW::Send_Recv_R()
////{
////    uiint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
////    CAssyModel *pAssy;
////    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
////        pAssy = mpGMGModel->getAssyModel(iLevel);
////
////        uiint imesh, nNumOfMesh = pAssy->getNumOfMesh();
////        CMesh *pMesh;
////        for(imesh=0; imesh < nNumOfMesh; imesh++){
////            pMesh = pAssy->getMesh(imesh);
////
////            uiint icmesh, nNumOfCommMesh = pMesh->getNumOfCommMesh();
////            CCommMesh2 *pCommMesh;
////            for(icmesh=0; icmesh < nNumOfCommMesh; icmesh++){
////                pCommMesh = pMesh->getCommMesh2IX(icmesh);
////
////                int transRank = pCommMesh->getTrasmitRank();
////
////                uiint icnode, nNumOfCommNode = pCommMesh->getCommNodeSize();
////                CCommNode *pCommNode;
////                for(icnode=0; icnode < nNumOfCommNode; icnode++){
////                    pCommNode = pCommMesh->getCommNodeIX(icnode);
////                    CNode* pNode = pCommNode->getNode();
////
////                    uiint dof_size = pNode->numOfTotalParam();
////                    double sendbuf[dof_size];
////                    double recvbuf[dof_size];
////
////                    uiint idof;
////                    // 送信bufにNodeの変数値を代入
////                    //
////                    for(idof=0; idof < dof_size; idof++){
////                        if(pNode->getType()==NodeType::Scalar) sendbuf[idof] = pNode->getScalar(idof);
////                        if(pNode->getType()==NodeType::Vector) sendbuf[idof] = pNode->getVector(idof);
////                        if(pNode->getType()==NodeType::ScalarVector){
////                            uiint nNumScalar = pNode->numOfScalarParam();
////                            if(idof < nNumScalar){
////                                sendbuf[idof] = pNode->getScalar(idof);
////                            }else{
////                                sendbuf[idof] = pNode->getVector(idof-nNumScalar);
////                            }
////                        }
////                    }
////
////                    // 送信
////                    mpMPI->Send(sendbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD);
////
////
////                    // 受信
////                    MPI_Status sta;
////                    mpMPI->Recv(recvbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD, &sta);
////
////                    // 受信値をNodeの変数に代入
////                    //
////                    for(idof=0; idof < dof_size; idof++){
////                        if(pNode->getType()==NodeType::Scalar) pNode->setScalar(recvbuf[idof], idof);
////                        if(pNode->getType()==NodeType::Vector) pNode->setVector(recvbuf[idof], idof);
////                        if(pNode->getType()==NodeType::ScalarVector){
////                            uiint nNumScalar = pNode->numOfScalarParam();
////                            if(idof < nNumScalar){
////                                pNode->setScalar(recvbuf[idof], idof);
////                            }else{
////                                pNode->setVector(recvbuf[idof-nNumScalar], idof);
////                            }
////                        }
////                    };// idof (Total DOF)
////                };// icnode (CommNode)
////            };// icmesh (CommMesh)
////        };// imesh (Mesh)
////    };// iLevel
////}


//--
// CommMesh2 (通信Mesh)  { select された AssyModel,Meshを対象 }
// 
// CommNode (通信Node) : Visualizer 用途
//--
uiint CMW::GetNumOfCommMesh()
{
    return mpMesh->getCommMesh2Size();//"CommMesh2"
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
// グループ  { select された AssyModel,Meshを対象 }
//--
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


//--
// Logger
//--
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
void CMW::LoggerInfo(const uiint& mode, char* message)
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
void CMW::LoggerInfo(const uiint& mode, const char* message)
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
// Logger_Parameter for Fortran
//--
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


