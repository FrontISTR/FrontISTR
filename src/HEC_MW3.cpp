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

// MPI, Logger の初期化
//
int CMW::Initialize(int argc, char** argv, const char* path)
{
    // MPI 初期化
    mpMPI->Initialize(argc, argv);
    // rank 取得
    uint rank= mpMPI->getRank();


    //-- このLogger設定はテスト用 --
    //
    mpLogger->setMode(Utility::LoggerMode::MWDebug);

    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error,   Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn,    Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info,    Utility::LoggerDevice::Display);
    //-- テスト用途 end --

    // Logger初期化
    mpLogger->initializeLogFile(rank);// ログファイル-オープン
    mpLogger->InfoDisplay();
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");

    //cntファイルへのパス指定
    mpFileIO->setPathName(path);

    //cntファイルから入力ファイルのベースネームを取得
    mpFileIO->ReadCntFile();

    // ファイルパス
    msInputFileName= mpFileIO->getPathName();
    msOutputFileName=mpFileIO->getPathName();

    // ファイルパスにファイルベース名を追加, アンダースコアを追加
    msInputFileName += mpFileIO->getMeshFileBaseName() + "_";
    msOutputFileName+= mpFileIO->getMeshFileBaseName() + "_";

    // 自分のrankを"ファイルベース名"に追加
    msInputFileName  += boost::lexical_cast<string>(rank);
    msOutputFileName += boost::lexical_cast<string>(rank);

    // 拡張子を追加
    msInputFileName  += ".msh";
    msOutputFileName += ".out";

    return 1;
}

//  MPI, Logger 後処理
//
int CMW::Finalize()
{
    // MPI
    mpMPI->Finalize();

    // Logger
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();//ログファイル-クローズ
    return 1;
}


// FileIO内からFactoryコール、
// 1.Factory内でAssyModel(マルチグリッド階層数分の生成)
// 2.Mesh(Level=0)を生成
//
int CMW::FileRead()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");

    mpFileIO->ReadFile(msInputFileName);
    mb_file = true;

    return 1;
}

int CMW::FileWrite()// const uint& nmgLevel
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileWrite");

    uint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile(msOutputFileName, nLevel+1);//階層数+1 <= 出力数

    return 1;
}


//int CMW::Initialize_Matrix()// const uint& nmgLevel
//{
//#ifdef ADVANCESOFT_DEBUG
//   	printf(" enter CMW::Initialize_Matrix \n");
//#endif
//
//	mpAssyMatrix = new CAssyMatrix(mpAssy);
//	uint level = mpAssy->getMGLevel();
//	mpAssyMatrix->setCoarseMatrix( NULL );
//	if( level > 0 ) mpAssyMatrix->setCoarseMatrix( mpGMGModel->getAssyModel(level-1)->getAssyMatrix() );
//
//#ifdef ADVANCESOFT_DEBUG
//   	printf(" exit CMW::Initialize_Matrix \n");
//#endif
//    return 1;
//}
//
//int CMW::Initialize_Vector()// const uint& nmgLevel
//{
//#ifdef ADVANCESOFT_DEBUG
//   	printf(" enter CMW::Initialize_Vector \n");
//#endif
//
//	mpAssyVector = new CAssyVector(mpAssy);
//	mpAssyVector2 = new CAssyVector(mpAssy);
//#ifdef ADVANCESOFT_DEBUG
//   	printf(" exit CMW::Initialize_Vector \n");
//#endif
//    return 1;
//}

//
// 全Levelに方程式を生成(AssyMatrixのコースグリッドも、ここで設定される)
//
void CMW::GeneLinearAlgebra(const uint& nNumOfAlgebra, uint* vDOF)
{
    vuint vecDOF;
    vecDOF.reserve(nNumOfAlgebra);
    for(uint idof=0; idof < nNumOfAlgebra; idof++){
        vecDOF.push_back(vDOF[idof]);
    };

    CAssyModel *pAssyModel, *pCoaAssyModel;
    uint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
        pAssyModel = mpGMGModel->getAssyModel(iLevel);

        if(iLevel > 0){
            pCoaAssyModel = mpGMGModel->getAssyModel(iLevel-1);
            pAssyModel->GeneLinearAlgebra(vecDOF, pCoaAssyModel);
        }else{
            pAssyModel->GeneLinearAlgebra(vecDOF, NULL);
        }
    };
}
//
// iequ番目の方程式をロード{ Levelは選択済み == AssyModelは選択済み }
//
void CMW::SelectAlgebra(const uint& iequ)
{
    // AssyMatrixのコースグリッドはGeneLinearAlgebraで既に設定済み
    //
    mpAssyMatrix = mpAssy->getAssyMatrix(iequ);
    mpRHSAssyVector = mpAssy->getRHSAssyVector(iequ);
    mpSolAssyVector = mpAssy->getSolutionAssyVector(iequ);
}


int CMW::Matrix_Add_Elem(const uint& iMesh, const uint& iElem, double *ElemMatrix)
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

int CMW::Matrix_Add_Node(const uint& iMesh, const uint& iNodeID, const uint& jNodeID, double* NodalMatrix)
{
    return mpAssyMatrix->Matrix_Add_Nodal(iMesh, iNodeID, jNodeID, NodalMatrix);
}

// Matrix 0 clear
//
void CMW::Matrix_Clear(const uint& iMesh)
{
    mpAssyMatrix->Matrix_Clear(iMesh);
}
void CMW::Vector_Clear(const uint& iMesh)
{
    mpRHSAssyVector->Vector_Clear(iMesh);
    mpSolAssyVector->Vector_Clear(iMesh);
}

void CMW::Sample_Set_BC(uint iMesh)
{
#ifdef ADVANCESOFT_DEBUG
	printf(" enter CMW::Sample_Set_BC \n");
#endif

	double X, Y, Z, value0, value1, value2;
	uint iDof0, iDof1, iDof2;
	uint iNodeMax = getNodeSize( iMesh );

	for( uint iNode = 0; iNode < iNodeMax; iNode++){
		X = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getX();
		Y = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getY();
		Z = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getZ();
		if(abs( Z - 4.0 ) < 1.0e-5 && abs( X - 1.0 ) < 1.0e-5 ) {
			value0 = 1.0e6;
			iDof0 = 0;
			Set_BC_RHS( iMesh, iNode, iDof0, value0);
		};
		if( (abs( Z ) < 1.0e-5) || (abs( Z - 8.0 ) < 1.0e-5) ) {
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
int CMW::Set_BC_RHS(uint& iMesh, uint& iNode, uint& iDof, double& value)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMW::Set_BC (RHS) %d %d %d %e \n", iMesh, iNode, iDof, value);
#endif

	mpRHSAssyVector->setValue(iMesh, iNode, iDof, value);

#ifdef ADVANCESOFT_DEBUG
	printf("exit CMW::Set_BC (RHS) \n");
#endif
	return 1;
}
//
// 行列の対角項、解ベクトルへ境界値をセット
//
int CMW::Set_BC_Mat_SolVec(uint& iMesh, uint& iNode, uint& iDof, double& value1, double& value2)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMW::Set_BC (Mat_SolVec) %d %d %d %e %e \n", iMesh, iNode, iDof, value1, value2);
#endif

	mpAssyMatrix->setValue(iMesh, iNode, iDof, value1);
        //mpRHSAssyVector->setValue(iMesh, iNode, iDof, value2);
	mpSolAssyVector->setValue(iMesh, iNode, iDof, value2);

#ifdef ADVANCESOFT_DEBUG
	printf("exit CMW::Set_BC (Mat_SolVec) \n");
#endif

	return 1;
}
//
// 行列の対角項、右辺ベクトルへ境界値をセット
//
int CMW::Set_BC_Mat_RHS(uint& iMesh, uint& iNode, uint& iDof, double& value1, double& value2)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMW::Set_BC (Mat_RHS) %d %d %d %e %e \n", iMesh, iNode, iDof, value1, value2);
#endif

	mpAssyMatrix->setValue(iMesh, iNode, iDof, value1);
        mpRHSAssyVector->setValue(iMesh, iNode, iDof, value2);

#ifdef ADVANCESOFT_DEBUG
	printf("exit CMW::Set_BC (Mat_RHS) \n");
#endif

	return 1;
}

int CMW::Solve(uint& iter_max, double& tolerance, uint& method, uint& precondition)
{
#ifdef ADVANCESOFT_DEBUG
  printf(" enter CMW::Solve %d %e \n", iter_max, tolerance);
#endif

  bool flag_iter_log = false;
  bool flag_time_log = false;
  char cfile[100];
  
   switch( method ){
        case( 1 ):
	   {CSolverCG *solver = new CSolverCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 2 ):
	   {CSolverBiCGSTAB *solver = new CSolverBiCGSTAB(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 3 ):
           {CSolverGPBiCG *solver = new CSolverGPBiCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        case( 4 ):
	   {CSolverGMRES *solver = new CSolverGMRES(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
	   solver->solve(mpAssyMatrix, mpRHSAssyVector, mpSolAssyVector);}
           break;
        default:
           break;
    }


    //
    // MicroAVS UCD format 出力
    //
    uint iMesh, iMeshMax;
    iMeshMax = GetNumOfMeshPart();
    for( iMesh = 0; iMesh < iMeshMax; iMesh++){
	
	FILE *fp1;
	sprintf(cfile, "outUCD%02d%02d_%05d.inp", iMesh, mpAssy->getMGLevel(), mpMPI->getRank()); 
	fp1 = fopen(cfile, "w");
	fprintf(fp1,"%d %d 3 0 0 \n",mpAssy->getMesh(iMesh)->getNumOfNode(), mpAssy->getMesh(iMesh)->getNumOfElement());
	for(int i=0; i < mpAssy->getMesh(iMesh)->getNumOfNode(); i++)
	{
            double x = mpAssy->getMesh(iMesh)->getNodeIX(i)->getX();
            double y = mpAssy->getMesh(iMesh)->getNodeIX(i)->getY();
            double z = mpAssy->getMesh(iMesh)->getNodeIX(i)->getZ();
            fprintf(fp1,"%d %e %e %e\n",i, x, y, z);

	}// Node ループ END

        vector<CNode*> vNode;
        
	for(int i=0; i < mpAssy->getMesh(iMesh)->getNumOfElement(); i++)
	{
            vNode.clear();
            CElement *pElem = mpAssy->getMesh(iMesh)->getElementIX(i);
            vNode = pElem->getNode();
            
            if(pElem->getType()==ElementType::Hexa){
                fprintf(fp1,"%d 0 hex %d %d %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID(),vNode[6]->getID(),vNode[7]->getID());
            }
            if(pElem->getType()==ElementType::Hexa2){
                fprintf(fp1,"%d 0 hex2 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID(),vNode[6]->getID(),vNode[7]->getID(),
                        vNode[8]->getID(),vNode[9]->getID(),vNode[10]->getID(),vNode[11]->getID(),
                        vNode[12]->getID(),vNode[13]->getID(),vNode[14]->getID(),vNode[15]->getID(),
                        vNode[16]->getID(),vNode[17]->getID(),vNode[18]->getID(),vNode[19]->getID());
            }
            if(pElem->getType()==ElementType::Tetra){
                fprintf(fp1,"%d 0 tet %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID());
            }
            if(pElem->getType()==ElementType::Tetra2){
                fprintf(fp1,"%d 0 tet2 %d %d %d %d %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID(),vNode[6]->getID(),vNode[7]->getID(),
                        vNode[8]->getID(),vNode[9]->getID());
            }
            if(pElem->getType()==ElementType::Prism){
                fprintf(fp1,"%d 0 prism %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID());
            }
            if(pElem->getType()==ElementType::Prism2){
                fprintf(fp1,"%d 0 prism2 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID(),vNode[6]->getID(),vNode[7]->getID(),
                        vNode[8]->getID(),vNode[9]->getID(),vNode[10]->getID(),vNode[11]->getID(),
                        vNode[12]->getID(),vNode[13]->getID(),vNode[14]->getID());
            }
            if(pElem->getType()==ElementType::Quad){
                fprintf(fp1,"%d 0 quad %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID());
            }
            if(pElem->getType()==ElementType::Quad2){
                fprintf(fp1,"%d 0 quad2 %d %d %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID(),vNode[6]->getID(),vNode[7]->getID());
            }
            if(pElem->getType()==ElementType::Triangle){
                fprintf(fp1,"%d 0 tri %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID());
            }
            if(pElem->getType()==ElementType::Triangle2){
                fprintf(fp1,"%d 0 tri2 %d %d %d %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID(),vNode[3]->getID(),
                        vNode[4]->getID(),vNode[5]->getID());
            }
            if(pElem->getType()==ElementType::Beam){
                fprintf(fp1,"%d 0 line %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID());
            }
            if(pElem->getType()==ElementType::Beam2){
                fprintf(fp1,"%d 0 line2 %d %d %d \n",i,
                        vNode[0]->getID(),vNode[1]->getID(),vNode[2]->getID());
            }
	}// Element ループ END
        
        uint nDOF=mpSolAssyVector->getDOF();// 解ベクトルのDOF
        
	fprintf(fp1,"1 %d \n", nDOF);
	fprintf(fp1,"vector, aaa \n");
        
	for(int i=0; i < mpAssy->getMesh(iMesh)->getNumOfNode(); i++)
	{
            double val[nDOF];
            for(uint idof=0; idof < nDOF; idof++){
                val[idof] = mpSolAssyVector->getValue(iMesh, i, idof);// 解ベクトルの内容
            }
            //for(uint idof=0; idof < nDOF; idof++){
            //    mpAssy->getMesh(iMesh)->getNodeIX(i)->setVector(val[idof], idof);
            //}
            if(nDOF==1) fprintf(fp1,"%d %e\n",i, val[0]);
            if(nDOF==2) fprintf(fp1,"%d %e %e\n",i, val[0], val[1]);
            if(nDOF==3) fprintf(fp1,"%d %e %e %e\n",i, val[0], val[1], val[2]);
            if(nDOF==4) fprintf(fp1,"%d %e %e %e %e\n",i, val[0], val[1], val[2], val[3]);
            if(nDOF==5) fprintf(fp1,"%d %e %e %e %e %e\n",i, val[0], val[1], val[2], val[3], val[4]);
            if(nDOF==6) fprintf(fp1,"%d %e %e %e %e %e %e\n",i, val[0], val[1], val[2], val[3], val[4], val[5]);
            if(nDOF==7) fprintf(fp1,"%d %e %e %e %e %e %e %e\n",i, val[0], val[1], val[2], val[3], val[4], val[5], val[6]);

	}// Node ループ END

	fclose(fp1);

    }//iMesh ループ END

    printf(" --- end of output to a UCD file. (%s) --- \n", cfile);

#ifdef ADVANCESOFT_DEBUG
    printf(" exit CMW::Solve \n");
#endif
    return 1;
}

//
// 1.マルチグリッドのための,Refineメッシュ生成
// 2.MPCのため,接合Meshのマスター面にスレーブ点をセット
//
int CMW::Refine()
{
    uint ilevel,mgLevel;
    CAssyModel *pAssy;
    uint icon,numOfConMesh;
    CContactMesh *pConMesh;
    
    if(mb_file){// ファイル読み込み後
        
        mpFactory->refineMesh();       //MultiGridデータの生成(通信Meshも含む)
        mpFactory->refineContactMesh();//接合Mesh(MPCメッシュ)のMultiGridデータ生成
        
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

        mpFactory->refineCommMesh2();//要素分割型(節点共有型)の通信Mesh(CommMesh2)のRefine
        mpFactory->refineBoundary(); //境界条件の階層化
        
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
    uint ilevel, mgLevel = mpFactory->getMGLevel();
    
    for(ilevel=0; ilevel < mgLevel+1; ilevel++){
        pAssy = mpGMGModel->getAssyModel(ilevel);

        CMesh *pMesh;
        uint imesh, nNumOfMesh = pAssy->getNumOfMesh();

        // Mesh      メッシュ本体
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            pMesh = pAssy->getMesh(imesh);
            
            pMesh->deleteProgData();// Elementの各Edge,Face の Node*,Element* 配列 を解放
            pMesh->deleteAggregate_on_Node();// Node::AggElemIDの解放

            // CommMesh2 通信メッシュ
            CCommMesh2 *pCommMesh;
            uint icom, nNumOfCommMesh = pMesh->getCommMesh2Size();
            for(icom=0; icom < nNumOfCommMesh; icom++){
                pCommMesh = pMesh->getCommMesh2IX(icom);
                
                pCommMesh->deleteProgData();
            };

            // BoundaryVolumeMesh 境界メッシュ(体積)
            CBoundaryVolumeMesh *pBVolMesh;
            uint ibvol, nNumOfBVolMesh = pMesh->getNumOfBoundaryVolumeMesh();
            for(ibvol=0; ibvol < nNumOfBVolMesh; ibvol++){
                pBVolMesh = pMesh->getBndVolumeMeshIX(ibvol);

                pBVolMesh->deleteProgData();
            };

            // BoundaryFaceMesh 境界メッシュ(面)
            CBoundaryFaceMesh *pBFaceMesh;
            uint ibface, nNumOfBFaceMesh = pMesh->getNumOfBoundaryFaceMesh();
            for(ibface=0; ibface < nNumOfBFaceMesh; ibface++){
                pBFaceMesh = pMesh->getBndFaceMeshIX(ibface);

                pBFaceMesh->deleteProgData();
            };

            // BoundaryEdgeMesh 境界メッシュ(辺)
            CBoundaryEdgeMesh *pBEdgeMesh;
            uint ibedge, nNumOfBEdgeMesh = pMesh->getNumOfBoundaryEdgeMesh();
            for(ibedge=0; ibedge < nNumOfBEdgeMesh; ibedge++){
                pBEdgeMesh = pMesh->getBndEdgeMeshIX(ibedge);

                pBEdgeMesh->deleteProgData();
            };
        };

        // ContactMesh 接合面メッシュ
        CContactMesh *pConMesh;
        uint icmesh, nNumOfConMesh = pAssy->getNumOfContactMesh();

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
uint CMW::GetNumOfAssembleModel()
{
    return mpGMGModel->getNumOfLevel();
}

// Assemble Modelの選択
void CMW::SelectAssembleModel(const uint& mgLevel)
{
    mpAssy= mpGMGModel->getAssyModel(mgLevel);
}

// Meshパーツの個数
uint CMW::GetNumOfMeshPart()
{
    return mpAssy->getNumOfMesh();
}

// Meshの選択 : IDによる選択
void CMW::SelectMeshPart_ID(const uint& mesh_id)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh_ID(mesh_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
// Meshの選択 : 配列Index番号による選択
void CMW::SelectMeshPart_IX(const uint& index)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
// Elementの選択 : IDによる選択
void CMW::SelectElement_ID(const uint& elem_id)
{
    if(mpMesh){
        mpElement= mpMesh->getElement(elem_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
// Elementの選択 : 配列Index番号による選択
void CMW::SelectElement_IX(const uint& index)
{
    if(mpMesh){
        mpElement= mpMesh->getElementIX(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
// Elementのタイプ取得
uint CMW::GetElementType()
{
    return mpElement->getType();
}

// Elementの頂点数
uint CMW::GetNumOfElementVert()
{
    return mpElement->getNumOfNode();
}

// Elementの頂点のノードID
void CMW::GetElementVertNodeID(int* vNodeID)
{
    uint numOfNode;
    numOfNode= mpElement->getNumOfNode();

    CNode* pNode;
    uint ivert;
    for(ivert=0; ivert< numOfNode; ivert++){
        pNode= mpElement->getNode(ivert);

        vNodeID[ivert]= pNode->getID();
    };
}

// Elementの辺の数
uint CMW::GetNumOfElementEdge()
{
    return mpElement->getNumOfEdge();
}

// Elementの辺のノードID
void CMW::GetElementEdgeNodeID(int* vNodeID)
{
    uint numOfEdge;
    numOfEdge= mpElement->getNumOfEdge();

    CNode *pNode;
    uint iedge;
    for(iedge=0; iedge< numOfEdge; iedge++){
        pNode= mpElement->getEdgeInterNode(iedge);

        vNodeID[iedge]= pNode->getID();
    };
}

// ID
uint& CMW::getElementID(const uint& index)
{
    CElement* pElem = mpMesh->getElementIX(index);

    return pElem->getID();
}
uint& CMW::getNodeID(const uint& index)
{
    CNode* pNode = mpMesh->getNodeIX(index);

    return pNode->getID();
}


//
// ノード
//
void CMW::GetNodeCoord(const uint& node_id, double& x, double& y, double& z)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    x= pNode->getX();
    y= pNode->getY();
    z= pNode->getZ();
}
uint CMW::GetNumOfDOF(const uint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    return pNode->numOfTotalParam();
}
uint CMW::GetNumOfScalar(const uint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    return pNode->numOfScalarParam();
}
uint CMW::GetNumOfVector(const uint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    return pNode->numOfVectorParam();
}
uint& CMW::GetNodeType(const uint& node_id)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    return pNode->getType();
}
void CMW::SetNodeValue(const uint& node_id, double value[])
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    uint numOfDOF = pNode->numOfTotalParam();

    uint idof;
    if(NodeType::Scalar==pNode->getType()){
        for(idof=0; idof < numOfDOF; idof++){
            pNode->setScalar(value[idof], idof);
        };
    }
    if(NodeType::Vector==pNode->getType()){
        for(idof=0; idof < numOfDOF; idof++){
            pNode->setVector(value[idof], idof);
        };
    }
    if(NodeType::ScalarVector==pNode->getType()){
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetNodeValue(const uint&, double[])");
    }
}
void CMW::SetNodeValue(const uint& node_id, const uint& idof, const double& value)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::Scalar==pNode->getType()){
        pNode->setScalar(value, idof);
    }
    if(NodeType::Vector==pNode->getType()){
        pNode->setVector(value, idof);
    }
    if(NodeType::ScalarVector==pNode->getType()){
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetNodeValue(const uint&, const uint&, double&)");
    }
}

void CMW::GetNodeValue(const uint& node_id, double value[])
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::Scalar==pNode->getType()){
        uint idof, numOfScalar=pNode->numOfScalarParam();
        for(idof=0; idof < numOfScalar; idof++)
            value[idof] = pNode->getScalar(idof);
    }
    if(NodeType::Vector==pNode->getType()){
        uint idof, numOfVector=pNode->numOfVectorParam();
        for(idof=0; idof < numOfVector; idof++)
            value[idof] = pNode->getVector(idof);
    }
    if(NodeType::ScalarVector==pNode->getType()){
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetNodeValue(uint&, double[])");
    }
}
double& CMW::GetNodeValue(const uint& node_id, const uint& idof)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::Scalar==pNode->getType()){
        return pNode->getScalar(idof);
    }
    if(NodeType::Vector==pNode->getType()){
        return pNode->getVector(idof);
    }
    if(NodeType::ScalarVector==pNode->getType()){
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetNodeValue(uint&, uint&)");
        return pLogger->getDDummyValue();
    }
}
//
// ScalarVector Node
//
void CMW::SetSVNodeValue(const uint& node_id, double v_value[], double s_value[])
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::ScalarVector==pNode->getType()){
        uint idof;

        uint numOfScalar=pNode->numOfScalarParam();
        for(idof=0; idof < numOfScalar; idof++) pNode->setScalar(s_value[idof], idof);

        uint numOfVector=pNode->numOfVectorParam();
        for(idof=0; idof < numOfVector; idof++) pNode->setVector(v_value[idof], idof);

    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetSVNodeValue(uint&, uint&, double&, uint&, double&)");
    }
}
void CMW::SetSVNodeValue(const uint& node_id,
        const uint& v_dof, const double& v_value, const uint& s_dof, const double& s_value)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::ScalarVector==pNode->getType()){
        pNode->setScalar(s_value, s_dof);
        pNode->setVector(v_value, v_dof);
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::SetSVNodeValue(uint&, uint&, double&, uint&, double&)");
    }
}
void CMW::GetSVNodeValue(const uint& node_id, double v_value[], double s_value[])
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::ScalarVector==pNode->getType()){
        uint idof;
        uint numOfScalar = pNode->numOfScalarParam();
        uint numOfVector = pNode->numOfVectorParam();

        for(idof=0; idof < numOfScalar; idof++)
            s_value[idof] = pNode->getScalar(idof);
        for(idof=0; idof < numOfVector; idof++)
            v_value[idof] = pNode->getVector(idof);
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetSVNodeValue(uint&, double [], double [])");
    }
}
void CMW::GetSVNodeValue(const uint& node_id,
        const uint& v_dof, double& v_value, const uint& s_dof, double& s_value)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);

    if(NodeType::ScalarVector==pNode->getType()){
        s_value = pNode->getScalar(s_dof);
        v_value = pNode->getVector(v_dof);
    }else{
        Utility::CLogger *pLogger = Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "NodeType mismatch, CMW::GetSVNodeValue(uint&, uint&, double&, uint&, double&)");
    }
}




//--
// 積分点 情報: ShapeType別の積分点数を返す.
//--
// numOfInteg := 積分点数
uint& CMW::NumOfIntegPoint(const uint& shapeType)
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
void CMW::ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N)
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
void CMW::ShapeFunc_on_pt(uint shapeType, uint igauss, double N[])
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

    uint i;
    for(i=0; i < vN.size(); i++) N[i]=vN[i];
}


// N(まるごと)
void CMW::ShapeFunc(const uint& shapeType, vvdouble& N)
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
double& CMW::ShapeFunc_Hexa81(int igauss, int ishape)
{
    return mpShapeHexa->N81(igauss,ishape);
}
double& CMW::ShapeFunc_Hexa82(int igauss, int ishape)
{
    return mpShapeHexa->N82(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa201(int igauss, int ishape)
{
    return mpShapeHexa->N201(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa202(int igauss, int ishape)
{
    return mpShapeHexa->N202(igauss, ishape);
}
double& CMW::ShapeFunc_Hexa203(int igauss, int ishape)
{
    return mpShapeHexa->N203(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra41(int igauss, int ishape)
{
    return mpShapeTetra->N41(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra101(int igauss, int ishape)
{
    return mpShapeTetra->N101(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra104(int igauss, int ishape)
{
    return mpShapeTetra->N104(igauss, ishape);
}
double& CMW::ShapeFunc_Tetra1015(int igauss, int ishape)
{
    return mpShapeTetra->N1015(igauss, ishape);
}
double& CMW::ShapeFunc_Prism62(int igauss, int ishape)
{
    return mpShapePrism->N62(igauss, ishape);
}
double& CMW::ShapeFunc_Prism156(int igauss, int ishape)
{
    return mpShapePrism->N156(igauss, ishape);
}
double& CMW::ShapeFunc_Prism159(int igauss, int ishape)
{
    return mpShapePrism->N159(igauss, ishape);
}
double& CMW::ShapeFunc_Prism1518(int igauss, int ishape)
{
    return mpShapePrism->N1518(igauss, ishape);
}
double& CMW::ShapeFunc_Quad41(int igauss, int ishape)
{
    return mpShapeQuad->N41(igauss, ishape);
}
double& CMW::ShapeFunc_Quad84(int igauss, int ishape)
{
    return mpShapeQuad->N84(igauss, ishape);
}
double& CMW::ShapeFunc_Quad89(int igauss, int ishape)
{
    return mpShapeQuad->N89(igauss, ishape);
}
double& CMW::ShapeFunc_Triangle31(int igauss, int ishape)
{
    return mpShapeTriangle->N31(igauss, ishape);
}
double& CMW::ShapeFunc_Triangle63(int igauss, int ishape)
{
    return mpShapeTriangle->N63(igauss, ishape);
}
double& CMW::ShapeFunc_Line21(int igauss, int ishape)
{
    return mpShapeLine->N21(igauss, ishape);
}
double& CMW::ShapeFunc_Line32(int igauss, int ishape)
{
    return mpShapeLine->N32(igauss, ishape);
}



// dN/dr(積分点ごと)
void CMW::dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr)
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
void CMW::dNdr(const uint& shapeType, vvvdouble& dNdr)
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
void CMW::dNdr(const uint& shapeType, double dNdr[])
{
    vvvdouble vdNdr;
    uint numOfInteg, numOfShape, numOfAxis;

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

    uint igauss,ishape,iaxis;
    uint pos=0;
    for(igauss=0; igauss < numOfInteg; igauss++){
        for(ishape=0; ishape < numOfShape; ishape++){
            for(iaxis=0; iaxis < numOfAxis; iaxis++){
                dNdr[pos] = vdNdr[igauss][ishape][iaxis];
                pos++;
            }
        }
    }
}
double& CMW::dNdr_Hexa81_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeHexa->dNdr81(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa82_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeHexa->dNdr82(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa201_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeHexa->dNdr201(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa202_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeHexa->dNdr202(igauss, ishape, iaxis);
}
double& CMW::dNdr_Hexa203_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeHexa->dNdr203(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra41_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTetra->dNdr41(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra101_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTetra->dNdr101(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra104_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTetra->dNdr104(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tetra1015_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTetra->dNdr1015(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism62_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapePrism->dNdr62(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism156_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapePrism->dNdr156(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism159_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapePrism->dNdr159(igauss, ishape, iaxis);
}
double& CMW::dNdr_Prism1518_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapePrism->dNdr1518(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad41_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeQuad->dNdr41(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad84_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeQuad->dNdr84(igauss, ishape, iaxis);
}
double& CMW::dNdr_Quad89_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeQuad->dNdr89(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tri31_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTriangle->dNdr31(igauss, ishape, iaxis);
}
double& CMW::dNdr_Tri63_on_pt_on_shape(int igauss, int ishape, int iaxis)
{
    return mpShapeTriangle->dNdr63(igauss, ishape, iaxis);
}
double& CMW::dNdr_Line21_on_pt_on_shape(int igauss, int ishape)
{
    return mpShapeLine->dNdr21(igauss, ishape);
}
double& CMW::dNdr_Line32_on_pt_on_shape(int igauss, int ishape)
{
    return mpShapeLine->dNdr32(igauss, ishape);
}


// dN/dx 計算
//
void CMW::Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index)
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
void CMW::dNdx_on_pt(const uint& igauss, vvdouble& dNdX)
{
    dNdX = mvdNdx[igauss];
}
// dN/dx(まるごと)
//
void CMW::dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdX)
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
void CMW::dNdx(const uint& elemType, const uint& numOfInteg, const uint& ielem, double dNdx[])
{
    Calculate_dNdx(elemType, numOfInteg, ielem);

    uint numOfShape;
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

    uint igauss, ishape, iaxis;
    uint pos=0;
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
void CMW::detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ)
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
void CMW::Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w)
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
uint CMW::nodetype_s(){ return NodeType::Scalar;}
uint CMW::nodetype_v(){ return NodeType::Vector;}
uint CMW::nodetype_sv(){ return NodeType::ScalarVector;}
//----
// 要素タイプ for Fortran
//----
uint CMW::elemtype_hexa(){ return ElementType::Hexa;}
uint CMW::elemtype_hexa2(){ return ElementType::Hexa2;}
uint CMW::elemtype_tetra(){ return ElementType::Tetra;}
uint CMW::elemtype_tetra2(){ return ElementType::Tetra2;}
uint CMW::elemtype_prism(){ return ElementType::Prism;}
uint CMW::elemtype_prism2(){ return ElementType::Prism2;}
uint CMW::elemtype_quad(){ return ElementType::Quad;}
uint CMW::elemtype_quad2(){ return ElementType::Quad2;}
uint CMW::elemtype_triangle(){ return ElementType::Triangle;}
uint CMW::elemtype_triangle2(){ return ElementType::Triangle2;}
uint CMW::elemtype_line(){ return ElementType::Beam;}
uint CMW::elemtype_line2(){ return ElementType::Beam2;}
//----
// FrontISTR 要素タイプ
//----
uint CMW::fistr_elemtype_hexa(){ return FistrElementType::Hexa;}
uint CMW::fistr_elemtype_hexa2(){ return FistrElementType::Hexa2;}
uint CMW::fistr_elemtype_tetra(){ return FistrElementType::Tetra;}
uint CMW::fistr_elemtype_tetra2(){ return FistrElementType::Tetra2;}
uint CMW::fistr_elemtype_prism(){ return FistrElementType::Prism;}
uint CMW::fistr_elemtype_prism2(){ return FistrElementType::Prism2;}
uint CMW::fistr_elemtype_quad(){ return FistrElementType::Quad;}
uint CMW::fistr_elemtype_quad2(){ return FistrElementType::Quad2;}
uint CMW::fistr_elemtype_triangle(){ return FistrElementType::Triangle;}
uint CMW::fistr_elemtype_triangle2(){ return FistrElementType::Triangle2;}
uint CMW::fistr_elemtype_line(){ return FistrElementType::Beam;}
uint CMW::fistr_elemtype_line2(){ return FistrElementType::Beam2;}
//----
// FrontISTR 要素タイプ　=> MW3 要素タイプ 変換
//----
uint CMW::fistr_elemtype_to_mw3_elemtype(const uint& fistr_elemtype)
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
uint CMW::mw3_elemtype_to_fistr_elemtype(const uint& mw3_elemtype)
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
uint CMW::shapetype_hexa81(){ return ShapeType::Hexa81;}
uint CMW::shapetype_hexa82(){ return ShapeType::Hexa82;}
uint CMW::shapetype_hexa201(){ return ShapeType::Hexa201;}
uint CMW::shapetype_hexa202(){ return ShapeType::Hexa202;}
uint CMW::shapetype_hexa203(){ return ShapeType::Hexa203;}
uint CMW::shapetype_tetra41(){ return ShapeType::Tetra41;}
uint CMW::shapetype_tetra101(){ return ShapeType::Tetra101;}
uint CMW::shapetype_tetra104(){ return ShapeType::Tetra104;}
uint CMW::shapetype_tetra1015(){ return ShapeType::Tetra1015;}
uint CMW::shapetype_prism62(){ return ShapeType::Prism62;}
uint CMW::shapetype_prism156(){ return ShapeType::Prism156;}
uint CMW::shapetype_prism159(){ return ShapeType::Prism159;}
uint CMW::shapetype_prism1518(){ return ShapeType::Prism1518;}
uint CMW::shapetype_quad41(){ return ShapeType::Quad41;}
uint CMW::shapetype_quad84(){ return ShapeType::Quad84;}
uint CMW::shapetype_quad89(){ return ShapeType::Quad89;}
uint CMW::shapetype_tri31(){ return ShapeType::Triangle31;}
uint CMW::shapetype_tri63(){ return ShapeType::Triangle63;}
uint CMW::shapetype_line21(){ return ShapeType::Line21;}
uint CMW::shapetype_line32(){ return ShapeType::Line32;}


//--
// Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
//--
uint CMW::GetNumOfBoundaryNodeMesh(){ return mpMesh->getNumOfBoundaryNodeMesh();}
uint CMW::GetNumOfBoundaryFaceMesh(){ return mpMesh->getNumOfBoundaryFaceMesh();}
uint CMW::GetNumOfBoundaryEdgeMesh(){ return mpMesh->getNumOfBoundaryEdgeMesh();}
uint CMW::GetNumOfBoundaryVolumeMesh(){ return mpMesh->getNumOfBoundaryVolumeMesh();}
// BNode
uint CMW::GetNumOfBNode_BNodeMesh(const uint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    
    return pBNodeMesh->getNumOfBNode();
}
uint CMW::GetNumOfBNode_BFaceMesh(const uint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    
    return pBFaceMesh->getNumOfBNode();
}
uint CMW::GetNumOfBNode_BEdgeMesh(const uint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    
    return pBEdgeMesh->getNumOfBNode();
}
uint CMW::GetNumOfBNode_BVolumeMesh(const uint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    
    return pBVolMesh->getNumOfBNode();
}
// DOF
uint CMW::GetNumOfDOF_BNodeMesh(const uint& ibmesh, const uint& ibnode)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    
    return pBNode->getNumOfDOF();
}
uint CMW::GetNumOfDOF_BFaceMesh(const uint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    
    return pBFaceMesh->getNumOfDOF();
}
uint CMW::GetNumOfDOF_BEdgeMesh(const uint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    
    return pBEdgeMesh->getNumOfDOF();
}
uint CMW::GetNumOfDOF_BVolumeMesh(const uint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    
    return pBVolMesh->getNumOfDOF();
}
//--
// BoundaryNode 境界値
//--
double& CMW::GetBNodeValue_BNodeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    
    return pBNode->getValue(idof);
}
double& CMW::GetBNodeValue_BFaceMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);
    
    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);
    
    return pBNode->getValue(idof, mgLevel);
}
double& CMW::GetBNodeValue_BEdgeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);
    
    CBoundaryNode *pBNode;
    pBNode = pBEdgeMesh->getBNodeIX(ibnode);
    
    return pBNode->getValue(idof, mgLevel);
}
double& CMW::GetBNodeValue_BVolumeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    
    CBoundaryNode *pBNode;
    pBNode = pBVolMesh->getBNodeIX(ibnode);
    
    return pBNode->getValue(idof, mgLevel);
}
//--
// Boundary Node の NodeID
//--
uint& CMW::GetNodeID_BNode_BNodeMesh(const uint& ibmesh, const uint& ibnode)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);
    
    CBoundarySBNode *pBNode;
    pBNode = pBNodeMesh->getBNodeIX(ibnode);
    
    CNode *pNode;
    pNode = pBNode->getNode();
    
    return pNode->getID();
}
uint& CMW::GetNodeID_BNode_BFaceMesh(const uint& ibmesh, const uint& ibnode)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);

    CBoundaryNode *pBNode;
    pBNode = pBFaceMesh->getBNodeIX(ibnode);

    CNode *pNode;
    pNode = pBNode->getNode();

    return pNode->getID();
}
uint& CMW::GetNodeID_BNode_BEdgeMesh(const uint& ibmesh, const uint& ibnode)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);

    CBoundaryNode *pBNode;
    pBNode = pBEdgeMesh->getBNodeIX(ibnode);

    CNode *pNode;
    pNode = pBNode->getNode();

    return pNode->getID();
}
uint& CMW::GetNodeID_BNode_BVolumeMesh(const uint& ibmesh, const uint& ibnode)
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
uint CMW::GetNumOfBFace(const uint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);

    return pBFaceMesh->getNumOfBFace();
}
double& CMW::GetBFaceValue(const uint& ibmesh, const uint& ibface, const uint& idof)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);

    CBoundaryFace *pBFace;
    pBFace = pBFaceMesh->getBFaceIX(ibface);

    return pBFace->getBndValue(idof);
}
uint CMW::GetNumOfBEdge(const uint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);

    return pBEdgeMesh->getNumOfEdge();
}
double& CMW::GetBEdgeValue(const uint& ibmesh, const uint& ibedge, const uint& idof)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);

    CBoundaryEdge *pBEdge;
    pBEdge = pBEdgeMesh->getBEdgeIX(ibedge);

    return pBEdge->getBndValue(idof);
}
uint CMW::GetNumOfBVolume(const uint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);

    return pBVolMesh->getNumOfVolume();
}
double& CMW::GetBVolumeValue(const uint& ibmesh, const uint& ibvol, const uint& idof)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);

    CBoundaryVolume *pBVol;
    pBVol = pBVolMesh->getBVolumeIX(ibvol);

    return pBVol->getBndValue(idof);
}
//--
// Boundaryの名称
//--
uint CMW::GetBNodeMesh_NameLength(const uint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);

    string sBndName = pBNodeMesh->getName();

    return sBndName.length();
}
string& CMW::GetBNodeMesh_Name(const uint& ibmesh)
{
    CBoundaryNodeMesh *pBNodeMesh;
    pBNodeMesh = mpMesh->getBndNodeMeshIX(ibmesh);

    return pBNodeMesh->getName();
}
uint CMW::GetBFaceMesh_NameLength(const uint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);

    string sBndName = pBFaceMesh->getName();

    return sBndName.length();
}
string& CMW::GetBFaceMesh_Name(const uint& ibmesh)
{
    CBoundaryFaceMesh *pBFaceMesh;
    pBFaceMesh = mpMesh->getBndFaceMeshIX(ibmesh);

    return pBFaceMesh->getName();
}
uint CMW::GetBVolumeMesh_NameLength(const uint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    
    string sBndName = pBVolMesh->getName();
    
    return sBndName.length();
}
string& CMW::GetBVolumeMesh_Name(const uint& ibmesh)
{
    CBoundaryVolumeMesh *pBVolMesh;
    pBVolMesh = mpMesh->getBndVolumeMeshIX(ibmesh);
    
    return pBVolMesh->getName();
}
uint CMW::GetBEdgeMesh_NameLength(const uint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);

    string sBndName = pBEdgeMesh->getName();

    return sBndName.length();
}
string& CMW::GetBEdgeMesh_Name(const uint& ibmesh)
{
    CBoundaryEdgeMesh *pBEdgeMesh;
    pBEdgeMesh = mpMesh->getBndEdgeMeshIX(ibmesh);

    return pBEdgeMesh->getName();
}





//--
// MPI (直接使用)
//--
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
int& CMW::GetRank()//自分のプロセス-ランクを取得
{
    return mpMPI->getRank();
}
int& CMW::GetNumOfProcess()
{
    return mpMPI->getNumOfProcess();
}

// bufの値を送信 => 受信した値をbufに代入、Nodeにセット
// # bufのサイズ == NumOfCommNode * dof_size
//
void CMW::Send_Recv_R(double* buf, int dof_size)
{
    uint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
    CAssyModel *pAssy;
    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
        pAssy = mpGMGModel->getAssyModel(iLevel);

        uint imesh, nNumOfMesh = pAssy->getNumOfMesh();
        CMesh *pMesh;
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            pMesh = pAssy->getMesh(imesh);

            uint icmesh, nNumOfCommMesh = pMesh->getNumOfCommMesh();
            CCommMesh2 *pCommMesh;
            for(icmesh=0; icmesh < nNumOfCommMesh; icmesh++){
                pCommMesh = pMesh->getCommMesh2IX(icmesh);

                int transRank = pCommMesh->getTrasmitRank();

                uint icnode, nNumOfCommNode = pCommMesh->getCommNodeSize();
                CCommNode *pCommNode;
                for(icnode=0; icnode < nNumOfCommNode; icnode++){
                    pCommNode = pCommMesh->getCommNodeIX(icnode);
                    
                    double sendbuf[dof_size];
                    double recvbuf[dof_size];

                    uint idof;
                    // 送信buf へ代入
                    for(idof=0; idof < dof_size; idof++) sendbuf[idof] = buf[icnode*dof_size + idof];

                    // 送信
                    mpMPI->Send(sendbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD);


                    // 受信
                    MPI_Status sta;
                    mpMPI->Recv(recvbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD, &sta);
                    
                    CNode* pNode = pCommNode->getNode();
                    //
                    // ・受信bufから引数bufへ代入
                    // ・受信値をNodeの変数に代入
                    //
                    for(idof=0; idof < dof_size; idof++){
                        
                        buf[icnode*dof_size + idof] = recvbuf[idof];

                        if(pNode->getType()==NodeType::Scalar) pNode->setScalar(recvbuf[idof], idof);
                        if(pNode->getType()==NodeType::Vector) pNode->setVector(recvbuf[idof], idof);
                        if(pNode->getType()==NodeType::ScalarVector){
                            uint nNumScalar = pNode->numOfScalarParam();
                            if(idof < nNumScalar){
                                pNode->setScalar(recvbuf[idof], idof);
                            }else{
                                pNode->setVector(recvbuf[idof-nNumScalar], idof);
                            }
                        }
                    };// idof (Total DOF)
                };// icnode (CommNode)
            };// icmesh (CommMesh)
        };// imesh (Mesh)
    };// iLevel
}

//
// 通信Nodeの値を入れ替えて更新
//
void CMW::Send_Recv_R()
{
    uint iLevel, nNumOfLevel = mpGMGModel->getNumOfLevel();
    CAssyModel *pAssy;
    for(iLevel=0; iLevel < nNumOfLevel; iLevel++){
        pAssy = mpGMGModel->getAssyModel(iLevel);

        uint imesh, nNumOfMesh = pAssy->getNumOfMesh();
        CMesh *pMesh;
        for(imesh=0; imesh < nNumOfMesh; imesh++){
            pMesh = pAssy->getMesh(imesh);

            uint icmesh, nNumOfCommMesh = pMesh->getNumOfCommMesh();
            CCommMesh2 *pCommMesh;
            for(icmesh=0; icmesh < nNumOfCommMesh; icmesh++){
                pCommMesh = pMesh->getCommMesh2IX(icmesh);

                int transRank = pCommMesh->getTrasmitRank();

                uint icnode, nNumOfCommNode = pCommMesh->getCommNodeSize();
                CCommNode *pCommNode;
                for(icnode=0; icnode < nNumOfCommNode; icnode++){
                    pCommNode = pCommMesh->getCommNodeIX(icnode);
                    CNode* pNode = pCommNode->getNode();

                    uint dof_size = pNode->numOfTotalParam();
                    double sendbuf[dof_size];
                    double recvbuf[dof_size];

                    uint idof;
                    // 送信bufにNodeの変数値を代入
                    //
                    for(idof=0; idof < dof_size; idof++){
                        if(pNode->getType()==NodeType::Scalar) sendbuf[idof] = pNode->getScalar(idof);
                        if(pNode->getType()==NodeType::Vector) sendbuf[idof] = pNode->getVector(idof);
                        if(pNode->getType()==NodeType::ScalarVector){
                            uint nNumScalar = pNode->numOfScalarParam();
                            if(idof < nNumScalar){
                                sendbuf[idof] = pNode->getScalar(idof);
                            }else{
                                sendbuf[idof] = pNode->getVector(idof-nNumScalar);
                            }
                        }
                    }

                    // 送信
                    mpMPI->Send(sendbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD);


                    // 受信
                    MPI_Status sta;
                    mpMPI->Recv(recvbuf, dof_size, MPI_DOUBLE, transRank, 0, MPI_COMM_WORLD, &sta);

                    // 受信値をNodeの変数に代入
                    //
                    for(idof=0; idof < dof_size; idof++){
                        if(pNode->getType()==NodeType::Scalar) pNode->setScalar(recvbuf[idof], idof);
                        if(pNode->getType()==NodeType::Vector) pNode->setVector(recvbuf[idof], idof);
                        if(pNode->getType()==NodeType::ScalarVector){
                            uint nNumScalar = pNode->numOfScalarParam();
                            if(idof < nNumScalar){
                                pNode->setScalar(recvbuf[idof], idof);
                            }else{
                                pNode->setVector(recvbuf[idof-nNumScalar], idof);
                            }
                        }
                    };// idof (Total DOF)
                };// icnode (CommNode)
            };// icmesh (CommMesh)
        };// imesh (Mesh)
    };// iLevel
}


//--
// グループ  { select された AssyModel,Meshを対象 }
//--
uint CMW::GetNumOfElementGroup()
{
    return mpMesh->getNumOfElemGrp();
}
uint CMW::GetNumOfElementID(const uint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);

    return pElemGrp->getNumOfElementID();
}
uint& CMW::GetElementID_with_ElementGroup(const uint& iGrp, const uint& index)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);

    return pElemGrp->getElementID(index);
}
uint CMW::GetElementGroupName_Length(const uint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);

    return pElemGrp->getNameLength();
}
string& CMW::GetElementGroupName(const uint& iGrp)
{
    CElementGroup *pElemGrp = mpMesh->getElemGrpIX(iGrp);

    return pElemGrp->getName();
}


//--
// Logger
//--
void CMW::LoggerMode(const uint& mode)
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
void CMW::LoggerDevice(const uint& mode, const uint& device)
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
void CMW::LoggerInfo(const uint& mode, char* message)
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
uint CMW::getErrorMode()
{
    return Utility::LoggerMode::Error;
}
uint CMW::getWarnMode()
{
    return Utility::LoggerMode::Warn;
}
uint CMW::getInfoMode()
{
    return Utility::LoggerMode::Info;
}
uint CMW::getDebugMode()
{
    return Utility::LoggerMode::Debug;
}
uint CMW::getDiskDevice()
{
    return Utility::LoggerDevice::Disk;
}
uint CMW::getDisplayDevice()
{
    return Utility::LoggerDevice::Display;
}


