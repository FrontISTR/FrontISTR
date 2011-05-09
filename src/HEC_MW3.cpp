/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   HEC_MW3.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
#ifdef WIN32
#include "HEC_MW3.hxx"
#else
#include "HEC_MW3.h"
#endif
using namespace pmw;
CMW::CMW(void)
{
    mpGMGModel= CGMGModel::Instance();
    mpFactory = CMeshFactory::Instance();
    mpFactory->setGMGModel(mpGMGModel);
    mpFileIO = FileIO::CFileIO::Instance();
    mb_file = false;
    mpFileIO->setFactory(mpFactory);
    mpMPI = CHecMPI::Instance();
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
}
CMW::~CMW(void)
{
}
int CMW::Initialize(int argc, char** argv, const char* path)
{
    mpMPI->Initialize(argc, argv);
    uint rank= mpMPI->getRank();
	mpLogger->setMode(Utility::LoggerMode::Info);
    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info, Utility::LoggerDevice::Display);
    mpLogger->initializeLogFile(rank);
    mpLogger->InfoDisplay();
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");
    mpFileIO->setPathName(path);
    mpFileIO->ReadCntFile();
    msInputFileName= mpFileIO->getPathName();
    msOutputFileName=mpFileIO->getPathName();
    msInputFileName += mpFileIO->getMeshFileBaseName() + "_";
    msOutputFileName+= mpFileIO->getMeshFileBaseName() + "_";
    msInputFileName  += boost::lexical_cast<string>(rank);
    msOutputFileName += boost::lexical_cast<string>(rank);
    msInputFileName  += ".msh";
    msOutputFileName += ".out";
    return 1;
}
int CMW::Finalize()
{
    mpMPI->Finalize();
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();
    return 1;
}
int CMW::FileRead()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");
    mpFileIO->ReadFile(msInputFileName);
    mb_file = true;
    return 1;
}
int CMW::FileWrite()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileWrite");
    uint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile(msOutputFileName, nLevel+1);
    return 1;
}
int CMW::Initialize_Matrix()
{
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CMW::Initialize_Matrix \n");
#endif
	mpAssyMatrix = new CAssyMatrix(mpAssy);
	uint level = mpAssy->getMGLevel();
	mpAssyMatrix->setCoarseMatrix( NULL );
	if( level > 0 ) mpAssyMatrix->setCoarseMatrix( mpGMGModel->getAssyModel(level-1)->getAssyMatrix() );
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CMW::Initialize_Matrix \n");
#endif
    return 1;
}
int CMW::Initialize_Vector()
{
#ifdef ADVANCESOFT_DEBUG
   	printf(" enter CMW::Initialize_Vector \n");
#endif
	mpAssyVector = new CAssyVector(mpAssy);
	mpAssyVector2 = new CAssyVector(mpAssy);
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CMW::Initialize_Vector \n");
#endif
    return 1;
}
int CMW::Matrix_Add_Elem(int iMesh, int iElem, double *ElemMatrix)
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
void CMW::Sample_Set_BC(int iMesh)
{
#ifdef ADVANCESOFT_DEBUG
	printf(" enter CMW::Sample_Set_BC \n");
#endif
	double X, Y, Z, value0, value1, value2;
	int iDof0, iDof1, iDof2;
	int iNodeMax = getNodeSize( iMesh );
	for( int iNode = 0; iNode < iNodeMax; iNode++){  
		X = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getX();
		Y = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getY();
		Z = mpAssy->getMesh(iMesh)->getNodeIX(iNode)->getZ();
		if(abs( Z - 4.0 ) < 1.0e-5 && abs( X - 1.0 ) < 1.0e-5 ) {
			value0 = 1.0e6;
			iDof0 = 0;
			Set_BC( iMesh, iNode, iDof0, value0);
		};
		if( (abs( Z ) < 1.0e-5) || (abs( Z - 8.0 ) < 1.0e-5) ) {
			value1 = 1.0e15;
			value2 = 0.0;
			iDof0 = 0; iDof1 = 1; iDof2 = 2;
			Set_BC( iMesh, iNode, iDof0, value1, value2);
			Set_BC( iMesh, iNode, iDof1, value1, value2);
			Set_BC( iMesh, iNode, iDof2, value1, value2);
		}
	};
#ifdef ADVANCESOFT_DEBUG
 	printf(" enter CMW::Sample_Set_BC \n");
#endif
}
int CMW::Set_BC(int iMesh, int iNode, int iDof, double value)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMW::Set_BC (1) %d %d %d %e \n", iMesh, iNode, iDof, value);
#endif
	mpAssyVector->setValue(iMesh, iNode, iDof, value);
#ifdef ADVANCESOFT_DEBUG
	printf("exit CMW::Set_BC (1) \n");
#endif
	return 1;
}
int CMW::Set_BC(int iMesh, int iNode, int iDof, double value1, double value2)
{
#ifdef ADVANCESOFT_DEBUG
	printf("enter CMW::Set_BC (2) %d %d %d %e %e \n", iMesh, iNode, iDof, value1, value2);
#endif
	mpAssyMatrix->setValue(iMesh, iNode, iDof, value1);
	mpAssyVector->setValue(iMesh, iNode, iDof, value2);
#ifdef ADVANCESOFT_DEBUG
	printf("exit CMW::Set_BC (2) \n");
#endif
	return 1;
}
int CMW::Solve(uint iter_max, double tolerance, uint method, uint precondition)
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
			solver->solve(mpAssyMatrix, mpAssyVector, mpAssyVector2);}
            break;
        case( 2 ):
			{CSolverBiCGSTAB *solver = new CSolverBiCGSTAB(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
			solver->solve(mpAssyMatrix, mpAssyVector, mpAssyVector2);}
            break;
        case( 3 ):
			{CSolverGPBiCG *solver = new CSolverGPBiCG(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
			solver->solve(mpAssyMatrix, mpAssyVector, mpAssyVector2);}
            break;
        case( 4 ):
			{CSolverGMRES *solver = new CSolverGMRES(iter_max, tolerance, method, precondition, flag_iter_log, flag_time_log);
			solver->solve(mpAssyMatrix, mpAssyVector, mpAssyVector2);}
            break;
        default:
            break;
    }
	uint iMesh, iMeshMax;
	GetNumOfMeshPart(iMeshMax);
	for( iMesh = 0; iMesh < iMeshMax; iMesh++){
	FILE *fp1;
	sprintf(cfile, "outUCD%02d%02d_%05d.inp", iMesh, mpAssy->getMGLevel(), mpMPI->getRank()); 
	fp1 = fopen(cfile, "w");
	fprintf(fp1,"%d %d 3 0 0 \n",mpAssy->getMesh(iMesh)->getNumOfNode(), mpAssy->getMesh(iMesh)->getNumOfElement());
	for(int i=0; i < mpAssy->getMesh(iMesh)->getNumOfNode(); i++)
	{
		fprintf(fp1,"%d %e %e %e\n",i, mpAssy->getMesh(iMesh)->getNodeIX(i)->getX(),
		mpAssy->getMesh(iMesh)->getNodeIX(i)->getY(), mpAssy->getMesh(iMesh)->getNodeIX(i)->getZ());
	}
	for(int i=0; i < mpAssy->getMesh(0)->getNumOfElement(); i++)
	{
		fprintf(fp1,"%d 0 hex %d %d %d %d %d %d %d %d \n",i,
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(0)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(1)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(2)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(3)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(4)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(5)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(6)->getID(),
		mpAssy->getMesh(iMesh)->getElementIX(i)->getNode(7)->getID());
	}
	fprintf(fp1,"1 3 \n");
	fprintf(fp1,"vector, aaa \n");
	for(int i=0; i < mpAssy->getMesh(iMesh)->getNumOfNode(); i++)
	{
		double val0 = mpAssyVector2->getValue(iMesh, i, 0);
		double val1 = mpAssyVector2->getValue(iMesh, i, 1);
		double val2 = mpAssyVector2->getValue(iMesh, i, 2);
		fprintf(fp1,"%d %e %e %e\n",i,val0, val1, val2);
		mpAssy->getMesh(iMesh)->getNodeIX(i)->setVector(val0, 0);
		mpAssy->getMesh(iMesh)->getNodeIX(i)->setVector(val1, 1);
		mpAssy->getMesh(iMesh)->getNodeIX(i)->setVector(val2, 2);
	}
	fclose(fp1);
	}
   	printf(" --- end of output to a UCD file. (%s) --- \n", cfile);
#ifdef ADVANCESOFT_DEBUG
   	printf(" exit CMW::Solve \n");
#endif
    return 1;
}
int CMW::Refine()
{
    uint ilevel,mgLevel;
    CAssyModel *pAssy;
    uint icon,numOfConMesh;
    CContactMesh *pConMesh;
    if(mb_file){
        mpFactory->refineMesh();       
        mpFactory->refineContactMesh();
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
        mpFactory->refineCommMesh2();
        return 1;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error," Not read file @MWMain::Refine()");
        return 0;
    }
    return 1;
}
void CMW::GetNumOfAssembleModel(uint& numOfAssembleModel)
{
    numOfAssembleModel= mpGMGModel->getNumOfLevel();
}
void CMW::SelectAssembleModel(const uint& mgLevel)
{
    mpAssy= mpGMGModel->getAssyModel(mgLevel);
}
void CMW::GetNumOfMeshPart(uint& numOfMeshPart)
{
    numOfMeshPart= mpAssy->getNumOfMesh();
}
void CMW::SelectMeshPart_ID(const uint& mesh_id)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh_ID(mesh_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
void CMW::SelectMeshPart_IX(const uint& index)
{
    if(mpAssy){
        mpMesh= mpAssy->getMesh(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"AssyModel pointer Null!, use SelectAssyModel");
    }
}
void CMW::SelectElement_ID(const uint& elem_id)
{
    if(mpMesh){
        mpElement= mpMesh->getElement(elem_id);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
void CMW::SelectElement_IX(const uint& index)
{
    if(mpMesh){
        mpElement= mpMesh->getElementIX(index);
    }else{
        mpLogger->Info(Utility::LoggerMode::Error,"Mesh Pointer Null!, use SelectMesh");
    }
}
void CMW::GetElementType(uint& elemType)
{
    elemType= mpElement->getType();
}
void CMW::GetNumOfElementVert(uint& numOfVert)
{
    numOfVert= mpElement->getNumOfNode();
}
void CMW::GetElementVertNodeID(vuint& vNodeID)
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
void CMW::GetNumOfElementEdge(uint& numOfEdge)
{
    numOfEdge= mpElement->getNumOfEdge();
}
void CMW::GetElementEdgeNodeID(vuint& vNodeID)
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
void CMW::GetNodeCoord(const uint& node_id, double& x, double& y, double& z)
{
    CNode *pNode;
    pNode= mpMesh->getNode(node_id);
    x= pNode->getX();
    y= pNode->getY();
    z= pNode->getZ();
}
void CMW::NumOfIntegPoint(const uint& shapeType, uint& numOfInteg)
{
    numOfInteg= mpShapeCatalog->NumOfIntegPoint(shapeType);
}
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
        case(ShapeType::HexaNic111):
            dNdr= mpShapeHexaNic->dNdr111();
            break;
        case(ShapeType::HexaNic118):
            dNdr= mpShapeHexaNic->dNdr118();
            break;
        case(ShapeType::HexaNic1127):
            dNdr= mpShapeHexaNic->dNdr1127();
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
void CMW::Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index)
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
void CMW::dNdx_on_pt(const uint& igauss, vvdouble& dNdX)
{
    dNdX = mvdNdx[igauss];
}
void CMW::dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdX)
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
void CMW::detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ)
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
void CMW::Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w)
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
/*
    switch(elemType){
        case(ElementType::Hexa):
            if(numOfInteg==1){
                ;
            }
            if(numOfInteg==8){
                ;
            }
            break;
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
        case(ElementType::Tetra):
            if(numOfInteg==1){
                ;
            }
            break;
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
