/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   HEC_MW3.hxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#define PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#include "TypeDef.h"
#include "HEC_MPI.h"
#include "FileIO.h"
#include "Logger.h"
#include "MeshFactory.h"
#include "GMGModel.h"
#include "AssyMatrix.h"
#include "ShapeHexa.h"
#include "ShapeHexaNic.h"
#include "ShapeTetra.h"
#include "ShapePrism.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeLine.h"
#include "ShapeFunctionCatalog.h"
#include "ISTR2Edge.h"
#include "Edge2ISTR.h"
#include "Jacobian.h"
#include <boost/lexical_cast.hpp>
typedef pair<uint,uint> integPair;
namespace pmw{
class CMW
{
public:
    __declspec(dllexport) static CMW*  Instance(){
	static CMW moMW;
	return &moMW;
    }
private:
    __declspec(dllexport) CMW(void);
public:
    __declspec(dllexport) virtual ~CMW(void);
protected:
    FileIO::CFileIO   *mpFileIO;
    Utility::CLogger  *mpLogger;    
    bool mb_file;
    string msInputFileName;
    string msOutputFileName;    
    pmw::CGMGModel    *mpGMGModel;
    pmw::CMeshFactory  *mpFactory;
    pmw::CHecMPI      *mpMPI;
    pmw::CAssyModel  *mpAssy;
    pmw::CAssyMatrix *mpAssyMatrix;
    pmw::CAssyVector *mpAssyVector;
    pmw::CAssyVector *mpAssyVector2;
    pmw::CMesh       *mpMesh;
    pmw::CElement    *mpElement;
    pmw::CShapeHexa     *mpShapeHexa;
    pmw::CShapeHexaNic  *mpShapeHexaNic;
    pmw::CShapeTetra    *mpShapeTetra;
    pmw::CShapePrism    *mpShapePrism;
    pmw::CShapeQuad     *mpShapeQuad;
    pmw::CShapeTriangle *mpShapeTriangle;
    pmw::CShapeLine     *mpShapeLine;
    pmw::CShapeFunctionCatalog *mpShapeCatalog;
    pmw::CISTR2Edge *mpISTR2Edge;
    pmw::CEdge2ISTR *mpEdge2ISTR;
    vvvdouble mvdNdx;
public:
    int __declspec(dllexport) Initialize(int argc, char** argv, const char* path);
    int __declspec(dllexport) Finalize();  
    int __declspec(dllexport) FileRead(); 
    int __declspec(dllexport) FileWrite();
    int __declspec(dllexport) Initialize_Matrix(); 
    int __declspec(dllexport) Initialize_Vector(); 
    int __declspec(dllexport) Matrix_Add_Elem(int iMesh, int iElem, double *ElemMatrix);
    int __declspec(dllexport) Set_BC(int iMesh, int iNode, int iDof, double value1, double value2);
    int __declspec(dllexport) Set_BC(int iMesh, int iNode, int iDof, double value);
	void __declspec(dllexport) Sample_Set_BC(int iMesh);
    int __declspec(dllexport) Solve(uint iter_max, double tolerance, uint method, uint precondition);
    int __declspec(dllexport) Refine();
    void __declspec(dllexport) SelectAssyModel(const uint& mgLevel);
    void __declspec(dllexport) SelectMesh_ID(const uint& mesh_id);
    void __declspec(dllexport) SelectMesh_IX(const uint& index);
    uint  __declspec(dllexport) getNodeSize(){ return mpMesh->getNodeSize();}
    uint  __declspec(dllexport) getElementSize(){ return mpMesh->getElementSize();}
    uint  __declspec(dllexport) getNodeSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    uint  __declspec(dllexport) getElementSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
	void __declspec(dllexport) StoreMatrix(){
		mpAssy->setAssyMatrix(mpAssyMatrix);
		mpAssy->setAssyVector(mpAssyVector);
		mpAssy->setAssyVector2(mpAssyVector2);
	}
	void __declspec(dllexport) LoadMatrix(){
		mpAssyMatrix = mpAssy->getAssyMatrix();
		mpAssyVector = mpAssy->getAssyVector();
		mpAssyVector2 = mpAssy->getAssyVector2();
	}
    void __declspec(dllexport) GetNumOfAssembleModel(uint& numOfAssembleModel);
    void __declspec(dllexport) SelectAssembleModel(const uint& mgLevel);
    void __declspec(dllexport) GetNumOfMeshPart(uint& numOfMeshPart);
    void __declspec(dllexport) SelectMeshPart_ID(const uint& mesh_id);
    void __declspec(dllexport) SelectMeshPart_IX(const uint& index);
    void __declspec(dllexport) SelectElement_ID(const uint& elem_id);
    void __declspec(dllexport) SelectElement_IX(const uint& index);
    void __declspec(dllexport) GetElementType(uint& elemType);
    void __declspec(dllexport) GetNumOfElementVert(uint& numOfVert);
    void __declspec(dllexport) GetElementVertNodeID(vuint& vNodeID);
    void __declspec(dllexport) GetNumOfElementEdge(uint& numOfEdge);
    void __declspec(dllexport) GetElementEdgeNodeID(vuint& vNodeID);
    void __declspec(dllexport) GetNodeCoord(const uint& node_id, double& x, double& y, double& z);
    void __declspec(dllexport) NumOfIntegPoint(const uint& shapeType, uint& numOfInteg);
    void __declspec(dllexport) ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N);
    void __declspec(dllexport) ShapeFunc(const uint& shapeType, vvdouble& N);
    void __declspec(dllexport) dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr);
    void __declspec(dllexport) dNdr(const uint& shapeType, vvvdouble& dNdr);
    void __declspec(dllexport) Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index);
    void __declspec(dllexport) dNdx_on_pt(const uint& igauss, vvdouble& dNdx);
    void __declspec(dllexport) dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdx);
    void __declspec(dllexport) detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ);
    void __declspec(dllexport) Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w);
};
}
#endif
