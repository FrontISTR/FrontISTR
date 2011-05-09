/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   HEC_MW3.h
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
    static CMW* Instance(){
	static CMW moMW;
	return &moMW;
    }
private:
    CMW(void);
public:
    virtual ~CMW(void);
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
    int Initialize(int argc, char** argv, const char* path);
    int Finalize();  
    int FileRead(); 
    int FileWrite();
    int Initialize_Matrix(); 
    int Initialize_Vector(); 
    int Matrix_Add_Elem(int iMesh, int iElem, double *ElemMatrix);
    int Set_BC(int iMesh, int iNode, int iDof, double value1, double value2);
    int Set_BC(int iMesh, int iNode, int iDof, double value);
	void Sample_Set_BC(int iMesh);
    int Solve(uint iter_max, double tolerance, uint method, uint precondition);
    int Refine();
    void SelectAssyModel(const uint& mgLevel);
    void SelectMesh_ID(const uint& mesh_id);
    void SelectMesh_IX(const uint& index);
    uint  getNodeSize(){ return mpMesh->getNodeSize();}
    uint  getElementSize(){ return mpMesh->getElementSize();}
    uint  getNodeSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    uint  getElementSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
	void StoreMatrix(){
		mpAssy->setAssyMatrix(mpAssyMatrix);
		mpAssy->setAssyVector(mpAssyVector);
		mpAssy->setAssyVector2(mpAssyVector2);
	}
	void LoadMatrix(){
		mpAssyMatrix = mpAssy->getAssyMatrix();
		mpAssyVector = mpAssy->getAssyVector();
		mpAssyVector2 = mpAssy->getAssyVector2();
	}
    void GetNumOfAssembleModel(uint& numOfAssembleModel);
    void SelectAssembleModel(const uint& mgLevel);
    void GetNumOfMeshPart(uint& numOfMeshPart);
    void SelectMeshPart_ID(const uint& mesh_id);
    void SelectMeshPart_IX(const uint& index);
    void SelectElement_ID(const uint& elem_id);
    void SelectElement_IX(const uint& index);
    void GetElementType(uint& elemType);
    void GetNumOfElementVert(uint& numOfVert);
    void GetElementVertNodeID(vuint& vNodeID);
    void GetNumOfElementEdge(uint& numOfEdge);
    void GetElementEdgeNodeID(vuint& vNodeID);
    void GetNodeCoord(const uint& node_id, double& x, double& y, double& z);
    void NumOfIntegPoint(const uint& shapeType, uint& numOfInteg);
    void ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N);
    void ShapeFunc(const uint& shapeType, vvdouble& N);
    void dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr);
    void dNdr(const uint& shapeType, vvvdouble& dNdr);
    void Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index);
    void dNdx_on_pt(const uint& igauss, vvdouble& dNdx);
    void dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdx);
    void detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ);
    void Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w);
};
}
#endif
