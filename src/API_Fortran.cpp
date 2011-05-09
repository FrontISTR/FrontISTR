//
//  API_Fortran.cpp
//
//                  2009.4.20
//                  2009.2.27
//                  k.Takeda
#ifdef MSVC
#include "API_Fortran.hxx"
#else
#include "API_Fortran.h"
#endif

#include "HEC_MW3.h"

pmw::CMW *pMW;
//----
// HEC_MW3 construct & destruct
//----
#include <string.h>
#include <math.h>
int mw_initialize_(int* argc, char** argv, char* path) // for C
{
    pMW = pmw::CMW::Instance();

    return pMW->Initialize(*argc, argv, path);// argc, argv <= MPI_Init 引数
}
int mw_initialize_1_(char* argv1, int* argv1_len, char* path, int* path_len)
{
    int argc = 1;// 1: exe program name
    char** argv;

    size_t nLength = (size_t)*argv1_len + 1;
    char carg[nLength];
    strncpy(carg, argv1, nLength);// 1: exe program name

    argv = new char*[1];
    argv[0] = carg;

    nLength = (size_t)*path_len + 1;
    char cpath[nLength];
    strncpy(cpath, path, nLength);// cnt file path

    return pMW->Initialize(argc, argv, cpath);
}
int mw_initialize_2_(char* argv1, int* argv1_len, char* argv2, int* argv2_len, char* path, int* path_len)
{
    int argc = 2;// 1: exe program name, 2: input file name
    char** argv;

    size_t nLength = (size_t)*argv1_len + 1;
    char carg1[nLength];
    strncpy(carg1, argv1, nLength);// 1: exe program name

    nLength = (size_t)*argv2_len + 1;
    char carg2[nLength];
    strncpy(carg2, argv2, nLength);// 2: input file name

    argv = new char*[2];
    argv[0] = carg1;
    argv[1] = carg2;

    nLength = (size_t)*path_len + 1;
    char cpath[nLength];
    strncpy(cpath, path, nLength);// cnt file path

    return pMW->Initialize(argc, argv, cpath);
}
int mw_finalize_()
{
    return pMW->Finalize();
}

//----
// file i/o API
//----
int mw_file_read_()
{
    return pMW->FileRead();
}
int mw_file_write_()
{
    return pMW->FileWrite();
}

//----
// linear solver API
//----
//int mw_initialize_matrix_()
//{
//    return pMW->Initialize_Matrix();
//}
//int mw_initialize_vector_()
//{
//    return pMW->Initialize_Vector();
//}
void mw_gene_linear_algebra_(int* num_of_algebra, int dof[])
{
    uint nNumOfAlgebra = (uint)*num_of_algebra;
    uint vDOF[nNumOfAlgebra];
    for(uint i=0; i < nNumOfAlgebra; i++) vDOF[i] = dof[i];

    pMW->GeneLinearAlgebra(nNumOfAlgebra, vDOF);
}
void mw_select_algebra_(int* ieq)
{
    uint iEqu = *ieq;
    pMW->SelectAlgebra(iEqu);
}

// matrix add elem
int mw_matrix_add_elem_(int* imesh,  int* ielem,  double elem_matrix[])//standard
{
    unsigned int iMesh = *imesh;
    unsigned int iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, elem_matrix);
}

int mw_matrix_add_node_(int* imesh, int* i_nid, int* j_nid, double nodal_matrix[])
{
    uint iMesh = *imesh;
    uint iNodeID = *i_nid;
    uint jNodeID = *j_nid;

    return pMW->Matrix_Add_Node(iMesh, iNodeID, jNodeID, nodal_matrix);
}

// matrix 0 clear
//
void mw_matrix_clear_(int* imesh)
{
    unsigned int iMesh = *imesh;
    pMW->Matrix_Clear(iMesh);
}
void mw_vector_clear_(int* imesh)
{
    unsigned int iMesh = *imesh;
    pMW->Vector_Clear(iMesh);
}

int mw_matrix_add_elem_24_(int* imesh, int* ielem, double elem_matrix[][24])//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
{
    uint nNumOfCol=24;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_60_(int* imesh, int* ielem, double elem_matrix[][60])//Hexa  20Node * 3DOF
{
    uint nNumOfCol=60;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_12_(int* imesh, int* ielem, double elem_matirx[][12])//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
{
    uint nNumOfCol=12;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_30_(int* imesh, int* ielem, double elem_matirx[][30])//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
{
    uint nNumOfCol=30;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_18_(int* imesh, int* ielem, double elem_matirx[][18])//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
{
    uint nNumOfCol=18;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_45_(int* imesh, int* ielem, double elem_matirx[][45])//Prism 15Node * 3DOF
{
    uint nNumOfCol=45;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_20_(int* imesh, int* ielem, double elem_matirx[][20])//Quad   4Node * 5DOF
{
    uint nNumOfCol=20;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_40_(int* imesh, int* ielem, double elem_matirx[][40])//Quad   8Node * 5DOF
{
    uint nNumOfCol=40;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matrix_add_elem_15_(int* imesh, int* ielem, double elem_matirx[][15])//Tri    3Node * 5DOF, Beam 3Node * 5DOF
{
    uint nNumOfCol=15;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_9_(int* imesh, int* ielem, double elem_matirx[][9])  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
{
    uint nNumOfCol=9;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_48_(int* imesh, int* ielem, double elem_matirx[][48])//Quad   8Node * 6DOF
{
    uint nNumOfCol=48;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_6_(int* imesh, int* ielem, double elem_matirx[][6])  //Beam   2Node * 3DOF
{
    uint nNumOfCol=6;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
int mw_matirx_add_elem_10_(int* imesh, int* ielem, double elem_matirx[][10])//Beam   2Node * 5DOF
{
    uint nNumOfCol=10;

    double mat[nNumOfCol*nNumOfCol];

    for(int i=0; i < nNumOfCol; i++)
        for(int ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uint iMesh = *imesh;
    uint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}

// set_bc
//
int mw_matrix_set_bc_(int* imesh, int* inode, int* idof, double* value1, double* value2)
{
    unsigned int iMesh = *imesh;
    unsigned int iNode = *inode;
    unsigned int iDOF  = *idof;
    double Value1 = *value1;
    double Value2 = *value2;

    // matrix-D and solution_vector
    return pMW->Set_BC_Mat_SolVec(iMesh, iNode, iDOF, Value1, Value2);
}
int mw_matrix_rhs_set_bc_(int* imesh, int* inode, int* idof, double* value1, double* value2)
{
    unsigned int iMesh = *imesh;
    unsigned int iNode = *inode;
    unsigned int iDOF  = *idof;
    double Value1 = *value1;
    double Value2 = *value2;

    // matrix-D and rhs_vector
    return pMW->Set_BC_Mat_RHS(iMesh, iNode, iDOF, Value1, Value2);
}
int mw_rhs_set_bc_(int* imesh, int* inode, int* idof, double* value)
{
    unsigned int iMesh = *imesh;
    unsigned int iNode = *inode;
    unsigned int iDOF = *idof;
    double Value = *value;

    // rhs_vector
    return pMW->Set_BC_RHS(iMesh, iNode, iDOF, Value);
}

// solver
int mw_solve_(int* iter_max, double* tolerance, int* method, int* pre_condition)
{
    unsigned int nIterMax = *iter_max;
    double dTolerance = *tolerance;
    unsigned int nMethod = *method;
    unsigned int nPreCondition = *pre_condition;

    return pMW->Solve(nIterMax, dTolerance, nMethod, nPreCondition);
}
//void mw_store_matrix_()
//{
//    pMW->StoreMatrix();
//}
//void mw_load_matrix_()
//{
//    pMW->LoadMatrix();
//}

//----
// MG construct (refine)
//----
int mw_refine_()
{
    return pMW->Refine();
}
int mw_mg_construct_()
{
    return pMW->Refine();
}

void mw_finalize_refine_()
{
    pMW->FinalizeRefine();// release memory (final proc for mesh construct)
}
void mw_finalize_mg_construct_()
{
    pMW->FinalizeRefine();// release memory (final proc for mesh construct)
}

//----
// model
//----
//
// assemble model
//
int mw_get_num_of_assemble_model_()
{
    return (int)pMW->GetNumOfAssembleModel();
}
void mw_select_assemble_model_(int* mglevel)
{
    uint nMGLevel = *mglevel;
    pMW->SelectAssembleModel(nMGLevel);
}
//
// mesh part
//
int mw_get_num_of_mesh_part_()
{
    return (int)pMW->GetNumOfMeshPart();
}
void mw_select_mesh_part_with_id_(int* mesh_id)
{
    unsigned int nMeshID = *mesh_id;

    pMW->SelectMeshPart_ID(nMeshID);
}
void mw_select_mesh_part_(int* index)
{
    unsigned int nIndex = *index;

    pMW->SelectMeshPart_IX(nIndex);
}
//
// element
//
void mw_select_element_with_id_(int* elem_id)
{
    unsigned int nElemID = *elem_id;

    pMW->SelectElement_ID(nElemID);
}
void mw_select_element_(int* index)
{
    unsigned int nIndex = *index;

    pMW->SelectElement_IX(nIndex);
}
int mw_get_element_type_()
{
    return (int)pMW->GetElementType();
}
int mw_get_num_of_element_vert_()
{
    return (int)pMW->GetNumOfElementVert();
}

void mw_get_element_vert_node_id_(int v_node_id[])
{
    pMW->GetElementVertNodeID(v_node_id);
}
int mw_get_num_of_element_edge_()
{
    return (int)pMW->GetNumOfElementEdge();
}
void mw_get_element_edge_node_id_(int v_node_id[])
{
    pMW->GetElementEdgeNodeID(v_node_id);
}
//
// node
//
void mw_get_node_coord_(int* node_id, double* x, double* y, double* z)
{
    unsigned int nNodeID = (unsigned int)*node_id;
    double X,Y,Z;

    pMW->GetNodeCoord(nNodeID,X,Y,Z);
    
    *x = X;
    *y = Y;
    *z = Z;
}
int mw_get_dof_(int* node_id)
{
    unsigned int nNodeID = (unsigned int)*node_id;

    return (int)pMW->GetNumOfDOF(nNodeID);
}
int mw_get_dof_scalar_(int* node_id)
{
    unsigned int nNodeID = (unsigned int)*node_id;

    return (int)pMW->GetNumOfScalar(nNodeID);
}
int mw_get_dof_vector_(int* node_id)
{
    unsigned int nNodeID = (unsigned int)*node_id;

    return (int)pMW->GetNumOfVector(nNodeID);
}

void mw_set_node_value_(int* node_id, double value[])
{
    unsigned int nNodeID = (unsigned int)*node_id;

    pMW->SetNodeValue(nNodeID, value);
}
void mw_set_node_value_with_dof_(int* node_id, int* idof, double* value)
{
    unsigned int nNodeID = (unsigned int)*node_id;
    unsigned int iDOF = *idof;

    pMW->SetNodeValue(nNodeID, iDOF, *value);
}
void mw_get_node_value_(int* node_id, double value[])
{
    unsigned int nNodeID = (unsigned int)*node_id;

    pMW->GetNodeValue(nNodeID, value);
}
void mw_get_node_value_with_dof_(int* node_id, int* idof, double* value)
{
    unsigned int nNodeID = (unsigned int)*node_id;

    *value = pMW->GetNodeValue(nNodeID, (unsigned int)*idof);
}
//
// scalar_vector node
//
void mw_set_sv_node_value_(int* node_id, double v_value[], double s_value[])
{
    pMW->SetSVNodeValue(*node_id, v_value, s_value);
}
void mw_set_sv_node_value_with_dof_(int* node_id, 
        int* v_dof, double* v_value, int* s_dof, double* s_value)
{
    pMW->SetSVNodeValue(*node_id, *v_dof, *v_value, *s_dof, *s_value);
}
void mw_get_sv_node_value_(int* node_id, double v_value[], double s_value[])
{
    pMW->GetSVNodeValue(*node_id, v_value, s_value);
}
void mw_get_sv_node_value_with_dof_(int* node_id, 
        int* v_dof, double* v_value, int* s_dof, double* s_value)
{
    pMW->GetSVNodeValue(*node_id, *v_dof, *v_value, *s_dof, *s_value);
}


// node size, element size
int mw_get_num_of_node_()
{
    return (int)pMW->getNodeSize();
}
int mw_get_num_of_node_with_mesh_(int* imesh)
{
    unsigned int iMesh = *imesh;

    return (int)pMW->getNodeSize(iMesh);
}
int mw_get_num_of_element_()
{
    return (int)pMW->getElementSize();
}
int mw_get_num_of_element_with_mesh_(int* imesh)
{
    unsigned int iMesh = *imesh;

    return (int)pMW->getElementSize(iMesh);
}

//----
// node type
//----
int mw_nodetype_s_()
{
    return pMW->nodetype_s();
}
int mw_nodetype_v_()
{
    return pMW->nodetype_v();
}
int mw_nodetype_sv_()
{
    return pMW->nodetype_sv();
}
//----
// element type
//----
int mw_elemtype_hexa_()
{
    return pMW->elemtype_hexa();
}
int mw_elemtype_hexa2_()
{
    return pMW->elemtype_hexa2();
}
int mw_elemtype_tetra_()
{
    return pMW->elemtype_tetra();
}
int mw_elemtype_tetra2_()
{
    return pMW->elemtype_tetra2();
}
int mw_elemtype_prism_()
{
    return pMW->elemtype_prism();
}
int mw_elemtype_prism2_()
{
    return pMW->elemtype_prism2();
}
int mw_elemtype_quad_()
{
    return pMW->elemtype_quad();
}
int mw_elemtype_quad2_()
{
    return pMW->elemtype_quad2();
}
int mw_elemtype_triangle_()
{
    return pMW->elemtype_triangle();
}
int mw_elemtype_triangle2_()
{
    return pMW->elemtype_triangle2();
}
int mw_elemtype_line_()
{
    return pMW->elemtype_line();
}
int mw_elemtype_line2_()
{
    return pMW->elemtype_line2();
}




//----
// shape function
//----
int mw_get_num_of_integ_point_(int* shape_type)
{
    unsigned int nShapeType = *shape_type;

    return (int)pMW->NumOfIntegPoint(nShapeType);
}
void mw_shape_function_on_pt_(int* shape_type, int* igauss, double N[])
{
    unsigned int nShapeType = *shape_type;
    unsigned int iGauss = *igauss;

    pMW->ShapeFunc_on_pt(nShapeType, iGauss, N);
}
void mw_shape_function_hexa81_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Hexa81(*igauss, *ishape);
}
void mw_shape_function_hexa82_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Hexa82(*igauss, *ishape);
}
void mw_shape_function_hexa201_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Hexa201(*igauss, *ishape);
}
void mw_shape_function_hexa202_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Hexa202(*igauss, *ishape);
}
void mw_shape_function_hexa203_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Hexa203(*igauss, *ishape);
}
void mw_shape_function_tetra41_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Tetra41(*igauss, *ishape);
}
void mw_shape_function_tetra101_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Tetra101(*igauss, *ishape);
}
void mw_shape_function_tetra104_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Tetra104(*igauss, *ishape);
}
void mw_shape_function_tetra1015_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Tetra1015(*igauss, *ishape);
}
void mw_shape_function_prism62_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Prism62(*igauss, *ishape);
}
void mw_shape_function_prism156_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Prism156(*igauss, *ishape);
}
void mw_shape_function_prism159_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Prism159(*igauss, *ishape);
}
void mw_shape_function_prism1518_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Prism1518(*igauss, *ishape);
}
void mw_shape_function_quad41_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Quad41(*igauss, *ishape);
}
void mw_shape_function_quad84_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Quad84(*igauss, *ishape);
}
void mw_shape_function_quad89_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Quad89(*igauss, *ishape);
}
void mw_shape_function_tri31_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Triangle31(*igauss, *ishape);
}
void mw_shape_function_tri63_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Triangle63(*igauss, *ishape);
}
void mw_shape_function_line21_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Line21(*igauss, *ishape);
}
void mw_shape_function_line32_(int* igauss, int* ishape, double* N)
{
    *N = pMW->ShapeFunc_Line32(*igauss, *ishape);
}


//----
// shape function deriv (rst coord)
//----
void mw_dndr_(int* shape_type, double dndr[])
{
    pMW->dNdr(*shape_type, dndr);
}
void mw_dndr_hexa81_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Hexa81_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_hexa82_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Hexa82_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_hexa201_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Hexa201_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_hexa202_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Hexa202_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_hexa203_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Hexa203_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tetra41_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tetra41_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tetra101_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tetra101_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tetra104_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tetra104_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tetra1015_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tetra1015_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_prism62_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Prism62_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_prism156_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Prism156_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_prism159_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Prism159_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_prism1518_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Prism1518_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_quad41_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Quad41_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_quad84_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Quad84_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_quad89_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Quad89_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tri31_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tri31_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_tri63_(int* igauss, int* ishape, int* iaxis, double* dndr)
{
    *dndr = pMW->dNdr_Tri63_on_pt_on_shape(*igauss, *ishape, *iaxis);
}
void mw_dndr_line21_(int* igauss, int* ishape, double* dndr)
{
    *dndr = pMW->dNdr_Line21_on_pt_on_shape(*igauss, *ishape);
}
void mw_dndr_line32_(int* igauss, int* ishape, double* dndr)
{
    *dndr = pMW->dNdr_Line32_on_pt_on_shape(*igauss, *ishape);
}

//----
// shape function deriv (xyz coord)
//----
void mw_dndx_(int* elem_type, int* num_of_integ, int* ielem, double dndx[])// [igauss][ishape][iaxis] : 積分点数 節点数 3(座標)
{
    pMW->dNdx(*elem_type, *num_of_integ, *ielem, dndx);
}
void mw_det_jacobian_(int* elem_type, int* num_of_integ, int* igauss, double* det_j)
{
    pMW->detJacobian(*elem_type, *num_of_integ, *igauss, *det_j);
}
void mw_weight_(int* elem_type, int* num_of_integ, int* igauss, double* w)
{
    pMW->Weight(*elem_type, *num_of_integ, *igauss, *w);
}

//----
// shape function type
//----
int mw_shapetype_hexa81_()
{
    return pMW->shapetype_hexa81();
}
int mw_shapetype_hexa82_()
{
    return pMW->shapetype_hexa82();
}
int mw_shapetype_hexa201_()
{
    return pMW->shapetype_hexa201();
}
int mw_shapetype_hexa202_()
{
    return pMW->shapetype_hexa202();
}
int mw_shapetype_hexa203_()
{
    return pMW->shapetype_hexa203();
}
int mw_shapetype_tetra41_()
{
    return pMW->shapetype_tetra41();
}
int mw_shapetype_tetra101_()
{
    return pMW->shapetype_tetra101();
}
int mw_shapetype_tetra104_()
{
    return pMW->shapetype_tetra104();
}
int mw_shapetype_tetra1015_()
{
    return pMW->shapetype_tetra1015();
}
int mw_shapetype_prism62_()
{
    return pMW->shapetype_prism62();
}
int mw_shapetype_prism156_()
{
    return pMW->shapetype_prism156();
}
int mw_shapetype_prism159_()
{
    return pMW->shapetype_prism159();
}
int mw_shapetype_prism1518_()
{
    return pMW->shapetype_prism1518();
}
int mw_shapetype_quad41_()
{
    return pMW->shapetype_quad41();
}
int mw_shapetype_quad84_()
{
    return pMW->shapetype_quad84();
}
int mw_shapetype_quad89_()
{
    return pMW->shapetype_quad89();
}
int mw_shapetype_tri31_()
{
    return pMW->shapetype_tri31();
}
int mw_shapetype_tri63_()
{
    return pMW->shapetype_tri63();
}
int mw_shapetype_line21_()
{
    return pMW->shapetype_line21();
}
int mw_shapetype_line32_()
{
    return pMW->shapetype_line32();
}


//--
// boundary mesh
//--
int mw_get_num_of_boundary_bnode_mesh_(){ return pMW->GetNumOfBoundaryNodeMesh();}
int mw_get_num_of_boundary_bface_mesh_(){ return pMW->GetNumOfBoundaryFaceMesh();}
int mw_get_num_of_boundary_bedge_mesh_(){ return pMW->GetNumOfBoundaryEdgeMesh();}
int mw_get_num_of_boundary_bvolume_mesh_(){ return pMW->GetNumOfBoundaryVolumeMesh();}

int mw_get_num_of_bnode_in_bnode_mesh_(int* ibmesh){ return pMW->GetNumOfBNode_BNodeMesh(*ibmesh);}
int mw_get_num_of_bnode_in_bface_mesh_(int* ibmesh){ return pMW->GetNumOfBNode_BFaceMesh(*ibmesh);}
int mw_get_num_of_bnode_in_bedge_mesh_(int* ibmesh){ return pMW->GetNumOfBNode_BEdgeMesh(*ibmesh);}
int mw_get_num_of_bnode_in_bvolume_mesh_(int* ibmesh){ return pMW->GetNumOfBNode_BVolumeMesh(*ibmesh);}

int mw_get_num_of_dof_in_bnode_mesh_(int* ibmesh, int* ibnode){ return pMW->GetNumOfDOF_BNodeMesh(*ibmesh, *ibnode);}
int mw_get_num_of_dof_in_bface_mesh_(int* ibmesh){ return pMW->GetNumOfDOF_BFaceMesh(*ibmesh);}
int mw_get_num_of_dof_in_bedge_mesh_(int* ibmesh){ return pMW->GetNumOfDOF_BEdgeMesh(*ibmesh);}
int mw_get_num_of_dof_in_bvolume_mesh_(int* ibmesh){ return pMW->GetNumOfDOF_BVolumeMesh(*ibmesh);}
//--
// value of boundary node
//--
double mw_get_bnode_value_in_bnode_mesh_(int* ibmesh, int* ibnode, int* idof)
{
    return pMW->GetBNodeValue_BNodeMesh(*ibmesh, *ibnode, *idof);
}
double mw_get_bnode_value_in_bface_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel)
{
    return pMW->GetBNodeValue_BFaceMesh(*ibmesh, *ibnode, *idof, *mglevel);
}
double mw_get_bnode_value_in_bedge_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel)
{
    return pMW->GetBNodeValue_BEdgeMesh(*ibmesh, *ibnode, *idof, *mglevel);
}
double mw_get_bnode_value_in_bvolume_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel)
{
    return pMW->GetBNodeValue_BVolumeMesh(*ibmesh, *ibnode, *idof, *mglevel);
}
// node id (in boundary node)
int mw_get_node_id_in_bnode_mesh_(int* ibmesh, int* ibnode)
{
    return pMW->GetNodeID_BNode_BNodeMesh(*ibmesh, *ibnode);
}
int mw_get_node_id_in_bface_mesh_(int* ibmesh, int* ibnode)
{
    return pMW->GetNodeID_BNode_BFaceMesh(*ibmesh, *ibnode);
}
int mw_get_node_id_in_bedge_mesh_(int* ibmesh, int* ibnode)
{
    return pMW->GetNodeID_BNode_BEdgeMesh(*ibmesh, *ibnode);
}
int mw_get_node_id_in_bvolume_mesh_(int* ibmesh, int* ibnode)
{
    return pMW->GetNodeID_BNode_BVolumeMesh(*ibmesh, *ibnode);
}
//--
// value of boundary face, edge, volume
//--
int mw_get_num_of_bface_(int* ibmesh)
{
    return pMW->GetNumOfBFace(*ibmesh);
}
double mw_get_bface_value_(int* ibmesh, int* ibface, int* idof)
{
    return pMW->GetBFaceValue(*ibmesh, *ibface, *idof);
}
int mw_get_num_of_bedge_(int* ibmesh)
{
    return pMW->GetNumOfBEdge(*ibmesh);
}
double mw_get_bedge_value_(int* ibmesh, int* ibedge, int* idof)
{
    return pMW->GetBEdgeValue(*ibmesh, *ibedge, *idof);
}
int mw_get_num_of_bvolume_(int* ibmesh)
{
    return pMW->GetNumOfBVolume(*ibmesh);
}
double mw_get_bvolume_value_(int* ibmesh, int* ibvol, int* idof)
{
    return pMW->GetBVolumeValue(*ibmesh, *ibvol, *idof);
}



//--
// mpi
//--
int mw_mpi_int_(){ return MPI_INT; }        // MPI_INT
int mw_mpi_double_(){ return MPI_DOUBLE; }  // MPI_DOUBLE
int mw_mpi_comm_(){ return MPI_COMM_WORLD; }// MPI_COMM_WORLD

int mw_mpi_sum_(){ return MPI_SUM;}// op  ,use allreduce_r 
int mw_mpi_max_(){ return MPI_MAX;}// op  ,use allreduce_r
int mw_mpi_min_(){ return MPI_MIN;}// op  ,use allreduce_r

void mw_allreduce_r_(double val[], int* val_size, int* op)
{
    double rval[*val_size];

    pMW->AllReduce(val, rval, *val_size, MPI_DOUBLE, *op, MPI_COMM_WORLD);
}
void mw_allreduce_i_(int val[], int* val_size, int* op)
{
    int rval[*val_size];

    pMW->AllReduce(val, rval, *val_size, MPI_INT, *op, MPI_COMM_WORLD);
}
int mw_barrier_()
{
    return pMW->Barrier(MPI_COMM_WORLD);
}
int mw_abort_(int* error)
{
    return pMW->Abort(MPI_COMM_WORLD, *error);
}
int mw_allgather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt)
{
    return pMW->AllGather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, MPI_COMM_WORLD);
}
int mw_allgather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt)
{
    return pMW->AllGather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, MPI_COMM_WORLD);
}
int mw_gather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root)
{
    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
int mw_gather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root)
{
    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, MPI_COMM_WORLD);
}
int mw_scatter_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root)
{
    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
int mw_scatter_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root)
{
    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, MPI_COMM_WORLD);
}

//int mw_allgather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* comm)
//{
//    return pMW->Allgather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *comm);
//}
//int mw_allgather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* comm)
//{
//    return pMW->Allgather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *comm);
//}
//int mw_gather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root, int* comm)
//{
//    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, *comm);
//}
//int mw_gather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root, int* comm)
//{
//    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, *comm);
//}
//int mw_scatter_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root, int* comm)
//{
//    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, *comm);
//}
//int mw_scatter_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root, int* comm)
//{
//    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, *comm);
//}
int mw_get_rank_()
{
    return pMW->GetRank();
}
int mw_get_num_of_process_()
{
    return pMW->GetNumOfProcess();
}
void mw_send_recv_r2_(double buf[], int* dof_size)// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
{
    pMW->Send_Recv_R(buf, *dof_size);
}
void mw_send_recv_r_()// 通信Nodeの値を入れ替えて更新
{
    pMW->Send_Recv_R();
}

//--
// Element_Group { select AssyModel, select Mesh }
//--
int mw_get_num_of_elementgroup_()
{
    return (int)pMW->GetNumOfElementGroup();
}
int mw_get_num_of_element_id_(int* igrp)
{
    uint iGrp = *igrp;

    return (int)pMW->GetNumOfElementID(iGrp);
}
int mw_get_element_id_with_elementgroup_(int* igrp, int* index)
{
    uint iGrp = *igrp;
    uint i = *index;

    return (int)pMW->GetElementID_with_ElementGroup(iGrp, i);
}
int mw_get_elementgroup_name_length_(int* igrp)
{
    uint iGrp = *igrp;
    int nLength = (int)pMW->GetElementGroupName_Length(iGrp) + 1;

    return nLength;
}
void mw_get_elementgroup_name_(int* igrp, char* name, int* name_len)
{
    uint iGrp = *igrp;
    string sName = pMW->GetElementGroupName(iGrp);
    
    uint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > (uint)*name_len){
        uint nLength = (uint)*name_len-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uint nLength = (uint)*name_len;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}

//----
// logger
//----
void mw_logger_set_mode_(int* mode)
{
    pMW->LoggerMode(*mode);
}
void mw_logger_set_device_(int* mode, int* device)
{
    pMW->LoggerDevice(*mode, *device);
}
void mw_logger_info_ (int* mode, char* message, int* str_len)
{
    size_t nLength = (size_t)*str_len + 1;
    
    char cmsg[nLength];
    strncpy(cmsg, message, nLength);
    
    pMW->LoggerInfo(*mode, cmsg);
}
//----
// logger parameter
//----
int mw_get_error_mode_()
{
    return pMW->getErrorMode();
}
int mw_get_warn_mode_()
{
    return pMW->getWarnMode();
}
int mw_get_info_mode_()
{
    return pMW->getInfoMode();
}
int mw_get_debug_mode_()
{
    return pMW->getDebugMode();
}
int mw_get_disk_device_()
{
    return pMW->getDiskDevice();
}
int mw_get_display_device_()
{
    return pMW->getDisplayDevice();
}








