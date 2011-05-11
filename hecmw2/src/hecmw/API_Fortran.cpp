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

// MW3 standard 
iint mw_initialize_(int* argc, char** argv)
{
    pMW = pmw::CMW::Instance();

    return pMW->Initialize(*argc, argv);
}
// fstr style
iint mw_initialize_fstr_(int* argc, char** argv, char* ctrlname)
{
    pMW = pmw::CMW::Instance();

    string sCtrlName= ctrlname;

    return pMW->Initialize_fstr(*argc, argv, sCtrlName);
}
iint mw_finalize_()
{
    return pMW->Finalize();
}

//----
// banner
//----
void mw_banner_()
{
    pMW->Banner_Display();
}

//----
// file i/o API
//----
iint mw_file_read_(char* basename)// mesh_file read:General
{
    string sBasename= basename;

    return pMW->FileRead(sBasename, false);
}
//iint mw_file_read_bin_(char* basename)// mesh_file(binary) read:General
//{
//    string sBasename= basename;
//
//    return pMW->FileRead(sBasename, true);
//}
iint mw_file_read_fstr_() // mesh_file read:FrontISTR
{
    return pMW->FileRead_fstr(false);
}
//iint mw_file_read_bin_fstr_()// mesh_file read (binary file):FrontISTR
//{
//    return pMW->FileRead_fstr(true);
//}
iint mw_file_debug_write_() //data_check file write
{
    return pMW->FileDebugWrite();
}
//----
// fstr file_name (hecmw_ctrl)
//----
iint mw_get_fstr_filename_length_mesh_()
{
    string sName = pMW->getFstr_FileName_Mesh();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr mesh filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_control_()
{
    string sName = pMW->getFstr_FileName_Control();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr control filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_result_()
{
    string sName = pMW->getFstr_FileName_Result();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr result filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_restart_()
{
    string sName = pMW->getFstr_FileName_Restart();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr restart filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_part_in_()
{
    string sName = pMW->getFstr_FileName_PartIn();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr part_in filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_part_out_()
{
    string sName = pMW->getFstr_FileName_PartOut();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr part_out filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_mesh_()
{
    string sName = pMW->getFstr_FileName_VisMesh();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_mesh filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_in_()
{
    string sName = pMW->getFstr_FileName_VisIn();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_in filename length");
        return ERROR;
    }else{
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_out_()
{
    string sName = pMW->getFstr_FileName_VisOut();
    uiint nLength= sName.length();

    if(nLength > IINT_MAX){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_out filename length");
        return ERROR;
    }else{
        return nLength;
    }
}

void mw_get_fstr_filename_mesh_(char name[], iint* len)
{
    uiint ncLen = strlen(name);
    
    string sName = pMW->getFstr_FileName_Mesh();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr mesh filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_control_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_Control();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr control filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_result_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_Result();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr result filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_restart_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_Restart();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr restart filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_part_in_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_PartIn();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr part_in filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_part_out_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_PartOut();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr part_out filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_mesh_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_VisMesh();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_mesh filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_in_(char name[], iint* len)
{
    uiint ncLen = strlen(name);

    string sName = pMW->getFstr_FileName_VisIn();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_in filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_out_(char name[], iint* len)
{
    uiint ncLen = strlen(name);
    
    string sName = pMW->getFstr_FileName_VisOut();
    uiint i, nLength = sName.length();

    if(ncLen != nLength){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_out filename length");
    }

    for(i=0; i < nLength; i++){
        name[i]=sName[i];
    }
}

//----
// result
//----
void mw_rlt_start_(iint* step)
{
    if(*step < 0) return;

    uiint nStep = *step;

    pMW->PrintRlt_Start(nStep, false);
}
void mw_rlt_start_bin_(iint* step)
{
    if(*step < 0) return;

    uiint nStep = *step;

    pMW->PrintRlt_Start(nStep, true);
}
/* 
 * format = %d:iint*, %f:double*(fixed), %e:double*(scientific), %s:const char*
 */
iint mw_rlt_print_(iint* width, char* format, ... )
{
    if(*width < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rlt_print_, minus value is invalid argument. width");
        
        return ERROR;
    }

    vector<void*> param;

    iint *pnVal;  double *pdVal;  string sVal; char *cVal;
    uiint nLength = strlen(format);
    uiint nWidth = *width;

    va_list list;

    va_start( list, format);
    for(uiint i=0; i < nLength; i++){
        if(format[i] == '%'){
            ++i;
            switch( format[i] ){
            case('d'):
                pnVal= va_arg( list, iint* );

                param.push_back(pnVal);
                break;
            case('f'):
                pdVal= va_arg( list, double* );

                param.push_back(pdVal);
                break;
            case('e'):
                pdVal= va_arg( list, double* );

                param.push_back(pdVal);
                break;
            case('s'):
                sVal= va_arg( list, char* );

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


    FileIO::CFileIO* pFileIO= FileIO::CFileIO::Instance();

    pFileIO->PrintResult(nWidth,format, param);
    
    return SUCCESS;
}
void mw_rlt_end_()
{
    pMW->PrintRlt_End();
}
// output *.inp 
void mw_print_avs_basis_()// basis variable output
{
    pMW->PrintMicroAVS_Basis();
}
void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof)
{
    if( *imesh < 0 ){
        return;
    }
    if( *ndof < 0 ){
        return;
    }
    uiint iMesh= *imesh;
    uiint nDOF = *ndof;

    pMW->recAVS_Label(iMesh, label, unit, nDOF);
}
void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value)
{
    if( *imesh < 0 ){
        return;
    }
    if( *num_of_node < 0 ){
        return;
    }
    uiint iMesh= *imesh;
    uiint nNumOfNode= *num_of_node;

    pMW->recAVS_Variable(iMesh, nNumOfNode, label, value);
}
void mw_print_avs_fem_()// record variable output
{
    pMW->PrintMicroAVS_FEM();
}


//----
// restart, linear_algebra_equation info
//----
iint mw_file_write_res_(iint* step) // restart file write
{
    if(*step < 0){
        return ERROR;
    }
    uiint nStep= *step;
    return pMW->FileWriteRes(nStep, false);
}
iint mw_set_restart_(iint* step)// restart file read, linear_algebra_equation info
{
    if(*step < 0){
        return ERROR;
    }
    uiint nStep= *step;
    return pMW->SetRestart(nStep, false);
}

iint mw_file_write_res_bin_(iint* step)// restart file write (binary file)
{
    if(*step < 0){
        return ERROR;
    }
    uiint nStep= *step;
    return pMW->FileWriteRes(nStep, true);
}
iint mw_set_restart_bin_(iint* step)// restart file read, linear_algebra_equation info (binary file)
{
    if(*step < 0){
        return ERROR;
    }
    uiint nStep= *step;
    return pMW->SetRestart(nStep, true);
}


//----
// linear solver API
//----
void mw_gene_linear_algebra_(iint* num_of_algebra, iint dof[])
{
    if(*num_of_algebra < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gene_linear_algebra_, minus value is invalid argument. num_of_algebra");
        
        return;
    }

    uiint nNumOfAlgebra = (uiint)*num_of_algebra;
    uiint vDOF[nNumOfAlgebra];
    for(uiint i=0; i < nNumOfAlgebra; i++) vDOF[i] = dof[i];

    pMW->GeneLinearAlgebra(nNumOfAlgebra, vDOF);
}
void mw_select_algebra_(iint* ieq)
{
    if(*ieq < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_algebra_, minus value is invalid argument. ieq");
        
        return;
    }

    uiint iEqu = *ieq;
    pMW->SelectAlgebra(iEqu);
}

// matrix add elem
iint mw_matrix_add_elem_(iint* imesh,  iint* ielem,  double elem_matrix[])//standard
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_elem_, minus value is invalid argument. imesh");
        
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_elem_, minus value is invalid argument. ielem");
        
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, elem_matrix);
}

iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[])
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. imesh");
        
        return ERROR;
    }
    if(*i_index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. i_index");
        
        return ERROR;
    }
    if(*j_index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. j_index");
        
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint inode = *i_index;
    uiint jnode = *j_index;

    return pMW->Matrix_Add_Node(iMesh, inode, jnode, nodal_matrix);
}

// matrix 0 clear
//
void mw_matrix_clear_(iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_clear_, minus value is invalid argument. imesh");
        
        return;
    }

    uiint iMesh = *imesh;
    pMW->Matrix_Clear(iMesh);
}
void mw_vector_clear_(iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_vector_clear_, minus value is invalid argument. imesh");
        
        return;
    }

    uiint iMesh = *imesh;
    pMW->Vector_Clear(iMesh);
}

iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24])//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_24_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_24_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=24;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60])//Hexa  20Node * 3DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_60_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_60_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=60;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12])//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_12_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_12_, minus value is invalid argument. ielem");
        return ERROR;
    }
    
    uiint nNumOfCol=12;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30])//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_30_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_30_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=30;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18])//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_18_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_18_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=18;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45])//Prism 15Node * 3DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_45_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_45_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=45;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20])//Quad   4Node * 5DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_20_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_20_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=20;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40])//Quad   8Node * 5DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_40_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_40_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=40;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15])//Tri    3Node * 5DOF, Beam 3Node * 5DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_15_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_15_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=15;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9])  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_9_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_9_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=9;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48])//Quad   8Node * 6DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_48_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_48_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=48;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6])  //Beam   2Node * 3DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_6_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_6_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=6;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}
iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10])//Beam   2Node * 5DOF
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_10_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ielem < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_10_, minus value is invalid argument. ielem");
        return ERROR;
    }

    uiint nNumOfCol=10;

    double mat[nNumOfCol*nNumOfCol];

    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];

    uiint iMesh = *imesh;
    uiint iElem = *ielem;

    return pMW->Matrix_Add_Elem(iMesh, iElem, mat);
}

// set boundary condition
//
iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* inode, iint* dof, double* diagval, double* solval)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. inode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint iNode = *inode;
    uiint nDOF  = *dof;
    double dDiagValue = *diagval;
    double dSolValue = *solval;

    // matirix_diag=diagval, off_diag=0  , X=solval
    return pMW->Set_BC_Mat_RHS2(iMesh, iNode, nDOF, dDiagValue, dSolValue);
}
// penalty method
iint mw_matrix_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* diagval, double* rhsval)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. inode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint iNode = *inode;
    uiint nDOF  = *dof;
    double dDiagValue = *diagval;
    double dRHSValue = *rhsval;

    // matrix_diag, rhs_vector
    return pMW->Set_BC_Mat_RHS(iMesh, iNode, nDOF, dDiagValue, dRHSValue);
}
iint mw_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* value)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. inode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint iNode = *inode;
    uiint nDOF = *dof;
    double Value = *value;

    // rhs_vector
    return pMW->Set_BC_RHS(iMesh, iNode, nDOF, Value);
}
iint mw_rhs_add_bc_(iint* imesh, iint* inode, iint* dof, double* value)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. inode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iMesh = *imesh;
    uiint iNode = *inode;
    uiint nDOF = *dof;
    double Value = *value;

    // add rhs_vector
    return pMW->Add_BC_RHS(iMesh, iNode, nDOF, Value);
}

// solver
iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition)
{
    if(*iter_max < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. iter_max");
        return ERROR;
    }
    if(*method < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. method");
        return ERROR;
    }
    if(*pre_condition < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. pre_condition");
        return ERROR;
    }

    iint nIterMax = *iter_max;
    double dTolerance = *tolerance;
    iint nMethod = *method;
    iint nPreCondition = *pre_condition;

    return pMW->Solve(nIterMax, dTolerance, nMethod, nPreCondition);
}
//--
// solution_vector copy,  at select MG-Level && select Equation
//--
void mw_get_solution_vector_(double buf[], iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_vector_, minus value is invalid argument. imesh");

        return;
    }

    uiint iMesh = *imesh;
    pMW->GetSolution_Vector(buf, iMesh);
}
void mw_get_solution_assy_vector_(double buf[])
{
    pMW->GetSolution_AssyVector(buf);
}
//--
// rhs_vector copy,  at select MG-Level && select Equation
//--
void mw_get_rhs_vector_(double buf[], iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_vector_, minus value is invalid argument. imesh");
        
        return;
    }
    
    uiint iMesh = *imesh;
    pMW->GetRHS_Vector(buf, iMesh);
}
void mw_get_rhs_assy_vector_(double buf[])
{
    pMW->GetRHS_AssyVector(buf);
}

//--
// solution vector value, rhs vector value
//--
double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. imesh");

        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. inode");

        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. idof");

        return ERROR;
    }

    return pMW->GetSolutionVector_Val(*imesh, *inode, *idof);
}
double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. imesh");
        
        return ERROR;
    }
    if(*inode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. inode");
        
        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. idof");
        
        return ERROR;
    }

    return pMW->GetRHSVector_Val(*imesh, *inode, *idof);
}
//--
// solution vector dof, rhs vector dof
//--
iint mw_get_solution_assy_vector_dof_()
{
    uiint nDOF= pMW->GetSolutionVector_DOF();

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_dof_, return value is over INT_MAX");

        return IINT_MAX;
    }
}
iint mw_get_rhs_assy_vector_dof_()
{
    uiint nDOF= pMW->GetRHSVector_DOF();

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_dof_, return value is over INT_MAX");

        return IINT_MAX;
    }
}

//--
// AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
//--
void mw_mult_vector_(double vX[], double vB[])
{
    pMW->multVector(vX, vB);
}




//----
// MG construct (refine)
//----
iint mw_refine_(int* num_of_refine)
{
    if(*num_of_refine < 0){
        return ERROR;
    }
    uiint nNumOfRefine = *num_of_refine;

    return pMW->Refine(nNumOfRefine);
}
iint mw_mg_construct_(int* num_of_refine)
{
    if(*num_of_refine < 0){
        return ERROR;
    }
    uiint nNumOfRefine = *num_of_refine;

    return pMW->Refine(nNumOfRefine);
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
// number of assemble model == number of mglevel
//
iint mw_get_num_of_assemble_model_()
{
    uiint nNumOfLevel= pMW->GetNumOfAssembleModel();

    if(nNumOfLevel < IINT_MAX){
        return nNumOfLevel;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_assemble_model_, return value is over INT_MAX");

        return IINT_MAX;
    }
}
void mw_select_assemble_model_(iint* mglevel)
{
    if(*mglevel < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_assemble_model_, minus value is invalid argument. mglevel");
        
        return;
    }
    
    uiint nMGLevel = *mglevel;
    pMW->SelectAssembleModel(nMGLevel);
}
//
// mesh part
//
iint mw_get_num_of_mesh_part_()
{
    uiint nNumOfMesh= pMW->GetNumOfMeshPart();

    if(nNumOfMesh < IINT_MAX){
        return nNumOfMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_mesh_part_, return value is over INT_MAX");
        
        return IINT_MAX;
    }
}
void mw_select_mesh_part_with_id_(iint* mesh_id)
{
    if(*mesh_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_mesh_part_with_id_, minus value is invalid argument. mesh_id");
        
        return;
    }

    uiint nMeshID = *mesh_id;

    pMW->SelectMeshPart_ID(nMeshID);
}
void mw_select_mesh_part_(iint* index)
{
    if(*index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_mesh_part_, minus value is invalid argument. index");
        
        return;
    }

    uiint nIndex = *index;

    pMW->SelectMeshPart_IX(nIndex);
}
//
// element
//
void mw_select_element_with_id_(iint* elem_id)
{
    if(*elem_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_element_with_id_, minus value is invalid argument. elem_id");
        
        return;
    }

    uiint nElemID = *elem_id;

    pMW->SelectElement_ID(nElemID);
}
void mw_select_element_(iint* index)
{
    if(*index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_element_, minus value is invalid argument. index");

        return;
    }

    uiint nIndex = *index;

    pMW->SelectElement_IX(nIndex);
}
iint mw_get_element_type_()
{
    return pMW->GetElementType();
}
iint mw_get_num_of_element_vert_()
{
    return pMW->GetNumOfElementVert();
}

void mw_get_element_vert_node_id_(iint v_node_id[])
{
    pMW->GetElementVertNodeID(v_node_id);
}
iint mw_get_num_of_element_edge_()
{
    return pMW->GetNumOfElementEdge();
}
void mw_get_element_edge_node_id_(iint v_node_id[])
{
    pMW->GetElementEdgeNodeID(v_node_id);
}

// id && index
iint mw_get_node_id_(iint* index)
{
    if(*index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_, minus value is invalid argument. index");
        return ERROR;
    }
    uiint nIndex= *index;

    uiint nNodeID= pMW->getNodeID(nIndex);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_element_id_(iint* index)
{
    if(*index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_, minus value is invalid argument. index");
        return ERROR;
    }
    uiint nIndex= *index;

    uiint nElemID= pMW->getElementID(nIndex);
    
    if(nElemID < IINT_MAX){
        return nElemID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_, return value is over INT_MAX");
        return IINT_MAX;
    } 
}
iint mw_get_node_index_(iint* id)
{
    if(*id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_index_, minus value is invalid argument. id");
        return ERROR;
    }
    uiint nID= *id;

    uiint nIndex= pMW->getNodeIndex(nID);

    if(nIndex < IINT_MAX){
        return nIndex;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_index_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_element_index_(iint* id)
{
    if(*id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_index_, minus value is invalid argument. id");
        return ERROR;
    }
    uiint nID= *id;

    uiint nIndex= pMW->getElementIndex(nID);

    if(nIndex < IINT_MAX){
        return nIndex;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_index_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//--
// node connectivity,  itemU[]:upper triangle, itemL[]:lower triangle
//--
void mw_construct_node_connect_fem_(iint* node_id)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_construct_node_connect_fem_, minus value is invalid argument. node_id");
        return ;
    }

    uiint nNodeID = *node_id;
    pMW->constructNodeConnectFEM(nNodeID);
}
void mw_get_node_connect_fem_size_(iint* num_itemu, iint* num_iteml)
{
    uiint nItemU, nItemL;
    pMW->getNodeConnectFEM_Size(nItemU, nItemL);

    *num_itemu = nItemU;
    *num_iteml = nItemL;
}
void mw_get_node_connect_fem_item_(iint itemU[], iint itemL[])
{
    pMW->getNodeConnectFEM_Item_F(itemU, itemL);
}
//--
// node around element_id
//--
iint mw_get_num_of_aggregate_element_(iint* node_id)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_coord_, minus value is invalid argument. node_id");
        return ERROR;
    }
    uiint nNodeID= *node_id;

    uiint nNumOfAggElem= pMW->getNumOfAggregateElement(nNodeID);

    if(nNumOfAggElem > IINT_MAX){
        return IINT_MAX;
    }else{
        return nNumOfAggElem;
    }

}
iint mw_get_aggregate_element_id_(iint* node_id, iint* index)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_coord_, minus value is invalid argument. node_id");
        return ERROR;
    }
    uiint nNodeID= *node_id;
    uiint nIndex = *index;

    uiint nElemID= pMW->getAggregateElementID(nNodeID, nIndex);

    if(nElemID > IINT_MAX){
        return IINT_MAX;
    }else{
        return nElemID;
    }
}

//--
// node
//--
void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_coord_, minus value is invalid argument. node_id");
        return;
    }

    uiint nNodeID = *node_id;
    double X,Y,Z;

    pMW->GetNodeCoord(nNodeID,X,Y,Z);
    
    *x = X;
    *y = Y;
    *z = Z;
}
iint mw_get_dof_(iint* node_id)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_, minus value is invalid argument. node_id");
        return ERROR;
    }

    uiint nNodeID = *node_id;
    uiint nNumOfDOF= pMW->GetNumOfDOF(nNodeID);

    if(nNumOfDOF < IINT_MAX){
        return nNumOfDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_scalar_(iint* node_id)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_scalar_, minus value is invalid argument. node_id");
        return ERROR;
    }

    uiint nNodeID = *node_id;
    uiint nNumOfScalar= pMW->GetNumOfScalar(nNodeID);
    
    if(nNumOfScalar < IINT_MAX){
        return nNumOfScalar;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_scalar_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_vector_(iint* node_id)
{
    if(*node_id < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_vector_, minus value is invalid argument. node_id");
        return ERROR;
    }

    uiint nNodeID = *node_id;
    uiint nNumOfVector= pMW->GetNumOfVector(nNodeID);

    if(nNumOfVector < IINT_MAX){
        return nNumOfVector;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_vector_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

////void mw_set_node_value_(iint* node_id, double value[])
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_node_value_, minus value is invalid argument. node_id");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////
////    pMW->SetNodeValue(nNodeID, value);
////}
////void mw_set_node_value_with_dof_(iint* node_id, iint* idof, double* value)
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_node_value_with_dof_, minus value is invalid argument. node_id");
////        return;
////    }
////    if(*idof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_node_value_with_dof_, minus value is invalid argument. idof");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////    uiint iDOF = *idof;
////
////    pMW->SetNodeValue(nNodeID, iDOF, *value);
////}
////void mw_get_node_value_(iint* node_id, double value[])
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_value_, minus value is invalid argument. node_id");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////
////    pMW->GetNodeValue(nNodeID, value);
////}
////void mw_get_node_value_with_dof_(iint* node_id, iint* idof, double* value)
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_value_with_dof_, minus value is invalid argument. node_id");
////        return;
////    }
////    if(*idof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_value_with_dof_, minus value is invalid argument. idof");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////    uiint iDOF = *idof;
////
////    *value = pMW->GetNodeValue(nNodeID, iDOF);
////}
//
// scalar_vector node
//
////void mw_set_sv_node_value_(iint* node_id, double v_value[], double s_value[])
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_sv_node_value_, minus value is invalid argument. node_id");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////
////    pMW->SetSVNodeValue(nNodeID, v_value, s_value);
////}
////void mw_set_sv_node_value_with_dof_(iint* node_id,
////        iint* v_dof, double* v_value, iint* s_dof, double* s_value)
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_sv_node_value_with_dof_, minus value is invalid argument. node_id");
////        return;
////    }
////    if(*v_dof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_sv_node_value_with_dof_, minus value is invalid argument. v_dof");
////        return;
////    }
////    if(*s_dof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_set_sv_node_value_with_dof_, minus value is invalid argument. s_dof");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////    uiint nvDOF = *v_dof;
////    uiint nsDOF = *s_dof;
////
////    pMW->SetSVNodeValue(nNodeID, nvDOF, *v_value, nsDOF, *s_value);
////}
////void mw_get_sv_node_value_(iint* node_id, double v_value[], double s_value[])
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_sv_node_value_, minus value is invalid argument. node_id");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////
////    pMW->GetSVNodeValue(nNodeID, v_value, s_value);
////}
////void mw_get_sv_node_value_with_dof_(iint* node_id,
////        iint* v_dof, double* v_value, iint* s_dof, double* s_value)
////{
////    if(*node_id < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_sv_node_value_with_dof_, minus value is invalid argument. node_id");
////        return;
////    }
////    if(*v_dof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_sv_node_value_with_dof_, minus value is invalid argument. v_dof");
////        return;
////    }
////    if(*s_dof < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_get_sv_node_value_with_dof_, minus value is invalid argument. s_dof");
////        return;
////    }
////
////    uiint nNodeID = *node_id;
////    uiint nvDOF = *v_dof;
////    uiint nsDOF = *s_dof;
////
////    pMW->GetSVNodeValue(nNodeID, nvDOF, *v_value, nsDOF, *s_value);
////}


// node size, element size
iint mw_get_num_of_node_()
{
    uiint nNumOfNode= pMW->getNodeSize();

    if(nNumOfNode < IINT_MAX){
        return nNumOfNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_with_mesh_(iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_with_mesh_, minus value is invalid argument. imesh");
        return ERROR;
    }
    uiint iMesh = *imesh;
    
    uiint nNumOfNode= pMW->getNodeSize(iMesh);

    if(nNumOfNode < IINT_MAX){
        return nNumOfNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_with_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_()
{
    uiint nNumOfElem= pMW->getElementSize();

    if(nNumOfElem < IINT_MAX){
        return nNumOfElem;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_with_mesh_(iint* imesh)
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_with_mesh_, minus value is invalid argument. imesh");
        return ERROR;
    }

    uiint iMesh = *imesh;

    uiint nNumOfElem= pMW->getElementSize(iMesh);

    if(nNumOfElem < IINT_MAX){
        return nNumOfElem;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_with_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//----
// node type
//----
iint mw_nodetype_s_()
{
    return pMW->nodetype_s();
}
iint mw_nodetype_v_()
{
    return pMW->nodetype_v();
}
iint mw_nodetype_sv_()
{
    return pMW->nodetype_sv();
}
//----
// element type
//----
iint mw_elemtype_hexa_(){  return pMW->elemtype_hexa();}
iint mw_elemtype_hexa2_(){  return pMW->elemtype_hexa2();}
iint mw_elemtype_tetra_(){  return pMW->elemtype_tetra();}
iint mw_elemtype_tetra2_(){ return pMW->elemtype_tetra2();}
iint mw_elemtype_prism_(){  return pMW->elemtype_prism();}
iint mw_elemtype_prism2_(){  return pMW->elemtype_prism2();}
iint mw_elemtype_quad_(){  return pMW->elemtype_quad();}
iint mw_elemtype_quad2_(){  return pMW->elemtype_quad2();}
iint mw_elemtype_triangle_(){  return pMW->elemtype_triangle();}
iint mw_elemtype_triangle2_(){  return pMW->elemtype_triangle2();}
iint mw_elemtype_line_(){  return pMW->elemtype_line();}
iint mw_elemtype_line2_(){  return pMW->elemtype_line2();}

//----
// frontISTR element type
//----
iint mw_fistr_elemtype_hexa_(){  return pMW->fistr_elemtype_hexa();}
iint mw_fistr_elemtype_hexa2_(){ return pMW->fistr_elemtype_hexa2();}
iint mw_fistr_elemtype_tetra_(){  return pMW->fistr_elemtype_tetra();}
iint mw_fistr_elemtype_tetra2_(){ return pMW->fistr_elemtype_tetra2();}
iint mw_fistr_elemtype_prism_(){  return pMW->fistr_elemtype_prism();}
iint mw_fistr_elemtype_prism2_(){ return pMW->fistr_elemtype_prism2();}
iint mw_fistr_elemtype_quad_(){   return pMW->fistr_elemtype_quad();}
iint mw_fistr_elemtype_quad2_(){  return pMW->fistr_elemtype_quad2();}
iint mw_fistr_elemtype_triangle_(){ return pMW->fistr_elemtype_triangle();}
iint mw_fistr_elemtype_triangle2_(){ return pMW->fistr_elemtype_triangle2();}
iint mw_fistr_elemtype_line_(){  return pMW->fistr_elemtype_line();}
iint mw_fistr_elemtype_line2_(){ return pMW->fistr_elemtype_line2();}
//----
// frontISTR element type => MW3 element type
//----
iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype)
{
    if(*fistr_elemtype < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_fistr_elemtype_to_mw3_elemtype_, minus value is invalid argument. fistr_elemtype");
        return ERROR;
    }
    uiint nElemType = *fistr_elemtype;

    return pMW->fistr_elemtype_to_mw3_elemtype(nElemType);
}
//----
// MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
//----
iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype)
{
    if(*mw3_elemtype < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_mw3_elemtype_to_fistr_elemtype_, minus value is invalid argument. mw3_elemtype");
        return ERROR;
    }
    uiint nElemType = *mw3_elemtype;

    return pMW->mw3_elemtype_to_fistr_elemtype(nElemType);
}

//----
// shape function
//----
iint mw_get_num_of_integ_point_(iint* shape_type)
{
    if(*shape_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_integ_point_, minus value is invalid argument. shape_type");
        return ERROR;
    }

    uiint nShapeType = *shape_type;

    return pMW->NumOfIntegPoint(nShapeType);
}
void mw_shape_function_on_pt_(iint* shape_type, iint* igauss, double N[])
{
    if(*shape_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_on_pt_, minus value is invalid argument. shape_type");
        return;
    }
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_on_pt_, minus value is invalid argument. igauss");
        return;
    }

    uiint nShapeType = *shape_type;
    uiint iGauss = *igauss;

    pMW->ShapeFunc_on_pt(nShapeType, iGauss, N);
}
void mw_shape_function_hexa81_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa81_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa81_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;
    
    *N = pMW->ShapeFunc_Hexa81(iGauss, nShape);
}
void mw_shape_function_hexa82_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa82_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa82_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Hexa82(iGauss, nShape);
}
void mw_shape_function_hexa201_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa201_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa201_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Hexa201(iGauss, nShape);
}
void mw_shape_function_hexa202_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa202_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa202_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Hexa202(iGauss, nShape);
}
void mw_shape_function_hexa203_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa203_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa203_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Hexa203(iGauss, nShape);
}
void mw_shape_function_tetra41_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra41_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Tetra41(iGauss, nShape);
}
void mw_shape_function_tetra101_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra101_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra101_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Tetra101(iGauss, nShape);
}
void mw_shape_function_tetra104_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra104_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra104_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Tetra104(iGauss, nShape);
}
void mw_shape_function_tetra1015_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra1015_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra1015_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Tetra1015(iGauss, nShape);
}
void mw_shape_function_prism62_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Prism62(iGauss, nShape);
}
void mw_shape_function_prism156_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Prism156(iGauss, nShape);
}
void mw_shape_function_prism159_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism159_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism159_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Prism159(iGauss, nShape);
}
void mw_shape_function_prism1518_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism1518_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism1518_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;
    
    *N = pMW->ShapeFunc_Prism1518(iGauss, nShape);
}
void mw_shape_function_quad41_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad41_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Quad41(iGauss, nShape);
}
void mw_shape_function_quad84_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad84_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad84_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Quad84(iGauss, nShape);
}
void mw_shape_function_quad89_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad89_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad89_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Quad89(iGauss, nShape);
}
void mw_shape_function_tri31_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri31_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri31_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Triangle31(iGauss, nShape);
}
void mw_shape_function_tri63_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri63_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri63_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;

    *N = pMW->ShapeFunc_Triangle63(iGauss, nShape);
}
void mw_shape_function_line21_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line21_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line21_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;
    
    *N = pMW->ShapeFunc_Line21(iGauss, nShape);
}
void mw_shape_function_line32_(iint* igauss, iint* ishape, double* N)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line32_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line32_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint nShape = *ishape;
    
    *N = pMW->ShapeFunc_Line32(iGauss, nShape);
}


//----
// shape function deriv (rst coord)
//----
void mw_dndr_(iint* shape_type, double dndr[])
{
    if(*shape_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_, minus value is invalid argument. shape_type");
        return;
    }

    uiint nShapeType = *shape_type;

    pMW->dNdr(nShapeType, dndr);
}
void mw_dndr_hexa81_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa81_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa81_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa81_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Hexa81_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_hexa82_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa82_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa82_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa82_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Hexa82_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_hexa201_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa201_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa201_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa201_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Hexa201_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_hexa202_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa202_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa202_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa202_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Hexa202_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_hexa203_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa203_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa203_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa203_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Hexa203_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tetra41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra41_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra41_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Tetra41_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tetra101_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra101_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra101_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra101_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;
    
    *dndr = pMW->dNdr_Tetra101_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tetra104_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra104_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra104_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra104_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Tetra104_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tetra1015_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra1015_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra1015_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra1015_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Tetra1015_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_prism62_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism62_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism62_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism62_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Prism62_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_prism156_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism156_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism156_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Prism156_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_prism159_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism159_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism159_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism159_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Prism159_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_prism1518_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism1518_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism1518_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism1518_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Prism1518_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_quad41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad41_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad41_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Quad41_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_quad84_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad84_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad84_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad84_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Quad84_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_quad89_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad89_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad89_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad89_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Quad89_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tri31_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri31_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri31_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri31_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;

    *dndr = pMW->dNdr_Tri31_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_tri63_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri63_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri63_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri63_, minus value is invalid argument. iaxis");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    uiint iAxis  = *iaxis;
    
    *dndr = pMW->dNdr_Tri63_on_pt_on_shape(iGauss, iShape, iAxis);
}
void mw_dndr_line21_(iint* igauss, iint* ishape, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line21_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line21_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;

    *dndr = pMW->dNdr_Line21_on_pt_on_shape(iGauss, iShape);
}
void mw_dndr_line32_(iint* igauss, iint* ishape, double* dndr)
{
    if(*igauss < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line32_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line32_, minus value is invalid argument. ishape");
        return;
    }

    uiint iGauss = *igauss;
    uiint iShape = *ishape;

    *dndr = pMW->dNdr_Line32_on_pt_on_shape(iGauss, iShape);
}

//----
// shape function deriv (xyz coord)
//----
void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[])// [igauss][ishape][iaxis] : 積分点数 節点数 3(座標)
{
    if(*elem_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndx_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndx_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*ielem  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndx_, minus value is invalid argument. ielem");
        return;
    }

    uiint nElemType = *elem_type;
    uiint nNumOfInteg = *num_of_integ;
    uiint iElem  = *ielem;

    pMW->dNdx(nElemType, nNumOfInteg, iElem, dndx);
}
void mw_det_jacobian_(iint* elem_type, iint* num_of_integ, iint* igauss, double* det_j)
{
    if(*elem_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_det_jacobian_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_det_jacobian_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*igauss  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_det_jacobian_, minus value is invalid argument. igauss");
        return;
    }

    uiint nElemType = *elem_type;
    uiint nNumOfInteg = *num_of_integ;
    uiint iGauss  = *igauss;

    pMW->detJacobian(nElemType, nNumOfInteg, iGauss, *det_j);
}
void mw_weight_(iint* elem_type, iint* num_of_integ, iint* igauss, double* w)
{
    if(*elem_type < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*igauss  < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. igauss");
        return;
    }

    uiint nElemType = *elem_type;
    uiint nNumOfInteg = *num_of_integ;
    uiint iGauss  = *igauss;

    pMW->Weight(nElemType, nNumOfInteg, iGauss, *w);
}

//----
// shape function type
//----
iint mw_shapetype_hexa81_()
{
    return pMW->shapetype_hexa81();
}
iint mw_shapetype_hexa82_()
{
    return pMW->shapetype_hexa82();
}
iint mw_shapetype_hexa201_()
{
    return pMW->shapetype_hexa201();
}
iint mw_shapetype_hexa202_()
{
    return pMW->shapetype_hexa202();
}
iint mw_shapetype_hexa203_()
{
    return pMW->shapetype_hexa203();
}
iint mw_shapetype_tetra41_()
{
    return pMW->shapetype_tetra41();
}
iint mw_shapetype_tetra101_()
{
    return pMW->shapetype_tetra101();
}
iint mw_shapetype_tetra104_()
{
    return pMW->shapetype_tetra104();
}
iint mw_shapetype_tetra1015_()
{
    return pMW->shapetype_tetra1015();
}
iint mw_shapetype_prism62_()
{
    return pMW->shapetype_prism62();
}
iint mw_shapetype_prism156_()
{
    return pMW->shapetype_prism156();
}
iint mw_shapetype_prism159_()
{
    return pMW->shapetype_prism159();
}
iint mw_shapetype_prism1518_()
{
    return pMW->shapetype_prism1518();
}
iint mw_shapetype_quad41_()
{
    return pMW->shapetype_quad41();
}
iint mw_shapetype_quad84_()
{
    return pMW->shapetype_quad84();
}
iint mw_shapetype_quad89_()
{
    return pMW->shapetype_quad89();
}
iint mw_shapetype_tri31_()
{
    return pMW->shapetype_tri31();
}
iint mw_shapetype_tri63_()
{
    return pMW->shapetype_tri63();
}
iint mw_shapetype_line21_()
{
    return pMW->shapetype_line21();
}
iint mw_shapetype_line32_()
{
    return pMW->shapetype_line32();
}


//--
// boundary mesh
//--
// number of boudary_mesh
iint mw_get_num_of_boundary_bnode_mesh_()
{ 
    uiint nNumOfBNodeMesh= pMW->GetNumOfBoundaryNodeMesh();

    if(nNumOfBNodeMesh < IINT_MAX){
        return nNumOfBNodeMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bface_mesh_()
{ 
    uiint nNumOfBFaceMesh= pMW->GetNumOfBoundaryFaceMesh();

    if(nNumOfBFaceMesh < IINT_MAX){
        return nNumOfBFaceMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bedge_mesh_()
{
    uiint nNumOfBEdgeMesh = pMW->GetNumOfBoundaryEdgeMesh();

    if(nNumOfBEdgeMesh < IINT_MAX){
        return nNumOfBEdgeMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bvolume_mesh_()
{ 
    uiint nNumOfBVolMesh= pMW->GetNumOfBoundaryVolumeMesh();

    if(nNumOfBVolMesh < IINT_MAX){
        return nNumOfBVolMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
// BND type for each boundary_mesh { Neumann || Dirichlet }
iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    return pMW->GetBNDType_BNodeMesh(iBMesh);
}
iint mw_get_bnd_type_bface_mesh_(iint*ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    return pMW->GetBNDType_BFaceMesh(iBMesh);
}
iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    return pMW->GetBNDType_BEdgeMesh(iBMesh);
}
iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh)
{ 
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    return pMW->GetBNDType_BVolumeMesh(iBMesh);
}
// BND type number
iint mw_get_neumann_type_(){ return pMW->getNeumannType();}
iint mw_get_dirichlet_type_(){ return pMW->getDirichletType();}

// number of bnode for each boundary_mesh
iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh)
{ 
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfBNode= pMW->GetNumOfBNode_BNodeMesh(iBMesh);

    if(nNumOfBNode < IINT_MAX){
        return nNumOfBNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh)
{ 
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfBNode= pMW->GetNumOfBNode_BFaceMesh(iBMesh);

    if(nNumOfBNode < IINT_MAX){
        return nNumOfBNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    
    uiint nNumOfBNode= pMW->GetNumOfBNode_BEdgeMesh(iBMesh);

    if(nNumOfBNode < IINT_MAX){
        return nNumOfBNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh)
{ 
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfBNode= pMW->GetNumOfBNode_BVolumeMesh(iBMesh);

    if(nNumOfBNode < IINT_MAX){
        return nNumOfBNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

// number of DOF for each boundary_mesh
iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;

    uiint nNumOfDOF= pMW->GetNumOfDOF_BNodeMesh(iBMesh, iBNode);

    if(nNumOfDOF < IINT_MAX){
        return nNumOfDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    
    uiint iBMesh= *ibmesh;

    uiint nNumOfDOF= pMW->GetNumOfDOF_BFaceMesh(iBMesh);

    if(nNumOfDOF < IINT_MAX){
        return nNumOfDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfDOF= pMW->GetNumOfDOF_BEdgeMesh(iBMesh);

    if(nNumOfDOF < IINT_MAX){
        return nNumOfDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    uiint nNumOfDOF= pMW->GetNumOfDOF_BVolumeMesh(iBMesh);

    if(nNumOfDOF < IINT_MAX){
        return nNumOfDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

// DOF number for each boundary_mesh ( DOF_index => DOF Number )
iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. idof");
        return ERROR;
    }
    
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint iDOF= *idof;

    uiint nDOF= pMW->GetDOF_BNodeMesh(iBMesh, iBNode, iDOF);

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, minus value is invalid argument. idof");
        return ERROR;
    }
    
    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;

    uiint nDOF= pMW->GetDOF_BFaceMesh(iBMesh, iDOF);

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, minus value is invalid argument. idof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;

    uiint nDOF= pMW->GetDOF_BEdgeMesh(iBMesh, iDOF);

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*idof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, minus value is invalid argument. idof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;
    
    uiint nDOF= pMW->GetDOF_BVolumeMesh(iBMesh, iDOF);

    if(nDOF < IINT_MAX){
        return nDOF;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//--
// value of boundary node
//--
double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;

    return pMW->GetBNodeValue_BNodeMesh(iBMesh, iBNode, nDOF);
}
double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;

    return pMW->GetBNodeValue_BFaceMesh(iBMesh, iBNode, nDOF, *mglevel);
}
double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. dof");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;

    return pMW->GetBNodeValue_BEdgeMesh(iBMesh, iBNode, nDOF, *mglevel);
}
double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. dof");
        return ERROR;
    }
    
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;
    
    return pMW->GetBNodeValue_BVolumeMesh(iBMesh, iBNode, nDOF, *mglevel);
}

// node id (in boundary node)
iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    
    uiint nNodeID= pMW->GetNodeID_BNode_BNodeMesh(iBMesh, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;

    uiint nNodeID= pMW->GetNodeID_BNode_BFaceMesh(iBMesh, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;

    uiint nNodeID= pMW->GetNodeID_BNode_BEdgeMesh(iBMesh, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    
    uiint nNodeID= pMW->GetNodeID_BNode_BVolumeMesh(iBMesh, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//--
// value of boundary face, edge, volume
//--
iint mw_get_num_of_bface_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bface_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfBFace= pMW->GetNumOfBFace(iBMesh);

    if(nNumOfBFace < IINT_MAX){
        return nNumOfBFace;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibface < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. ibface");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. dof");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;
    uiint nDOF= *dof;
    
    return pMW->GetBFaceValue(iBMesh, iBFace, nDOF);
}
iint mw_get_num_of_bedge_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bedge_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    uiint nNumOfBEdge= pMW->GetNumOfBEdge(iBMesh);

    if(nNumOfBEdge < IINT_MAX){
        return nNumOfBEdge;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibedge < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. ibedge");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBEdge= *ibedge;
    uiint nDOF= *dof;

    return pMW->GetBEdgeValue(iBMesh, iBEdge, nDOF);
}
iint mw_get_num_of_bvolume_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bvolume_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    uiint nNumOfBVol= pMW->GetNumOfBVolume(iBMesh);

    if(nNumOfBVol < IINT_MAX){
        return nNumOfBVol;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibvol < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. ibvol");
        return ERROR;
    }
    if(*dof < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. dof");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;
    uiint nDOF= *dof;

    return pMW->GetBVolumeValue(iBMesh, iBVol, nDOF);
}
//--
// node_id, face, edge, volume
//--
iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface)//BFaceのノード数
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bface_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibface < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bface_, minus value is invalid argument. ibface");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;

    return pMW->GetNumOfNode_BFace(iBMesh, iBFace);
}
iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibface < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibface");
        return ERROR;
    }
    if(*ibnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibnode");
        return ERROR;
    }
    
    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;
    uiint iBNode= *ibnode;

    uiint nNodeID= pMW->GetNodeID_BFace(iBMesh, iBFace, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge)//BEdgeのノード数
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bedge_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibedge < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bedge_, minus value is invalid argument. ibedge");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBEdge= *ibedge;

    return pMW->GetNumOfNode_BEdge(iBMesh, iBEdge);
}
iint mw_get_node_id_bedge_(iint* ibmesh, iint* ibedge, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, minus value is invalid argument. ibedge");
        return ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, minus value is invalid argument. ibnode");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBEdge= *ibedge;
    uiint iBNode= *ibnode;
    
    uiint nNodeID= pMW->GetNodeID_BEdge(iBMesh, iBEdge, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol)//BVolumeのノード数
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bvolume_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibvol < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bvolume_, minus value is invalid argument. ibvol");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;

    return pMW->GetNumOfNode_BVolume(iBMesh, iBVol);
}
iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibvol < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibvol");
        return ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibnode");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;
    uiint iBNode= *ibnode;

    uiint nNodeID= pMW->GetNodeID_BVolume(iBMesh, iBVol, iBNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//--
// boundary_mesh name
//--
iint mw_get_bnode_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_namelength_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    return pMW->GetBNodeMesh_NameLength(iBMesh);
}
void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_name_, minus value is invalid argument. name_len");
        return;
    }

    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;

    string sName = pMW->GetBNodeMesh_Name(iBMesh);
    
    uiint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > nNameLen){
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}
iint mw_get_bface_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_namelength_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;
    
    return pMW->GetBFaceMesh_NameLength(iBMesh);
}
void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_name_, minus value is invalid argument. name_len");
        return;
    }

    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;
    
    string sName = pMW->GetBFaceMesh_Name(iBMesh);
    
    uiint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > nNameLen){
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}
iint mw_get_bvolume_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_namelength_, minus value is invalid argument. ibmesh");
        return ERROR;
    }

    uiint iBMesh= *ibmesh;

    return pMW->GetBVolumeMesh_NameLength(iBMesh);
}
void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_name_, minus value is invalid argument. name_len");
        return;
    }

    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;

    string sName = pMW->GetBVolumeMesh_Name(iBMesh);
    
    uiint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > nNameLen){
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}
iint mw_get_bedge_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_namelength_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    uiint iBMesh= *ibmesh;

    return pMW->GetBEdgeMesh_NameLength(iBMesh);
}
void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_name_, minus value is invalid argument. name_len");
        return;
    }

    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;

    string sName = pMW->GetBEdgeMesh_Name(iBMesh);
    
    uiint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > nNameLen){
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}
//--
// entity_id of boundary_mesh (for FrontISTR)
//--
iint mw_get_edge_id_bedge_(iint* ibmesh, iint* ibedge)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibedge < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, minus value is invalid argument. ibedge");
        return ERROR;
    }

    uiint iBMesh = *ibmesh;
    uiint iBEdge = *ibedge;

    uiint nEdgeID= pMW->GetEdgeID_BEdge(iBMesh, iBEdge);

    if(nEdgeID < IINT_MAX){
        return nEdgeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bedge_(iint*ibmesh, iint* ibedge)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibedge < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, minus value is invalid argument. ibedge");
        return ERROR;
    }

    uiint iBMesh = *ibmesh;
    uiint iBEdge = *ibedge;

    uiint nElemID= pMW->GetElemID_BEdge(iBMesh, iBEdge);

    if(nElemID < IINT_MAX){
        return nElemID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibface < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, minus value is invalid argument. ibface");
        return ERROR;
    }

    uiint iBMesh = *ibmesh;
    uiint iBFace = *ibface;

    uiint nFaceID= pMW->GetFaceID_BFace(iBMesh, iBFace);

    if(nFaceID < IINT_MAX){
        return nFaceID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibface < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, minus value is invalid argument. ibface");
        return ERROR;
    }

    uiint iBMesh = *ibmesh;
    uiint iBFace = *ibface;

    uiint nElemID= pMW->GetElemID_BFace(iBMesh, iBFace);

    if(nElemID < IINT_MAX){
        return nElemID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol)
{
    if(*ibmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, minus value is invalid argument. ibmesh");
        return ERROR;
    }
    if(*ibvol < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, minus value is invalid argument. ibvol");
        return ERROR;
    }

    uiint iBMesh = *ibmesh;
    uiint iBVol = *ibvol;

    uiint nElemID= pMW->GetElemID_BVolume(iBMesh, iBVol);

    if(nElemID < IINT_MAX){
        return nElemID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
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

int mw_get_rank_()
{
    return pMW->GetRank();
}
int mw_get_num_of_process_()
{
    return pMW->GetNumOfProcess();
}

void mw_allreduce_r_(double val[], int* val_size, int* op)
{
    if(*val_size < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allreduce_r_, minus value is invalid argument. val_size");
        return;
    }

    double rval[*val_size];

    pMW->AllReduce(val, rval, *val_size, MPI_DOUBLE, *op, MPI_COMM_WORLD);
}
void mw_allreduce_i_(iint val[], int* val_size, int* op)
{
    if(*val_size < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allreduce_i_, minus value is invalid argument. val_size");
        return;
    }

    iint rval[*val_size];

    pMW->AllReduce(val, rval, *val_size, MPI_INT, *op, MPI_COMM_WORLD);
}
int mw_barrier_()
{
    return pMW->Barrier(MPI_COMM_WORLD);
}
int mw_abort_(iint* error)
{
    return pMW->Abort(MPI_COMM_WORLD, *error);
}
int mw_allgather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    
    return pMW->AllGather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, MPI_COMM_WORLD);
}
int mw_allgather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_i_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_i_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    
    return pMW->AllGather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, MPI_COMM_WORLD);
}
int mw_gather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }

    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
int mw_gather_i_(int sendbuf[], iint* sendcnt, int recvbuf[], int* recvcnt, int* root)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_i_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_i_, minus value is invalid argument. sendcnt");
        return ERROR;
    }

    return pMW->Gather((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, MPI_COMM_WORLD);
}
int mw_scatter_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_r_, minus value is invalid argument. sendcnt");
        return ERROR;
    }

    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
int mw_scatter_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root)
{
    if(*sendcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_i_, minus value is invalid argument. sendcnt");
        return ERROR;
    }
    if(*recvcnt < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_i_, minus value is invalid argument. recvcnt");
        return ERROR;
    }
    
    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_INT, (void*)recvbuf, *recvcnt, MPI_INT, *root, MPI_COMM_WORLD);
}
int mw_bcast_i_(int buf[], int* cnt, int* root)
{
    return pMW->Bcast(buf, *cnt, MPI_INT, *root, MPI_COMM_WORLD);
}
int mw_bcast_r_(double buf[], int* cnt, int* root)
{
    return pMW->Bcast(buf, *cnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
int mw_bcast_s_(char buf[], int* cnt, int* root)
{
    return pMW->Bcast(buf, *cnt, MPI_CHAR, *root, MPI_COMM_WORLD);
}

//
// 以下の4メソッドは、ペア NumOfNeibPE--getTransRank--Send_Recv_R, _I
//
int mw_get_num_of_neibpe_(int* imesh)//Meshパーツが通信する相手の数
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_neibpe_, minus value is invalid argument. imesh");
        return ERROR;
    }
    
    uiint iMesh = *imesh;
    
    uiint nNumOfNeibPE= pMW->GetNumOfNeibPE(iMesh);

    if(nNumOfNeibPE < IINT32_MAX){
        return nNumOfNeibPE;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_neibpe_, return value is over INT_MAX");
        return IINT32_MAX;
    }
}
int mw_get_transrank_(int* imesh, int* ipe)//通信Mesh毎のランク番号
{
    if(*imesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, minus value is invalid argument. imesh");
        return ERROR;
    }
    if(*ipe < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, minus value is invalid argument. ipe");
        return ERROR;
    }

    int iMesh = *imesh;
    int iPE = *ipe;

    uiint nTransRank= pMW->GetTransRank(iMesh, iPE);

    if(nTransRank < IINT32_MAX){
        return nTransRank;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, return value is over INT_MAX");
        return IINT32_MAX;
    }
}
void mw_send_recv_r_(double buf[], int* num_of_node, int* dof_size, int* trans_rank)//bufの値を送信-受信
{
    if(*num_of_node < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. num_of_node");
        return;
    }
    if(*dof_size < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. dof_size");
        return;
    }
    if(*trans_rank < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. trans_rank");
        return;
    }

    int nNumOfNode = *num_of_node;
    int nNumOfDOF = *dof_size;
    int transRank = *trans_rank;

    pMW->Send_Recv_R(buf, nNumOfNode, nNumOfDOF, transRank);
}
void mw_send_recv_i_(int buf[], int* num_of_node, int* dof_size, int* trans_rank)//bufの値を送信-受信
{
    if(*num_of_node < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. num_of_node");
        return;
    }
    if(*dof_size < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. dof_size");
        return;
    }
    if(*trans_rank < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. trans_rank");
        return;
    }

    int nNumOfNode = *num_of_node;
    int nNumOfDOF = *dof_size;
    int transRank = *trans_rank;

    pMW->Send_Recv_I(buf, nNumOfNode, nNumOfDOF, transRank);
}

////void mw_send_recv_r2_(double buf[], iint* dof_size)// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
////{
////    if(*dof_size < 0){
////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
////        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_r2_, minus value is invalid argument. dof_size");
////        return;
////    }
////    uiint nNumOfDOF= *dof_size;
////
////    pMW->Send_Recv_R(buf, nNumOfDOF);
////}
////void mw_send_recv_r_()// 通信Nodeの値を入れ替えて更新
////{
////    pMW->Send_Recv_R();
////}



//--
// CommMesh2 (通信Mesh) { select された AssyModel,Meshを対象 }
// CommNode (通信Node) :  Visualizer 用途
//--
iint mw_get_num_of_comm_mesh_()
{
    uiint nNumOfCommMesh= pMW->GetNumOfCommMesh();

    if(nNumOfCommMesh < IINT_MAX){
        return nNumOfCommMesh;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_comm_node_(iint* icmesh)
{
    if(*icmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_node_, minus value is invalid argument. icmesh");
        return ERROR;
    }
    uiint iComMesh= *icmesh;

    uiint nNumOfCommNode= pMW->GetNumOfCommNode(iComMesh);

    if(nNumOfCommNode < IINT_MAX){
        return nNumOfCommNode;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode)//MeshのNodeID
{
    if(*icmesh < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, minus value is invalid argument. icmesh");
        return ERROR;
    }
    if(*icnode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, minus value is invalid argument. icnode");
        return ERROR;
    }
    uiint iComMesh= *icmesh;
    uiint iComNode= *icnode;

    uiint nNodeID= pMW->GetNodeID_CommNode(iComMesh, iComNode);

    if(nNodeID < IINT_MAX){
        return nNodeID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}



//--
// Element_Group { select AssyModel, select Mesh }
//--
iint mw_get_num_of_elementgroup_()
{
    uiint nNumOfElemGrp= pMW->GetNumOfElementGroup();

    if(nNumOfElemGrp < IINT_MAX){
        return nNumOfElemGrp;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_elementgroup_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_id_(iint* igrp)
{
    if(*igrp < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_id_, minus value is invalid argument. igrp");
        return ERROR;
    }
    uiint iGrp = *igrp;

    uiint nNumOfElem= pMW->GetNumOfElementID(iGrp);

    if(nNumOfElem < IINT_MAX){
        return nNumOfElem;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_id_, INT_MAX over value");
        return IINT_MAX;
    }
}
iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index)
{
    if(*igrp < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, minus value is invalid argument. igrp");
        return ERROR;
    }
    if(*index < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, minus value is invalid argument. index");
        return ERROR;
    }
    uiint iGrp = *igrp;
    uiint i = *index;

    uiint nElemID= pMW->GetElementID_with_ElementGroup(iGrp, i);

    if(nElemID < IINT_MAX){
        return nElemID;
    }else{
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elementgroup_name_length_(iint* igrp)
{
    if(*igrp < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_length_, minus value is invalid argument. igrp");
        return ERROR;
    }

    uiint iGrp = *igrp;
    iint nLength = pMW->GetElementGroupName_Length(iGrp) + 1;

    return nLength;
}
void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len)
{
    if(*igrp < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_, minus value is invalid argument. igrp");
        return;
    }
    if(*name_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_, minus value is invalid argument. name_len");
        return;
    }

    uiint iGrp = *igrp;
    uiint nNameLen= *name_len;

    string sName = pMW->GetElementGroupName(iGrp);
    
    uiint i, nNumOfChar = sName.length();
    
    if(nNumOfChar > nNameLen){
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++){
            name[i] = sName[i];
        };
        name[nLength]='\0';
    }else{
        for(i=0; i < nNumOfChar; i++){
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++){
            name[i] = '\0';
        };
    }
}

//----
// logger
//----
void mw_logger_set_mode_(iint* mode)
{
    if(*mode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_mode_, minus value is invalid argument. mode");
        return;
    }
    uiint nMode= *mode;

    pMW->LoggerMode(nMode);
}
void mw_logger_set_device_(iint* mode, iint* device)
{
    if(*mode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_device_, minus value is invalid argument. mode");
        return;
    }
    if(*device < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_device_, minus value is invalid argument. device");
        return;
    }
    uiint nMode= *mode;
    uiint nDevice= *device;

    pMW->LoggerDevice(nMode, nDevice);
}
void mw_logger_info_ (iint* mode, char* message, iint* str_len)
{
    if(*mode < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_info_, minus value is invalid argument. mode");
        return;
    }
    if(*str_len < 0){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_info_, minus value is invalid argument. str_len");
        return;
    }
    uiint nMode= *mode;
    uiint nStrLen= *str_len;

    uiint nLength = nStrLen + 1;
    
    char cmsg[nLength];
    strncpy(cmsg, message, nLength);//Null追加
    
    pMW->LoggerInfo(nMode, cmsg);
}
//----
// logger parameter
//----
iint mw_get_error_mode_()
{
    return pMW->getErrorMode();
}
iint mw_get_warn_mode_()
{
    return pMW->getWarnMode();
}
iint mw_get_info_mode_()
{
    return pMW->getInfoMode();
}
iint mw_get_debug_mode_()
{
    return pMW->getDebugMode();
}
iint mw_get_disk_device_()
{
    return pMW->getDiskDevice();
}
iint mw_get_display_device_()
{
    return pMW->getDisplayDevice();
}








