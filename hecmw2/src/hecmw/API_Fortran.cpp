/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/API_Fortran.cpp
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
#ifdef MSVC
#include "HEC_MW3.hxx"
#include "API_Fortran.hxx"
#else
#include "HEC_MW3.h"
#include "API_Fortran.h"
#endif

pmw::CMW *pMW;
iint mw_initialize_(int* argc, char** argv)
{
    pMW = pmw::CMW::Instance();

    return pMW->Initialize(*argc, argv);
}
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
void mw_banner_()
{
    pMW->Banner_Display();
}
iint mw_revocap_refine_(char* filename, iint* n_refine)
{
#ifdef REVOCAP_REFINE //=========================================~REVOCAP_REFINE
    string sFileName= filename;

    if(*n_refine < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_revocap_refine, second arg is incorrect");
        return MW_ERROR;
    }
    uiint nNumRefine = *n_refine;

    return pMW->RevocapRefine(sFileName, nNumRefine);
#endif //=========================================================REVOCAP_REFINE

    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn, "mw_revocap_refine_ is not available.");

    return MW_ERROR;
}

iint mw_file_read_(char* basename)
{
    string sBasename= basename;
    return pMW->FileRead(sBasename, false);
}
iint mw_file_read_fstr_()
{
    return pMW->FileRead_fstr(false);
}
iint mw_file_debug_write_()
{
    return pMW->FileDebugWrite();
}
iint mw_get_fstr_filename_length_mesh_()
{
    string sName = pMW->getFstr_FileName_Mesh();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr mesh filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_control_()
{
    string sName = pMW->getFstr_FileName_Control();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr control filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_result_()
{
    string sName = pMW->getFstr_FileName_Result();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr result filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_restart_()
{
    string sName = pMW->getFstr_FileName_Restart();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr restart filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_part_in_()
{
    string sName = pMW->getFstr_FileName_PartIn();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr part_in filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_part_out_()
{
    string sName = pMW->getFstr_FileName_PartOut();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr part_out filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_mesh_()
{
    string sName = pMW->getFstr_FileName_VisMesh();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_mesh filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_in_()
{
    string sName = pMW->getFstr_FileName_VisIn();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_in filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_vis_out_()
{
    string sName = pMW->getFstr_FileName_VisOut();
    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr vis_out filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}
iint mw_get_fstr_filename_length_cadfit_()
{
    string sName= pMW->getFstr_RefineCADFitName();

    uiint nLength= sName.length();
    if(nLength > IINT_MAX) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "int size over : fstr cad_fit filename length");
        return MW_ERROR;
    } else {
        return nLength;
    }
}

void mw_get_fstr_filename_mesh_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_mesh, second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_Mesh();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_mesh, name length mismatch : fstr mesh filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_control_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_control, second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_Control();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_control, name length mismatch : fstr control filename length");
        return;
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_result_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_result, second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_Result();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_fstr_filename_result,  name length mismatch : fstr result filename length");
        return;
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_restart_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_Restart();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr restart filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_part_in_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_PartIn();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr part_in filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_part_out_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_PartOut();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr part_out filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_mesh_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_VisMesh();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_mesh filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_in_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_VisIn();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_in filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_vis_out_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }
    string sName = pMW->getFstr_FileName_VisOut();
    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr vis_out filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
void mw_get_fstr_filename_cadfit_(char name[], iint* len)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    if(*len < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "second arg is incorrect");
        return;
    }

    string sName= pMW->getFstr_RefineCADFitName();

    uiint i, nLength = sName.length();
    if(*len < nLength) {
        pLogger->Info(Utility::LoggerMode::Error, "name length mismatch : fstr cad_fit filename length");
    }

    for(i=0; i < nLength; i++) {
        name[i]=sName[i];
    }
}
iint mw_get_fstr_refine_num_()
{
    return pMW->getFstr_RefineNum();
}
iint mw_get_fstr_refine_type_()
{
    return pMW->getFstr_RefineType();
}



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
//
// format= %d:int32, %f:double(fixed), %e:double(scientific), %s:const char*
// --
// # Fortran引数は、ポインター
// --
iint mw_rlt_print_(iint* width, const char* format, ... )
{
    if(*width < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rlt_print_, minus value is invalid argument. width");
        return MW_ERROR;
    }
    vector<void*> param;
    vint    vnVal;
    vuint   vuVal;
    vdouble vdVal, veVal;
    vstring vsVal;

    uiint nLength= strlen(format);
    uiint nWidth = *width;

    va_list list;
    va_start( list, format);
    for(uiint i=0; i < nLength; i++) {
        if(format[i] == '%') {
            ++i;
            switch( format[i] ) {
            case('d'): {
                iint* nVal = va_arg(list, iint*);     //Fortran引数は、ポインター
                vnVal.push_back(*nVal);
            }
            break;
            case('u'): {
                uiint* uVal = va_arg(list, uiint*);     //Fortran引数は、ポインター
                vuVal.push_back(*uVal);
            }
            break;
            case('f'): {
                double* dVal= va_arg(list, double*);    //Fortran引数は、ポインター
                vdVal.push_back(*dVal);
            }
            break;
            case('e'): {
                double* dVal= va_arg(list, double*);    //Fortran引数は、ポインター
                veVal.push_back(*dVal);
            }
            break;
            case('s'): {
                string sVal= va_arg(list, const char*);
                vsVal.push_back(sVal);
            }
            break;
            default:
                break;
            }
        }
    };
    va_end( list );

    uiint icase_d(0), icase_u(0), icase_f(0), icase_e(0), icase_s(0);
    for(uiint i=0; i < nLength; i++) {
        if(format[i] == '%') {
            ++i;
            switch( format[i] ) {
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
                ++icase_s;
                break;
            }
        }
    };

    FileIO::CFileIO* pFileIO= FileIO::CFileIO::Instance();

    pFileIO->PrintResult(nWidth, format, param);

    return MW_SUCCESS;
}
void mw_rlt_end_()
{
    pMW->PrintRlt_End();
}
void mw_print_avs_basis_(iint* ieq)
{
    uiint iEq = *ieq;
    pMW->PrintMicroAVS_Basis(iEq);
}
void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *ndof < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nDOF = *ndof;
    pMW->recAVS_Label(iMesh, label, unit, nDOF);
}
void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *num_of_node < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nNumOfNode= *num_of_node;
    pMW->recAVS_Variable(iMesh, nNumOfNode, label, value);
}
void mw_print_avs_fem_()
{
    pMW->PrintMicroAVS_FEM();
}
//--
// vtk
//--
void mw_rec_vtk_label_(iint* imesh, char* label, char* unit, iint* ndof)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *ndof < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nDOF = *ndof;
    pMW->recVTK_Label(iMesh, label, unit, nDOF);
}
void mw_rec_vtk_variable_(iint* imesh, iint* num_of_node, char* label, double* value)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *num_of_node < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nNumOfNode= *num_of_node;
    pMW->recVTK_Variable(iMesh, nNumOfNode, label, value);
}
void mw_print_vtk_fem_()
{
    pMW->PrintVTK_FEM();
}
//--
// field view
//--
void mw_rec_uns_label_(iint* imesh, char* label, char* unit, iint* ndof)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *ndof < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nDOF = *ndof;
    pMW->recUNS_Label(iMesh, label, unit, nDOF);
}
void mw_rec_uns_variable_(iint* imesh, iint* num_of_node, char* label, double* value)
{
    if( *imesh < 0 ) {
        return;
    }
    if( *num_of_node < 0 ) {
        return;
    }
    uiint iMesh= *imesh;
    uiint nNumOfNode= *num_of_node;
    pMW->recUNS_Variable(iMesh, nNumOfNode, label, value);
}
void mw_print_uns_fem_()
{
    pMW->PrintUNS_FEM();
}
//--
// res
//--
iint mw_file_write_res_(iint* step)
{
    if(*step < 0) {
        return MW_ERROR;
    }
    uiint nStep= *step;
    return pMW->FileWriteRes(nStep, false);
}
iint mw_set_restart_(iint* step, iint* num_of_level, iint* num_of_algebra)
{
    if(*step < 0) {
        return MW_ERROR;
    }
    if(*num_of_level < 0) {
        return MW_ERROR;
    }
    if(*num_of_algebra < 0) {
        return MW_ERROR;
    }
    uiint nStep= *step, nNumOfLevel= *num_of_level, nNumEquation= *num_of_algebra;


    uiint nRetVal = pMW->SetRestart(nStep, nNumOfLevel, nNumEquation, false);

    return nRetVal;
}
iint mw_file_write_res_bin_(iint* step)
{
    if(*step < 0) {
        return MW_ERROR;
    }
    uiint nStep= *step;
    return pMW->FileWriteRes(nStep, true);
}
iint mw_set_restart_bin_(iint* step, iint* num_of_level, iint* num_of_algebra)
{
    if(*step < 0) {
        return MW_ERROR;
    }
    if(*num_of_level < 0) {
        return MW_ERROR;
    }
    if(*num_of_algebra < 0) {
        return MW_ERROR;
    }
    uiint nStep= *step, nNumOfLevel= *num_of_level, nNumEquation= *num_of_algebra;


    uiint nRetVal= pMW->SetRestart(nStep, nNumOfLevel, nNumEquation, true);

    return nRetVal;
}
//--
// 線形方程式生成:接合面なし
//--
void mw_gene_linear_algebra_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof)
{
    if(*num_of_algebra<0 || *global_num_of_mesh<0 ) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gene_linear_algebra_, minus value is invalid argument. num_of_algebra or num_of_mesh");
        return;
    }
    uiint nNumOfAlgebra = (uiint)*num_of_algebra;
    uiint nGlobalNumOfMesh = (uiint)*global_num_of_mesh;
    //
    // DOF:方程式-Mesh(パーツ) の順に配列
    //
    uiint** vvDOF= new uiint*[nNumOfAlgebra];
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        vvDOF[ieq]= new uiint[nGlobalNumOfMesh];
    };
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        for(uiint imesh=0; imesh < nGlobalNumOfMesh; imesh++) {
            vvDOF[ieq][imesh]= dof[ieq*nGlobalNumOfMesh + imesh];
        };
    };

    //--
    // Ax=b 生成
    //--
    pMW->GeneLinearAlgebra(nNumOfAlgebra, nGlobalNumOfMesh, vvDOF);

    // temp 配列 削除
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        delete[] vvDOF[ieq];
    };
    delete[] vvDOF;
}
//--
// 線形方程式生成:接合面あり
//--
void mw_gene_linear_algebra_assymodel_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof, iint* num_of_contact, double* transCoeff)
{
    if(*num_of_algebra<0 || *global_num_of_mesh<0 ) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gene_linear_algebra_assymodel_, minus value is invalid argument. num_of_algebra or num_of_mesh");
        return;
    }
    uiint nNumOfAlgebra = (uiint)*num_of_algebra;
    uiint nGlobalNumOfMesh = (uiint)*global_num_of_mesh;
    //
    // DOF:方程式-Mesh(パーツ) の順に配列
    //
    uiint** vvDOF= new uiint*[nNumOfAlgebra];
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        vvDOF[ieq]= new uiint[nGlobalNumOfMesh];
    };
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        for(uiint imesh=0; imesh < nGlobalNumOfMesh; imesh++) {
            vvDOF[ieq][imesh]= dof[ieq*nGlobalNumOfMesh + imesh];
        };
    };
    //
    // 伝達率:接合面-方程式 の順に配列
    //
    uiint nNumOfContact= (uiint)*num_of_contact;
    double** vvTransCoeff= new double*[nNumOfContact];
    for(uiint icont=0; icont < nNumOfContact; icont++) {
        vvTransCoeff[icont]= new double[nNumOfAlgebra];
    };
    for(uiint icont=0; icont < nNumOfContact; icont++) {
        for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
            vvTransCoeff[icont][ieq]= transCoeff[icont*nNumOfAlgebra + ieq];
        }
    }

    //--
    // Ax=b 生成
    //--
    pMW->GeneLinearAlgebra(nNumOfAlgebra, nGlobalNumOfMesh, vvDOF, vvTransCoeff);

    // temp 配列 削除
    for(uiint ieq=0; ieq < nNumOfAlgebra; ieq++) {
        delete[] vvDOF[ieq];
    };
    delete[] vvDOF;

    for(uiint icont=0; icont < nNumOfContact; icont++) {
        delete[] vvTransCoeff[icont];
    };
    delete[] vvTransCoeff;
}
void mw_select_algebra_(iint* ieq)
{
    if(*ieq < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_algebra_, minus value is invalid argument. ieq");
        return;
    }
    uiint iEqu = *ieq;
    pMW->SelectAlgebra(iEqu);
}
iint mw_matrix_add_elem_(iint* imesh,  iint* ielem,  double elem_matrix[])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_elem_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_elem_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    return pMW->Matrix_Add_Elem(iMesh, iElem, elem_matrix);
}
iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*i_index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. i_index");
        return MW_ERROR;
    }
    if(*j_index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_add_node_, minus value is invalid argument. j_index");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint inode = *i_index;
    uiint jnode = *j_index;
    return pMW->Matrix_Add_Node(iMesh, inode, jnode, nodal_matrix);
}
void mw_matrix_clear_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_clear_, minus value is invalid argument. imesh");
        return;
    }
    uiint iMesh = *imesh;
    pMW->Matrix_Clear(iMesh);
}
void mw_vector_clear_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_vector_clear_, minus value is invalid argument. imesh");
        return;
    }
    uiint iMesh = *imesh;
    pMW->Vector_Clear(iMesh);
}
void mw_assy_matrix_clear_()
{
    pMW->AssyMatrix_Clear();
}
void mw_assy_vector_clear_()
{
    pMW->AssyVector_Clear();
}



iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_24_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_24_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=24;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_60_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_60_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=60;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matrix[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_12_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_12_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=12;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_30_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_30_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=30;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_18_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_18_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=18;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_45_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_45_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=45;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_20_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_20_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=20;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_40_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_40_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=40;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_15_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_15_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=15;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_9_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_9_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=9;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_48_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_48_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=48;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_6_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_6_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=6;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10])
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_10_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ielem < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matirx_add_elem_10_, minus value is invalid argument. ielem");
        return MW_ERROR;
    }
    uiint nNumOfCol=10;
    double* mat = new double[nNumOfCol*nNumOfCol];
    for(uiint i=0; i < nNumOfCol; i++)
        for(uiint ii=0; ii < nNumOfCol; ii++)
            mat[nNumOfCol*i+ii] = elem_matirx[i][ii];
    uiint iMesh = *imesh;
    uiint iElem = *ielem;
    uiint result = pMW->Matrix_Add_Elem(iMesh, iElem, mat);
    delete [] mat;
    return result;
}
iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* solval)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc2_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF  = *dof;
    double dDiagValue = *diagval;
    double dSolValue = *solval;

    return pMW->Set_BC_Mat_RHS2(iMesh, iNodeID, nDOF, dDiagValue, dSolValue);
}
iint mw_matrix_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* rhsval)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_matrix_rhs_set_bc_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF  = *dof;
    double dDiagValue = *diagval;
    double dRHSValue = *rhsval;

    return pMW->Set_BC_Mat_RHS(iMesh, iNodeID, nDOF, dDiagValue, dRHSValue);
}
iint mw_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_set_bc_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF = *dof;
    double Value = *value;

    return pMW->Set_BC_RHS(iMesh, iNodeID, nDOF, Value);
}
iint mw_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_rhs_add_bc_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF = *dof;
    double Value = *value;

    return pMW->Add_BC_RHS(iMesh, iNodeID, nDOF, Value);
}

//--
// 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.
//--
iint mw_nl_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_set_bc_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_set_bc_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_set_bc_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF = *dof;
    double Value = *value;

    return pMW->Set_BC_NL_RHS( iMesh, iNodeID, nDOF, Value);
}
//--
// 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.
//--
iint mw_nl_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_add_bc_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_add_bc_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_nl_rhs_add_bc_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint iNodeID = *node_id;
    uiint nDOF = *dof;
    double Value = *value;

    return pMW->Add_BC_NL_RHS( iMesh, iNodeID, nDOF, Value);
}

//--
// solver
//--
iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition)
{
    if(*iter_max < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. iter_max");
        return MW_ERROR;
    }
    if(*method < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. method");
        return MW_ERROR;
    }
    if(*pre_condition < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_solve_, minus value is invalid argument. pre_condition");
        return MW_ERROR;
    }
    iint nIterMax = *iter_max;
    double dTolerance = *tolerance;
    iint nMethod = *method;
    iint nPreCondition = *pre_condition;

    return pMW->Solve(nIterMax, dTolerance, nMethod, nPreCondition);
}
void mw_get_solution_vector_(double buf[], iint* imesh)
{
    if(*imesh < 0) {
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
void mw_get_rhs_vector_(double buf[], iint* imesh)
{
    if(*imesh < 0) {
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

void mw_get_rhs_load_(double buf[], iint* imesh)// RHSベクトルをsumupしたベクトル
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_vector_, minus value is invalid argument. imesh");
        return;
    }
    uiint iMesh = *imesh;

    pMW->GetRHS_Load(buf, iMesh);
}
void mw_get_rhs_assy_load_(double buf[])// RHSアセンブルベクトルをsumupしたベクトル
{
    pMW->GetRHS_AssyLoad(buf);
}

double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*inode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_val_, minus value is invalid argument. idof");
        return MW_ERROR;
    }

    return pMW->GetSolutionVector_Val(*imesh, *inode, *idof);
}
double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*inode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. inode");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_val_, minus value is invalid argument. idof");
        return MW_ERROR;
    }

    return pMW->GetRHSVector_Val(*imesh, *inode, *idof);
}
iint mw_get_solution_assy_vector_dof_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_dof_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    uiint iMesh= (uiint)*imesh;

    uiint nDOF= pMW->GetSolutionVector_DOF(iMesh);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_solution_assy_vector_dof_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_rhs_assy_vector_dof_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_dof_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    uiint iMesh= (uiint)*imesh;

    uiint nDOF= pMW->GetRHSVector_DOF(iMesh);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_rhs_assy_vector_dof_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
//--
// Debug : AssyMatrix dump
//--
void mw_dump_assy_matrix_()
{
    pMW->dump_AssyMatrix();
}

void mw_dump_rhs_assy_vector_()
{
    pMW->dump_RHSAssyVector();
}


//////--
////// update, sumup in vec[],
//////--
////void mw_update_add_(iint* mglevel, iint* imesh, double vec[], iint* num_of_dof)
////{
////    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0){
////        return;
////    }
////    uiint iLevel = *mglevel;
////    uiint iMesh = *imesh;
////    uiint nNumDOF = *num_of_dof;
////
////    pMW->Update_Add(iLevel, iMesh, vec, nNumDOF);
////}
////void mw_update_(iint* mglevel, iint* imesh, double vec[], iint* num_of_dof)
////{
////    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0){
////        return;
////    }
////    uiint iLevel = *mglevel;
////    uiint iMesh = *imesh;
////    uiint nNumDOF = *num_of_dof;
////
////    pMW->Update(iLevel, iMesh, vec, nNumDOF);
////}
void mw_sumup_(iint* mglevel, iint* imesh, double vec[], iint* num_of_dof)
{
    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0){
        return;
    }
    uiint iLevel = *mglevel;
    uiint iMesh = *imesh;
    uiint nNumDOF = *num_of_dof;

    pMW->Sumup(iLevel, iMesh, vec, nNumDOF);
}
////void mw_update_average_(iint* mglevel, iint* imesh, double vec[], iint* num_of_dof)
////{
////    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0){
////        return;
////    }
////    uiint iLevel = *mglevel;
////    uiint iMesh = *imesh;
////    uiint nNumDOF = *num_of_dof;
////
////    pMW->Update_Average(iLevel, iMesh, vec, nNumDOF);
////}
////void mw_update_add_rhs_(iint* mglevel, iint* imesh, iint* ieq, iint* num_of_dof)
////{
////    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0 || *ieq < 0 ){
////        return;
////    }
////    uiint iLevel = *mglevel;
////    uiint iMesh = *imesh;
////    uiint nNumDOF = *num_of_dof;
////    uiint iEq = *ieq;
////
////    pMW->Update_Add(iLevel, iMesh, iEq, nNumDOF, RHS_VEC);
////
////}
////void mw_update_average_rhs_(iint* mglevel, iint* imesh, iint* ieq, iint* num_of_dof)
////{
////    if(*mglevel < 0 || *imesh < 0 || *num_of_dof < 0 || *ieq < 0 ){
////        return;
////    }
////    uiint iLevel = *mglevel;
////    uiint iMesh = *imesh;
////    uiint nNumDOF = *num_of_dof;
////    uiint iEq = *ieq;
////
////    pMW->Update_Average(iLevel, iMesh, iEq, nNumDOF, RHS_VEC);
////}
////void mw_update_assy_add_(iint* mglevel, double assy_vec[], iint vdof[], iint* dof_size)
////{
////    if(*mglevel < 0 || *dof_size < 0 ){
////        return;
////    }
////    uiint iLevel= *mglevel;
////    uiint nNumDOF = *dof_size;
////    uiint* vDOF = (uiint*)calloc(nNumDOF, sizeof(uiint));
////
////    for(uiint i=0; i < nNumDOF; i++) vDOF[i]=vdof[i];
////
////    pMW->Update_Assy_Add(iLevel, assy_vec, vDOF);
////
////    free( vDOF );
////}
////void mw_update_assy_(iint* mglevel, double assy_vec[], iint vdof[], iint* dof_size)
////{
////    if(*mglevel < 0 || *dof_size < 0 ){
////        return;
////    }
////    uiint iLevel= *mglevel;
////    uiint nNumDOF = *dof_size;
////    uiint* vDOF = (uiint*)calloc(nNumDOF, sizeof(uiint));
////
////    for(uiint i=0; i < nNumDOF; i++) vDOF[i]=vdof[i];
////
////    pMW->Update_Assy(iLevel, assy_vec, vDOF);
////
////    free( vDOF );
////}
////void mw_sumup_assy_(iint* mglevel, double assy_vec[], iint vdof[], iint* dof_size)
////{
////    if(*mglevel < 0 || *dof_size < 0 ){
////        return;
////    }
////    uiint iLevel= *mglevel;
////    uiint nNumDOF = *dof_size;
////    uiint* vDOF = (uiint*)calloc(nNumDOF, sizeof(uiint));
////
////    pMW->Sumup_Assy(iLevel, assy_vec, vDOF);
////
////    free( vDOF );
////}
////void mw_update_assy_average_(iint* mglevel, double assy_vec[], iint vdof[], iint* dof_size)
////{
////    if(*mglevel < 0 || *dof_size < 0 ){
////        return;
////    }
////    uiint iLevel= *mglevel;
////    uiint nNumDOF = *dof_size;
////    uiint* vDOF = (uiint*)calloc(nNumDOF, sizeof(uiint));
////
////    for(uiint i=0; i < nNumDOF; i++) vDOF[i]=vdof[i];
////
////    pMW->Update_Assy_Average(iLevel, assy_vec, vDOF);
////
////    free( vDOF );
////}

//--
// y=Ax : x update, y sumup
//--
void mw_matvec_assy_(iint* mglevel, iint* ieq, double assy_x[], double assy_y[])
{
    if(*mglevel < 0 || *ieq < 0) {
        return ;
    }
    uiint iLevel= *mglevel;
    uiint iEq= *ieq;

    pMW->MatVec_Assy(iLevel, iEq, assy_x, assy_y);
}
void mw_matvec_(iint* mglevel, iint* imesh, iint* ieq, double x[], double y[])
{
    if(*mglevel < 0 || *imesh < 0 || *ieq < 0 ) {
        return ;
    }
    uiint iLevel= *mglevel;
    uiint iMesh = *imesh;
    uiint iEq= *ieq;

    pMW->MatVec(iLevel, iMesh, iEq, x, y);
}

//--
// mesh data construct
//--
iint mw_refine_(iint* num_of_refine)
{
    if(*num_of_refine < 0) {
        return MW_ERROR;
    }
    uiint nNumOfRefine = *num_of_refine;
    return pMW->Refine(nNumOfRefine);
}
iint mw_mg_construct_(iint* num_of_refine)
{
    if(*num_of_refine < 0) {
        return MW_ERROR;
    }
    uiint nNumOfRefine = *num_of_refine;
    return pMW->Refine(nNumOfRefine);
}
void mw_finalize_refine_()
{
    pMW->FinalizeRefine();
}
void mw_finalize_mg_construct_()
{
    pMW->FinalizeRefine();
}


iint mw_get_num_of_assemble_model_()
{
    uiint nNumOfLevel= pMW->GetNumOfAssembleModel();
    if(nNumOfLevel < IINT_MAX) {
        return nNumOfLevel;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_assemble_model_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
void mw_select_assemble_model_(iint* mglevel)
{
    if(*mglevel < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_assemble_model_, minus value is invalid argument. mglevel");
        return;
    }
    uiint nMGLevel = *mglevel;
    pMW->SelectAssembleModel(nMGLevel);
}
iint mw_get_num_of_mesh_part_()
{
    uiint nNumOfMesh= pMW->GetNumOfMeshPart();
    if(nNumOfMesh < IINT_MAX) {
        return nNumOfMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_mesh_part_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
void mw_select_mesh_part_with_id_(iint* mesh_id)
{
    if(*mesh_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_mesh_part_with_id_, minus value is invalid argument. mesh_id");
        return;
    }
    uiint nMeshID = *mesh_id;
    pMW->SelectMeshPart_ID(nMeshID);
}
void mw_select_mesh_part_(iint* index)
{
    if(*index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_mesh_part_, minus value is invalid argument. index");
        return;
    }
    uiint nIndex = *index;
    pMW->SelectMeshPart_IX(nIndex);
}
iint mw_get_mesh_part_id_(iint* mesh_index)//-- mesh_index to  mesh_id
{
    if(*mesh_index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_mesh_part_id_, minus value is invalid argument. mesh_index");
        return MW_ERROR;
    }
    uiint nIndex= *mesh_index;

    return pMW->GetMeshID_Num(nIndex);
}
iint mw_get_mesh_part_index_(iint* mesh_id)//-- mesh_id    to  mesh_index
{
    if(*mesh_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_mesh_part_index_, minus value is invalid argument. mesh_id");
        return MW_ERROR;
    }
    uiint nID= *mesh_id;

    return pMW->GetMeshIndex_Num(nID);
}

void mw_select_element_with_id_(iint* elem_id)
{
    if(*elem_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_select_element_with_id_, minus value is invalid argument. elem_id");
        return;
    }
    uiint nElemID = *elem_id;
    pMW->SelectElement_ID(nElemID);
}
void mw_select_element_(iint* index)
{
    if(*index < 0) {
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
void mw_get_element_vert_node_index_(iint v_node_index[])
{
    pMW->GetElementVertNodeIndex(v_node_index);
}
iint mw_get_num_of_element_edge_()
{
    return pMW->GetNumOfElementEdge();
}
iint mw_get_num_of_element_face_()
{
    return pMW->GetNumOfElementFace();
}
void mw_get_element_edge_node_id_(iint v_node_id[])
{
    pMW->GetElementEdgeNodeID(v_node_id);
}

//--
// パーティショナー
//--
iint mw_get_element_face_element_id_(iint* iface)
{
    if(*iface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_face_element_id_, minus value is invalid argument. index");
        return IINT_MIN;
    }

    uiint nFaceIndex=(uiint)*iface;

    uiint nElemID= pMW->GetElementFaceElementID(nFaceIndex);

    if(nElemID < IINT_MAX) {
        return (iint)nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_face_element_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_Of_element_edge_element_(iint* iedge)
{
    if(*iedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_edge_element_, minus value is invalid argument. index");
        return IINT_MIN;
    }

    uiint nEdgeIndex=(uiint)*iedge;

    uiint nNumOfEdgeElem= pMW->GetNumOfElementEdgeElement(nEdgeIndex);

    if(nNumOfEdgeElem < IINT_MAX) {
        return (iint)nNumOfEdgeElem;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_edge_element_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_element_edge_element_id_(iint* iedge, iint* i)
{
    if(*iedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_edge_element_id_, minus value is invalid argument. index");
        return IINT_MIN;
    }
    if(*i < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_edge_element_id_, minus value is invalid argument. index");
        return IINT_MIN;
    }

    uiint nEdgeIndex= (uiint)*iedge;
    uiint nIndex= (uiint)*i;

    uiint nElemID= pMW->GetElementEdgeElementID( nEdgeIndex, nIndex );

    if(nElemID < IINT_MAX) {
        return (iint)nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_edge_element_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
void mw_setup_neighbors_()
{
    pMW->setupNeighbors();
}
//------------------------


iint mw_get_node_id_(iint* index)
{
    if(*index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_, minus value is invalid argument. index");
        return MW_ERROR;
    }
    uiint nIndex= *index;
    uiint nNodeID= pMW->getNodeID(nIndex);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_element_id_(iint* index)
{
    if(*index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_, minus value is invalid argument. index");
        return MW_ERROR;
    }
    uiint nIndex= *index;
    uiint nElemID= pMW->getElementID(nIndex);
    if(nElemID < IINT_MAX) {
        return nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_index_(iint* id)
{
    if(*id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_index_, minus value is invalid argument. id");
        return MW_ERROR;
    }
    uiint nID= *id;
    uiint nIndex= pMW->getNodeIndex(nID);
    if(nIndex < IINT_MAX) {
        return nIndex;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_index_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_element_index_(iint* id)
{
    if(*id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_index_, minus value is invalid argument. id");
        return MW_ERROR;
    }
    uiint nID= *id;
    uiint nIndex= pMW->getElementIndex(nID);
    if(nIndex < IINT_MAX) {
        return nIndex;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_index_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
//--
// parent node
//--
iint mw_get_num_of_parent_node_(iint* id, iint* mglevel)
{
    if(*id < 0 || *mglevel < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_parent_node_, minus value is invalid argument. id | mglevel");
        return MW_ERROR;
    }

    uiint nID= *id;
    uiint nLevel= *mglevel;

    iint nNumOfParent= pMW->getNumOfParentNode( nID, nLevel);

    return nNumOfParent;
}
iint mw_get_parent_node_id_(iint* id, iint* mglevel, iint* index)
{
    if(*id < 0 || *mglevel < 0 || *index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_parent_node_id_, minus value is invalid argument. id | mglevel | index");
        return MW_ERROR;
    }

    uiint nID= *id;
    uiint nLevel= *mglevel;
    uiint nIndex= *index;

    iint nParentID= pMW->getParentNodeID(nID, nLevel, nIndex);

    return nParentID;
}

void mw_construct_node_connect_fem_(iint* node_id)
{
    if(*node_id < 0) {
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
iint mw_get_num_of_aggregate_element_(iint* node_id)
{
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_coord_, minus value is invalid argument. node_id");
        return MW_ERROR;
    }
    uiint nNodeID= *node_id;
    uiint nNumOfAggElem= pMW->getNumOfAggregateElement(nNodeID);
    if(nNumOfAggElem > IINT_MAX) {
        return IINT_MAX;
    } else {
        return nNumOfAggElem;
    }
}
iint mw_get_aggregate_element_id_(iint* node_id, iint* index)
{
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_coord_, minus value is invalid argument. node_id");
        return MW_ERROR;
    }
    uiint nNodeID= *node_id;
    uiint nIndex = *index;
    uiint nElemID= pMW->getAggregateElementID(nNodeID, nIndex);
    if(nElemID > IINT_MAX) {
        return IINT_MAX;
    } else {
        return nElemID;
    }
}
void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z)
{
    if(*node_id < 0) {
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
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_, minus value is invalid argument. node_id");
        return MW_ERROR;
    }
    uiint nNodeID = *node_id;
    uiint nNumOfDOF= pMW->GetNumOfDOF(nNodeID);
    if(nNumOfDOF < IINT_MAX) {
        return nNumOfDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_scalar_(iint* node_id)
{
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_scalar_, minus value is invalid argument. node_id");
        return MW_ERROR;
    }
    uiint nNodeID = *node_id;
    uiint nNumOfScalar= pMW->GetNumOfScalar(nNodeID);
    if(nNumOfScalar < IINT_MAX) {
        return nNumOfScalar;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_scalar_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_vector_(iint* node_id)
{
    if(*node_id < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_vector_, minus value is invalid argument. node_id");
        return MW_ERROR;
    }
    uiint nNodeID = *node_id;
    uiint nNumOfVector= pMW->GetNumOfVector(nNodeID);
    if(nNumOfVector < IINT_MAX) {
        return nNumOfVector;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_vector_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_()
{
    uiint nNumOfNode= pMW->getNodeSize();
    if(nNumOfNode < IINT_MAX) {
        return nNumOfNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_with_mesh_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_with_mesh_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint nNumOfNode= pMW->getNodeSize(iMesh);
    if(nNumOfNode < IINT_MAX) {
        return nNumOfNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_with_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_()
{
    uiint nNumOfElem= pMW->getElementSize();
    if(nNumOfElem < IINT_MAX) {
        return nNumOfElem;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_with_mesh_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_with_mesh_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint nNumOfElem= pMW->getElementSize(iMesh);
    if(nNumOfElem < IINT_MAX) {
        return nNumOfElem;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_with_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
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
iint mw_elemtype_hexa_()
{
    return pMW->elemtype_hexa();
}
iint mw_elemtype_hexa2_()
{
    return pMW->elemtype_hexa2();
}
iint mw_elemtype_tetra_()
{
    return pMW->elemtype_tetra();
}
iint mw_elemtype_tetra2_()
{
    return pMW->elemtype_tetra2();
}
iint mw_elemtype_prism_()
{
    return pMW->elemtype_prism();
}
iint mw_elemtype_prism2_()
{
    return pMW->elemtype_prism2();
}
iint mw_elemtype_quad_()
{
    return pMW->elemtype_quad();
}
iint mw_elemtype_quad2_()
{
    return pMW->elemtype_quad2();
}
iint mw_elemtype_triangle_()
{
    return pMW->elemtype_triangle();
}
iint mw_elemtype_triangle2_()
{
    return pMW->elemtype_triangle2();
}
iint mw_elemtype_line_()
{
    return pMW->elemtype_line();
}
iint mw_elemtype_line2_()
{
    return pMW->elemtype_line2();
}
iint mw_fistr_elemtype_hexa_()
{
    return pMW->fistr_elemtype_hexa();
}
iint mw_fistr_elemtype_hexa2_()
{
    return pMW->fistr_elemtype_hexa2();
}
iint mw_fistr_elemtype_tetra_()
{
    return pMW->fistr_elemtype_tetra();
}
iint mw_fistr_elemtype_tetra2_()
{
    return pMW->fistr_elemtype_tetra2();
}
iint mw_fistr_elemtype_prism_()
{
    return pMW->fistr_elemtype_prism();
}
iint mw_fistr_elemtype_prism2_()
{
    return pMW->fistr_elemtype_prism2();
}
iint mw_fistr_elemtype_quad_()
{
    return pMW->fistr_elemtype_quad();
}
iint mw_fistr_elemtype_quad2_()
{
    return pMW->fistr_elemtype_quad2();
}
iint mw_fistr_elemtype_triangle_()
{
    return pMW->fistr_elemtype_triangle();
}
iint mw_fistr_elemtype_triangle2_()
{
    return pMW->fistr_elemtype_triangle2();
}
iint mw_fistr_elemtype_line_()
{
    return pMW->fistr_elemtype_line();
}
iint mw_fistr_elemtype_line2_()
{
    return pMW->fistr_elemtype_line2();
}
iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype)
{
    if(*fistr_elemtype < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_fistr_elemtype_to_mw3_elemtype_, minus value is invalid argument. fistr_elemtype");
        return MW_ERROR;
    }
    uiint nElemType = *fistr_elemtype;
    return pMW->fistr_elemtype_to_mw3_elemtype(nElemType);
}
iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype)
{
    if(*mw3_elemtype < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_mw3_elemtype_to_fistr_elemtype_, minus value is invalid argument. mw3_elemtype");
        return MW_ERROR;
    }
    uiint nElemType = *mw3_elemtype;
    return pMW->mw3_elemtype_to_fistr_elemtype(nElemType);
}
iint mw_get_num_of_integ_point_(iint* shape_type)
{
    if(*shape_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_integ_point_, minus value is invalid argument. shape_type");
        return MW_ERROR;
    }
    uiint nShapeType = *shape_type;
    return pMW->NumOfIntegPoint(nShapeType);
}
void mw_shape_function_on_pt_(iint* shape_type, iint* igauss, double N[])
{
    if(*shape_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_on_pt_, minus value is invalid argument. shape_type");
        return;
    }
    if(*igauss < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa81_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa82_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa201_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa202_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_hexa203_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra101_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra104_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tetra1015_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism159_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_prism1518_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad84_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_quad89_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri31_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_tri63_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line21_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line32_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_shape_function_line32_, minus value is invalid argument. ishape");
        return;
    }
    uiint iGauss = *igauss;
    uiint nShape = *ishape;
    *N = pMW->ShapeFunc_Line32(iGauss, nShape);
}
void mw_dndr_(iint* shape_type, double dndr[])
{
    if(*shape_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_, minus value is invalid argument. shape_type");
        return;
    }
    uiint nShapeType = *shape_type;
    pMW->dNdr(nShapeType, dndr);
}
void mw_dndr_hexa81_(iint* igauss, iint* ishape, iint* iaxis, double* dndr)
{
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa81_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa81_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa82_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa82_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa201_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa201_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa202_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa202_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa203_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_hexa203_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra41_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra101_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra101_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra104_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra104_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra1015_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tetra1015_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism62_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism62_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism156_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism156_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism159_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism159_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism1518_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_prism1518_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad41_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad41_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad84_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad84_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad89_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_quad89_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri31_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri31_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri63_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_tri63_, minus value is invalid argument. ishape");
        return;
    }
    if(*iaxis  < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line21_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
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
    if(*igauss < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line32_, minus value is invalid argument. igauss");
        return;
    }
    if(*ishape < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndr_line32_, minus value is invalid argument. ishape");
        return;
    }
    uiint iGauss = *igauss;
    uiint iShape = *ishape;
    *dndr = pMW->dNdr_Line32_on_pt_on_shape(iGauss, iShape);
}
void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[])
{
    if(*elem_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndx_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_dndx_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*ielem  < 0) {
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
    if(*elem_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_det_jacobian_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_det_jacobian_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*igauss  < 0) {
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
    if(*elem_type < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. elem_type");
        return;
    }
    if(*num_of_integ < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. num_of_integ");
        return;
    }
    if(*igauss  < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_weight_, minus value is invalid argument. igauss");
        return;
    }
    uiint nElemType = *elem_type;
    uiint nNumOfInteg = *num_of_integ;
    uiint iGauss  = *igauss;
    pMW->Weight(nElemType, nNumOfInteg, iGauss, *w);
}
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
iint mw_get_num_of_boundary_bnode_mesh_()
{
    uiint nNumOfBNodeMesh= pMW->GetNumOfBoundaryNodeMesh();
    if(nNumOfBNodeMesh < IINT_MAX) {
        return nNumOfBNodeMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bface_mesh_()
{
    uiint nNumOfBFaceMesh= pMW->GetNumOfBoundaryFaceMesh();
    if(nNumOfBFaceMesh < IINT_MAX) {
        return nNumOfBFaceMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bedge_mesh_()
{
    uiint nNumOfBEdgeMesh = pMW->GetNumOfBoundaryEdgeMesh();
    if(nNumOfBEdgeMesh < IINT_MAX) {
        return nNumOfBEdgeMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_boundary_bvolume_mesh_()
{
    uiint nNumOfBVolMesh= pMW->GetNumOfBoundaryVolumeMesh();
    if(nNumOfBVolMesh < IINT_MAX) {
        return nNumOfBVolMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_boundary_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBNDType_BNodeMesh(iBMesh);
}
iint mw_get_bnd_type_bface_mesh_(iint*ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBNDType_BFaceMesh(iBMesh);
}
iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBNDType_BEdgeMesh(iBMesh);
}
iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnd_type_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBNDType_BVolumeMesh(iBMesh);
}
iint mw_get_neumann_type_()
{
    return pMW->getNeumannType();
}
iint mw_get_dirichlet_type_()
{
    return pMW->getDirichletType();
}
iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBNode= pMW->GetNumOfBNode_BNodeMesh(iBMesh);
    if(nNumOfBNode < IINT_MAX) {
        return nNumOfBNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBNode= pMW->GetNumOfBNode_BFaceMesh(iBMesh);
    if(nNumOfBNode < IINT_MAX) {
        return nNumOfBNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBNode= pMW->GetNumOfBNode_BEdgeMesh(iBMesh);
    if(nNumOfBNode < IINT_MAX) {
        return nNumOfBNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBNode= pMW->GetNumOfBNode_BVolumeMesh(iBMesh);
    if(nNumOfBNode < IINT_MAX) {
        return nNumOfBNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bnode_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nNumOfDOF= pMW->GetNumOfDOF_BNodeMesh(iBMesh, iBNode);
    if(nNumOfDOF < IINT_MAX) {
        return nNumOfDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfDOF= pMW->GetNumOfDOF_BFaceMesh(iBMesh);
    if(nNumOfDOF < IINT_MAX) {
        return nNumOfDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfDOF= pMW->GetNumOfDOF_BEdgeMesh(iBMesh);
    if(nNumOfDOF < IINT_MAX) {
        return nNumOfDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfDOF= pMW->GetNumOfDOF_BVolumeMesh(iBMesh);
    if(nNumOfDOF < IINT_MAX) {
        return nNumOfDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_dof_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, minus value is invalid argument. idof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint iDOF= *idof;
    uiint nDOF= pMW->GetDOF_BNodeMesh(iBMesh, iBNode, iDOF);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, minus value is invalid argument. idof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;
    uiint nDOF= pMW->GetDOF_BFaceMesh(iBMesh, iDOF);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, minus value is invalid argument. idof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;
    uiint nDOF= pMW->GetDOF_BEdgeMesh(iBMesh, iDOF);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*idof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, minus value is invalid argument. idof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iDOF= *idof;
    uiint nDOF= pMW->GetDOF_BVolumeMesh(iBMesh, iDOF);
    if(nDOF < IINT_MAX) {
        return nDOF;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_dof_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bnode_mesh_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;
    return pMW->GetBNodeValue_BNodeMesh(iBMesh, iBNode, nDOF);
}
double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bface_mesh_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;
    return pMW->GetBNodeValue_BFaceMesh(iBMesh, iBNode, nDOF, *mglevel);
}
double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bedge_mesh_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;
    return pMW->GetBNodeValue_BEdgeMesh(iBMesh, iBNode, nDOF, *mglevel);
}
double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_value_in_bvolume_mesh_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nDOF= *dof;
    return pMW->GetBNodeValue_BVolumeMesh(iBMesh, iBNode, nDOF, *mglevel);
}
iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BNode_BNodeMesh(iBMesh, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bnode_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BNode_BFaceMesh(iBMesh, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bface_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BNode_BEdgeMesh(iBMesh, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bedge_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BNode_BVolumeMesh(iBMesh, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_in_bvolume_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_bface_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bface_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBFace= pMW->GetNumOfBFace(iBMesh);
    if(nNumOfBFace < IINT_MAX) {
        return nNumOfBFace;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. ibface");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_value_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;
    uiint nDOF= *dof;
    return pMW->GetBFaceValue(iBMesh, iBFace, nDOF);
}
iint mw_get_num_of_bedge_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bedge_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBEdge= pMW->GetNumOfBEdge(iBMesh);
    if(nNumOfBEdge < IINT_MAX) {
        return nNumOfBEdge;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. ibedge");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_value_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBEdge= *ibedge;
    uiint nDOF= *dof;
    return pMW->GetBEdgeValue(iBMesh, iBEdge, nDOF);
}
iint mw_get_num_of_bvolume_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bvolume_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint nNumOfBVol= pMW->GetNumOfBVolume(iBMesh);
    if(nNumOfBVol < IINT_MAX) {
        return nNumOfBVol;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibvol < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. ibvol");
        return MW_ERROR;
    }
    if(*dof < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_value_, minus value is invalid argument. dof");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;
    uiint nDOF= *dof;
    return pMW->GetBVolumeValue(iBMesh, iBVol, nDOF);
}
iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bface_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bface_, minus value is invalid argument. ibface");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;
    return pMW->GetNumOfNode_BFace(iBMesh, iBFace);
}
iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibface");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBFace= *ibface;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BFace(iBMesh, iBFace, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bedge_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bedge_, minus value is invalid argument. ibedge");
        return MW_ERROR;
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
        return MW_ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, minus value is invalid argument. ibedge");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBEdge= *ibedge;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BEdge(iBMesh, iBEdge, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bvolume_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibvol < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_node_bvolume_, minus value is invalid argument. ibvol");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;
    return pMW->GetNumOfNode_BVolume(iBMesh, iBVol);
}
iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibvol < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibvol");
        return MW_ERROR;
    }
    if(*ibnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, minus value is invalid argument. ibnode");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    uiint iBVol= *ibvol;
    uiint iBNode= *ibnode;
    uiint nNodeID= pMW->GetNodeID_BVolume(iBMesh, iBVol, iBNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_bnode_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_namelength_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBNodeMesh_NameLength(iBMesh);
}
void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bnode_mesh_name_, minus value is invalid argument. name_len");
        return;
    }
    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;
    string sName = pMW->GetBNodeMesh_Name(iBMesh);
    uiint i, nNumOfChar = sName.length();
    if(nNumOfChar > nNameLen) {
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++) {
            name[i] = sName[i];
        };
        name[nLength]='\0';
    } else {
        for(i=0; i < nNumOfChar; i++) {
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++) {
            name[i] = '\0';
        };
    }
}
iint mw_get_bface_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_namelength_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBFaceMesh_NameLength(iBMesh);
}
void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bface_mesh_name_, minus value is invalid argument. name_len");
        return;
    }
    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;
    string sName = pMW->GetBFaceMesh_Name(iBMesh);
    uiint i, nNumOfChar = sName.length();
    if(nNumOfChar > nNameLen) {
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++) {
            name[i] = sName[i];
        };
        name[nLength]='\0';
    } else {
        for(i=0; i < nNumOfChar; i++) {
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++) {
            name[i] = '\0';
        };
    }
}
iint mw_get_bvolume_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_namelength_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBVolumeMesh_NameLength(iBMesh);
}
void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bvolume_mesh_name_, minus value is invalid argument. name_len");
        return;
    }
    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;
    string sName = pMW->GetBVolumeMesh_Name(iBMesh);
    uiint i, nNumOfChar = sName.length();
    if(nNumOfChar > nNameLen) {
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++) {
            name[i] = sName[i];
        };
        name[nLength]='\0';
    } else {
        for(i=0; i < nNumOfChar; i++) {
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++) {
            name[i] = '\0';
        };
    }
}
iint mw_get_bedge_mesh_namelength_(iint* ibmesh)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_namelength_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    uiint iBMesh= *ibmesh;
    return pMW->GetBEdgeMesh_NameLength(iBMesh);
}
void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_name_, minus value is invalid argument. ibmesh");
        return;
    }
    if(*name_len < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_bedge_mesh_name_, minus value is invalid argument. name_len");
        return;
    }
    uiint iBMesh= *ibmesh;
    uiint nNameLen= *name_len;
    string sName = pMW->GetBEdgeMesh_Name(iBMesh);
    uiint i, nNumOfChar = sName.length();
    if(nNumOfChar > nNameLen) {
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++) {
            name[i] = sName[i];
        };
        name[nLength]='\0';
    } else {
        for(i=0; i < nNumOfChar; i++) {
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++) {
            name[i] = '\0';
        };
    }
}
iint mw_get_edge_id_bedge_(iint* ibmesh, iint* ibedge)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, minus value is invalid argument. ibedge");
        return MW_ERROR;
    }
    uiint iBMesh = *ibmesh;
    uiint iBEdge = *ibedge;
    uiint nEdgeID= pMW->GetEdgeID_BEdge(iBMesh, iBEdge);
    if(nEdgeID < IINT_MAX) {
        return nEdgeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_edge_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bedge_(iint*ibmesh, iint* ibedge)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibedge < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, minus value is invalid argument. ibedge");
        return MW_ERROR;
    }
    uiint iBMesh = *ibmesh;
    uiint iBEdge = *ibedge;
    uiint nElemID= pMW->GetElemID_BEdge(iBMesh, iBEdge);
    if(nElemID < IINT_MAX) {
        return nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bedge_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, minus value is invalid argument. ibface");
        return MW_ERROR;
    }
    uiint iBMesh = *ibmesh;
    uiint iBFace = *ibface;
    uiint nFaceID= pMW->GetFaceID_BFace(iBMesh, iBFace);
    if(nFaceID < IINT_MAX) {
        return nFaceID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_face_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibface < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, minus value is invalid argument. ibface");
        return MW_ERROR;
    }
    uiint iBMesh = *ibmesh;
    uiint iBFace = *ibface;
    uiint nElemID= pMW->GetElemID_BFace(iBMesh, iBFace);
    if(nElemID < IINT_MAX) {
        return nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bface_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol)
{
    if(*ibmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, minus value is invalid argument. ibmesh");
        return MW_ERROR;
    }
    if(*ibvol < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, minus value is invalid argument. ibvol");
        return MW_ERROR;
    }
    uiint iBMesh = *ibmesh;
    uiint iBVol = *ibvol;
    uiint nElemID= pMW->GetElemID_BVolume(iBMesh, iBVol);
    if(nElemID < IINT_MAX) {
        return nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elem_id_bvolume_, return value is over INT_MAX");
        return IINT_MAX;
    }
}


MPI_Datatype mw_mpi_int_()
{
    return MPI_INT;
}
MPI_Datatype mw_mpi_iint_()
{
    return pMW->MPI_IINT();
}
MPI_Datatype mw_mpi_uiint_()
{
    return pMW->MPI_UIINT();
}
MPI_Datatype mw_mpi_double_()
{
    return MPI_DOUBLE;
}
MPI_Comm mw_mpi_comm_()
{
    return MPI_COMM_WORLD;
}
MPI_Op mw_mpi_sum_()
{
    return MPI_SUM;
}
MPI_Op mw_mpi_max_()
{
    return MPI_MAX;
}
MPI_Op mw_mpi_min_()
{
    return MPI_MIN;
}

int mw_get_rank_()
{
    return pMW->GetRank();
}
int mw_get_num_of_process_()
{
    return pMW->GetNumOfProcess();
}
void mw_allreduce_r_(double val[], iint* val_size, MPI_Op* op)
{
    if(*val_size < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allreduce_r_, minus value is invalid argument. val_size");
        return;
    }
    double* rval = new double[*val_size];

    pMW->AllReduce(val, rval, *val_size, MPI_DOUBLE, *op, MPI_COMM_WORLD);

    for(iint i=0; i < *val_size; i++) val[i]= rval[i];

    delete [] rval;
}
void mw_allreduce_i_(iint val[], iint* val_size, MPI_Op* op)
{
    if(*val_size < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allreduce_i_, minus value is invalid argument. val_size");
        return;
    }
    iint* rval = new iint[*val_size];
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    pMW->AllReduce(val, rval, *val_size, MPI_IINT, *op, MPI_COMM_WORLD);

    for(iint i=0; i < *val_size; i++) val[i]= rval[i];

    delete [] rval;
}
iint mw_barrier_()
{
    return pMW->Barrier(MPI_COMM_WORLD);
}
iint mw_abort_(iint* error)
{
    return pMW->Abort(MPI_COMM_WORLD, *error);
}
iint mw_allgather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }

    return pMW->AllGather((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, MPI_COMM_WORLD);
}
iint mw_allgather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_i_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_allgather_i_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->AllGather( (void*)sendbuf, *sendcnt, MPI_IINT, (void*)recvbuf, *recvcnt, MPI_IINT, MPI_COMM_WORLD );
}
iint mw_gather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }

    return pMW->Gather( (void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
iint mw_gather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_i_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_gather_i_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->Gather( (void*)sendbuf, *sendcnt, MPI_IINT, (void*)recvbuf, *recvcnt, MPI_IINT, *root, MPI_COMM_WORLD);
}
iint mw_scatter_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_r_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }

    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_DOUBLE, (void*)recvbuf, *recvcnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
iint mw_scatter_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root)
{
    if(*sendcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_i_, minus value is invalid argument. sendcnt");
        return MW_ERROR;
    }
    if(*recvcnt < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_scatter_i_, minus value is invalid argument. recvcnt");
        return MW_ERROR;
    }
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->Scatter((void*)sendbuf, *sendcnt, MPI_IINT, (void*)recvbuf, *recvcnt, MPI_IINT, *root, MPI_COMM_WORLD);
}
iint mw_bcast_i_(iint buf[], iint* cnt, int* root)
{
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->Bcast(buf, *cnt, MPI_IINT, *root, MPI_COMM_WORLD);
}
iint mw_bcast_r_(double buf[], iint* cnt, int* root)
{
    return pMW->Bcast(buf, *cnt, MPI_DOUBLE, *root, MPI_COMM_WORLD);
}
iint mw_bcast_s_(char buf[], iint* cnt, int* root)
{
    return pMW->Bcast(buf, *cnt, MPI_CHAR, *root, MPI_COMM_WORLD);
}
iint mw_get_num_of_neibpe_(iint* imesh)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_neibpe_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    uiint iMesh = *imesh;
    uiint nNumOfNeibPE= pMW->GetNumOfNeibPE(iMesh);

    if(nNumOfNeibPE < IINT_MAX) {
        return nNumOfNeibPE;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_neibpe_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_transrank_(iint* imesh, iint* ipe)
{
    if(*imesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, minus value is invalid argument. imesh");
        return MW_ERROR;
    }
    if(*ipe < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, minus value is invalid argument. ipe");
        return MW_ERROR;
    }
    iint iMesh = *imesh;
    iint iPE = *ipe;
    uiint nTransRank= pMW->GetTransRank(iMesh, iPE);

    if(nTransRank < IINT_MAX) {
        return nTransRank;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_transrank_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

iint mw_send_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank)
{
    iint nNumOfNode = *num_of_node;
    iint nNumOfDOF = *dof_size;
    int tag=0;

    iint nCount = nNumOfNode*nNumOfDOF;

    return pMW->Send((void*)buf, nCount, MPI_DOUBLE, *trans_rank, tag, MPI_COMM_WORLD);
}
iint mw_send_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank)
{
    iint nNumOfNode = *num_of_node;
    iint nNumOfDOF = *dof_size;
    int tag=0;

    iint nCount = nNumOfNode*nNumOfDOF;
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->Send((void*)buf, nCount, MPI_IINT, *trans_rank, tag, MPI_COMM_WORLD);
}
iint mw_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank)
{
    iint nNumOfNode = *num_of_node;
    iint nNumOfDOF = *dof_size;
    int tag=0;
    MPI_Status sta;

    iint nCount = nNumOfNode*nNumOfDOF;

    return pMW->Recv((void*)buf, nCount, MPI_DOUBLE, *trans_rank, tag, MPI_COMM_WORLD, &sta);
}
iint mw_recv_i_(iint buf[], iint* num_of_node, int* dof_size, int* trans_rank)
{
    iint nNumOfNode = *num_of_node;
    int nNumOfDOF = *dof_size;
    int tag=0;
    MPI_Status sta;

    iint nCount = nNumOfNode*nNumOfDOF;
    MPI_Datatype MPI_IINT= pMW->MPI_IINT();

    return pMW->Recv((void*)buf, nCount, MPI_IINT, *trans_rank, tag, MPI_COMM_WORLD, &sta);
}

void mw_send_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank)
{
    if(*num_of_node < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. num_of_node");
        return;
    }
    if(*dof_size < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. dof_size");
        return;
    }
    if(*trans_rank < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. trans_rank");
        return;
    }
    iint nNumOfNode = *num_of_node;
    iint nNumOfDOF = *dof_size;
    int transRank = *trans_rank;

    pMW->Send_Recv_R(buf, nNumOfNode, nNumOfDOF, transRank);
    //// pMW->Sendrecv_R(buf, nNumOfNode, nNumOfDOF, transRank);
}
void mw_send_recv_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank)
{
    if(*num_of_node < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. num_of_node");
        return;
    }
    if(*dof_size < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. dof_size");
        return;
    }
    if(*trans_rank < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_send_recv_, minus value is invalid argument. trans_rank");
        return;
    }
    iint nNumOfNode = *num_of_node;
    iint nNumOfDOF = *dof_size;
    int transRank = *trans_rank;

    pMW->Send_Recv_I(buf, nNumOfNode, nNumOfDOF, transRank);
    //// pMW->Sendrecv_I(buf, nNumOfNode, nNumOfDOF, transRank);
}
//--
// MPI End
//--


//--
// 通信テーブル:CommMesh2
//--
iint mw_get_num_of_comm_mesh_()
{
    uiint nNumOfCommMesh= pMW->GetNumOfCommMesh();
    if(nNumOfCommMesh < IINT_MAX) {
        return nNumOfCommMesh;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_mesh_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_comm_node_(iint* icmesh)
{
    if(*icmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_node_, minus value is invalid argument. icmesh");
        return MW_ERROR;
    }
    uiint iComMesh= *icmesh;
    uiint nNumOfCommNode= pMW->GetNumOfCommNode(iComMesh);
    if(nNumOfCommNode < IINT_MAX) {
        return nNumOfCommNode;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_comm_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode)
{
    if(*icmesh < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, minus value is invalid argument. icmesh");
        return MW_ERROR;
    }
    if(*icnode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, minus value is invalid argument. icnode");
        return MW_ERROR;
    }
    uiint iComMesh= *icmesh;
    uiint iComNode= *icnode;
    uiint nNodeID= pMW->GetNodeID_CommNode(iComMesh, iComNode);
    if(nNodeID < IINT_MAX) {
        return nNodeID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_node_id_comm_node_, return value is over INT_MAX");
        return IINT_MAX;
    }
}

//--
// 接合面:ContactMesh
//--
iint mw_get_num_of_contact_()
{
    uiint nNumOfContact= pMW->GetNumOfContactMesh();

    if(nNumOfContact < IINT_MAX) {
        return nNumOfContact;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_contact_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_contact_id_(iint* icont)
{
    uiint nID= pMW->GetContactMeshID(*icont);

    if(nID < IINT_MAX) {
        return nID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_contact_id_, return value is over INT_MAX");
        return IINT_MAX;
    }
}






iint mw_get_num_of_elementgroup_()
{
    uiint nNumOfElemGrp= pMW->GetNumOfElementGroup();
    if(nNumOfElemGrp < IINT_MAX) {
        return nNumOfElemGrp;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_elementgroup_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_num_of_element_id_(iint* igrp)
{
    if(*igrp < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_id_, minus value is invalid argument. igrp");
        return MW_ERROR;
    }
    uiint iGrp = *igrp;
    uiint nNumOfElem= pMW->GetNumOfElementID(iGrp);
    if(nNumOfElem < IINT_MAX) {
        return nNumOfElem;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_num_of_element_id_, INT_MAX over value");
        return IINT_MAX;
    }
}
iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index)
{
    if(*igrp < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, minus value is invalid argument. igrp");
        return MW_ERROR;
    }
    if(*index < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, minus value is invalid argument. index");
        return MW_ERROR;
    }
    uiint iGrp = *igrp;
    uiint i = *index;
    uiint nElemID= pMW->GetElementID_with_ElementGroup(iGrp, i);
    if(nElemID < IINT_MAX) {
        return nElemID;
    } else {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_element_id_with_elementgroup_, return value is over INT_MAX");
        return IINT_MAX;
    }
}
iint mw_get_elementgroup_name_length_(iint* igrp)
{
    if(*igrp < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_length_, minus value is invalid argument. igrp");
        return MW_ERROR;
    }
    uiint iGrp = *igrp;
    iint nLength = pMW->GetElementGroupName_Length(iGrp);
    return nLength;
}
void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len)
{
    if(*igrp < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_, minus value is invalid argument. igrp");
        return;
    }
    if(*name_len < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_get_elementgroup_name_, minus value is invalid argument. name_len");
        return;
    }
    uiint iGrp = *igrp;
    uiint nNameLen= *name_len;
    string sName = pMW->GetElementGroupName(iGrp);
    uiint i, nNumOfChar = sName.length();
    if(nNumOfChar > nNameLen) {
        uiint nLength = nNameLen-1;
        for(i=0; i < nLength; i++) {
            name[i] = sName[i];
        };
        name[nLength]='\0';
    } else {
        for(i=0; i < nNumOfChar; i++) {
            name[i] = sName[i];
        };
        uiint nLength = nNameLen;
        for(i=nNumOfChar; i < nLength; i++) {
            name[i] = '\0';
        };
    }
}
void mw_logger_set_mode_(iint* mode)
{
    if(*mode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_mode_, minus value is invalid argument. mode");
        return;
    }
    uiint nMode= *mode;
    pMW->LoggerMode(nMode);
}
void mw_logger_set_device_(iint* mode, iint* device)
{
    if(*mode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_device_, minus value is invalid argument. mode");
        return;
    }
    if(*device < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_set_device_, minus value is invalid argument. device");
        return;
    }
    uiint nMode= *mode;
    uiint nDevice= *device;
    pMW->LoggerDevice(nMode, nDevice);
}
void mw_logger_info_mssg_ (iint* mode, const char* message)
{
    if(*mode < 0) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_info_mssg_, minus value is invalid argument. mode");
        return;
    }
    uiint nMode= *mode;

    pMW->LoggerInfoMssg(nMode, message);

}
//--
// format= %d:iint, %f:double(fixed), %e:double(scientific), %s:const char*
//--
void mw_logger_info_(iint* mode, const char* format, ...)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();

    if(*mode < 0) {
        pLogger->Info(Utility::LoggerMode::Error, "mw_logger_info_, minus value is invalid argument. mode");
        return;
    }
    uiint nMode= *mode;


    vector<void*> param;
    vint    vnVal;
    vuint   vuVal;
    vdouble vdVal, veVal;
    vstring vsVal;

    if(nMode==Utility::LoggerMode::Debug || nMode==Utility::LoggerMode::Error ||
       nMode==Utility::LoggerMode::Warn  || nMode==Utility::LoggerMode::Info) {

        uiint nLength = strlen(format);

        va_list list;
        va_start( list, format);
        for(uiint i=0; i < nLength; i++) {
            if(format[i] == '%') {
                ++i;
                switch( format[i] ) {
                case('d'): {
                    iint* nVal = va_arg(list, iint*);     //Fortran引数はポインター
                    vnVal.push_back(*nVal);
                }
                break;
                case('u'): {
                    uiint* nVal = va_arg(list, uiint*);     //Fortran引数はポインター
                    vuVal.push_back(*nVal);
                }
                break;
                case('f'): {
                    double* dVal= va_arg( list, double* );    //Fortran引数はポインター
                    vdVal.push_back(*dVal);
                }
                break;
                case('e'): {
                    double* dVal= va_arg( list, double* );    //Fortran引数はポインター
                    veVal.push_back(*dVal);
                }
                break;
                case('s'): {
                    string sVal= va_arg( list, const char* );
                    vsVal.push_back(sVal);
                }
                break;
                default:
                    break;
                }
            }
        };
        va_end( list );

        uiint icase_d(0), icase_u(0), icase_f(0), icase_e(0), icase_s(0);
        for(uiint i=0; i < nLength; i++) {
            if(format[i] == '%') {
                ++i;
                switch( format[i] ) {
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

        pLogger->Info(nMode, format, param);//--------------- Logger出力

    } else {
        pLogger->Info(Utility::LoggerMode::Error, "invalid logger mode, CMW::LoggerInfo");
    }
}
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

