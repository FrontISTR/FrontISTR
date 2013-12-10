/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : FrontSTR Input File Converter                     *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

/*
	CNFDB_507 Ver.1.0
*/

// 507 Loads

#include "CNFData.h"
#include "CNFDB_507.h"

using namespace std;

void CNFDB_507::cstructural_load_rec::copy_from (
	const CNFDB_507::cstructural_load_rec& r)
{
	int i;

	loadID     = r.loadID;
        loadtype   = r.loadtype;
	color      = r.color;
	layer      = r.layer;
	define_sys = r.define_sys;
	subtype    = r.subtype;
	is_expanded= r.is_expanded;
	for( i=0; i<3; i++) dof_face[i]  = r.dof_face[i];
	for( i=0; i<5; i++) value[i]     = r.value[i];
	for( i=0; i<5; i++) functions[i] = r.functions[i];
	Enclosure = r.Enclosure;
	can_shade = r.can_shade;
	can_be_shaded = r.can_be_shaded;
	for( i=0; i<2; i++) add1_id[i]   = r.add1_id[i];
	for( i=0; i<3; i++) dir_func[i]  = r.dir_func[i];
	for( i=0; i<3; i++) direction[i] = r.direction[i];
}

CNFDB_507::cstructural_load_rec::cstructural_load_rec(
	const CNFDB_507::cstructural_load_rec& r)
{
	copy_from( r );
}


CNFDB_507::cstructural_load_rec& CNFDB_507::cstructural_load_rec::operator = (
	const CNFDB_507::cstructural_load_rec& r)
{
	copy_from(r);
	return *this;
}


//=============================================================================

void CNFDB_507::cgeometric_load_rec::copy_from( const CNFDB_507::cgeometric_load_rec& r)
{
	int i, j;
	// ##1
	loadID     = r.loadID;
	loadtype   = r.loadtype;
	color      = r.color;
	layer      = r.layer;
	define_sys = r.define_sys;
	subtype    = r.subtype;
	is_expanded= r.is_expanded;
	// ##2
	for( i=0; i<3; i++) dof_face[i]  = r.dof_face[i];
	// ##3
	for( i=0; i<5; i++) value[i]  = r.value[i];
	// ##4
	for( i=0; i<5; i++) functions[i]  = r.functions[i];
	// ##5
	Enclosure = r.Enclosure;
	can_shade = r.can_shade;
	can_be_shaded = r.can_be_shaded;
	for( i=0; i<2; i++) addl_id[i]  = r.addl_id[i];
	// ##6
	for( i=0; i<3; i++) dir_func[i] = r.dir_func[i];
	// ##7
	for( i=0; i<3; i++) direction[i] = r.direction[i];
	// ##8
	dir_mode = r.dir_mode;
	dir_id   = r.dir_id;
	// ##9
	for( i=0; i<3; i++) dir_base[i] = r.dir_base[i];
	// ##10
	for( i=0; i<3; i++) dir_vector[i] = r.dir_vector[i];
	// ##11
	var_mode = r.var_mode;
	var_funcID = r.var_funcID;
	// ##12
	strcpy( var_name, r.var_name );
	// ##13
	strcpy( var_equation, r.var_equation );
	// ##14-17
       for( i=0; i<4; i++) {
	      for( j=0;j<3; j++) var_locate[i][j] = r.var_locate[i][j];
	}
	// ##18
	for( i=0; i<4; i++) var_value[i] = r.var_value[i];
	// ##19
	adjust_midside = r.adjust_midside;
	is_expanded2 = r.is_expanded2;
}


CNFDB_507::cgeometric_load_rec::cgeometric_load_rec( const CNFDB_507::cgeometric_load_rec& r)
{
	copy_from( r );
}

CNFDB_507::cgeometric_load_rec&
CNFDB_507::cgeometric_load_rec::operator=( const CNFDB_507::cgeometric_load_rec& r)
{
	copy_from( r );
	return *this;
}


//=============================================================================



CNFDB_507::CNFDB_507()
 : CNFDataBlock(507)
{}



void CNFDB_507::read_structural_load( CNFData* nfd, std::vector<cstructural_load_rec>& list )
{
	char buff[256];
	cstructural_load_rec rec;

	list.clear();
	while(1) {
		// ##1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIB", 
			&rec.loadID,
			&rec.loadtype,
			&rec.color,
			&rec.layer,
			&rec.define_sys,
			&rec.subtype,
			&rec.is_expanded );
		if( rec.loadID == -1 ) break;
		// ##2
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III", 
			&rec.dof_face[0], &rec.dof_face[1], &rec.dof_face[2] );
		// ##3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFFF", 
			&rec.value[0], &rec.value[1], &rec.value[2], &rec.value[3], &rec.value[4] );
		// ##4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIII", 
			&rec.functions[0], &rec.functions[1], &rec.functions[2], &rec.functions[3], &rec.functions[4] );
		// ##5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "BBBII",
			&rec.Enclosure,
			&rec.can_shade,
			&rec.can_be_shaded,
			&rec.add1_id[0],
			&rec.add1_id[1]);
		// ##6
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III", 
			&rec.dir_func[0], &rec.dir_func[1], &rec.dir_func[2] );
		// ##7
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", 
			&rec.direction[0], &rec.direction[1], &rec.direction[2] );
		// ----------------------------------------------------------
		list.push_back( rec );
	}
}


void CNFDB_507::read_geometric_load( CNFData* nfd, std::vector<cgeometric_load_rec>& list )
{
	char buff[256];
	int i;
	cgeometric_load_rec rec;
	list.clear();
	while(1) {
		// ##1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIB", 
			&rec.loadID,
			&rec.loadtype,
			&rec.color,
			&rec.layer,
			&rec.define_sys,
			&rec.subtype,
			&rec.is_expanded );
		if( rec.loadID == -1 ) break;
		// ##2
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III", 
			&rec.dof_face[0], &rec.dof_face[1], &rec.dof_face[2] );
		// ##3
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFFF", 
			&rec.value[0], &rec.value[1], &rec.value[2], &rec.value[3], &rec.value[4] );
		// ##4
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIII", 
			&rec.functions[0], &rec.functions[1], &rec.functions[2], &rec.functions[3], &rec.functions[4] );
		// ##5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "BBBII",
			&rec.Enclosure,
			&rec.can_shade,
			&rec.can_be_shaded,
			&rec.addl_id[0],
			&rec.addl_id[1] );
		// ##6
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III", 
			&rec.dir_func[0], &rec.dir_func[1], &rec.dir_func[2] );
		// ##7
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", 
			&rec.direction[0], &rec.direction[1], &rec.direction[2] );
		// ##8
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "II", 
			&rec.dir_mode,
			&rec.dir_id );
		// ##9
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", 
			&rec.dir_base[0], &rec.dir_base[1], &rec.dir_base[2] );
		// ##10
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", 
			&rec.dir_vector[0], &rec.dir_vector[1], &rec.dir_vector[2] );
		// ##11
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "II", 
			&rec.var_mode,
			&rec.var_funcID );
		// ##12
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, rec.var_name, sizeof(rec.var_name) );
		// ##13
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff,  rec.var_equation, sizeof(rec.var_equation) );
		// ##14-17
		for(i=0; i<4; i++){
			nfd->ReadLineEx( buff );
			nfd->ReadRecord( buff, "FFF", 
				&rec.var_locate[i][0], &rec.var_locate[i][1], &rec.var_locate[i][2] );
		}
		// ##18
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFF", 
			&rec.var_value[0], &rec.var_value[1], &rec.var_value[2], &rec.var_value[3] );
		// ##19
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "BB", 
			&rec.adjust_midside,
			&rec.is_expanded2 );
		//-------------------------------------------------------------------
		list.push_back( rec );
	}
}


void CNFDB_507::read_temp_load( CNFData* nfd, std::vector<ctemp_load_rec>& list )
{
	char buff[256];
	ctemp_load_rec rec;

	list.clear();
	while(1) {
		// ##1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIFFIB", 
			&rec.ID,
			&rec.color,
			&rec.layer,
			&rec.temp,
			&rec.temp_co,
			&rec.funcID,
			&rec.is_extended );
			if( rec.ID == -1 ) break;
		//-------------------------------------
		list.push_back( rec );
	}
}




void CNFDB_507::Read( CNFData* nfd )
{
	char buff[256];

	// #1
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "I", &setID);
	// #2
		nfd->ReadLineEx( buff );
		nfd->ReadStr( buff, title, sizeof(title));

	// #3
		nfd->ReadLineEx( buff );
		if( nfd->version >= 8.0 ){
			nfd->ReadRecord( buff, "IFBBBBF",
				&CSys,
				&Def_temp,
				&temp_on,
				&grav_on,
				&omega_on,
				&Ref_temp_on,
				&Ref_temp );
		} else {
			nfd->ReadRecord( buff, "IFBBB",
				&CSys,
				&Def_temp,
				&temp_on,
				&grav_on,
				&omega_on );
		}
	// #4,5
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &grav[0], &grav[1], &grav[2]);
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &grav[3], &grav[4], &grav[5]);
	// #6
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &origin[0], &origin[1], &origin[2]);
	// #7
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &omega[0], &omega[1], &omega[2]);
	// #8
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFI",
			&stef_boltz,
			&abs_temp,
			&free_cnv_exp,
			&rad_space_element );
	// #9
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFF",
			&fc_flu_cond,
			&fc_flu_cp,
			&fc_flu_vis,
			&fc_flu_dens );
	// #10
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFF",
			&fc_cons_coeff,
			&fc_reynolds,
			&fc_pran_in,
			&fc_pran_out);
	// #11
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "III",
			&tfc_flu_cond,
			&tfc_flu_cp,
			&tfc_flu_vis);
	// #12
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "BBB",
			&alt_free_conv,
			&fc_flu_flag,
			&fc_conv_flow );
	// #13
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFF",
			&nl_arc_scale,
			&nl_arcmaxadj,
			&nl_arcminadj,
			&nl_bounds_rb );
	// #14
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFF", &nl_conv[0], &nl_conv[1], &nl_conv[2]);
	// #15
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "FFFFFF",
			&nl_fstress,
			&nl_lseach_tol,
			&nl_mxadj_init,
			&nl_max_rot,
			&nl_stab_tol,
			&nl_time_inc );
	// #16
		nfd->ReadLineEx( buff );
		if( nfd->version >= 6.0 ){
			nfd->ReadRecord( buff, "FFFFFFFFF",
				&dyn_damp_ov,
				&dyn_dampW3,
				&dyn_dampW4,
				&dyn_keep_freq[0],
				&dyn_keep_freq[1],
				&dyn_trans_dt,
				&dyn_min_freq,
				&dyn_max_freq,
				&dyn_cluster_freq );
		} else {
			nfd->ReadRecord( buff, "FFFFFF",
				&dyn_damp_ov,
				&dyn_dampW3,
				&dyn_dampW4,
				&dyn_keep_freq[0],
				&dyn_keep_freq[1],
				&dyn_trans_dt );
		}
	// #17
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIIII",
			&nl_arc_const,
			&nl_arc_iter,
			&nl_arc_maxst,
			&nl_div_limit,
			&nl_dom_pdstp,
			&nl_increment,
			&nl_inter_out,
			&nl_kstep,
			&nl_mx_bisect );
	// #18
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIIII",
			&nl_max_iter,
			&nl_max_lsrch,
			&nl_out_iter,
			&nl_quasi_newt,
			&nl_sol_strat,
			&nl_stiff_meth,
			&nl_skip_adj,
			&nl_sol_over );
	// #19
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "IIIIIII",
			&dyn_freq_tab,
			&dyn_damptab,
			&dyn_keep_md,
			&dyn_tran_ts,
			&dyn_out_int,
			&dyn_rand_psd,
			&dyn_on_freq );
	// #20
		nfd->ReadLineEx( buff );
		nfd->ReadRecord( buff, "BBBBBBB",
			&nl_on,
			&nl_conv_flag[0],
			&nl_conv_flag[1],
			&nl_conv_flag[2],
			&nl_mnewt_ls,
			&nl_mnewt_qn,
			&nl_mnewt_bs );
	// #21
		nfd->ReadLineEx( buff );
		if( nfd->version >= 7.0 ){
			nfd->ReadRecord( buff, "BBBBBBIII",
				&dyn_on,
				&dyn_type,
				&dyn_damp_method,
				&dyn_massfrm,
				&dyn_datarec,
				&dyn_log_inter,
				&dyn_freq_type,
				&dyn_psd_type,
				&dyn_psd_interpol );
		} else if( nfd->version >= 6.0 ){
			nfd->ReadRecord( buff, "BBBBBBII",
				&dyn_on,
				&dyn_type,
				&dyn_damp_method,
				&dyn_massfrm,
				&dyn_datarec,
				&dyn_log_inter,
				&dyn_freq_type,
				&dyn_psd_type );
		} else {
			nfd->ReadRecord( buff, "BBBBB",
				&dyn_on,
				&dyn_type,
				&dyn_damp_method,
				&dyn_massfrm,
				&dyn_datarec);
		}
	//------------------------------------
		read_structural_load( nfd, structural_load_list );
	//------------------------------------
		read_geometric_load( nfd, geometric_load_list );
	//------------------------------------
		read_temp_load( nfd, ndtemp_load_list );
	//------------------------------------
		read_temp_load( nfd, eltemp_load_list );
}



//*****************************************************************************


void CNFDB_507::write_structural_load( class CNFData* nfd, FILE* fp, std::vector<cstructural_load_rec>& list )
{
	vector<cstructural_load_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter ++ ){
		// ##1
		nfd->WriteData( fp, "IIIIIIBn", 
			iter->loadID,
			iter->loadtype,
			iter->color,
			iter->layer,
			iter->define_sys,
			iter->subtype,
			iter->is_expanded );
		// ##2
		nfd->WriteData( fp, "IIIn", 
			iter->dof_face[0], iter->dof_face[1], iter->dof_face[2] );
		// ##3
		nfd->WriteData( fp, "FFFFFn", 
			iter->value[0], iter->value[1], iter->value[2], iter->value[3], iter->value[4] );
		// ##4
		nfd->WriteData( fp, "IIIIIn", 
			iter->functions[0], iter->functions[1], iter->functions[2], iter->functions[3], iter->functions[4] );
		// ##5
		nfd->WriteData( fp, "BBBIIn",
			iter->Enclosure,
			iter->can_shade,
			iter->can_be_shaded,
			iter->add1_id[0],
			iter->add1_id[1]);
		// ##6
		nfd->WriteData( fp, "IIIn", 
			iter->dir_func[0], iter->dir_func[1], iter->dir_func[2] );
		// ##7
		nfd->WriteData( fp, "FFFn", 
			iter->direction[0], iter->direction[1], iter->direction[2] );
	}
	nfd->WriteData( fp, "IIIIIIBn", -1,-1,-1,-1,-1,-1,0);
}



void CNFDB_507::write_geometric_load( class CNFData* nfd, FILE* fp, std::vector<cgeometric_load_rec>& list )
{
	vector<cgeometric_load_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter ++ ){
		// ##1
		nfd->WriteData( fp, "IIIIIIBn", 
			iter->loadID,
			iter->loadtype,
			iter->color,
			iter->layer,
			iter->define_sys,
			iter->subtype,
			iter->is_expanded );
		// ##2
		nfd->WriteData( fp, "IIIn", 
			iter->dof_face[0], iter->dof_face[1], iter->dof_face[2] );
		// ##3
		nfd->WriteData( fp, "FFFFFn", 
			iter->value[0], iter->value[1], iter->value[2], iter->value[3], iter->value[4] );
		// ##4
		nfd->WriteData( fp, "IIIIIn", 
			iter->functions[0], iter->functions[1], iter->functions[2], iter->functions[3], iter->functions[4] );
		// ##5
		nfd->WriteData( fp, "BBBIIn",
			iter->Enclosure,
			iter->can_shade,
			iter->can_be_shaded,
			iter->addl_id[0],
			iter->addl_id[1] );
		// ##6
		nfd->WriteData( fp, "IIIn", 
			iter->dir_func[0], iter->dir_func[1], iter->dir_func[2] );
		// ##7
		nfd->WriteData( fp, "FFFn", 
			iter->direction[0], iter->direction[1], iter->direction[2] );
		// ##8
		nfd->WriteData( fp, "IIn", 
			iter->dir_mode,
			iter->dir_id );
		// ##9
		nfd->WriteData( fp, "FFFn", 
			iter->dir_base[0], iter->dir_base[1], iter->dir_base[2] );
		// ##10
		nfd->WriteData( fp, "FFFn", 
			iter->dir_vector[0], iter->dir_vector[1], iter->dir_vector[2] );
		// ##11
		nfd->WriteData( fp, "IIn", 
			iter->var_mode,
			iter->var_funcID );
		// ##12
		nfd->WriteStr( fp, iter->var_name );
		// ##13
		nfd->WriteStr( fp, iter->var_equation );
		// ##14-17
		for(int i=0; i<4; i++){
			nfd->WriteData( fp, "FFFn", 
				iter->var_locate[i][0], iter->var_locate[i][1], iter->var_locate[i][2] );
		}
		// ##18
		nfd->WriteData( fp, "FFFFn", 
			iter->var_value[0], iter->var_value[1], iter->var_value[2], iter->var_value[3] );
		// ##19
		nfd->WriteData( fp, "BBn", 
			iter->adjust_midside,
			iter->is_expanded2 );
	}
	nfd->WriteData( fp, "IIIIIIBn", -1,-1,-1,-1,-1,-1,0);
}


void CNFDB_507::write_temp_load( class CNFData* nfd, FILE* fp, std::vector<ctemp_load_rec>& list )
{
	vector<ctemp_load_rec>::iterator iter;
	for(iter = list.begin(); iter != list.end(); iter ++ ){
		// ##1
		nfd->WriteData( fp, "IIIFFIBn", 
			iter->ID,
			iter->color,
			iter->layer,
			iter->temp,
			iter->temp_co,
			iter->funcID,
			iter->is_extended );
	}
	nfd->WriteData( fp, "IIIFFIBn",-1,-1,-1, 0.0, 0.0, -1, 0);
}



void CNFDB_507::WriteData( class CNFData* nfd, FILE* fp )
{
	// #1
		nfd->WriteData( fp, "In", setID);
	// #2
		nfd->WriteStr( fp, title );

	// #3
		if( nfd->version >= 8.0 ){
			nfd->WriteData( fp, "IFBBBBFn",
				CSys,
				Def_temp,
				temp_on,
				grav_on,
				omega_on,
				Ref_temp_on,
				Ref_temp );
		} else {
			nfd->WriteData( fp, "IFBBBn",
				CSys,
				Def_temp,
				temp_on,
				grav_on,
				omega_on );
		}
	// #4,5
		nfd->WriteData( fp, "FFFn", grav[0], grav[1], grav[2]);
		nfd->WriteData( fp, "FFFn", grav[3], grav[4], grav[5]);
	// #6
		nfd->WriteData( fp, "FFFn", origin[0], origin[1], origin[2]);
	// #7
		nfd->WriteData( fp, "FFFn", omega[0], omega[1], omega[2]);
	// #8
		nfd->WriteData( fp, "FFFIn",
			stef_boltz,
			abs_temp,
			free_cnv_exp,
			rad_space_element );
	// #9
		nfd->WriteData( fp, "FFFFn",
			fc_flu_cond,
			fc_flu_cp,
			fc_flu_vis,
			fc_flu_dens );
	// #10
		nfd->WriteData( fp, "FFFFn",
			fc_cons_coeff,
			fc_reynolds,
			fc_pran_in,
			fc_pran_out);
	// #11
		nfd->WriteData( fp, "IIIn",
			tfc_flu_cond,
			tfc_flu_cp,
			tfc_flu_vis);
	// #12
		nfd->WriteData( fp, "BBBn",
			alt_free_conv,
			fc_flu_flag,
			fc_conv_flow );
	// #13
		nfd->WriteData( fp, "FFFFn",
			nl_arc_scale,
			nl_arcmaxadj,
			nl_arcminadj,
			nl_bounds_rb );
	// #14
		nfd->WriteData( fp, "FFFn", nl_conv[0], nl_conv[1], nl_conv[2]);
	// #15
		nfd->WriteData( fp, "FFFFFFn",
			nl_fstress,
			nl_lseach_tol,
			nl_mxadj_init,
			nl_max_rot,
			nl_stab_tol,
			nl_time_inc );
	// #16
		if( nfd->version >= 6.0 ){
			nfd->WriteData( fp, "FFFFFFFFFn",
				dyn_damp_ov,
				dyn_dampW3,
				dyn_dampW4,
				dyn_keep_freq[0],
				dyn_keep_freq[1],
				dyn_trans_dt,
				dyn_min_freq,
				dyn_max_freq,
				dyn_cluster_freq );
		} else {
			nfd->WriteData( fp, "FFFFFFn",
				dyn_damp_ov,
				dyn_dampW3,
				dyn_dampW4,
				dyn_keep_freq[0],
				dyn_keep_freq[1],
				dyn_trans_dt );
		}
	// #17
		nfd->WriteData( fp, "IIIIIIIIIn",
			nl_arc_const,
			nl_arc_iter,
			nl_arc_maxst,
			nl_div_limit,
			nl_dom_pdstp,
			nl_increment,
			nl_inter_out,
			nl_kstep,
			nl_mx_bisect );
	// #18
		nfd->WriteData( fp, "IIIIIIIIn",
			nl_max_iter,
			nl_max_lsrch,
			nl_out_iter,
			nl_quasi_newt,
			nl_sol_strat,
			nl_stiff_meth,
			nl_skip_adj,
			nl_sol_over );
	// #19
		nfd->WriteData( fp, "IIIIIIIn",
			dyn_freq_tab,
			dyn_damptab,
			dyn_keep_md,
			dyn_tran_ts,
			dyn_out_int,
			dyn_rand_psd,
			dyn_on_freq );
	// #20
		nfd->WriteData( fp, "BBBBBBBn",
			nl_on,
			nl_conv_flag[0],
			nl_conv_flag[1],
			nl_conv_flag[2],
			nl_mnewt_ls,
			nl_mnewt_qn,
			nl_mnewt_bs );
	// #21
		if( nfd->version >= 7.0 ){
			nfd->WriteData( fp, "BBBBBBIIIn",
				dyn_on,
				dyn_type,
				dyn_damp_method,
				dyn_massfrm,
				dyn_datarec,
				dyn_log_inter,
				dyn_freq_type,
				dyn_psd_type,
				dyn_psd_interpol );
		} else if( nfd->version >= 6.0 ){
			nfd->WriteData( fp, "BBBBBBIIn",
				dyn_on,
				dyn_type,
				dyn_damp_method,
				dyn_massfrm,
				dyn_datarec,
				dyn_log_inter,
				dyn_freq_type,
				dyn_psd_type );
		} else {
			nfd->WriteData( fp, "BBBBBn",
				dyn_on,
				dyn_type,
				dyn_damp_method,
				dyn_massfrm,
				dyn_datarec);
		}
	//------------------------------------
		write_structural_load( nfd, fp, structural_load_list );
	//------------------------------------
		write_geometric_load( nfd, fp, geometric_load_list );
	//------------------------------------
		write_temp_load( nfd, fp, ndtemp_load_list );
	//------------------------------------
		write_temp_load( nfd, fp, eltemp_load_list );
}


