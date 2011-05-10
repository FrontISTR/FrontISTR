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

#ifndef CNFDB_507H
#define CNFDB_507H


#include <vector>
#include "CNFDataBlock.h"

// 507 Loads



class CNFDB_507 : public CNFDataBlock {
public:
	CNFDB_507();
	virtual ~CNFDB_507() {}
	
	virtual void Read( CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	class cstructural_load_rec {
	public:
		// ##1
		nf_int loadID;
		nf_int loadtype;
		nf_int color;
		nf_int layer;
		nf_int define_sys;
		nf_int subtype;
		nf_bool is_expanded;
		// ##2
		nf_int dof_face[3];
		// ##3
		nf_float value[5];
		// ##4
		nf_int functions[5];
		// ##5
		nf_bool Enclosure;
		nf_bool can_shade;
		nf_bool can_be_shaded;
		nf_int  add1_id[2];
		// ##6
		nf_int dir_func[3];
		// ##7
		nf_float direction[3];
		cstructural_load_rec() {}
		cstructural_load_rec( const cstructural_load_rec& r);
		cstructural_load_rec& operator=( const cstructural_load_rec& r);
		void copy_from( const cstructural_load_rec& r );
	};	
	// -----------------------------
	class cgeometric_load_rec {
	public:
		// ##1
		nf_int loadID;
		nf_int loadtype;
		nf_int color;
		nf_int layer;
		nf_int define_sys;
		nf_int subtype;
		nf_bool is_expanded;
		// ##2
		nf_int dof_face[3];
		// ##3
		nf_float value[5];
		// ##4
		nf_int functions[5];
		// ##5
		nf_bool Enclosure;
		nf_bool can_shade;
		nf_bool can_be_shaded;
		nf_int addl_id[2];
		// ##6
		nf_int dir_func[3];
		// ##7
		nf_float direction[3];
		// ##8
		nf_int dir_mode;
		nf_int dir_id;
		// ##9
		nf_float dir_base[3];
		// ##10
		nf_float dir_vector[3];
		// ##11
		nf_int var_mode;
		nf_int var_funcID;
		// ##12
		nf_char var_name[256];
		// ##13
		nf_char var_equation[256];
		// ##14-17
		nf_float var_locate[4][3];
		// ##18
		nf_float var_value[4];
		// ##19
		nf_bool adjust_midside;
		nf_bool is_expanded2;
		cgeometric_load_rec() {}
		cgeometric_load_rec( const cgeometric_load_rec& r);
		cgeometric_load_rec& operator=( const cgeometric_load_rec& r);
		void copy_from( const cgeometric_load_rec& r);
	};	
	// -----------------------------
	class ctemp_load_rec {
	public:
		// ##1
		nf_int ID;
		nf_int color;
		nf_int layer;
		nf_float temp;
		nf_float temp_co;
		nf_int funcID;
		nf_bool is_extended;
	};
	// -----------------------------

	// #1
		nf_int setID;
	// #2
		nf_char title[26];
	// #3
		nf_int CSys;
		nf_float Def_temp;
		nf_bool temp_on;
		nf_bool grav_on;
		nf_bool omega_on;
		nf_bool Ref_temp_on;
		nf_float Ref_temp;
	// #4,5
		nf_float grav[6];
	// #6
		nf_float origin[3];
	// #7
		nf_float omega[3];
	// #8
		nf_float stef_boltz;
		nf_float abs_temp;
		nf_float free_cnv_exp;
		nf_int rad_space_element;
	// #9
		nf_float fc_flu_cond;
		nf_float fc_flu_cp;
		nf_float fc_flu_vis;
		nf_float fc_flu_dens;
	// #10
		nf_float fc_cons_coeff;
		nf_float fc_reynolds;
		nf_float fc_pran_in;
		nf_float fc_pran_out;
	// #11
		nf_int tfc_flu_cond;
		nf_int tfc_flu_cp;
		nf_int tfc_flu_vis;
	// #12
		nf_bool alt_free_conv;
		nf_bool fc_flu_flag;
		nf_bool fc_conv_flow;
	// #13
		nf_float nl_arc_scale;
		nf_float nl_arcmaxadj;
		nf_float nl_arcminadj;
		nf_float nl_bounds_rb;
	// #14
		nf_float nl_conv[3];
	// #15
		nf_float nl_fstress;
		nf_float nl_lseach_tol;
		nf_float nl_mxadj_init;
		nf_float nl_max_rot;
		nf_float nl_stab_tol;
		nf_float nl_time_inc;
	// #16
		nf_float dyn_damp_ov;
		nf_float dyn_dampW3;
		nf_float dyn_dampW4;
		nf_float dyn_keep_freq[2];
		nf_float dyn_trans_dt;
		nf_float dyn_min_freq;
		nf_float dyn_max_freq;
		nf_float dyn_cluster_freq;
	// #17
		nf_int nl_arc_const;
		nf_int nl_arc_iter;
		nf_int nl_arc_maxst;
		nf_int nl_div_limit;
		nf_int nl_dom_pdstp;
		nf_int nl_increment;
		nf_int nl_inter_out;
		nf_int nl_kstep;
		nf_int nl_mx_bisect;
	// #18
		nf_int nl_max_iter;
		nf_int nl_max_lsrch;
		nf_int nl_out_iter;
		nf_int nl_quasi_newt;
		nf_int nl_sol_strat;
		nf_int nl_stiff_meth;
		nf_int nl_skip_adj;
		nf_int nl_sol_over;
	// #19
		nf_int dyn_freq_tab;
		nf_int dyn_damptab;
		nf_int dyn_keep_md;
		nf_int dyn_tran_ts;
		nf_int dyn_out_int;
		nf_int dyn_rand_psd;
		nf_int dyn_on_freq;
	// #20
		nf_bool nl_on;
		nf_bool nl_conv_flag[3];
		nf_bool nl_mnewt_ls;
		nf_bool nl_mnewt_qn;
		nf_bool nl_mnewt_bs;
	// #21
		nf_bool dyn_on;
		nf_bool dyn_type;
		nf_bool dyn_damp_method;
		nf_bool dyn_massfrm;
		nf_bool dyn_datarec;
		nf_bool dyn_log_inter;
		nf_int dyn_freq_type;
		nf_int dyn_psd_type;
		nf_int dyn_psd_interpol;
	// -----------------------------------------
	std::vector<cstructural_load_rec> structural_load_list;
	std::vector<cgeometric_load_rec> geometric_load_list;
	std::vector<ctemp_load_rec> ndtemp_load_list;
	std::vector<ctemp_load_rec> eltemp_load_list;
protected:
	void read_structural_load( class CNFData* nfd, std::vector<cstructural_load_rec>& list );
	void read_geometric_load( class CNFData* nfd, std::vector<cgeometric_load_rec>& list );
	void read_temp_load( class CNFData* nfd, std::vector<ctemp_load_rec>& list );
	void write_structural_load( class CNFData* nfd, FILE* fp, std::vector<cstructural_load_rec>& list );
	void write_geometric_load( class CNFData* nfd, FILE* fp, std::vector<cgeometric_load_rec>& list );
	void write_temp_load( class CNFData* nfd, FILE* fp, std::vector<ctemp_load_rec>& list );
};






#endif
