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
	CNFDB_408 Ver.1.0
*/

#ifndef CNFDB_408H
#define CNFDB_408H


#include <vector>
#include "CNFDataBlock.h"


// 408 Group

class CNFDB_408 : public CNFDataBlock {
public:
	CNFDB_408();

	virtual ~CNFDB_408() {}
	virtual void Read( class CNFData* nfd );
	virtual void WriteData( class CNFData* nfd, FILE* fp );

public:
	class cplclip_rec {
	public:
		// ##1
		nf_int plclip_on;
		nf_int plclip_neg;
		// ##2
		nf_float plclip_base[3];
		// ##3
		nf_float plclip_norm[3];
		cplclip_rec() {}
		cplclip_rec( const cplclip_rec& c )
		: plclip_on( c.plclip_on), plclip_neg(c.plclip_neg) {
			plclip_base[0] = c.plclip_base[0];
			plclip_base[1] = c.plclip_base[1];
			plclip_base[2] = c.plclip_base[2];
			plclip_norm[0] = c.plclip_norm[0];
			plclip_norm[1] = c.plclip_norm[1];
			plclip_norm[2] = c.plclip_norm[2];
		}
		cplclip_rec& operator=( const cplclip_rec& c ) {
			plclip_on = c.plclip_on;
			plclip_neg = c.plclip_neg;
			plclip_base[0] = c.plclip_base[0];
			plclip_base[1] = c.plclip_base[1];
			plclip_base[2] = c.plclip_base[2];
			plclip_norm[0] = c.plclip_norm[0];
			plclip_norm[1] = c.plclip_norm[1];
			plclip_norm[2] = c.plclip_norm[2];
			return *this;
		}
	};
	class crule_rec {
	public:
		// ##1
		nf_int rule_type;

		class centry_rec {
		public:
			nf_int startID;
			nf_int stopID;
			nf_int incID;
			nf_int include;
		};
		std::vector<centry_rec> entry_set;
		crule_rec& operator = (const crule_rec& rec) {
			rule_type = rec.rule_type;
			entry_set = rec.entry_set;
			return *this;
		}
	};
	class clist_rec {
	public:
		nf_int list_type;
		std::vector<nf_int> entityID;
	};

	// #1
		nf_int ID;
		nf_int need_eval;
		nf_bool prev_enum;
	// #2
		nf_char title[26];
	// #3
		nf_int layer[2];
		nf_int layer_method;
	// #4
		nf_int coclip_on;
		nf_int coclip_dof;
		nf_int coclip_meth;
		nf_int coclip_csys;
		nf_float coclip_min;
		nf_float coclip_max;
	// #5
		nf_int plclip_meth;
		nf_int plclip_in;
	// #--------------------------------
		cplclip_rec plclip_set[6];
	// #--------------------------------
		nf_int max_rules;
	// #--------------------------------
		std::vector<crule_rec> rule_set;
	// #--------------------------------
		nf_int max_lists;
	// #--------------------------------
		std::vector<clist_rec> list_set;
protected:
	void read_rule( class CNFData* nfd, std::vector<crule_rec>& list );
	void read_list( class CNFData* nfd, std::vector<clist_rec>& list );
	void write_rule( class CNFData* nfd, FILE* fp, std::vector<crule_rec>& list );
	void write_list( class CNFData* nfd, FILE* fp, std::vector<clist_rec>& list );
};

#endif
