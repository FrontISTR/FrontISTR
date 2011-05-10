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

/* conv_neu2fstr_static ver.1.0 */


#include <vector>
#include <stdlib.h>

#include "conv_neu2fstr_heat.h"
#include "cconv_mat.h"
#include "conv_util.h"

using namespace std;
using namespace n2h_util;

static
void ItoA( int i, char* s )
{
	sprintf(s, "%d", i);
}


//-----------------------------------------------------------------------------
// FIXTEMP
//-----------------------------------------------------------------------------

namespace {

class cntemp {
public:
	int nid;
	double value;
	cntemp( int id=-1, double v=0 ) :nid(id), value(v) {}
};

inline bool operator == (const cntemp& a, const cntemp& b ) { return a.nid == b.nid; }
inline bool operator <  (const cntemp& a, const cntemp& b ) { return a.nid <  b.nid; }
inline bool operator >  (const cntemp& a, const cntemp& b ) { return a.nid >  b.nid; }


} // of namespace



static
void SetFixtemp( CNFData& neu, CHECData& hec )
{
	if( neu.DB_507.size() == 0 ) return;

	set<cntemp> ntlist;

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		// node temperature -------------------------------------
		vector<CNFDB_507::ctemp_load_rec>::iterator siter;
		for( siter = p->ndtemp_load_list.begin();
		     siter != p->ndtemp_load_list.end(); siter++) {
			ntlist.insert( cntemp( siter->ID, siter->temp ));
		}
		// element temperature ----------------------------------
		for( siter = p->eltemp_load_list.begin();
		     siter != p->eltemp_load_list.end(); siter++) {
			CHECDB_Element::CElemItem* e = hec.GetElemItem( siter->ID );
			if(!e) {
				printf("##Warning : Not exist element %d (FIXTEMP)\n", (int)siter->ID );
				continue; 
			}
			for(int i=0; i<e->node_n; i++ ){
				ntlist.insert( cntemp( e->node[i], siter->temp ));
			}
		}
	}

	if( ntlist.size() > 0 ){
		CFSTRDB_Fixtemp* dblock = new CFSTRDB_Fixtemp();
		set<cntemp>::iterator it;
		for(it = ntlist.begin(); it != ntlist.end(); it++) {
			char name[256];
			ItoA( it->nid, name);
			CFSTRDB_Fixtemp::CItem item( name, it->value );
			dblock->ItemList.push_back( item );
		}
		hec.DB.push_back( dblock );
	}
}


//-----------------------------------------------------------------------------
// CFLUX
//-----------------------------------------------------------------------------


static
void SetCFlux( CNFData& neu, CHECData& hec )
{
	const int FEA_load_nHeatFlux = 10;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_CFlux* dblock = new CFSTRDB_CFlux();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_nHeatFlux ) continue;
			char name[256];
			ItoA( siter->loadID, name);
			CFSTRDB_CFlux::CItem item( name, siter->value[0] );
			dblock->ItemList.push_back( item );
		}
	}

	if( dblock->ItemList.size()>0 ){
		hec.DB.push_back( dblock );
	} else {
		delete dblock;
	}
}


//-----------------------------------------------------------------------------
// DFLUX
//-----------------------------------------------------------------------------


static
void SetDFlux( CNFData& neu, CHECData& hec )
{
	const int FEA_load_eHeatFlux = 44;
	const int FEA_load_eHeatGen = 47;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_DFlux* dblock = new CFSTRDB_DFlux();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_eHeatFlux 
			   && siter->loadtype != FEA_load_eHeatGen ) continue;
 			int eid = siter->loadID;
			int surf_no = siter->dof_face[0];
			int hec_e_type = hec.GetElemType( eid );
			if( hec_e_type == 0 ) continue;
			int load_type;
			double load;
			if( siter->loadtype == FEA_load_eHeatFlux ) {
				int fg_front;
				load_type = CFSTRDB_DFlux::TYPE_S0 + hec_face_no( hec_e_type, surf_no, fg_front );
				load = fg_front ? siter->value[0]: -siter->value[0];
			} else {
				load_type = CFSTRDB_DFlux::TYPE_BF;
				load = siter->value[0];
			}
			char name[256];
			ItoA( siter->loadID, name);
			CFSTRDB_DFlux::CItem item( name, load_type, load );
			dblock->ItemList.push_back( item );
		}
	}

	if( dblock->ItemList.size()>0 ){
		hec.DB.push_back( dblock );
	} else {
		delete dblock;
	}
}


//-----------------------------------------------------------------------------
// FILM
//-----------------------------------------------------------------------------


static
void SetFilm( CNFData& neu, CHECData& hec )
{
	const int FEA_load_eConvection = 45;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_Film* dblock = new CFSTRDB_Film();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_eConvection ) continue;
			int eid = siter->loadID;
			int surf_no = siter->dof_face[0];
			int hec_e_type = hec.GetElemType( eid );
			if( hec_e_type == 0 ) continue;

			double value = siter->value[0]; // coeff of heat transfer
			double sink  = siter->value[2]; // emvironmental temperature

			int fg_front;
			int load_type = CFSTRDB_Film::TYPE_F0 + hec_face_no( hec_e_type, surf_no, fg_front );

			char name[256];
			ItoA( siter->loadID, name);
			CFSTRDB_Film::CItem item( name, load_type, value, sink);
			dblock->ItemList.push_back( item );
		}
	}

	if( dblock->ItemList.size()>0 ){
		hec.DB.push_back( dblock );
	} else {
		delete dblock;
	}
}



//-----------------------------------------------------------------------------
// RADIATE
//-----------------------------------------------------------------------------


static
void SetRadiate( CNFData& neu, CHECData& hec )
{
	const int FEA_load_eRadiation = 46;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_Radiate* dblock = new CFSTRDB_Radiate();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_eRadiation ) continue;
			int eid = siter->loadID;
			int surf_no = siter->dof_face[0];
			int hec_e_type = hec.GetElemType( eid );
			if( hec_e_type == 0 ) continue;

			double value = siter->value[0];
			double sink  = siter->value[2]; // emvironmental temperature

			int fg_front;
			int load_type = CFSTRDB_Radiate::TYPE_R0 + hec_face_no( hec_e_type, surf_no, fg_front );

			char name[256];
			ItoA( siter->loadID, name);
			CFSTRDB_Radiate::CItem item( name, load_type, value, sink);
			dblock->ItemList.push_back( item );
		}
	}

	if( dblock->ItemList.size()>0 ){
		hec.DB.push_back( dblock );
	} else {
		delete dblock;
	}
}


//=============================================================================
// conv_neu2fstr_heat
//=============================================================================

void conv_neu2fstr_heat( CNFData& neu, CHECData& hec )
{
	SetFixtemp( neu, hec );
	SetCFlux  ( neu, hec );
	SetDFlux  ( neu, hec );
	SetFilm   ( neu, hec );
	SetRadiate( neu, hec );
}

