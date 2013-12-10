/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/07/20                                        *
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
	conv_neu2fstr_static ver.1.01
*/


#include <vector>
#include <set>
#include <stdlib.h>

#include "conv_neu2fstr_static.h"
#include "cconv_mat.h"
#include "conv_util.h"

using namespace std;
using namespace n2h_util;


static
void ItoA( int i, char* s )
{
	sprintf(s, "%d", i);
}



static
int mesh_dof = 3;

static
bool check_mesh_for_fstr( CHECData& hec )
{
	bool fg_dof[4] = { false, false, false, false };
	int dof_tbl[4] = { 1,2,3,6 };

	vector<CHECDataBlock*>::iterator iter;
	for(iter = hec.DB.begin(); iter != hec.DB.end(); iter++) {
		if((*iter)->data_type == HECDB_ELEMENT){
			CHECDB_Element* e = (CHECDB_Element*)*iter;
			switch( e->DOFNumber( e->type )){
			case 1: fg_dof[0] = true; break;
			case 2: fg_dof[1] = true; break;
			case 3: fg_dof[2] = true; break;
			case 6: fg_dof[3] = true; break;
			default:
				assert(0);
			}
		}
	}

	int c = 0;
	for( int i=0; i<4 ; i++) {
		if( fg_dof[i] ) {
			c++;
			mesh_dof = dof_tbl[i];
		}
	}
	if( c != 1 ) {
		printf("### FrontSTR does not permit plural DOF\n");
		fflush(stdout);
		return false;
	}
	return true;
}



//-----------------------------------------------------------------------------
// BOUNADARY
//-----------------------------------------------------------------------------



class cbitem {
public:
	int nid; //node id
	double value;
	cbitem( int id=-1, double v=0 ): nid(id), value(v) {}
	cbitem( const cbitem& i ) : nid(i.nid), value(i.value) {}
	cbitem& operator=( const cbitem& i ) {
		nid = i.nid;
		value = i.value;
		return *this;
	}
	void set_value( double v ) { value = v; }
};


inline bool operator == (const cbitem& a, const cbitem& b ) { return a.nid == b.nid; }
inline bool operator <  (const cbitem& a, const cbitem& b ) { return a.nid <  b.nid; }
inline bool operator >  (const cbitem& a, const cbitem& b ) { return a.nid >  b.nid; }

/*
static
CHECDB_NGroup* create_new_ngrp( CHECData& hec, const char* name )
{
	CHECDB_NGroup* ngrp = new CHECDB_NGroup();
	strcpy( ngrp->name, name );
	return ngrp;
}
*/

static
void set_boundary_node_by_506( set<cbitem>* bitem, int dof_n,
	vector<CNFDB_506::cconst_item>& item, CHECData& hec )
{
	vector<CNFDB_506::cconst_item>::iterator iter;
	for(iter = item.begin(); iter != item.end(); iter++) {
		for( int i=0; i<dof_n ; i++) {
			if( iter->DOF[i] ){
				bitem[i].insert( cbitem(iter->ID));
			}
		}
	}
}



// execute after all executing set_boundary_node_by_506

static
void set_boundary_node_by_507( set<cbitem>* bitem,
	vector<CNFDB_507::cstructural_load_rec>& item, CHECData& hec )
{
	const int load_dof = 3;
	const int FEA_load_nDisplacement = 3;
	vector<cbitem> new_bitem[load_dof];

	vector<CNFDB_507::cstructural_load_rec>::iterator iter;
	for( iter = item.begin(); iter != item.end(); iter++) {
		if( iter->loadtype != FEA_load_nDisplacement ) continue;
		set<cbitem>::iterator biter;
		for(int i=0; i<load_dof; i++) { // not dof_n !!
//			biter = find( bitem[i].begin(), bitem[i].end(), cbitem( iter->loadID ));
			for( biter = bitem[i].begin(); biter != bitem[i].end(); biter ++ ){
				if( biter->nid == iter->loadID ) break;
			}
			if( biter == bitem[i].end()) continue;
			bitem[i].erase( biter );
			new_bitem[i].push_back( cbitem(iter->loadID, iter->value[i]));
		}
	}

	for(int i=0; i<load_dof; i++) {
		vector<cbitem>::iterator biter;
		for( biter = new_bitem[i].begin(); biter != new_bitem[i].end(); biter ++ ){
			bitem[i].insert( *biter );
		}
	}
}


// CAUTION) NOT SUPPORTED CURVE AND SURFACE NOW!!
static
void SetBoundary( CNFData& neu, CHECData& hec )
{
	const int max_dof_n = 6;
	set<cbitem> bitem[ max_dof_n ];
	int dof_n = mesh_dof;

	// search 506
	{
		vector<CNFDB_506*>::iterator iter;
		for(iter = neu.DB_506.begin(); iter != neu.DB_506.end(); iter++){
			CNFDB_506* p = *iter;
			set_boundary_node_by_506( bitem, dof_n, p->const_nodes, hec );
		}
	}
	// search 507
	{
		vector<CNFDB_507*>::iterator iter;
		for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
			CNFDB_507* p = *iter;
			set_boundary_node_by_507( bitem, p->structural_load_list, hec );
		}
	}
	// registration
	{
		CFSTRDB_Boundary* bc = new CFSTRDB_Boundary();
		for(int i=0; i<dof_n; i++) {
			set<cbitem>::iterator iter;
			for( iter = bitem[i].begin(); iter != bitem[i].end(); iter++){
				char name[256];
				ItoA( iter->nid, name);
				int dof = i+1;
				CFSTRDB_Boundary::CItem item( name, dof, dof, iter->value );
				bc->ItemList.push_back( item );
			}
		}
		if( bc->ItemList.size() > 0 ){
			hec.DB.push_back(bc);
		} else {
			delete bc;
		}
	}
}


//-----------------------------------------------------------------------------
// CLOAD
//-----------------------------------------------------------------------------


static
void SetCLoad( CNFData& neu, CHECData& hec )
{
	const int load_dof = 3;
	const int FEA_load_nForce = 1;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_CLoad* cload = new CFSTRDB_CLoad();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_nForce ) continue;
			for( int i=0; i<load_dof; i++) {
				if( siter->dof_face[i] == 0) continue;
				char name[256];
				ItoA( siter->loadID, name);
				int dof = i+1;
				CFSTRDB_CLoad::CItem item( name, dof, siter->value[i] );
				cload->ItemList.push_back( item );
			}
		}
	}

	if( cload->ItemList.size()>0 ){
		hec.DB.push_back( cload );
	} else {
		delete cload;
	}
}


//-----------------------------------------------------------------------------
// DLOAD
//-----------------------------------------------------------------------------


static
void set_dload_grav( CNFDB_507* block, CFSTRDB_DLoad* dload )
{
	if( !block->grav_on ) return;

	double dir[3];
	double g;

	dir[0] = block->grav[0];
	dir[1] = block->grav[1];
	dir[2] = block->grav[2];
	g = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
	if( g == 0 ) return;
	g = sqrt(g);
	dir[0] /= g;
	dir[1] /= g;
	dir[2] /= g;

	CFSTRDB_DLoad::CItem item( "ALL", CFSTRDB_DLoad::TYPE_GRAV );
	item.param[0] = g;
	item.param[1] = dir[0];
	item.param[2] = dir[1];
	item.param[3] = dir[2];
	dload->ItemList.push_back( item );
}


static
void set_dload_cent( CNFDB_507* block, CFSTRDB_DLoad* dload )
{
	if( !block->omega_on ) return;

	double dir[3];
	double omega;

	dir[0] = block->omega[0];
	dir[1] = block->omega[1];
	dir[2] = block->omega[2];
	omega = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
	if( omega == 0 ) return;
	omega = sqrt(omega);
	dir[0] /= omega;
	dir[1] /= omega;
	dir[2] /= omega;

	CFSTRDB_DLoad::CItem item( "ALL", CFSTRDB_DLoad::TYPE_CENT );
	item.param[0] = omega;
	item.param[1] = block->origin[0];
	item.param[2] = block->origin[1];
	item.param[3] = block->origin[2];
	item.param[4] = dir[0];
	item.param[5] = dir[1];
	item.param[6] = dir[2];
	dload->ItemList.push_back( item );
}


static
void SetDLoad( CNFData& neu, CHECData& hec )
{
	const int FEA_load_ePressure = 42;

	if( neu.DB_507.size() == 0 ) return;

	CFSTRDB_DLoad* dload = new CFSTRDB_DLoad();

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;

		set_dload_grav( p, dload );
		set_dload_cent( p, dload );

		vector<CNFDB_507::cstructural_load_rec>::iterator siter;
		for( siter = p->structural_load_list.begin();
		     siter != p->structural_load_list.end(); siter++) {
			if( siter->loadtype != FEA_load_ePressure ) continue;
			int eid = siter->loadID;
			int surf_no = siter->dof_face[0];
			int hec_e_type = hec.GetElemType( eid );
			if( hec_e_type == 0 ) continue;
			int fg_front;
			int load_type = CFSTRDB_DLoad::TYPE_P0 + hec_face_no( hec_e_type, surf_no, fg_front );
			double load = fg_front ? siter->value[0]: -siter->value[0];
			char name[256];
			ItoA( siter->loadID, name);
			CFSTRDB_DLoad::CItem item( name, load_type );
			item.param[0] = load;
			dload->ItemList.push_back( item );
		}
	}

	if( dload->ItemList.size()>0 ){
		hec.DB.push_back( dload );
	} else {
		delete dload;
	}
}


//-----------------------------------------------------------------------------
// TEMPERATURE
//-----------------------------------------------------------------------------


class cnode_temp {
public:
	int nid;
	double value;
	cnode_temp( int id=-1, double v=0 ) : nid(id), value(v) {}
};

inline bool operator==( const cnode_temp& a, const cnode_temp& b ) { return a.nid == b.nid; }
inline bool operator< ( const cnode_temp& a, const cnode_temp& b ) { return a.nid <  b.nid; }
inline bool operator> ( const cnode_temp& a, const cnode_temp& b ) { return a.nid >  b.nid; }


static
void set_temp_node( vector<CNFDB_507::ctemp_load_rec>& tload, set<cnode_temp>& nt )
{
	vector<CNFDB_507::ctemp_load_rec>::iterator iter;
	for(iter = tload.begin(); iter != tload.end(); iter++) {
		cnode_temp item( iter->ID, iter->temp );
		nt.insert( item );
	}
}

static
void set_temp_elem( CHECData& hec, vector<CNFDB_507::ctemp_load_rec>& tload, set<cnode_temp>& nt )
{
	vector<CNFDB_507::ctemp_load_rec>::iterator iter;
	for(iter = tload.begin(); iter != tload.end(); iter++) {
		CHECDB_Element::CElemItem* elem = hec.GetElemItem( iter->ID );
		if( !elem ) continue;
		for(int i=0; i<elem->node_n; i++) {
			cnode_temp item( elem->node[i], iter->temp );
			nt.insert( item );
		}
	}
}



static
void SetTemperature( CNFData& neu, CHECData& hec )
{
	if( neu.DB_507.size() == 0 ) return;

	 set<cnode_temp> nt;

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		set_temp_node( p->ndtemp_load_list, nt );
		set_temp_elem( hec, p->eltemp_load_list, nt );
	}

	if( nt.size() > 0 ){
		CFSTRDB_Temperature* temp = new CFSTRDB_Temperature();
		hec.DB.push_back( temp );
	}
}


//-----------------------------------------------------------------------------
// REFTEMP
//-----------------------------------------------------------------------------

static
void SetReftemp( CNFData& neu, CHECData& hec )
{
	if( neu.DB_507.size() == 0 ) return;
	if( neu.version < 8.0 ) return;

	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		if( p->Ref_temp_on ){
			CFSTRDB_Reftemp* reftemp = new CFSTRDB_Reftemp();
			reftemp->value = p->Ref_temp;
			hec.DB.push_back( reftemp );
		}
	}
}


//=============================================================================
// conv_neu2fstr_static
//=============================================================================


void conv_neu2fstr_static( CNFData& neu, CHECData& hec )
{
	if(!check_mesh_for_fstr( hec )) throw;
	SetBoundary   ( neu, hec );
	SetCLoad      ( neu, hec );
	SetDLoad      ( neu, hec );
	SetTemperature( neu, hec );
	SetReftemp    ( neu, hec );
}

