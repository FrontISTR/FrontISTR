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
	conv_neu2hec ver.2.0
	-------------------------------------
	Rule for name of materials
	MAT<ID>   <ID> : ID number of material
*/



#include <vector>
#include "conv_neu2hec.h"
#include "cconv_mat.h"
#include "conv_util.h"
#include "CConvMessage.h"

using namespace std;
using namespace n2h_util;

static int solution = sol_static;
static bool fg_element_error_stop = true;

//=============================================================================
// Coordinate Converting
//=============================================================================

static
void set_conv_mat( cconv_mat& m, double o[], double rot[] )
{
	cconv_mat a, rx, ry, rz;
	a.transfer( -o[0], -o[1], -o[2] );
	rx.rotate( 'x', -rot[0] );
	ry.rotate( 'y', -rot[1] );
	rz.rotate( 'z', -rot[2] );

	m = a * rx;
	m *= ry;
	m *= rz;
}


static
CNFDB_405* search_405( CNFData& neu, int id )
{
	if(id<=0) return 0;

	vector<CNFDB_405*>::iterator iter;
	for(iter = neu.DB_405.begin(); iter != neu.DB_405.end(); iter++){
		if( (*iter)->ID == id ) return *iter;
	}
	return 0;
}


static
bool set_conv_mat(  CNFData& neu, int id, cconv_mat& m )
{
	m.unit();

	CNFDB_405* coord = search_405( neu, id );
	if(!coord) {
		return (id==0);
	}

	while(1) {
		cconv_mat a;
		set_conv_mat( a, coord->origin, coord->rot );
		switch( coord->type ){
		case 0:
			a.type = coord_t_cartesian;
			break;
		case 1:
			a.type = coord_t_cylinder;
			break;
		case 2:
			a.type = coord_t_sphere;
			break;
		default:
			assert(0);
		}
		m = m * a;
		if( coord->define_sys == 0 ) break;
		coord = search_405( neu, coord->define_sys );
		if(!coord) {
			return false;
		}
	}
	return true;
}


//=============================================================================
// Converting Block Data
//=============================================================================

//-----------------------------------------------------------------------------
// HEADER
//-----------------------------------------------------------------------------

static
void SetHeader( CNFData& neu, CHECData& hec )
{
	CHECDB_Header* header = new CHECDB_Header();
	strcpy( header->title, neu.title );
	hec.DB.push_back( header );
}

//-----------------------------------------------------------------------------
// NODE
//-----------------------------------------------------------------------------

static
void SetNode( CNFData& neu, CHECData& hec )
{
	if( neu.DB_403.size() == 0 ) return;

	cconv_mat m;
	int coord_id = 0;
	m.unit();

	CHECDB_Node* node = new CHECDB_Node();
	CHECDB_Node::CNodeItem nitem;

	int node_n = 0;
	vector<CNFDB_403*>::iterator iter;
	for(iter = neu.DB_403.begin(); iter != neu.DB_403.end(); iter++){
		CNFDB_403* n = *iter;
		if( n->define_sys != coord_id ) {
			coord_id = n->define_sys;
			if(!set_conv_mat( neu, coord_id, m ))
				throw CConvMessage( CONV_COORDINATE_ERROR, "NODE" );
		}
		nitem.ID = n->ID;
		m.convert( n->x, n->y, n->z, nitem.x, nitem.y, nitem.z );
		node->NodeList.insert(nitem);
		node_n++;
	}

	hec.DB.push_back(node);
	fprintf(stdout,"%d nodes converted\n", node_n);
	fflush(stdout);
}


//-----------------------------------------------------------------------------
// ELEMENT
//-----------------------------------------------------------------------------


static
CHECDB_Element* search_element( CHECData& hec, int type, int sec_id=0, int option=0 )
{
	vector<CHECDataBlock*>::iterator iter;
	for( iter = hec.DB.begin(); iter != hec.DB.end(); iter++) {
		if( (*iter)->data_type == HECDB_ELEMENT ){
			CHECDB_Element* e = (CHECDB_Element*)(*iter);
			if( e->type == type
			 && e->sec_id == sec_id
			 && e->option == option ) return e;
		}
	}
	return 0;
}



static int line_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Line2 );

	switch( e->type ){ // material property 
	case NEU_ELEM_PROP_ROD:
	case NEU_ELEM_PROP_LINK:
		if( fg ) return 111;
		else     return 112;
	case NEU_ELEM_PROP_BEAM:
	case NEU_ELEM_PROP_BEAM2:
	case NEU_ELEM_PROP_CURVEBEAM:
	case NEU_ELEM_PROP_BAR:
		if( fg ) return 611;
		else     return 612;
		break;
	default:
		throw CConvMessage( CONV_INVALID_ELEMENT_PROPERTY, "Line%d, type:%d", (fg?2:3), e->type );
	}
}


static int tri_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Tri3 );

	switch( e->type ){
	case NEU_ELEM_PROP_PLANESTRAIN:
	case NEU_ELEM_PROP_PLANESTRAIN2:
		if( fg ) return 231;
		else     return 232;
	case NEU_ELEM_PROP_PLATE:
	case NEU_ELEM_PROP_PLATE2:
		if( fg ) return 731;
		else     return 732;
	default:
		throw CConvMessage( CONV_INVALID_ELEMENT_PROPERTY, "Tri%d, type:%d", (fg?3:6), e->type );
	}
}


static int quad_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Quad4 );
	switch( e->type ){
	case NEU_ELEM_PROP_PLANESTRAIN:
	case NEU_ELEM_PROP_PLANESTRAIN2:
		if( fg ) return 241;
		else     return 242;
	case NEU_ELEM_PROP_PLATE:
	case NEU_ELEM_PROP_PLATE2:
		if( fg ) return 741;
		else     return 742;
	default:
		throw CConvMessage( CONV_INVALID_ELEMENT_PROPERTY, "Quad%d, type:%d", (fg?4:8), e->type );
	}
}


static int tetra_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Tetra4 );
	if( fg ) return 341;
	else     return 342;
}

static int wedge_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Wedge6 );
	if( fg ) return 351;
	else     return 352;
}


static int brick_elem_type( CNFDB_404* e )
{
	bool fg = ( e->topology == CNFDB_404::top_Brick8 );
	if( fg ) return 361;
	else     return 362;
}


enum {
	opt_plane_stress = 0,
	opt_plane_strain = 1
};


static int get_hec_elem_option( CNFDB_404* e, int hec_type )
{
	switch( hec_type ){
	case 231: case 232: case 241: case 242:
		if( e->formulation2 == 0 )
			return opt_plane_stress;
		else
			return opt_plane_strain;
	default:
		return 0;
	}
}


static
void SetElement( CNFData& neu, CHECData& hec )
{
	const int con_table[12][20] = {
		// line
		{ 0,1 },	// 111, 611
		{ 0,1,2 },	// 112, 612

		// tri
		{ 0,1,2 },	// 231, 731
		//{ 0,1,2,4,5,3 },// 232, 732
		  { 0,1,2,5,6,4 },// 232, 732
		//  0 1 2 3 4 5
		//  0 1 2 4 5 6

		// quad
		{ 0,1,2,3 },		// 241, 741
		{ 0,1,2,3,4,5,6,7 },	// 242, 742

		// tetra
		{ 0,1,2,4 },// 341
		// 0 1 2 3 5  6  4 7 8  9
		{  0,1,2,4,9,10, 8,12,13,14 },// 342

		// Wedge
		{ 0,1,2,4,5,6 },			// 351
		//  1 2 3 4 5 6  7 8 9   10 11 12 13 14 15 -- hec
		//{ 0,1,2,3,4,5, 7,8,6,  13,14,12, 9,10,11},	// 352
		{   0,1,2,4,5,6, 9,10,8, 17,18,16,12,13,14},	// 352
		//  0 1 2 3 4 5 6 7  8  9  10 11 12 13 14
		//  0 1 2 4 5 6 8 9 10 12  13 14 16 17 18

		// brick
		{ 0,1,2,3,4,5,6,7, 8,9,10,11, 16,17,18,19, 12,13,14,15 },// 361
		{ 0,1,2,3,4,5,6,7, 8,9,10,11, 16,17,18,19, 12,13,14,15 } // 362
	};

	int elem_count = 0;
	int ignored_count = 0;

	if( neu.DB_404.size() == 0 ) return;

	CHECDB_Element* hec_e = 0;

	vector<CNFDB_404*>::iterator iter;
	for(iter = neu.DB_404.begin(); iter != neu.DB_404.end(); iter++){
		CNFDB_404* e = *iter;
		int con;
		int nn;
		int hec_type;
		int option = 0;
		int sec_id = e->propID;

		elem_count++;

		try {
			switch( e->topology ){
			case CNFDB_404::top_Line2:
				con = 0; nn = 2; hec_type = line_elem_type( e );
				break;
			case CNFDB_404::top_Line3:
				con = 1; nn = 3; hec_type = line_elem_type( e );
				break;
			case CNFDB_404::top_Tri3:
				con = 2; nn = 3; hec_type = tri_elem_type( e );
				break;
			case CNFDB_404::top_Tri6:
				con = 3; nn = 6; hec_type = tri_elem_type( e );
				break;
			case CNFDB_404::top_Quad4:
				con = 4; nn = 4; hec_type = quad_elem_type( e );
				break;
			case CNFDB_404::top_Quad8:
				con = 5; nn = 8; hec_type = quad_elem_type( e );
				break;
			case CNFDB_404::top_Tetra4:
				con = 6; nn = 4; hec_type = tetra_elem_type( e );
				break;
			case CNFDB_404::top_Tetra10:
				con = 7; nn =10; hec_type = tetra_elem_type( e );
				break;
			case CNFDB_404::top_Wedge6:
				con = 8; nn =6; hec_type = wedge_elem_type( e );
				break;
			case CNFDB_404::top_Wedge15:
				con = 9; nn =15; hec_type = wedge_elem_type( e );
				break;
			case CNFDB_404::top_Brick8:
				con = 10; nn = 8; hec_type = brick_elem_type( e );
				break;
			case CNFDB_404::top_Brick20:
				con = 11; nn =20; hec_type = brick_elem_type( e );
				break;
			default:
				//throw CConvMessage( CONV_NO_SUPPORTED_ELEMENT, "ELEMENT" );
                                fprintf(stderr, "##Warning: Non suported element type is found.\n");
                                continue;
			}

			option = get_hec_elem_option( e, hec_type );

		} catch( CConvMessage& emsg) {
			if( fg_element_error_stop ) throw emsg;
			fprintf( stdout, "%s --- ignored!!\n", emsg.Msg());
			fflush(stdout);
			ignored_count++;
			continue;
		}

		if( !hec_e
		|| hec_e->type   != hec_type
		|| hec_e->sec_id != sec_id
		|| hec_e->option != option ) {
			hec_e = search_element( hec, hec_type, sec_id, option );
			if(!hec_e) {
				hec_e = new CHECDB_Element();
				hec_e->type = hec_type;
				hec_e->option = option;
				hec_e->sec_id = sec_id;
				hec.DB.push_back(hec_e);
			}
		}

		CHECDB_Element::CElemItem eitem( hec_type, e->ID );

		for(int i=0; i<nn; i++) {
			eitem.node[i] = e->node[ con_table[con][i]];
		}
		hec_e->ElemList.insert( eitem );
	}

	fprintf( stdout, "%d/%d elements converted\n",
		elem_count-ignored_count, elem_count);
	fflush(stdout);
}


//-----------------------------------------------------------------------------
// MATERIAL
//-----------------------------------------------------------------------------


static
void set_material_static( CHECDB_Material* hec_m, CNFDB_601* m )
{
	create_mat_name( m->ID, hec_m->name );

	// ITEM 1 -- E & poisson
	{
		CHECDB_Material::CItem item;
		item.ID = 1;

		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->E(0) );
		rec.params.push_back( m->NU(0) );
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
	// ITEM 2 -- density
	{
		CHECDB_Material::CItem item;
		item.ID = 2;

		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->DENSITY() );
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
	// ITEM 3 -- thermal expansion
	{
		CHECDB_Material::CItem item;
		item.ID = 3;

		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->THERMAL_EXPANSION(0));
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
}


static
void set_material_heat( CHECDB_Material* hec_m, CNFDB_601* m )
{
	create_mat_name( m->ID, hec_m->name );

	// ITEM 1 - density
	{
		CHECDB_Material::CItem item;
		item.ID = 1;
		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->DENSITY() );
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
	// ITEM 2 -- specific heat
	{
		CHECDB_Material::CItem item;
		item.ID = 2;

		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->THERMAL_CAPACITY() );
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
	// ITEM 3 -- heat conductivity
	{
		CHECDB_Material::CItem item;
		item.ID = 3;

		CHECDB_Material::CItem::CItemRec rec;
		rec.params.push_back( m->THERMAL_CONDUCTIVITY(0) );
		item.RecList.push_back( rec );
		hec_m->ItemList.push_back( item );
	}
}



static
void SetMaterial( CNFData& neu, CHECData& hec )
{
	if( neu.DB_601.size() == 0 ) return;

	int mat_n = 0;
	vector<CNFDB_601*>::iterator iter;
	for(iter = neu.DB_601.begin(); iter != neu.DB_601.end(); iter++){
		CNFDB_601* m = *iter;
		CHECDB_Material* hec_m = new CHECDB_Material();
		switch( solution ){
		case sol_static:
		case sol_eigen:
			set_material_static( hec_m, m );
			break;
		case sol_heat:
			set_material_heat( hec_m, m );
			break;
		default:
			assert(0);
		}
		hec.DB.push_back(hec_m);
		mat_n++;
	}
	fprintf(stdout,"%d materials converted\n", mat_n);
	fflush(stdout);
}


//-----------------------------------------------------------------------------
// SECTION
// CAUTIION) Execute after SetElement
//-----------------------------------------------------------------------------


static
bool exist_elem_sec( CHECData& hec, int sec_id )
{
	vector<CHECDataBlock*>::iterator iter;
	for( iter = hec.DB.begin(); iter != hec.DB.end(); iter++) {
		if( (*iter)->data_type == HECDB_ELEMENT ){
			CHECDB_Element* e = (CHECDB_Element*)(*iter);
			if( e->sec_id == sec_id ) return true;
		}
	}
	return false;
}



enum {
	EESO_NOTHING = 0,
	EESO_SEC,
	EESO_SEC_OPT
};

static
int exist_elem_sec_opt( CHECData& hec, int sec_id, int option )
{
	bool fg = false;
	vector<CHECDataBlock*>::iterator iter;
	for( iter = hec.DB.begin(); iter != hec.DB.end(); iter++) {
		if( (*iter)->data_type == HECDB_ELEMENT ){
			CHECDB_Element* e = (CHECDB_Element*)(*iter);
			if( e->sec_id == sec_id ) {
				if( e->option == option ) return EESO_SEC_OPT;
				fg = true;
			}
		}
	}
	return fg ? EESO_SEC : EESO_NOTHING;
}


// execute after execution of SetElement
static
void SetSection( CNFData& neu, CHECData& hec )
{
	if( neu.DB_402.size() == 0 ) return;

	vector<CNFDB_402*>::iterator iter;
	for(iter = neu.DB_402.begin(); iter != neu.DB_402.end(); iter++){
		CNFDB_402* p = *iter;
		CHECDB_Section* sec = new CHECDB_Section();
		create_egrp_name_for_sec( p->ID, sec->egrp );
		create_mat_name( p->matID, sec->material );
		bool fg_nothing = false;
		switch( p->type ){
		case NEU_ELEM_PROP_PLANESTRAIN:
		case NEU_ELEM_PROP_PLANESTRAIN2:
			sec->type = CHECDB_Section::TYPE_SOLID;
			assert(!(p->num_val <= 0));
			sec->thickness = p->Value[0];
			switch( exist_elem_sec_opt( hec, p->ID, opt_plane_strain )){
			case EESO_SEC_OPT:
				sec->secopt = 1;
				break;
			case EESO_NOTHING:
				fg_nothing = true;
				break;
			}
			break;
		case NEU_ELEM_PROP_PLATE:
		case NEU_ELEM_PROP_PLATE2:
			if( !exist_elem_sec( hec, p->ID )) {
				fg_nothing = true;
				break;
			}
			sec->type = CHECDB_Section::TYPE_SHELL;
			assert(!(p->num_val <= 0));
			sec->thickness = p->Value[0];
			sec->integpoints = 3;
			break;
		default:
			sec->type = CHECDB_Section::TYPE_SOLID;
			if( !exist_elem_sec( hec, p->ID )) {
				fg_nothing = true;
			}
			break;
		}
		if( fg_nothing ){
			delete sec; sec=0;
			fprintf(stdout, "##Warning: Non used material-property ID:%d is eliminated\n", (int)p->ID);
			fflush(stdout);
		} else {
			hec.DB.push_back(sec);
		}
	}
}

//-----------------------------------------------------------------------------
// Generate Element Group for Section
//-----------------------------------------------------------------------------

class eg_rec_t {
public:
	int sec_id;
	CHECDB_EGroup* eg;
	eg_rec_t( int id=-1, CHECDB_EGroup* e=0):sec_id(id),eg(e) {}
	~eg_rec_t() {}
};

inline bool operator == (const eg_rec_t& a, const eg_rec_t& b ) { return a.sec_id == b.sec_id; }
inline bool operator <  (const eg_rec_t& a, const eg_rec_t& b ) { return a.sec_id <  b.sec_id; }
inline bool operator >  (const eg_rec_t& a, const eg_rec_t& b ) { return a.sec_id >  b.sec_id; }
inline bool operator <= (const eg_rec_t& a, const eg_rec_t& b ) { return a.sec_id <= b.sec_id; }
inline bool operator >= (const eg_rec_t& a, const eg_rec_t& b ) { return a.sec_id >= b.sec_id; }

static
void GenerateEGroupForSection( CHECData& hec )
{
	vector<eg_rec_t> eg_list;
	vector<eg_rec_t>::iterator egi;
	CHECDB_EGroup* eg;

	vector<CHECDataBlock*>::iterator iter;
	for(iter = hec.DB.begin(); iter != hec.DB.end(); iter++){
		if( (*iter)->data_type != HECDB_ELEMENT ) continue;
		CHECDB_Element* e = (CHECDB_Element*)*iter;
		if( e->ElemList.size() == 0 ) continue;
		for( egi = eg_list.begin(); egi != eg_list.end(); egi++) {
			if( egi->sec_id == e->sec_id ) break;
		}
		if(egi != eg_list.end()) {
			eg = egi->eg;
		} else {
			eg = new CHECDB_EGroup();
			create_egrp_name_for_sec( e->sec_id, eg->name );
			eg_list.push_back( eg_rec_t( e->sec_id, eg ));
		}
		set<CHECDB_Element::CElemItem>::iterator ei;
		for( ei=e->ElemList.begin(); ei!=e->ElemList.end(); ei++) {
			eg->ElemList.insert( ei->ID );
		}
	}

	for( egi=eg_list.begin(); egi != eg_list.end(); egi++) {
		hec.DB.push_back( egi->eg );
	}
}



//-----------------------------------------------------------------------------
// ZERO
//-----------------------------------------------------------------------------

static
void SetZero( CNFData& neu, CHECData& hec )
{
	if( neu.DB_507.size() == 0 ) return;

	bool fg_on = false;
	vector<CNFDB_507*>::iterator iter;
	for(iter = neu.DB_507.begin(); iter != neu.DB_507.end(); iter++){
		CNFDB_507* p = *iter;
		if( p->temp_on ){
			fg_on = true;
			CHECDB_Zero* zero = new CHECDB_Zero();
			zero->zero = p->Def_temp;
			hec.DB.push_back( zero );
		}
	}
	if( !fg_on ){
		CHECDB_Zero* zero = new CHECDB_Zero();
		zero->zero =  -273.16;
		hec.DB.push_back( zero );
	}
}

//=============================================================================
// conv_neu2hec
//=============================================================================


void conv_neu2hec( CNFData& neu, CHECData& hec, int sol )
{
	solution = sol;
	SetHeader  (neu, hec);
	SetNode    (neu, hec);
	SetElement (neu, hec);
	SetMaterial(neu, hec);
	SetSection (neu, hec);
	SetZero    (neu, hec);
	GenerateEGroupForSection( hec );
}

