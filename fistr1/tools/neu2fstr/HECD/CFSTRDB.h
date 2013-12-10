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
	CFSTRDB Ver.1.0
	--------------------------------------------
	Control Data of FRONT-STR( pSAN )
*/

#ifndef CFSTRDBH
#define CFSTRDBH

#include "CHECDataBlock.h"
#include "CHECDB.h"
#include <string.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef YES
#define YES 1
#define NO 0
#endif

enum {
	FSTRDB_SOLUTION = 100,
	FSTRDB_SOLVER,
	FSTRDB_WRITE,
	FSTRDB_ECHO,
	FSTRDB_STEP,

	FSTRDB_STATIC,
	FSTRDB_BOUNDARY,
	FSTRDB_CLOAD,
	FSTRDB_DLOAD,
	FSTRDB_TEMPERATURE,
	FSTRDB_REFTEMP,

	FSTRDB_EIGEN,

	FSTRDB_HEAT,
	FSTRDB_FIXTEMP,
	FSTRDB_CFLUX,
	FSTRDB_DFLUX,
	FSTRDB_SFLUX,
	FSTRDB_FILM,
	FSTRDB_SFILM,
	FSTRDB_RADIATE,
	FSTRDB_SRADIATE
};


CHECDataBlock* CreateFSTRDataBlock( const char* header_name );
bool IsFSTRDataBlockName( const char* name );


class CFSTRDataBlock : public CHECDataBlock {
public:
	CFSTRDataBlock(int dtype) : CHECDataBlock( dtype ) {}
	virtual bool IsMesh() { return false; }
};


//-----------------------------------------------
// For Common Used
//-----------------------------------------------


class CFSTRDB_Solution : public CFSTRDataBlock {
public:
	int type;
	enum {
		TYPE_UNKNOWN,
		TYPE_STATIC,
		TYPE_HEAT,
		TYPE_EIGEN
	};

	CFSTRDB_Solution();
	virtual ~CFSTRDB_Solution();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Solver : public CFSTRDataBlock {
public:
	// header line ---------------------

	char method[30];
	int precond;	// 1,2,3,10,11,12 or 21
	int nset;	// 0,-1 or 1
	int iterlog;	// 1:Yes, 0:No
	int timelog;	// 1:Yes, 0:No

	// 2nd line ------------------------

	int nier, iterPREmax, nrest;

	// 3rd line ------------------------

	double resid, fsigma_diag, sigma;

	// 4th line ------------------------

	double thresh, filter;

	CFSTRDB_Solver();
	virtual ~CFSTRDB_Solver();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Write : public CFSTRDataBlock {
public:
	int result;	// 1:Yes, 0:No
	int visual;	// 1:Yes, 0:No

	CFSTRDB_Write();
	virtual ~CFSTRDB_Write();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Echo : public CFSTRDataBlock {
public:
	int echo;

	CFSTRDB_Echo();
	virtual ~CFSTRDB_Echo();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Step : public CFSTRDataBlock {
public:
	int type;
	enum {
		TYPE_UNKNOWN,
		TYPE_STANDARD,
		TYPE_NLGEOM
	};
	int incmax;

	CFSTRDB_Step();
	virtual ~CFSTRDB_Step();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};

//-----------------------------------------------
// For Elastic Static Structual Analysis
//-----------------------------------------------

class CFSTRDB_Static : public CFSTRDataBlock {
public:
	// 2nd line -----------------
	double dtime, etime, itmax, eps;

	CFSTRDB_Static();
	virtual ~CFSTRDB_Static();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Boundary : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char ngrp[hec_name_size]; // ngrp name or node id
		int dof_ids;
		int dof_ide;
		double value;
		CItem(): dof_ids(0), dof_ide(0), value(0) {
			ngrp[0] = 0;
		}
		CItem( const char* name, int s, int e, double v=0 )
		 : dof_ids(s), dof_ide(e), value(v) {
			strcpy( ngrp, name );
		}
		CItem( const CItem& item )
		 : dof_ids(item.dof_ids), dof_ide(item.dof_ide), value(item.value) {
			strcpy( ngrp, item.ngrp );
		}
	};

	// data line -----------------

	std::vector<CItem> ItemList;

	CFSTRDB_Boundary();
	virtual ~CFSTRDB_Boundary();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_CLoad : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char ngrp[hec_name_size];
		int dof_id;
		double value;
		CItem() : dof_id(-1), value(0) { ngrp[0] = 0; }
		CItem(const char* name, int id, double v=0) : dof_id(id), value(v) {
			strcpy(ngrp,name);
		}
		CItem( const CItem& item ) : dof_id(item.dof_id), value(item.value) {
			strcpy( ngrp, item.ngrp );
		}
	};

	// data line -----------------

	std::vector<CItem> ItemList;

	CFSTRDB_CLoad();
	virtual ~CFSTRDB_CLoad();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_DLoad : public CFSTRDataBlock {
public:
	enum {
		TYPE_P0,
		TYPE_P1,
		TYPE_P2,
		TYPE_P3,
		TYPE_P4,
		TYPE_P5,
		TYPE_P6,
		TYPE_BX,
		TYPE_BY,
		TYPE_BZ,
		TYPE_GRAV,
		TYPE_CENT,
		TYPE_UNKNOWN
	};
	static int TypeNumber() { return 12; }
	static int ParamNumber( int type );
	static const char* LoadTypeName( int type );

	class CItem {
	public:
		char egrp[hec_name_size];
		int type;
		double param[7];
		CItem() : type( TYPE_UNKNOWN ) { egrp[0] = 0; init_param(); }
		CItem(const char* name, int t) : type(t) { strcpy(egrp,name); init_param();}
		CItem(const CItem& item) : type(item.type) {
			strcpy( egrp, item.egrp );
			param[0] = item.param[0];
			param[1] = item.param[1];
			param[2] = item.param[2];
			param[3] = item.param[3];
			param[4] = item.param[4];
			param[5] = item.param[5];
			param[6] = item.param[6];
		}
		void init_param() { param[0]=param[1]=param[2]=param[3]=param[4]=param[5]=param[6]=0.0; }
	};

	// data line -----------------

	std::vector<CItem> ItemList;

	CFSTRDB_DLoad();
	virtual ~CFSTRDB_DLoad();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Temperature : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char ngrp[hec_name_size];
		int value;
		CItem() : value(0) { ngrp[0] = 0; }
		CItem( const CItem& item ) : value(item.value) { strcpy( ngrp, item.ngrp ); }
		CItem( const char* name, int v=0 ) : value(v) { strcpy( ngrp, name ); }
	};

	// data line -----------------

	std::vector<CItem> ItemList;

	CFSTRDB_Temperature();
	virtual ~CFSTRDB_Temperature();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};



class CFSTRDB_Reftemp : public CFSTRDataBlock {
public:
	// data line -----------------

	double value;

	CFSTRDB_Reftemp();
	virtual ~CFSTRDB_Reftemp();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


//-----------------------------------------------
// For Eigen Analysis
//-----------------------------------------------


class CFSTRDB_Eigen : public CFSTRDataBlock {
public:
	// data line -----------------

	int nset;
	double lcztol;
	int lczmax;

	CFSTRDB_Eigen();
	virtual ~CFSTRDB_Eigen();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


//-----------------------------------------------
// For Heat Transfer Analysis
//-----------------------------------------------


class CFSTRDB_Heat : public CFSTRDataBlock {
public:
	// header line -----------------

	int restart; // TRUE/FALSE

	// 2nd line -------------------

	double dt, etime, dtime, deltmx;
	int itmax;
	double eps;

	CFSTRDB_Heat();
	virtual ~CFSTRDB_Heat();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Fixtemp : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char ngrp[hec_name_size];
		double value;
		CItem() : value(0) { ngrp[0] = 0; }
		CItem( const CItem& item ):value(item.value) { strcpy( ngrp, item.ngrp ); }
		CItem( const char* name, double v):value(v) { strcpy( ngrp, name ); }
	};

	// header line --------------------

	char amp[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_Fixtemp();
	virtual ~CFSTRDB_Fixtemp();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_CFlux : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char ngrp[hec_name_size];
		double value;
		CItem() : value(0) { ngrp[0] = 0; }
		CItem( const CItem& item ):value(item.value) { strcpy( ngrp, item.ngrp ); }
		CItem( const char* name, double v):value(v) { strcpy( ngrp, name ); }
	};

	// header line --------------------

	char amp[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_CFlux();
	virtual ~CFSTRDB_CFlux();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_DFlux : public CFSTRDataBlock {
public:
	enum {
		TYPE_S0,
		TYPE_S1,
		TYPE_S2,
		TYPE_S3,
		TYPE_S4,
		TYPE_S5,
		TYPE_S6,
		TYPE_BF,
		TYPE_UNKNOWN
	};

	static int TypeNumber() { return 8; }
	static const char* LoadTypeName( int type );

	class CItem {
	public:
		char egrp[hec_name_size];
		int type;
		double value;
		CItem() : type(TYPE_UNKNOWN), value(0) { egrp[0] = 0; }
		CItem( const CItem& item )
		:type(item.type),value(item.value) { strcpy( egrp, item.egrp ); }
		CItem( const char* name, int t, double v):type(t), value(v) { strcpy( egrp, name ); }
	};

	// header line --------------------

	char amp[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_DFlux();
	virtual ~CFSTRDB_DFlux();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_SFlux : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char sgrp[hec_name_size];
		double value;
		CItem() : value(0) { sgrp[0] = 0; }
		CItem( const CItem& item ):value(item.value) { strcpy( sgrp, item.sgrp ); }
		CItem( const char* name, double v): value(v) { strcpy( sgrp, name ); }
	};

	// header line --------------------

	char amp[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_SFlux();
	virtual ~CFSTRDB_SFlux();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Film : public CFSTRDataBlock {
public:
	enum {
		TYPE_F0,
		TYPE_F1,
		TYPE_F2,
		TYPE_F3,
		TYPE_F4,
		TYPE_F5,
		TYPE_F6,
		TYPE_UNKNOWN
	};

	static int TypeNumber() { return 7; }
	static const char* LoadTypeName( int type );

	class CItem {
	public:
		char egrp[hec_name_size];
		int type;
		double value;
		double sink;
		CItem() : type(TYPE_UNKNOWN), value(0), sink(0) { egrp[0] = 0; }
		CItem( const CItem& item )
		  :type(item.type),value(item.value),sink(item.sink) { strcpy( egrp, item.egrp ); }
		CItem( const char* name, int t, double v, double s)
		  :type(t), value(v), sink(s) { strcpy( egrp, name ); }
	};

	// header line --------------------

	char amp1[hec_name_size];
	char amp2[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_Film();
	virtual ~CFSTRDB_Film();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_SFilm : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char sgrp[hec_name_size];
		double value;
		double sink;
		CItem() : value(0),sink(0) { sgrp[0] = 0; }
		CItem( const CItem& item )
		  : value(item.value),sink(item.sink) { strcpy( sgrp, item.sgrp ); }
		CItem( const char* name, double v, double s)
		  : value(v), sink(s) { strcpy( sgrp, name ); }
	};

	// header line --------------------

	char amp1[hec_name_size];
	char amp2[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_SFilm();
	virtual ~CFSTRDB_SFilm();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_Radiate : public CFSTRDataBlock {
public:
	enum {
		TYPE_R0,
		TYPE_R1,
		TYPE_R2,
		TYPE_R3,
		TYPE_R4,
		TYPE_R5,
		TYPE_R6,
		TYPE_UNKNOWN
	};

	static int TypeNumber() { return 7; }
	static const char* LoadTypeName( int type );

	class CItem {
	public:
		char egrp[hec_name_size];
		int type;
		double value;
		double sink;
		CItem() : type(TYPE_UNKNOWN), value(0), sink(0) { egrp[0] = 0; }
		CItem( const CItem& item )
		  :type(item.type),value(item.value),sink(item.sink) { strcpy( egrp, item.egrp ); }
		CItem( const char* name, int t, double v, double s)
		  :type(t), value(v), sink(s) { strcpy( egrp, name ); }
	};

	// header line --------------------

	char amp1[hec_name_size];
	char amp2[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_Radiate();
	virtual ~CFSTRDB_Radiate();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CFSTRDB_SRadiate : public CFSTRDataBlock {
public:
	class CItem {
	public:
		char sgrp[hec_name_size];
		double value;
		double sink;
		CItem() : value(0),sink(0) { sgrp[0] = 0; }
		CItem( const CItem& item )
		  : value(item.value),sink(item.sink) { strcpy( sgrp, item.sgrp ); }
		CItem( const char* name, double v, double s)
		  : value(v), sink(s) { strcpy( sgrp, name ); }
	};

	// header line --------------------

	char amp1[hec_name_size];
	char amp2[hec_name_size];

	// data line ----------------------

	std::vector<CItem> ItemList;

	CFSTRDB_SRadiate();
	virtual ~CFSTRDB_SRadiate();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


#endif


