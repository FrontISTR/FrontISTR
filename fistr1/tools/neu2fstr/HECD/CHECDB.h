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
	CHECDB Ver.1.0
*/

#ifndef CHECDBH
#define CHECDBH

#include <set>
#include <vector>
#include "CHECDataBlock.h"
#include "CHECDB_Visual.h"

enum {
	HECDB_HEADER = 1,
	HECDB_NODE,
	HECDB_ELEMENT,
	HECDB_MATERIAL,
	HECDB_SECTION,
	HECDB_NGROUP,
	HECDB_EGROUP,
	HECDB_SGROUP,
	HECDB_AMPLITUDE,
	HECDB_ZERO,
	HECDB_VISUAL
};


class checdb_id_class {
public:
	int ID;
	checdb_id_class(int id=-1): ID(id) {}
	virtual ~checdb_id_class() {}
};


inline bool operator==( const checdb_id_class& a, const checdb_id_class& b) { return a.ID == b.ID; }
inline bool operator< ( const checdb_id_class& a, const checdb_id_class& b) { return a.ID <  b.ID; }
inline bool operator> ( const checdb_id_class& a, const checdb_id_class& b) { return a.ID >  b.ID; }


CHECDataBlock* CreateHECDataBlock( const char* header_name );
bool IsHECDataBlockName( const char* name );


class CHECDB_Header : public CHECDataBlock {
public:
	char title[hec_str_size];

	CHECDB_Header();
	virtual ~CHECDB_Header();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};



class CHECDB_Element : public CHECDataBlock {
public:
	static int NodeNumber(int type);
	static int FaceNumber(int type);
	static int DOFNumber(int type);
	static const int* Connectivity(int face_id );
	static bool CheckType( int type );
public:
	int type;
	// for generate SECTION
	int sec_id;
	int option;

	class CElemItem :public checdb_id_class {
	public:
		int* node;
		int node_n;
		CElemItem():checdb_id_class(),node(0), node_n(0) {}
		CElemItem(int type, int id=-1 )
		 :checdb_id_class(id),node(0), node_n(0) { Init(type);}
		CElemItem( const CElemItem& e)
		  :checdb_id_class(e.ID),node(0),node_n(0){
			InitNode( e.node_n );
			for(int i=0; i<node_n; i++) node[i] = e.node[i];
		}
		virtual ~CElemItem() { delete[] node; }
		void Init(int type ){
			InitNode( NodeNumber(type));
		}
		void InitNode(int n ){
			delete[] node;
			node_n = n;
			node = new int[node_n];
		}
	};
	std::set<CElemItem> ElemList;

	CHECDB_Element();
	virtual ~CHECDB_Element();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );

	int ElemNumber() { return ElemList.size(); }
	CElemItem* GetElem( int id );
};



class CHECDB_Node : public CHECDataBlock {
public:
	class CNodeItem:public checdb_id_class {
	public:
		// int ID;
		double x, y, z;
		CNodeItem( int id=-1, double X=0, double Y=0, double Z=0 )
		 : checdb_id_class(id), x(X), y(Y), z(Z) {}
	};
	std::set<CNodeItem> NodeList;

	CHECDB_Node();
	virtual ~CHECDB_Node();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );

	int NodeNumber() { return NodeList.size(); }
	CNodeItem* GetNode( int id );
};



class CHECDB_Material : public CHECDataBlock {
public:
	char name[hec_name_size];

	class CItem {
	public:
		int ID;
		CItem() : ID(-1) {}
		class CItemRec {
		public:
			bool last_is_temp;
			std::vector<double> params;
			CItemRec(): last_is_temp(false), params(){}
		};
		std::vector<CItemRec> RecList;
		int RecNumber() { return RecList.size(); }
		int SubItemNumber() {
			if( RecNumber()>0 ) {
				std::vector<CItemRec>::iterator i = RecList.begin();
				if( i->last_is_temp )
					return i->params.size()-1;
				else
					return i->params.size();
			} else
				return 0;
		}
		bool LastIsTemp() {
			if( RecNumber()>0 ) {
				std::vector<CItemRec>::iterator i = RecList.begin();
				return i->last_is_temp;
			} else
				return false;
		}
		void Clear() { RecList.clear(); }
		void Write( class CHECData* hecd );
		bool Read( class CHECData* hecd );
	};
	std::vector<CItem> ItemList;

	CHECDB_Material();
	virtual ~CHECDB_Material();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};



class CHECDB_Section : public CHECDataBlock {
public:
	// type
	enum {
		TYPE_UNKNOWN = 0,
		TYPE_SOLID,
		TYPE_SHELL,
		TYPE_BEAM,
		TYPE_INTERFACE
	};

	int type;
	char egrp[hec_name_size];
	char material[hec_name_size];
	int n_comp;
	int secopt;

	// type == TYPE_SOLID, TYPE_SHELL or TYPE_INTERFACE
	double thickness;

	// type == TYPE_SHELL
	int integpoints;

	// type == TYPE_INTERFACE
	double gapcon;
	double gaprad1;
	double gaprad2;

	CHECDB_Section();
	virtual ~CHECDB_Section();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CHECDB_NGroup : public CHECDataBlock {
public:
	char name[hec_name_size];
	std::set<int> NodeList;

	CHECDB_NGroup();
	virtual ~CHECDB_NGroup();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CHECDB_EGroup : public CHECDataBlock {
public:
	char name[hec_name_size];
	std::set<int> ElemList;

	CHECDB_EGroup();
	virtual ~CHECDB_EGroup();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CHECDB_SGroup : public CHECDataBlock {
public:
	char name[hec_name_size];

	class CItem {
	public:
		int elem, surf;
		CItem( int eid=-1, int sid=0 ) : elem(eid), surf(sid) {}
	};
	std::vector<CItem> ItemList;

	CHECDB_SGroup();
	virtual ~CHECDB_SGroup();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CHECDB_Amplitude : public CHECDataBlock {
public:
	char name[hec_name_size];
	char definition[hec_name_size];
	char time[hec_name_size];
	char value[hec_name_size];

	class CItem {
	public:
		double val, t;
		CItem( double v=0, double T=0):val(v),t(T) {}
	};
	std::vector<CItem> ItemList;

	CHECDB_Amplitude();
	virtual ~CHECDB_Amplitude();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};


class CHECDB_Zero : public CHECDataBlock {
public:
	double zero;

	CHECDB_Zero();
	virtual ~CHECDB_Zero();
	virtual void Clear();
	virtual void Write( class CHECData* hecd );
	virtual bool Read( class CHECData* hecd, char* header_line );
};




#endif
