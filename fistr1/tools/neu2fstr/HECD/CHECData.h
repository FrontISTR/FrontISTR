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
	CHECData Ver.1.0
*/


#ifndef CHECDataH
#define CHECDataH


#include <stdio.h>
#include <vector>
#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include "hecd_util.h"
#include "CHECDataBlock.h"

const int mw_fname_size = 256;
const int mw_comment_size = 256;

#include "CHECDB.h"


class CHECData {
public:
	FILE* fp;
	char fname[mw_fname_size];
	int line_count;
	std::vector<CHECDataBlock*> DB;

	CHECData();
	virtual ~CHECData();
	virtual void Clear();
	virtual void StoreDataBlock( CHECDataBlock* block );

	// ============ etc. ====================
	virtual bool IsDataBlockName( const char* name ) {
		return IsHECDataBlockName(name);
	}

	// ============ Utilities for Save ================

	virtual bool Save( const char* file_name );
	virtual void WriteLine( const char* s );

	// Header Line -----------------------------------------------

	// fmt : format of parameters.
	//       Each character specify the parameter's type
	//         'I': integer, 'F':float, 'S':string
	// ... : piars of parameter name and value
	virtual void WriteHeader( const char* name, const char* fmt="", ... );

	// Parameter Line -----------------------------------------------
	virtual  void WriteParameter(const char* fmt="", ... );

	// Data Line ---------------------------------------------------

	// fmt : formt of data like one of WriteHeader
	// ... : values to write the file.
	virtual void WriteData( const char* fmt, ...);

	virtual void ClearDataLineBuffer();
	virtual void AddDataLineItems( const char* fmt, ...);
	virtual void WriteDataLine();

	// ============ Utilities for Load ================

	virtual bool Load( const char* file_name );
	virtual bool AddLoad( const char* file_name ); // append to DB (not cleared)
	virtual bool ReadLine( char* s, int size=256 );
	virtual CHECDataBlock* CreateDataBlock( const char* header_name );
	virtual bool GetHeaderName( const char* header_line, char* header_name );
	virtual void PushReadLine( const char* s );

	// rcode[i] : 1 -- set, 0 -- not set, -1 -- error
	// fmt : string composed by 'I'(int), 'F'(double), 'E'(int) or 'S'(char*) characters
	// ... : pairs of a parameter name and a pointer of the parameter
	virtual bool ParseHeader( char* header_line, int* rcode, const char* fmt, ... );
	virtual bool ReadParameter( int* rcode, const char* fmt, ... );
	virtual bool vParseParameter( char* line, int* rcode, const char* fmt, va_list va );
	virtual bool ParseParameter( char* line, int* rcode, const char* fmt, ... );
	// fmt : string composed by 'I', 'F' or 'S'
	// ... : pointers of parameter
	virtual bool ReadData( int* rcode, const char* fmt, ... );
	// return : num. of data or  -(error position+1)
	virtual int ParseDoubleDataArray( char* line, double* data );
	virtual int ParseIntDataArray( char* line, int* data );

	// ============ Utilities for DataBlock ============
	virtual class CHECDB_Material* GetMaterial( const char* name );
	virtual class CHECDB_NGroup* GetNGroup( const char* name );
	virtual class CHECDB_EGroup* GetEGroup( const char* name );
	virtual class CHECDB_SGroup* GetSGroup( const char* name );

	virtual class CHECDB_Node::CNodeItem* GetNodeItem( int id );
	virtual class CHECDB_Element::CElemItem* GetElemItem( int id );
	virtual int GetElemType( int id ); // 0 : not existed

protected:
	char data_line_buffer[256];
	char header_line_buffer[256];
	bool fg_header_pushed;
};



#endif
