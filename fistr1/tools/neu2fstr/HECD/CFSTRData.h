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
	CFSTRData Ver.1.0
	----------------------------------------
	[method for loading mesh and control file]
	1) Clear()
	2) AddLoad( mesh_file )
	3) AddLoad( ctrl_file )
*/

#ifndef CFSTRDataH
#define CFSTRDataH

#include <stdio.h>
#include "CHECData.h"
#include "CHECDB.h"
#include "CFSTRDB.h"


class CFSTRData : public CHECData {
public:
	CFSTRData();
	virtual bool SaveMesh( const char* file_name, const char* comment="");
	virtual bool SaveCtrl( const char* file_name, const char* comment="");
	virtual CHECDataBlock* CreateDataBlock( const char* header_name );
	virtual bool IsDataBlockName( const char* name ) {
		return IsHECDataBlockName(name) || IsFSTRDataBlockName(name);
	}
protected:
	virtual void WriteComment( FILE* fp, const char* comment);
};



#endif
