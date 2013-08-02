/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.1                                               *
 *                                                                     *
 *     Last Update : 2006/09/27                                        *
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
	CNFDataBlock Ver.1.0
*/

#ifndef CNFDataBlockH
#define CNFDataBlockH

#include "CNFMessage.h"


#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


typedef int           nf_int;   // 4 byte
typedef unsigned char nf_bool;  // 1 byte
typedef double        nf_float; // 8 byte
typedef char          nf_char;  // 1 byte


class CNFDataBlock {
public:
	int DataBlockID;

	CNFDataBlock(int id) : DataBlockID(id) {}
	virtual ~CNFDataBlock() {}

	int Type() { return DataBlockID;}

	// read block data from file without DataBlockID member
	virtual void Read( class CNFData* nfd ) = 0;
	virtual void WriteData( class CNFData* nfd, FILE* fp=0 ) = 0;
};


// type of element property (not element type)

#define NEU_ELEM_PROP_UNKNOWN	0
#define NEU_ELEM_PROP_ROD	1	//used
#define NEU_ELEM_PROP_BAR	2	//used
#define NEU_ELEM_PROP_TUBE	3
#define NEU_ELEM_PROP_LINK	4	//used
#define NEU_ELEM_PROP_BEAM	5	//used
#define NEU_ELEM_PROP_BEAM2	37	//used
#define NEU_ELEM_PROP_SPRING	6
#define NEU_ELEM_PROP_DOFSPRING	7
#define NEU_ELEM_PROP_CURVEBEAM	8	//used
#define NEU_ELEM_PROP_GAP	9
#define NEU_ELEM_PROP_SHEAR	11
#define NEU_ELEM_PROP_SHEAR2	12
#define NEU_ELEM_PROP_MEMBRANE	13	//used
#define NEU_ELEM_PROP_MEMBRANE2	14	//used
#define NEU_ELEM_PROP_BENDING	15	//used
#define NEU_ELEM_PROP_BENDING2	16	//used
#define NEU_ELEM_PROP_PLATE	17	//used
#define NEU_ELEM_PROP_PLATE2	18	//used
#define NEU_ELEM_PROP_PLANESTRAIN 19	//used
#define NEU_ELEM_PROP_PLANESTRAIN2 20	//used
#define NEU_ELEM_PROP_LAMINATE	21	//used
#define NEU_ELEM_PROP_LAMINATE2	22	//used
#define NEU_ELEM_PROP_SOLID	25	//used
#define NEU_ELEM_PROP_SOLID2	26	//used
#define NEU_ELEM_PROP_MASS	27
#define NEU_ELEM_PROP_MASSMAT	28
#define NEU_ELEM_PROP_STIFFMAT	30
#define NEU_ELEM_PROP_SYMSHELL	35	//used
#define NEU_ELEM_PROP_SYMSHELL2	36	//used


#endif









