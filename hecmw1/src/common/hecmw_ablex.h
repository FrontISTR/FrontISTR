/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/12/03                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef HECMW_ABLEX_INCLUDED
#define HECMW_ABLEX_INCLUDED

#include <stdio.h>


enum {
	HECMW_ABLEX_NL = 1000,			
	HECMW_ABLEX_INT,				
	HECMW_ABLEX_DOUBLE,				
	HECMW_ABLEX_NAME,				
	HECMW_ABLEX_FILENAME,			
	HECMW_ABLEX_HEADER,				

	HECMW_ABLEX_H_AMPLITUDE = 2000,	
	HECMW_ABLEX_H_CONDUCTIVITY,		
	HECMW_ABLEX_H_DENSITY,			
	HECMW_ABLEX_H_ELASTIC,			
	HECMW_ABLEX_H_ELEMENT,			
	HECMW_ABLEX_H_ELSET,			
	HECMW_ABLEX_H_EQUATION,			
	HECMW_ABLEX_H_HEADING,			
	HECMW_ABLEX_H_INCLUDE,			
	HECMW_ABLEX_H_INITIAL,			
	HECMW_ABLEX_H_MATERIAL,			
	HECMW_ABLEX_H_NODE,				
	HECMW_ABLEX_H_NSET,				
	HECMW_ABLEX_H_SHELL_SECTION,	
	HECMW_ABLEX_H_SOLID_SECTION,	
	HECMW_ABLEX_H_SPECIFIC_HEAT,	
	HECMW_ABLEX_H_SYSTEM,			

	HECMW_ABLEX_K_ABSOLUTE = 3000,	
	HECMW_ABLEX_K_ANISOTROPIC,		
	HECMW_ABLEX_K_ELSET,			
	HECMW_ABLEX_K_ENGINEERING_CONSTANTS,	
	HECMW_ABLEX_K_DEFINITION,		
	HECMW_ABLEX_K_DEPENDENCIES,		
	HECMW_ABLEX_K_GENERATE,			
	HECMW_ABLEX_K_INPUT,			
	HECMW_ABLEX_K_ISOTROPIC,		
	HECMW_ABLEX_K_LAMINA,			
	HECMW_ABLEX_K_MATERIAL,			
	HECMW_ABLEX_K_NAME,				
	HECMW_ABLEX_K_NSET,				
	HECMW_ABLEX_K_ORTHOTROPIC,		
	HECMW_ABLEX_K_RELATIVE,			
	HECMW_ABLEX_K_STEP_TIME,		
	HECMW_ABLEX_K_SYSTEM,			
	HECMW_ABLEX_K_TABULAR,			
	HECMW_ABLEX_K_TEMPERATURE,		
	HECMW_ABLEX_K_TIME,				
	HECMW_ABLEX_K_TYPE,				
	HECMW_ABLEX_K_UNSORTED,			
	HECMW_ABLEX_K_VALUE,			

	/* element type */
	HECMW_ABLEX_E_B31 = 4000,		
	HECMW_ABLEX_E_B32,				
	HECMW_ABLEX_E_C3D4,				
	HECMW_ABLEX_E_C3D6,				
	HECMW_ABLEX_E_C3D8,				
	HECMW_ABLEX_E_C3D8I,			
	HECMW_ABLEX_E_C3D10,			
	HECMW_ABLEX_E_C3D15,			
	HECMW_ABLEX_E_C3D20,			
	HECMW_ABLEX_E_CAX3,				
	HECMW_ABLEX_E_CAX4,				
	HECMW_ABLEX_E_CAX4I,			
	HECMW_ABLEX_E_CAX4R,			
	HECMW_ABLEX_E_CAX6,				
	HECMW_ABLEX_E_CAX8,				
	HECMW_ABLEX_E_CAX8R,			
	HECMW_ABLEX_E_CPE3,				
	HECMW_ABLEX_E_CPE4,				
	HECMW_ABLEX_E_CPE4I,			
	HECMW_ABLEX_E_CPE4R,			
	HECMW_ABLEX_E_CPE6,				
	HECMW_ABLEX_E_CPE8,				
	HECMW_ABLEX_E_CPE8R,			
	HECMW_ABLEX_E_CPS3,				
	HECMW_ABLEX_E_CPS4,				
	HECMW_ABLEX_E_CPS4I,			
	HECMW_ABLEX_E_CPS4R,			
	HECMW_ABLEX_E_CPS6,				
	HECMW_ABLEX_E_CPS8,				
	HECMW_ABLEX_E_CPS8R,			
	HECMW_ABLEX_E_DC1D2,			
	HECMW_ABLEX_E_DC1D3,			
	HECMW_ABLEX_E_DC2D3,			
	HECMW_ABLEX_E_DC2D4,			
	HECMW_ABLEX_E_DC2D6,			
	HECMW_ABLEX_E_DC2D8,			
	HECMW_ABLEX_E_DC3D4,			
	HECMW_ABLEX_E_DC3D6,			
	HECMW_ABLEX_E_DC3D8,			
	HECMW_ABLEX_E_DC3D10,			
	HECMW_ABLEX_E_DC3D15,			
	HECMW_ABLEX_E_DC3D20,			
	HECMW_ABLEX_E_DCAX3,			
	HECMW_ABLEX_E_DCAX4,			
	HECMW_ABLEX_E_DCAX6,			
	HECMW_ABLEX_E_DCAX8,			
	HECMW_ABLEX_E_DINTER4,			
	HECMW_ABLEX_E_DINTER8,			
	HECMW_ABLEX_E_DS4,				
	HECMW_ABLEX_E_DS8,				
	HECMW_ABLEX_E_INTER4,			
	HECMW_ABLEX_E_INTER8,			
	HECMW_ABLEX_E_S3R,				
	HECMW_ABLEX_E_S4R,				
	HECMW_ABLEX_E_S8R,				
	HECMW_ABLEX_E_T3D2,				
	HECMW_ABLEX_E_T3D3				
};


extern double HECMW_ablex_get_number(void);


extern char *HECMW_ablex_get_text(void);


extern int HECMW_ablex_get_lineno(void);


extern int HECMW_ablex_next_token(void);


extern int HECMW_ablex_next_token_skip(int skip_token);


extern int HECMW_ablex_set_input(FILE *fp);


extern int HECMW_ablex_skip_line(void);


extern int HECMW_ablex_switch_to_include(const char *filename);


extern int HECMW_ablex_unput_token(void);


extern int HECMW_ablex_is_including(void);

#endif
