/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#ifndef INC_HECMW_COUPLE_DEFINE
#define INC_HECMW_COUPLE_DEFINE


#define HECMW_COUPLE_TRUE 1

#define HECMW_COUPLE_FALSE 0


#define HECMW_COUPLE_TYPE_UNDEF (-100)


#define HECMW_COUPLE_TYPE_MXN     101


#define HECMW_COUPLE_TYPE_MAXMN   102


#define HECMW_COUPLE_TYPE_MANUAL  103


#define HECMW_COUPLE_UNIT_UNDEF (-200)


#define HECMW_COUPLE_UNIT1        201


#define HECMW_COUPLE_UNIT2        202


#define HECMW_COUPLE_DIRECTION_UNDEF (-300)


#define HECMW_COUPLE_UNIT1_TO_UNIT2    301


#define HECMW_COUPLE_UNIT2_TO_UNIT1    302


#define HECMW_COUPLE_GROUP_UNDEF (-400)


#define HECMW_COUPLE_NODE_GROUP    401


#define HECMW_COUPLE_ELEMENT_GROUP 402


#define HECMW_COUPLE_SURFACE_GROUP 403


#define HECMW_COUPLE_IP_UNDEF      (-500)


#define HECMW_COUPLE_IP_NODE_TO_NODE 501


#define HECMW_COUPLE_IP_NODE_TO_ELEM 502


#define HECMW_COUPLE_IP_NODE_TO_SURF 503


#define HECMW_COUPLE_IP_ELEM_TO_NODE 511


#define HECMW_COUPLE_IP_ELEM_TO_ELEM 512


#define HECMW_COUPLE_IP_ELEM_TO_SURF 513


#define HECMW_COUPLE_IP_SURF_TO_NODE 521


#define HECMW_COUPLE_IP_SURF_TO_ELEM 522


#define HECMW_COUPLE_IP_SURF_TO_SURF 523


#define HECMW_COUPLE_MAP_UNDEF      (-600)


#define HECMW_COUPLE_MAP_NODE_TO_NODE 601


#define HECMW_COUPLE_MAP_NODE_TO_ELEM 602


#define HECMW_COUPLE_MAP_NODE_TO_SURF 603


#define HECMW_COUPLE_MAP_ELEM_TO_NODE 611


#define HECMW_COUPLE_MAP_ELEM_TO_ELEM 612


#define HECMW_COUPLE_MAP_ELEM_TO_SURF 613


#define HECMW_COUPLE_MAP_SURF_TO_NODE 621


#define HECMW_COUPLE_MAP_SURF_TO_ELEM 622


#define HECMW_COUPLE_MAP_SURF_TO_SURF 623


#define HECMW_COUPLE_TOLERANCE_DEFAULT (1.0E-04)


#define HECMW_COUPLE_BBCOEF_DEFAULT (1.05)


#define HECMW_COUPLE_BGCOEF_DEFAULT (1.05)

/*================================================================================================*/

#define HECMWCPL_E                      HECMW_COUPLE_E9999


#define HECMWCPL_E_INVALID_ARG          HECMW_COUPLE_E1001


#define HECMWCPL_W_INVALID_ARG          HECMW_COUPLE_W1001


#define HECMWCPL_E_INVALID_NULL_PTR     HECMW_COUPLE_E1002


#define HECMWCPL_E_INV_ETYPE            HECMW_COUPLE_E1101


#define HECMWCPL_E_UNSUP_ETYPE          HECMW_COUPLE_E1102


#define HECMWCPL_E_CTRL_INVALID_TOKEN   HECMW_COUPLE_E2001


#define HECMWCPL_E_CTRL_LONG_NAME       HECMW_COUPLE_E2002


#define HECMWCPL_E_CPLU                 HECMW_COUPLE_E2101


#define HECMWCPL_E_CPLU_NO_NAME         HECMW_COUPLE_E2111


#define HECMWCPL_E_CPLU_NO_NPROC        HECMW_COUPLE_E2112


#define HECMWCPL_E_CPLU_UNMATCH_RANKS   HECMW_COUPLE_E2113


#define HECMWCPL_E_CPL                  HECMW_COUPLE_E2201


#define HECMWCPL_E_CPL_NO_NAME          HECMW_COUPLE_E2211


#define HECMWCPL_E_CPL_NO_TYPE          HECMW_COUPLE_E2212


#define HECMWCPL_E_CPL_NO_UNIT1         HECMW_COUPLE_E2213


#define HECMWCPL_E_CPL_NO_UNIT2         HECMW_COUPLE_E2214


#define HECMWCPL_E_CPLB                 HECMW_COUPLE_E2301


#define HECMWCPL_E_CPLB_NO_NAME         HECMW_COUPLE_E2311


#define HECMWCPL_E_CPLB_NO_COUPLE       HECMW_COUPLE_E2312


#define HECMWCPL_E_CPLB_NO_DIRECTION    HECMW_COUPLE_E2313


#define HECMWCPL_E_CPLB_NO_UNIT1        HECMW_COUPLE_E2321


#define HECMWCPL_E_CPLB_NO_UNIT2        HECMW_COUPLE_E2322


#define HECMWCPL_E_CPLB_NO_GEOM         HECMW_COUPLE_E2323


#define HECMWCPL_E_CPLB_NO_DATA         HECMW_COUPLE_E2324


#define HECMWCPL_E_CPLB_NO_GRPNAME      HECMW_COUPLE_E2325


#define HECMWCPL_E_CPLB_UNMATCH_GRPTYPE HECMW_COUPLE_E2326


#define HECMWCPL_E_UNDEF_COUPLE_ID      HECMW_COUPLE_E3001


#define HECMWCPL_E_UNDEF_UNIT_ID        HECMW_COUPLE_E3002


#define HECMWCPL_E_UNDEF_BOUNDARY_ID    HECMW_COUPLE_E3003


#define HECMWCPL_E_UNDEF_GRPNAME        HECMW_COUPLE_E3011


#define HECMWCPL_E_INVALID_NPROC        HECMW_COUPLE_E3101


#define HECMWCPL_E_UNMATCH_PSIZE        HECMW_COUPLE_E3102


#define HECMWCPL_E_INVALID_RANKS        HECMW_COUPLE_E3103


#define HECMWCPL_E_DISCONTINUOUS_RANKS  HECMW_COUPLE_E3104


#define HECMWCPL_E_INVALID_CPLTYPE      HECMW_COUPLE_E3111


#define HECMWCPL_E_MULTIPLE_CPLTYPE     HECMW_COUPLE_E3112


#define HECMWCPL_E_INVALID_DIRECTION    HECMW_COUPLE_E3121


#define HECMWCPL_E_INVALID_GRPTYPE      HECMW_COUPLE_E3122


#define HECMWCPL_E_INVALID_GEOMTYPE     HECMW_COUPLE_E3122


#define HECMWCPL_E_INVALID_DATATYPE     HECMW_COUPLE_E3122


#define HECMWCPL_E_INVALID_UNITTYPE     HECMW_COUPLE_E3123


#define HECMWCPL_E_INVALID_MAPTYPE      HECMW_COUPLE_E3124


#define HECMWCPL_E_INVALID_IPTYPE       HECMW_COUPLE_E3124


#define HECMWCPL_E_NONSUPPORT_ETYPE     HECMW_COUPLE_E3201


#define HECMWCPL_E_NONSUPPORT_GEOMTYPE  HECMW_COUPLE_E3202


#define HECMWCPL_E_MPI                  HECMW_COUPLE_E8001


#define HECMWCPL_E_MPI_DATATYPE         HECMW_COUPLE_E8011


#define HECMWCPL_E_MPI_OP               HECMW_COUPLE_E8012

#endif	/* INC_HECMW_COUPLE_DEFINE */
