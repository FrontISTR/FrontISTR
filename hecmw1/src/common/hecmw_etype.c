/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdlib.h>
#include "hecmw_common_define.h"
#include "hecmw_etype.h"

extern int HECMW_get_etype_UTIL2HECMW(int etype) {
  switch (etype) {
    case HECMW_MESH_ETYPE_ROD1:
      return HECMW_ETYPE_ROD1;
    case HECMW_MESH_ETYPE_ROD2:
      return HECMW_ETYPE_ROD2;
    case HECMW_MESH_ETYPE_TRI1:
      return HECMW_ETYPE_TRI1;
    case HECMW_MESH_ETYPE_TRI2:
      return HECMW_ETYPE_TRI2;
    case HECMW_MESH_ETYPE_TRI22:
      return HECMW_ETYPE_TRI22;
    case HECMW_MESH_ETYPE_QUA1:
      return HECMW_ETYPE_QUA1;
    case HECMW_MESH_ETYPE_QUA2:
      return HECMW_ETYPE_QUA2;
    case HECMW_MESH_ETYPE_ROD31:
      return HECMW_ETYPE_ROD31;
    case HECMW_MESH_ETYPE_TET1:
      return HECMW_ETYPE_TET1;
    case HECMW_MESH_ETYPE_TET1_4:
      return HECMW_ETYPE_TET1_4;
    case HECMW_MESH_ETYPE_TET2:
      return HECMW_ETYPE_TET2;
    case HECMW_MESH_ETYPE_TET22:
      return HECMW_ETYPE_TET22;
    case HECMW_MESH_ETYPE_PRI1:
      return HECMW_ETYPE_PRI1;
    case HECMW_MESH_ETYPE_PRI2:
      return HECMW_ETYPE_PRI2;
    case HECMW_MESH_ETYPE_HEX1:
      return HECMW_ETYPE_HEX1;
    case HECMW_MESH_ETYPE_HEX1_4:
      return HECMW_ETYPE_HEX1_4;
    case HECMW_MESH_ETYPE_HEX2:
      return HECMW_ETYPE_HEX2;
    case HECMW_MESH_ETYPE_PYR1:
      return HECMW_ETYPE_PYR1;
    case HECMW_MESH_ETYPE_PYR2:
      return HECMW_ETYPE_PYR2;
    case HECMW_MESH_ETYPE_MST1:
      return HECMW_ETYPE_MST1;
    case HECMW_MESH_ETYPE_MST2:
      return HECMW_ETYPE_MST2;
    case HECMW_MESH_ETYPE_MSQ1:
      return HECMW_ETYPE_MSQ1;
    case HECMW_MESH_ETYPE_MSQ2:
      return HECMW_ETYPE_MSQ2;
    case HECMW_MESH_ETYPE_JTB1:
      return HECMW_ETYPE_JTB1;
    case HECMW_MESH_ETYPE_SPGDPT1:
        return HECMW_ETYPE_SPGDPT1;
    case HECMW_MESH_ETYPE_JTT1:
      return HECMW_ETYPE_JTT1;
    case HECMW_MESH_ETYPE_JTT2:
      return HECMW_ETYPE_JTT2;
    case HECMW_MESH_ETYPE_JTQ1:
      return HECMW_ETYPE_JTQ1;
    case HECMW_MESH_ETYPE_JTQ2:
      return HECMW_ETYPE_JTQ2;
    case HECMW_MESH_ETYPE_BEM1:
      return HECMW_ETYPE_BEM1;
    case HECMW_MESH_ETYPE_BEM2:
      return HECMW_ETYPE_BEM2;
    case HECMW_MESH_ETYPE_BEM3:
      return HECMW_ETYPE_BEM3; /* mixed beam-341*/
    case HECMW_MESH_ETYPE_SHT1:
      return HECMW_ETYPE_SHT1;
    case HECMW_MESH_ETYPE_SHT2:
      return HECMW_ETYPE_SHT2;
    case HECMW_MESH_ETYPE_SHQ1:
      return HECMW_ETYPE_SHQ1;
    case HECMW_MESH_ETYPE_SHQ2:
      return HECMW_ETYPE_SHQ2;
    case HECMW_MESH_ETYPE_SHQ3:
      return HECMW_ETYPE_SHQ3;
    case HECMW_MESH_ETYPE_SHT6:
      return HECMW_ETYPE_SHT6; /* mixed shell-solid */
    case HECMW_MESH_ETYPE_SHQ8:
      return HECMW_ETYPE_SHQ8;
    case HECMW_MESH_ETYPE_LN11:
      return HECMW_ETYPE_LN11;
    case HECMW_MESH_ETYPE_LN12:
      return HECMW_ETYPE_LN12;
    case HECMW_MESH_ETYPE_LN13:
      return HECMW_ETYPE_LN13;
    case HECMW_MESH_ETYPE_LN14:
      return HECMW_ETYPE_LN14;
    case HECMW_MESH_ETYPE_LN15:
      return HECMW_ETYPE_LN15;
    case HECMW_MESH_ETYPE_LN16:
      return HECMW_ETYPE_LN16;
    case HECMW_MESH_ETYPE_LN21:
      return HECMW_ETYPE_LN21;
    case HECMW_MESH_ETYPE_LN22:
      return HECMW_ETYPE_LN22;
    case HECMW_MESH_ETYPE_LN23:
      return HECMW_ETYPE_LN23;
    case HECMW_MESH_ETYPE_LN24:
      return HECMW_ETYPE_LN24;
    case HECMW_MESH_ETYPE_LN25:
      return HECMW_ETYPE_LN25;
    case HECMW_MESH_ETYPE_LN26:
      return HECMW_ETYPE_LN26;
    case HECMW_MESH_ETYPE_LN31:
      return HECMW_ETYPE_LN31;
    case HECMW_MESH_ETYPE_LN32:
      return HECMW_ETYPE_LN32;
    case HECMW_MESH_ETYPE_LN33:
      return HECMW_ETYPE_LN33;
    case HECMW_MESH_ETYPE_LN34:
      return HECMW_ETYPE_LN34;
    case HECMW_MESH_ETYPE_LN35:
      return HECMW_ETYPE_LN35;
    case HECMW_MESH_ETYPE_LN36:
      return HECMW_ETYPE_LN36;
    case HECMW_MESH_ETYPE_LN41:
      return HECMW_ETYPE_LN41;
    case HECMW_MESH_ETYPE_LN42:
      return HECMW_ETYPE_LN41;
    case HECMW_MESH_ETYPE_LN43:
      return HECMW_ETYPE_LN43;
    case HECMW_MESH_ETYPE_LN44:
      return HECMW_ETYPE_LN44;
    case HECMW_MESH_ETYPE_LN45:
      return HECMW_ETYPE_LN45;
    case HECMW_MESH_ETYPE_LN46:
      return HECMW_ETYPE_LN46;
    case HECMW_MESH_ETYPE_LN51:
      return HECMW_ETYPE_LN51;
    case HECMW_MESH_ETYPE_LN52:
      return HECMW_ETYPE_LN52;
    case HECMW_MESH_ETYPE_LN53:
      return HECMW_ETYPE_LN53;
    case HECMW_MESH_ETYPE_LN54:
      return HECMW_ETYPE_LN54;
    case HECMW_MESH_ETYPE_LN55:
      return HECMW_ETYPE_LN55;
    case HECMW_MESH_ETYPE_LN56:
      return HECMW_ETYPE_LN56;
    case HECMW_MESH_ETYPE_LN61:
      return HECMW_ETYPE_LN61;
    case HECMW_MESH_ETYPE_LN62:
      return HECMW_ETYPE_LN62;
    case HECMW_MESH_ETYPE_LN63:
      return HECMW_ETYPE_LN63;
    case HECMW_MESH_ETYPE_LN64:
      return HECMW_ETYPE_LN64;
    case HECMW_MESH_ETYPE_LN65:
      return HECMW_ETYPE_LN65;
    case HECMW_MESH_ETYPE_LN66:
      return HECMW_ETYPE_LN66;
    case HECMW_MESH_ETYPE_PTT1:
      return HECMW_ETYPE_PTT1;
    case HECMW_MESH_ETYPE_PTT2:
      return HECMW_ETYPE_PTT2;
    case HECMW_MESH_ETYPE_PTQ1:
      return HECMW_ETYPE_PTQ1;
    case HECMW_MESH_ETYPE_PTQ2:
      return HECMW_ETYPE_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_etype_HECMW2UTIL(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MESH_ETYPE_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MESH_ETYPE_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MESH_ETYPE_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MESH_ETYPE_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MESH_ETYPE_TRI22;
    case HECMW_ETYPE_QUA1:
      return HECMW_MESH_ETYPE_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MESH_ETYPE_QUA2;
    case HECMW_ETYPE_ROD31:
      return HECMW_MESH_ETYPE_ROD31;
    case HECMW_ETYPE_TET1:
      return HECMW_MESH_ETYPE_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MESH_ETYPE_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MESH_ETYPE_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_MESH_ETYPE_TET22;
    case HECMW_ETYPE_PRI1:
      return HECMW_MESH_ETYPE_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MESH_ETYPE_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MESH_ETYPE_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MESH_ETYPE_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MESH_ETYPE_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MESH_ETYPE_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MESH_ETYPE_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MESH_ETYPE_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MESH_ETYPE_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MESH_ETYPE_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MESH_ETYPE_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MESH_ETYPE_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MESH_ETYPE_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MESH_ETYPE_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MESH_ETYPE_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MESH_ETYPE_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MESH_ETYPE_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MESH_ETYPE_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MESH_ETYPE_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MESH_ETYPE_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MESH_ETYPE_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MESH_ETYPE_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MESH_ETYPE_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MESH_ETYPE_SHQ2;
    case HECMW_ETYPE_SHQ3:
      return HECMW_MESH_ETYPE_SHQ3;
    case HECMW_ETYPE_SHT6:
      return HECMW_MESH_ETYPE_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MESH_ETYPE_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MESH_ETYPE_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MESH_ETYPE_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MESH_ETYPE_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MESH_ETYPE_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MESH_ETYPE_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MESH_ETYPE_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MESH_ETYPE_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MESH_ETYPE_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MESH_ETYPE_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MESH_ETYPE_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MESH_ETYPE_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MESH_ETYPE_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MESH_ETYPE_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MESH_ETYPE_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MESH_ETYPE_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MESH_ETYPE_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MESH_ETYPE_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MESH_ETYPE_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MESH_ETYPE_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MESH_ETYPE_LN41;
    case HECMW_ETYPE_LN43:
      return HECMW_MESH_ETYPE_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MESH_ETYPE_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MESH_ETYPE_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MESH_ETYPE_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MESH_ETYPE_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MESH_ETYPE_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MESH_ETYPE_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MESH_ETYPE_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MESH_ETYPE_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MESH_ETYPE_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MESH_ETYPE_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MESH_ETYPE_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MESH_ETYPE_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MESH_ETYPE_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MESH_ETYPE_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MESH_ETYPE_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MESH_ETYPE_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MESH_ETYPE_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MESH_ETYPE_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MESH_ETYPE_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_etype_GeoFEM2HECMW(int etype) {
  switch (etype) {
    case HECMW_GEOFEM_ETYPE_ROD1:
      return HECMW_ETYPE_ROD1;
    case HECMW_GEOFEM_ETYPE_ROD2:
      return HECMW_ETYPE_ROD2;
    case HECMW_GEOFEM_ETYPE_TRI1:
      return HECMW_ETYPE_TRI1;
    case HECMW_GEOFEM_ETYPE_TRI2:
      return HECMW_ETYPE_TRI2;
    case HECMW_GEOFEM_ETYPE_QUA1:
      return HECMW_ETYPE_QUA1;
    case HECMW_GEOFEM_ETYPE_QUA2:
      return HECMW_ETYPE_QUA2;
    case HECMW_GEOFEM_ETYPE_TET1:
      return HECMW_ETYPE_TET1;
    case HECMW_GEOFEM_ETYPE_TET1_4:
      return HECMW_ETYPE_TET1_4;
    case HECMW_GEOFEM_ETYPE_TET2:
      return HECMW_ETYPE_TET2;
    case HECMW_GEOFEM_ETYPE_PRI1:
      return HECMW_ETYPE_PRI1;
    case HECMW_GEOFEM_ETYPE_PRI2:
      return HECMW_ETYPE_PRI2;
    case HECMW_GEOFEM_ETYPE_HEX1:
      return HECMW_ETYPE_HEX1;
    case HECMW_GEOFEM_ETYPE_HEX1_4:
      return HECMW_ETYPE_HEX1_4;
    case HECMW_GEOFEM_ETYPE_HEX2:
      return HECMW_ETYPE_HEX2;
    case HECMW_GEOFEM_ETYPE_MST1:
      return HECMW_ETYPE_MST1;
    case HECMW_GEOFEM_ETYPE_MST2:
      return HECMW_ETYPE_MST2;
    case HECMW_GEOFEM_ETYPE_MSQ1:
      return HECMW_ETYPE_MSQ1;
    case HECMW_GEOFEM_ETYPE_MSQ2:
      return HECMW_ETYPE_MSQ2;
    case HECMW_GEOFEM_ETYPE_JTB1:
      return HECMW_ETYPE_JTB1;
    case HECMW_GEOFEM_ETYPE_JTT1:
      return HECMW_ETYPE_JTT1;
    case HECMW_GEOFEM_ETYPE_JTT2:
      return HECMW_ETYPE_JTT2;
    case HECMW_GEOFEM_ETYPE_JTQ1:
      return HECMW_ETYPE_JTQ1;
    case HECMW_GEOFEM_ETYPE_JTQ2:
      return HECMW_ETYPE_JTQ2;
    case HECMW_GEOFEM_ETYPE_BEM1:
      return HECMW_ETYPE_BEM1;
    case HECMW_GEOFEM_ETYPE_BEM2:
      return HECMW_ETYPE_BEM2;
    case HECMW_GEOFEM_ETYPE_BEM3:
      return HECMW_ETYPE_BEM3;
    case HECMW_GEOFEM_ETYPE_SHT1:
      return HECMW_ETYPE_SHT1;
    case HECMW_GEOFEM_ETYPE_SHT2:
      return HECMW_ETYPE_SHT2;
    case HECMW_GEOFEM_ETYPE_SHQ1:
      return HECMW_ETYPE_SHQ1;
    case HECMW_GEOFEM_ETYPE_SHQ2:
      return HECMW_ETYPE_SHQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_max_node(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MAX_NODE_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MAX_NODE_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MAX_NODE_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MAX_NODE_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MAX_NODE_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_MAX_NODE_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MAX_NODE_QUA2;
    case HECMW_ETYPE_ROD31:
      return HECMW_MAX_NODE_ROD31;
    case HECMW_ETYPE_TET1:
      return HECMW_MAX_NODE_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MAX_NODE_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MAX_NODE_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_MAX_NODE_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_MAX_NODE_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MAX_NODE_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MAX_NODE_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MAX_NODE_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MAX_NODE_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MAX_NODE_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MAX_NODE_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MAX_NODE_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MAX_NODE_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MAX_NODE_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MAX_NODE_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MAX_NODE_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MAX_NODE_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MAX_NODE_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MAX_NODE_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MAX_NODE_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MAX_NODE_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MAX_NODE_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MAX_NODE_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MAX_NODE_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MAX_NODE_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MAX_NODE_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MAX_NODE_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MAX_NODE_SHQ2;
    case HECMW_ETYPE_SHQ3:
      return HECMW_MAX_NODE_SHQ3;
    case HECMW_ETYPE_SHT6:
      return HECMW_MAX_NODE_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MAX_NODE_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MAX_NODE_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MAX_NODE_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MAX_NODE_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MAX_NODE_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MAX_NODE_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MAX_NODE_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MAX_NODE_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MAX_NODE_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MAX_NODE_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MAX_NODE_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MAX_NODE_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MAX_NODE_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MAX_NODE_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MAX_NODE_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MAX_NODE_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MAX_NODE_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MAX_NODE_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MAX_NODE_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MAX_NODE_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MAX_NODE_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_MAX_NODE_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MAX_NODE_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MAX_NODE_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MAX_NODE_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MAX_NODE_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MAX_NODE_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MAX_NODE_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MAX_NODE_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MAX_NODE_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MAX_NODE_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MAX_NODE_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MAX_NODE_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MAX_NODE_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MAX_NODE_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MAX_NODE_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MAX_NODE_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MAX_NODE_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MAX_NODE_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MAX_NODE_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MAX_NODE_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_max_edge(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MAX_EDGE_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MAX_EDGE_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MAX_EDGE_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MAX_EDGE_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MAX_EDGE_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_MAX_EDGE_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MAX_EDGE_QUA2;
    case HECMW_ETYPE_TET1:
      return HECMW_MAX_EDGE_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MAX_EDGE_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MAX_EDGE_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_MAX_EDGE_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_MAX_EDGE_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MAX_EDGE_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MAX_EDGE_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MAX_EDGE_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MAX_EDGE_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MAX_EDGE_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MAX_EDGE_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MAX_EDGE_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MAX_EDGE_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MAX_EDGE_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MAX_EDGE_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MAX_EDGE_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MAX_EDGE_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MAX_EDGE_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MAX_EDGE_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MAX_EDGE_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MAX_EDGE_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MAX_EDGE_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MAX_EDGE_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MAX_EDGE_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MAX_EDGE_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MAX_EDGE_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MAX_EDGE_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MAX_EDGE_SHQ2;
    case HECMW_ETYPE_SHT6:
      return HECMW_MAX_EDGE_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MAX_EDGE_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MAX_EDGE_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MAX_EDGE_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MAX_EDGE_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MAX_EDGE_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MAX_EDGE_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MAX_EDGE_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MAX_EDGE_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MAX_EDGE_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MAX_EDGE_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MAX_EDGE_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MAX_EDGE_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MAX_EDGE_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MAX_EDGE_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MAX_EDGE_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MAX_EDGE_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MAX_EDGE_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MAX_EDGE_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MAX_EDGE_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MAX_EDGE_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MAX_EDGE_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_MAX_EDGE_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MAX_EDGE_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MAX_EDGE_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MAX_EDGE_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MAX_EDGE_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MAX_EDGE_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MAX_EDGE_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MAX_EDGE_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MAX_EDGE_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MAX_EDGE_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MAX_EDGE_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MAX_EDGE_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MAX_EDGE_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MAX_EDGE_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MAX_EDGE_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MAX_EDGE_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MAX_EDGE_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MAX_EDGE_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MAX_EDGE_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MAX_EDGE_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_max_surf(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MAX_SURF_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MAX_SURF_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MAX_SURF_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MAX_SURF_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MAX_SURF_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_MAX_SURF_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MAX_SURF_QUA2;
    case HECMW_ETYPE_TET1:
      return HECMW_MAX_SURF_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MAX_SURF_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MAX_SURF_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_MAX_SURF_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_MAX_SURF_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MAX_SURF_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MAX_SURF_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MAX_SURF_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MAX_SURF_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MAX_SURF_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MAX_SURF_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MAX_SURF_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MAX_SURF_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MAX_SURF_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MAX_SURF_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MAX_SURF_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MAX_SURF_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MAX_SURF_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MAX_SURF_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MAX_SURF_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MAX_SURF_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MAX_SURF_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MAX_SURF_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MAX_SURF_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MAX_SURF_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MAX_SURF_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MAX_SURF_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MAX_SURF_SHQ2;
    case HECMW_ETYPE_SHT6:
      return HECMW_MAX_SURF_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MAX_SURF_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MAX_SURF_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MAX_SURF_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MAX_SURF_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MAX_SURF_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MAX_SURF_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MAX_SURF_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MAX_SURF_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MAX_SURF_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MAX_SURF_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MAX_SURF_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MAX_SURF_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MAX_SURF_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MAX_SURF_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MAX_SURF_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MAX_SURF_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MAX_SURF_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MAX_SURF_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MAX_SURF_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MAX_SURF_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MAX_SURF_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_MAX_SURF_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MAX_SURF_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MAX_SURF_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MAX_SURF_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MAX_SURF_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MAX_SURF_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MAX_SURF_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MAX_SURF_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MAX_SURF_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MAX_SURF_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MAX_SURF_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MAX_SURF_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MAX_SURF_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MAX_SURF_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MAX_SURF_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MAX_SURF_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MAX_SURF_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MAX_SURF_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MAX_SURF_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MAX_SURF_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_max_tsuf(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MAX_TSUF_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MAX_TSUF_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MAX_TSUF_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MAX_TSUF_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MAX_TSUF_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_MAX_TSUF_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MAX_TSUF_QUA2;
    case HECMW_ETYPE_TET1:
      return HECMW_MAX_TSUF_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MAX_TSUF_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MAX_TSUF_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_MAX_TSUF_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MAX_TSUF_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MAX_TSUF_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MAX_TSUF_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MAX_TSUF_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MAX_TSUF_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MAX_TSUF_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MAX_TSUF_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MAX_TSUF_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MAX_TSUF_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MAX_TSUF_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MAX_TSUF_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MAX_TSUF_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MAX_TSUF_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MAX_TSUF_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MAX_TSUF_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MAX_TSUF_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MAX_TSUF_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MAX_TSUF_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MAX_TSUF_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MAX_TSUF_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MAX_TSUF_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MAX_TSUF_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MAX_TSUF_SHQ2;
    case HECMW_ETYPE_SHT6:
      return HECMW_MAX_TSUF_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MAX_TSUF_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MAX_TSUF_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MAX_TSUF_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MAX_TSUF_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MAX_TSUF_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MAX_TSUF_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MAX_TSUF_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MAX_TSUF_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MAX_TSUF_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MAX_TSUF_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MAX_TSUF_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MAX_TSUF_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MAX_TSUF_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MAX_TSUF_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MAX_TSUF_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MAX_TSUF_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MAX_TSUF_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MAX_TSUF_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MAX_TSUF_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MAX_TSUF_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MAX_TSUF_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_MAX_TSUF_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MAX_TSUF_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MAX_TSUF_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MAX_TSUF_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MAX_TSUF_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MAX_TSUF_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MAX_TSUF_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MAX_TSUF_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MAX_TSUF_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MAX_TSUF_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MAX_TSUF_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MAX_TSUF_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MAX_TSUF_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MAX_TSUF_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MAX_TSUF_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MAX_TSUF_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MAX_TSUF_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MAX_TSUF_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MAX_TSUF_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MAX_TSUF_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_max_qsuf(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_MAX_QSUF_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_MAX_QSUF_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_MAX_QSUF_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_MAX_QSUF_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_MAX_QSUF_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_MAX_QSUF_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_MAX_QSUF_QUA2;
    case HECMW_ETYPE_TET1:
      return HECMW_MAX_QSUF_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_MAX_QSUF_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_MAX_QSUF_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_MAX_QSUF_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_MAX_QSUF_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_MAX_QSUF_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_MAX_QSUF_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_MAX_QSUF_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_MAX_QSUF_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_MAX_QSUF_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_MAX_QSUF_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_MAX_QSUF_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_MAX_QSUF_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_MAX_QSUF_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_MAX_QSUF_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_MAX_QSUF_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_MAX_QSUF_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_MAX_QSUF_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_MAX_QSUF_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_MAX_QSUF_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_MAX_QSUF_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_MAX_QSUF_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_MAX_QSUF_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_MAX_QSUF_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_MAX_QSUF_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_MAX_QSUF_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_MAX_QSUF_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_MAX_QSUF_SHQ2;
    case HECMW_ETYPE_SHT6:
      return HECMW_MAX_QSUF_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_MAX_QSUF_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_MAX_QSUF_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_MAX_QSUF_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_MAX_QSUF_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_MAX_QSUF_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_MAX_QSUF_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_MAX_QSUF_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_MAX_QSUF_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_MAX_QSUF_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_MAX_QSUF_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_MAX_QSUF_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_MAX_QSUF_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_MAX_QSUF_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_MAX_QSUF_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_MAX_QSUF_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_MAX_QSUF_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_MAX_QSUF_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_MAX_QSUF_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_MAX_QSUF_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_MAX_QSUF_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_MAX_QSUF_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_MAX_QSUF_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_MAX_QSUF_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_MAX_QSUF_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_MAX_QSUF_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_MAX_QSUF_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_MAX_QSUF_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_MAX_QSUF_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_MAX_QSUF_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_MAX_QSUF_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_MAX_QSUF_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_MAX_QSUF_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_MAX_QSUF_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_MAX_QSUF_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_MAX_QSUF_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_MAX_QSUF_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_MAX_QSUF_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_MAX_QSUF_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_MAX_QSUF_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_MAX_QSUF_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_MAX_QSUF_PTQ2;
    default:
      return -1;
  }

  return -1;
}

extern char *HECMW_get_ucd_label(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_UCD_LABEL_ROD1;
    case HECMW_ETYPE_ROD2:
      return HECMW_UCD_LABEL_ROD2;
    case HECMW_ETYPE_TRI1:
      return HECMW_UCD_LABEL_TRI1;
    case HECMW_ETYPE_TRI2:
      return HECMW_UCD_LABEL_TRI2;
    case HECMW_ETYPE_TRI22:
      return HECMW_UCD_LABEL_TRI2;
    case HECMW_ETYPE_QUA1:
      return HECMW_UCD_LABEL_QUA1;
    case HECMW_ETYPE_QUA2:
      return HECMW_UCD_LABEL_QUA2;
    case HECMW_ETYPE_ROD31:
      return HECMW_UCD_LABEL_ROD31;
    case HECMW_ETYPE_TET1:
      return HECMW_UCD_LABEL_TET1;
    case HECMW_ETYPE_TET1_4:
      return HECMW_UCD_LABEL_TET1_4;
    case HECMW_ETYPE_TET2:
      return HECMW_UCD_LABEL_TET2;
    case HECMW_ETYPE_TET22:
      return HECMW_UCD_LABEL_TET2;
    case HECMW_ETYPE_PRI1:
      return HECMW_UCD_LABEL_PRI1;
    case HECMW_ETYPE_PRI2:
      return HECMW_UCD_LABEL_PRI2;
    case HECMW_ETYPE_HEX1:
      return HECMW_UCD_LABEL_HEX1;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_UCD_LABEL_HEX1_4;
    case HECMW_ETYPE_HEX2:
      return HECMW_UCD_LABEL_HEX2;
    case HECMW_ETYPE_PYR1:
      return HECMW_UCD_LABEL_PYR1;
    case HECMW_ETYPE_PYR2:
      return HECMW_UCD_LABEL_PYR2;
    case HECMW_ETYPE_MST1:
      return HECMW_UCD_LABEL_MST1;
    case HECMW_ETYPE_MST2:
      return HECMW_UCD_LABEL_MST2;
    case HECMW_ETYPE_MSQ1:
      return HECMW_UCD_LABEL_MSQ1;
    case HECMW_ETYPE_MSQ2:
      return HECMW_UCD_LABEL_MSQ2;
    case HECMW_ETYPE_JTB1:
      return HECMW_UCD_LABEL_JTB1;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_UCD_LABEL_SPGDPT1;
    case HECMW_ETYPE_JTT1:
      return HECMW_UCD_LABEL_JTT1;
    case HECMW_ETYPE_JTT2:
      return HECMW_UCD_LABEL_JTT2;
    case HECMW_ETYPE_JTQ1:
      return HECMW_UCD_LABEL_JTQ1;
    case HECMW_ETYPE_JTQ2:
      return HECMW_UCD_LABEL_JTQ2;
    case HECMW_ETYPE_BEM1:
      return HECMW_UCD_LABEL_BEM1;
    case HECMW_ETYPE_BEM2:
      return HECMW_UCD_LABEL_BEM2;
    case HECMW_ETYPE_BEM3:
      return HECMW_UCD_LABEL_BEM3; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_UCD_LABEL_SHT1;
    case HECMW_ETYPE_SHT2:
      return HECMW_UCD_LABEL_SHT2;
    case HECMW_ETYPE_SHQ1:
      return HECMW_UCD_LABEL_SHQ1;
    case HECMW_ETYPE_SHQ2:
      return HECMW_UCD_LABEL_SHQ2;
    case HECMW_ETYPE_SHQ3:
      return HECMW_UCD_LABEL_SHQ2;
    case HECMW_ETYPE_SHT6:
      return HECMW_UCD_LABEL_SHT6; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_UCD_LABEL_SHQ8;
    case HECMW_ETYPE_LN11:
      return HECMW_UCD_LABEL_LN11;
    case HECMW_ETYPE_LN12:
      return HECMW_UCD_LABEL_LN12;
    case HECMW_ETYPE_LN13:
      return HECMW_UCD_LABEL_LN13;
    case HECMW_ETYPE_LN14:
      return HECMW_UCD_LABEL_LN14;
    case HECMW_ETYPE_LN15:
      return HECMW_UCD_LABEL_LN15;
    case HECMW_ETYPE_LN16:
      return HECMW_UCD_LABEL_LN16;
    case HECMW_ETYPE_LN21:
      return HECMW_UCD_LABEL_LN21;
    case HECMW_ETYPE_LN22:
      return HECMW_UCD_LABEL_LN22;
    case HECMW_ETYPE_LN23:
      return HECMW_UCD_LABEL_LN23;
    case HECMW_ETYPE_LN24:
      return HECMW_UCD_LABEL_LN24;
    case HECMW_ETYPE_LN25:
      return HECMW_UCD_LABEL_LN25;
    case HECMW_ETYPE_LN26:
      return HECMW_UCD_LABEL_LN26;
    case HECMW_ETYPE_LN31:
      return HECMW_UCD_LABEL_LN31;
    case HECMW_ETYPE_LN32:
      return HECMW_UCD_LABEL_LN32;
    case HECMW_ETYPE_LN33:
      return HECMW_UCD_LABEL_LN33;
    case HECMW_ETYPE_LN34:
      return HECMW_UCD_LABEL_LN34;
    case HECMW_ETYPE_LN35:
      return HECMW_UCD_LABEL_LN35;
    case HECMW_ETYPE_LN36:
      return HECMW_UCD_LABEL_LN36;
    case HECMW_ETYPE_LN41:
      return HECMW_UCD_LABEL_LN41;
    case HECMW_ETYPE_LN42:
      return HECMW_UCD_LABEL_LN42;
    case HECMW_ETYPE_LN43:
      return HECMW_UCD_LABEL_LN43;
    case HECMW_ETYPE_LN44:
      return HECMW_UCD_LABEL_LN44;
    case HECMW_ETYPE_LN45:
      return HECMW_UCD_LABEL_LN45;
    case HECMW_ETYPE_LN46:
      return HECMW_UCD_LABEL_LN46;
    case HECMW_ETYPE_LN51:
      return HECMW_UCD_LABEL_LN51;
    case HECMW_ETYPE_LN52:
      return HECMW_UCD_LABEL_LN52;
    case HECMW_ETYPE_LN53:
      return HECMW_UCD_LABEL_LN53;
    case HECMW_ETYPE_LN54:
      return HECMW_UCD_LABEL_LN54;
    case HECMW_ETYPE_LN55:
      return HECMW_UCD_LABEL_LN55;
    case HECMW_ETYPE_LN56:
      return HECMW_UCD_LABEL_LN56;
    case HECMW_ETYPE_LN61:
      return HECMW_UCD_LABEL_LN61;
    case HECMW_ETYPE_LN62:
      return HECMW_UCD_LABEL_LN62;
    case HECMW_ETYPE_LN63:
      return HECMW_UCD_LABEL_LN63;
    case HECMW_ETYPE_LN64:
      return HECMW_UCD_LABEL_LN64;
    case HECMW_ETYPE_LN65:
      return HECMW_UCD_LABEL_LN65;
    case HECMW_ETYPE_LN66:
      return HECMW_UCD_LABEL_LN66;
    case HECMW_ETYPE_PTT1:
      return HECMW_UCD_LABEL_PTT1;
    case HECMW_ETYPE_PTT2:
      return HECMW_UCD_LABEL_PTT2;
    case HECMW_ETYPE_PTQ1:
      return HECMW_UCD_LABEL_PTQ1;
    case HECMW_ETYPE_PTQ2:
      return HECMW_UCD_LABEL_PTQ2;
    default:
      return NULL;
  }

  return NULL;
}

extern int HECMW_get_etype_class(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_CLASS_LINE;
    case HECMW_ETYPE_ROD2:
      return HECMW_CLASS_LINE;
    case HECMW_ETYPE_TRI1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_TRI2:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_TRI22:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_QUA1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_QUA2:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_ROD31:
      return HECMW_CLASS_LINE;
    case HECMW_ETYPE_TET1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_TET1_4:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_TET2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_TET22:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_PRI1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_PRI2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_HEX1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_HEX2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_PYR1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_PYR2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_MST1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_MST2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_MSQ1:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_MSQ2:
      return HECMW_CLASS_SOLID;
    case HECMW_ETYPE_JTB1:
      return HECMW_CLASS_JOINT;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_CLASS_LINE;
    case HECMW_ETYPE_JTT1:
      return HECMW_CLASS_JOINT;
    case HECMW_ETYPE_JTT2:
      return HECMW_CLASS_JOINT;
    case HECMW_ETYPE_JTQ1:
      return HECMW_CLASS_JOINT;
    case HECMW_ETYPE_JTQ2:
      return HECMW_CLASS_JOINT;
    case HECMW_ETYPE_BEM1:
      return HECMW_CLASS_LINE;
    case HECMW_ETYPE_BEM2:
      return HECMW_CLASS_LINE;
    case HECMW_ETYPE_BEM3:
      return HECMW_CLASS_LINE; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_SHT2:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_SHQ1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_SHQ2:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_SHQ3:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_SHT6:
      return HECMW_CLASS_SURF; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_LN11:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN12:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN13:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN14:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN15:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN16:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN21:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN22:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN23:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN24:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN25:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN26:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN31:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN32:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN33:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN34:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN35:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN36:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN41:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN42:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN43:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN44:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN45:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN46:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN51:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN52:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN53:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN54:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN55:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN56:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN61:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN62:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN63:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN64:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN65:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_LN66:
      return HECMW_CLASS_LINK;
    case HECMW_ETYPE_PTT1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_PTT2:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_PTQ1:
      return HECMW_CLASS_SURF;
    case HECMW_ETYPE_PTQ2:
      return HECMW_CLASS_SURF;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_etype_shape(int etype) {
  switch (etype) {
    case HECMW_ETYPE_ROD1:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_ROD2:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_TRI1:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_TRI2:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_TRI22:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_QUA1:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_QUA2:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_ROD31:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_TET1:
      return HECMW_SHAPE_TETRA;
    case HECMW_ETYPE_TET1_4:
      return HECMW_SHAPE_TETRA;
    case HECMW_ETYPE_TET2:
      return HECMW_SHAPE_TETRA;
    case HECMW_ETYPE_TET22:
      return HECMW_SHAPE_TETRA;
    case HECMW_ETYPE_PRI1:
      return HECMW_SHAPE_PRISM;
    case HECMW_ETYPE_PRI2:
      return HECMW_SHAPE_PRISM;
    case HECMW_ETYPE_HEX1:
      return HECMW_SHAPE_HEXA;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_SHAPE_HEXA;
    case HECMW_ETYPE_HEX2:
      return HECMW_SHAPE_HEXA;
    case HECMW_ETYPE_PYR1:
      return HECMW_SHAPE_PYRAM;
    case HECMW_ETYPE_PYR2:
      return HECMW_SHAPE_PYRAM;
    case HECMW_ETYPE_MST1:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_MST2:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_MSQ1:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_MSQ2:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_JTB1:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_JTT1:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_JTT2:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_JTQ1:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_JTQ2:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_BEM1:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_BEM2:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_BEM3:
      return HECMW_SHAPE_LINE; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_SHT2:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_SHQ1:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_SHQ2:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_SHQ3:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_SHT6:
      return HECMW_SHAPE_TRI; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_LN11:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN12:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN13:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN14:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN15:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN16:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN21:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN22:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN23:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN24:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN25:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN26:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN31:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN32:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN33:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN34:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN35:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN36:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN41:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN42:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN43:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN44:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN45:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN46:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN51:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN52:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN53:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN54:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN55:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN56:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN61:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN62:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN63:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN64:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN65:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_LN66:
      return HECMW_SHAPE_LINE;
    case HECMW_ETYPE_PTT1:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_PTT2:
      return HECMW_SHAPE_TRI;
    case HECMW_ETYPE_PTQ1:
      return HECMW_SHAPE_QUAD;
    case HECMW_ETYPE_PTQ2:
      return HECMW_SHAPE_QUAD;
    default:
      return -1;
  }

  return -1;
}

extern int HECMW_get_etype_vtk_shape(int etype) {
  switch (etype) {

    case HECMW_ETYPE_ROD1:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_ROD2:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_TRI1:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_TRI2:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_TRI22:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_QUA1:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_QUA2:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_ROD31:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_TET1:
      return HECMW_VTK_SHAPE_TETRA;
    case HECMW_ETYPE_TET1_4:
      return HECMW_VTK_SHAPE_TETRA;
    case HECMW_ETYPE_TET2:
      return HECMW_VTK_SHAPE_TETRA2;
    case HECMW_ETYPE_TET22:
      return HECMW_VTK_SHAPE_TETRA;
    case HECMW_ETYPE_PRI1:
      return HECMW_VTK_SHAPE_PRISM;
    case HECMW_ETYPE_PRI2:
      return HECMW_VTK_SHAPE_PRISM;
    case HECMW_ETYPE_HEX1:
      return HECMW_VTK_SHAPE_HEXA;
    case HECMW_ETYPE_HEX1_4:
      return HECMW_VTK_SHAPE_HEXA;
    case HECMW_ETYPE_HEX2:
      return HECMW_VTK_SHAPE_HEXA;
    case HECMW_ETYPE_PYR1:
      return HECMW_VTK_SHAPE_PYRAM;
    case HECMW_ETYPE_PYR2:
      return HECMW_VTK_SHAPE_PYRAM;
    case HECMW_ETYPE_MST1:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_MST2:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_MSQ1:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_MSQ2:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_JTB1:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_SPGDPT1:
        return HECMW_VTK_SHAPE_LINE;
        /* case HECMW_ETYPE_JTR1:
      return HECMW_VTK_SHAPE_LINE; */
    case HECMW_ETYPE_JTT1:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_JTT2:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_JTQ1:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_JTQ2:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_BEM1:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_BEM2:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_BEM3:
      return HECMW_VTK_SHAPE_LINE; /* mixed beam-341*/
    case HECMW_ETYPE_SHT1:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_SHT2:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_SHQ1:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_SHQ2:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_SHQ3:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_SHT6:
      return HECMW_VTK_SHAPE_TRI; /* mixed shell-solid */
    case HECMW_ETYPE_SHQ8:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_LN11:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN12:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN13:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN14:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN15:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN16:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN21:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN22:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN23:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN24:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN25:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN26:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN31:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN32:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN33:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN34:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN35:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN36:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN41:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN42:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN43:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN44:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN45:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN46:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN51:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN52:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN53:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN54:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN55:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN56:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN61:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN62:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN63:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN64:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN65:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_LN66:
      return HECMW_VTK_SHAPE_LINE;
    case HECMW_ETYPE_PTT1:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_PTT2:
      return HECMW_VTK_SHAPE_TRI;
    case HECMW_ETYPE_PTQ1:
      return HECMW_VTK_SHAPE_QUAD;
    case HECMW_ETYPE_PTQ2:
      return HECMW_VTK_SHAPE_QUAD;
    default:
      return -1;
  }

  return -1;
}


/*============================================================================*/

extern int HECMW_is_etype_rod(int etype) {
  switch (etype) {
    case 111: /* fall through */
    case 112:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_surface(int etype) {
  switch (etype) {
    case 231: /* fall through */
    case 232: /* fall through */
    case 241: /* fall through */
    case 242:
    case 2322:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_solid(int etype) {
  switch (etype) {
    case 301:
    case 341:  /* fall through */
    case 342:  /* fall through */
    case 3414: /* fall through */
    case 3422:
    case 351:  /* fall through */
    case 352:  /* fall through */
    case 361:  /* fall through */
    case 362:  /* fall through */
    case 3614: /* fall through */
    case 371:  /* fall through */
    case 372:
      return 1;
  }
  if (HECMW_is_etype_rod(etype)) return 1;
  if (HECMW_is_etype_surface(etype)) return 1;
  return 0;
}

extern int HECMW_is_etype_interface(int etype) {
  switch (etype) {
    case 431: /* fall through */
    case 432: /* fall through */
    case 441: /* fall through */
    case 442: /* fall through */
    case 501: /* fall through */
    case 511: /* fall through */
    case 531: /* fall through */
    case 532: /* fall through */
    case 541: /* fall through */
    case 542:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_beam(int etype) {
  switch (etype) {
    case 611: /* fall through */
    case 612:
    case 641:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_shell(int etype) {
  switch (etype) {
    case 731: /* fall through */
    case 732: /* fall through */
    case 741: /* fall through */
    case 743:
    case 761:
    case 781:
    case 742:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_link(int etype) {
  switch (etype) {
    case 911: /* fall through */
    case 912: /* fall through */
    case 913: /* fall through */
    case 914: /* fall through */
    case 915: /* fall through */
    case 916: /* fall through */
    case 921: /* fall through */
    case 922: /* fall through */
    case 923: /* fall through */
    case 924: /* fall through */
    case 925: /* fall through */
    case 926: /* fall through */
    case 931: /* fall through */
    case 932: /* fall through */
    case 933: /* fall through */
    case 934: /* fall through */
    case 935: /* fall through */
    case 936: /* fall through */
    case 941: /* fall through */
    case 942: /* fall through */
    case 943: /* fall through */
    case 944: /* fall through */
    case 945: /* fall through */
    case 946: /* fall through */
    case 951: /* fall through */
    case 952: /* fall through */
    case 953: /* fall through */
    case 954: /* fall through */
    case 955: /* fall through */
    case 956: /* fall through */
    case 961: /* fall through */
    case 962: /* fall through */
    case 963: /* fall through */
    case 964: /* fall through */
    case 965: /* fall through */
    case 966:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_33struct(int etype) {
  switch (etype) {
    case HECMW_ETYPE_BEM3: /* fall through */
    case HECMW_ETYPE_SHT6: /* fall through */
    case HECMW_ETYPE_SHQ8:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_truss(int etype) {
  switch (etype) {
    case 301:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_patch(int etype) {
  switch (etype) {
    case HECMW_ETYPE_PTT1: /* fall through */
    case HECMW_ETYPE_PTT2: /* fall through */
    case HECMW_ETYPE_PTQ1: /* fall through */
    case HECMW_ETYPE_PTQ2:
      return 1;
  }
  return 0;
}

extern int HECMW_is_etype_smoothing(int etype) {
  switch (etype) {
    case 881:
    case 891:
      return 1;
  }
  return 0;
}

extern const int *HECMW_get_surf_nodes(int etype, int sid, int *surf_etype)
{
  static const int snodes_tet[4][6] = {
    {1, 2, 3,  5,  6, 7},
    {1, 2, 4,  9,  8, 7},
    {2, 3, 4, 10,  9, 5},
    {3, 1, 4,  8, 10, 6} };
  static const int snodes_pri[5][8] = {
    {1, 2, 3,     7,  8,  9},
    {4, 5, 6,    10, 11, 12},
    {1, 2, 5, 4,  9, 14, 12, 13},
    {2, 3, 6, 5,  7, 15, 10, 14},
    {3, 1, 4, 6,  8, 13, 11, 15} };
  static const int snodes_hex[6][8] = {
    {1, 2, 3, 4,  9, 10, 11, 12},
    {5, 6, 7, 8, 13, 14, 15, 16},
    {1, 2, 6, 5,  9, 18, 13, 17},
    {2, 3, 7, 6, 10, 19, 14, 18},
    {3, 4, 8, 7, 11, 20, 15, 19},
    {4, 1, 5, 8, 12, 17, 16, 20} };
  static const int snodes_pyr[5][8] = {
    {4, 1, 5,    10, 13, 9},
    {2, 3, 5,    12, 11, 7},
    {1, 2, 5,    11, 10, 6},
    {3, 4, 5,    13, 12, 8},
    {4, 3, 2, 1,  8,  7, 6, 9} };

  switch (etype) {
  case HECMW_ETYPE_TET1:
    HECMW_assert( 0 < sid && sid <= 4 );
    *surf_etype = HECMW_ETYPE_PTT1;
    return snodes_tet[sid-1];
  case HECMW_ETYPE_TET2:
    HECMW_assert( 0 < sid && sid <= 4 );
    *surf_etype = HECMW_ETYPE_PTT2;
    return snodes_tet[sid-1];
  case HECMW_ETYPE_PRI1:
    HECMW_assert( 0 < sid && sid <= 5 );
    if (sid <= 2) {
      *surf_etype = HECMW_ETYPE_PTT1;
    } else {
      *surf_etype = HECMW_ETYPE_PTQ1;
    }
    return snodes_pri[sid-1];
  case HECMW_ETYPE_PRI2:
    HECMW_assert( 0 < sid && sid <= 5 );
    if (sid <= 2) {
      *surf_etype = HECMW_ETYPE_PTT2;
    } else {
      *surf_etype = HECMW_ETYPE_PTQ2;
    }
    return snodes_pri[sid-1];
  case HECMW_ETYPE_HEX1:
    HECMW_assert( 0 < sid && sid <= 6 );
    *surf_etype = HECMW_ETYPE_PTQ1;
    return snodes_hex[sid-1];
  case HECMW_ETYPE_HEX2:
    HECMW_assert( 0 < sid && sid <= 6 );
    *surf_etype = HECMW_ETYPE_PTQ2;
    return snodes_hex[sid-1];
  case HECMW_ETYPE_PYR1:
    HECMW_assert( 0 < sid && sid <= 5 );
    if (sid <= 4) {
      *surf_etype = HECMW_ETYPE_PTT1;
    } else {
      *surf_etype = HECMW_ETYPE_PTQ1;
    }
    return snodes_pyr[sid-1];
  case HECMW_ETYPE_PYR2:
    HECMW_assert( 0 < sid && sid <= 5 );
    if (sid <= 4) {
      *surf_etype = HECMW_ETYPE_PTT2;
    } else {
      *surf_etype = HECMW_ETYPE_PTQ2;
    }
    return snodes_pyr[sid-1];
  default:
    return NULL;
  }
  return NULL;
}

/* interface for fortran -- added by Kazuya Goto (AdvanceSoft) */

int hecmw_get_max_node_if(int *etype) { return HECMW_get_max_node(*etype); }
int hecmw_get_max_node_if_(int *etype) { return HECMW_get_max_node(*etype); }
int hecmw_get_max_node_if__(int *etype) { return HECMW_get_max_node(*etype); }
int HECMW_GET_MAX_NODE_IF(int *etype) { return HECMW_get_max_node(*etype); }

int hecmw_is_etype_rod_if(int *etype) { return HECMW_is_etype_rod(*etype); }
int hecmw_is_etype_rod_if_(int *etype) { return HECMW_is_etype_rod(*etype); }
int hecmw_is_etype_rod_if__(int *etype) { return HECMW_is_etype_rod(*etype); }
int HECMW_IS_ETYPE_ROD_IF(int *etype) { return HECMW_is_etype_rod(*etype); }

int hecmw_is_etype_surface_if(int *etype) {
  return HECMW_is_etype_surface(*etype);
}
int hecmw_is_etype_surface_if_(int *etype) {
  return HECMW_is_etype_surface(*etype);
}
int hecmw_is_etype_surface_if__(int *etype) {
  return HECMW_is_etype_surface(*etype);
}
int HECMW_IS_ETYPE_SURFACE_IF(int *etype) {
  return HECMW_is_etype_surface(*etype);
}

int hecmw_is_etype_solid_if(int *etype) { return HECMW_is_etype_solid(*etype); }
int hecmw_is_etype_solid_if_(int *etype) {
  return HECMW_is_etype_solid(*etype);
}
int hecmw_is_etype_solid_if__(int *etype) {
  return HECMW_is_etype_solid(*etype);
}
int HECMW_IS_ETYPE_SOLID_IF(int *etype) { return HECMW_is_etype_solid(*etype); }

int hecmw_is_etype_interface_if(int *etype) {
  return HECMW_is_etype_interface(*etype);
}
int hecmw_is_etype_interface_if_(int *etype) {
  return HECMW_is_etype_interface(*etype);
}
int hecmw_is_etype_interface_if__(int *etype) {
  return HECMW_is_etype_interface(*etype);
}
int HECMW_IS_ETYPE_INTERFACE_IF(int *etype) {
  return HECMW_is_etype_interface(*etype);
}

int hecmw_is_etype_beam_if(int *etype) { return HECMW_is_etype_beam(*etype); }
int hecmw_is_etype_beam_if_(int *etype) { return HECMW_is_etype_beam(*etype); }
int hecmw_is_etype_beam_if__(int *etype) { return HECMW_is_etype_beam(*etype); }
int HECMW_IS_ETYPE_BEAM_IF(int *etype) { return HECMW_is_etype_beam(*etype); }

int hecmw_is_etype_shell_if(int *etype) { return HECMW_is_etype_shell(*etype); }
int hecmw_is_etype_shell_if_(int *etype) {
  return HECMW_is_etype_shell(*etype);
}
int hecmw_is_etype_shell_if__(int *etype) {
  return HECMW_is_etype_shell(*etype);
}
int HECMW_IS_ETYPE_SHELL_IF(int *etype) { return HECMW_is_etype_shell(*etype); }

int hecmw_is_etype_link_if(int *etype) { return HECMW_is_etype_link(*etype); }
int hecmw_is_etype_link_if_(int *etype) { return HECMW_is_etype_link(*etype); }
int hecmw_is_etype_link_if__(int *etype) { return HECMW_is_etype_link(*etype); }
int HECMW_IS_ETYPE_LINK_IF(int *etype) { return HECMW_is_etype_link(*etype); }

int hecmw_is_etype_33struct_if(int *etype) {
  return HECMW_is_etype_33struct(*etype);
}
int hecmw_is_etype_33struct_if_(int *etype) {
  return HECMW_is_etype_33struct(*etype);
}
int hecmw_is_etype_33struct_if__(int *etype) {
  return HECMW_is_etype_33struct(*etype);
}
int HECMW_IS_ETYPE_33STRUCT_IF(int *etype) {
  return HECMW_is_etype_33struct(*etype);
}

int hecmw_is_etype_truss_if(int *etype) { return HECMW_is_etype_truss(*etype); }
int hecmw_is_etype_truss_if_(int *etype) {
  return HECMW_is_etype_truss(*etype);
}
int hecmw_is_etype_truss_if__(int *etype) {
  return HECMW_is_etype_truss(*etype);
}
int HECMW_IS_ETYPE_TRUSS_IF(int *etype) { return HECMW_is_etype_truss(*etype); }

int hecmw_is_etype_patch_if(int *etype) { return HECMW_is_etype_patch(*etype); }
int hecmw_is_etype_patch_if_(int *etype) { return HECMW_is_etype_patch(*etype); }
int hecmw_is_etype_patch_if__(int *etype) { return HECMW_is_etype_patch(*etype); }
int HECMW_IS_ETYPE_PATCH_IF(int *etype) { return HECMW_is_etype_patch(*etype); }

int hecmw_is_etype_smoothing_if(int *etype) { return HECMW_is_etype_smoothing(*etype); }
int hecmw_is_etype_smoothing_if_(int *etype) { return HECMW_is_etype_smoothing(*etype); }
int hecmw_is_etype_smoothing_if__(int *etype) { return HECMW_is_etype_smoothing(*etype); }
int HECMW_IS_ETYPE_SMOOTHING_IF(int *etype) { return HECMW_is_etype_smoothing(*etype); }
