/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdint.h>

#include "hecmw_fstr_endian.h"

/**
 * @fn void HECMW_is_big_endian(void)
 * @brief function for investigating endian of running CPU
 * @param[in] void
 * @return bool if running CPU is big endian, return true,
 * @details if running CPU is big endian, return true.
 *          if it is littel endian, return false.
 */
bool HECMW_is_big_endian(void) {
  union {
    uint32_t i;
    char c[4];
  } endian = {0x01020304};

  return endian.c[0] == true;
}

/**
 * @fn void HECMW_endian_str(void)
 * @brief function for investigating endian of running CPU
 * @param[in] void
 * @return const char* if running CPU is big endian, return "BigEndian",
 * @details if running CPU is big endian, return "BigEndian" as const char*.
 *          if it is littel endian, return "LittleEndian" as const char*.
 */
const char* HECMW_endian_str(void) {
  if (HECMW_is_big_endian() == true) {
    return "BigEndian";
  } else {
    return "LittleEndian";
  }
}
