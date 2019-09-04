/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef hecmw_is_big_endian
#define hecmw_is_big_endian

#include <stdbool.h>

/**
 * @fn void HECMW_is_big_endian(void)
 * @brief function for investigating endian of running CPU
 * @param[in] void
 * @return bool if running CPU is big endian, return true,
 * @details if running CPU is big endian, return true.
 *          if it is littel endian, return false.
 */
extern bool HECMW_is_big_endian(void);

/**
 * @fn void HECMW_is_big_endian(void)
 * @brief function for investigating endian of running CPU
 * @param[in] void
 * @return const char* if running CPU is big endian, return "BigEndian",
 * @details if running CPU is big endian, return "BigEndian" as const char*.
 *          if it is littel endian, return "LittleEndian" as const char*.
 */
extern const char* HECMW_endian_str(void);

#endif /* hecmw_is_big_endian */
