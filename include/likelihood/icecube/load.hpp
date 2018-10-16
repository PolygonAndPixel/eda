#ifndef _LOAD_H
#define _LOAD_H
/* The code here is based (or mostly taken) from I3PhotoSplineTable.cxx
 * SetupTable
 * 
 */ 
#include <iostream>
#include "types.hpp"
#include "helper/photospline/geo_type.h"
#include "helper/photospline/splinetable.h"

// index_t readsplinefitstable(const char *path, struct splinetable *table);
// index_t splinetable_read_key(const struct splinetable *table, 
//     value_t type, const char *key, void *result);

bool load_splines(splinetable *& tablestruct_, char * filename);

#endif 

// weak scaling linear (problem size per processor is fixed)