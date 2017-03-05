//Local
#include "common.h"

/****************************************************
 **  __constant__ memory used for storing filters  **
 ***************************************************/

__constant__ DTYPE c_kern_L[MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_H[MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_IL[MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_IH[MAX_FILTER_WIDTH];

__constant__ DTYPE c_kern_LL[MAX_FILTER_WIDTH * MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_LH[MAX_FILTER_WIDTH * MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_HL[MAX_FILTER_WIDTH * MAX_FILTER_WIDTH];
__constant__ DTYPE c_kern_HH[MAX_FILTER_WIDTH * MAX_FILTER_WIDTH];
