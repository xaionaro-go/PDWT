// Local
#include "filters.cu.h"

/****************************************************
 **  __constant__ memory used for storing filters  **
 ***************************************************/

__constant__ float float_c_kern_L[MAX_FILTER_WIDTH];
__constant__ float float_c_kern_H[MAX_FILTER_WIDTH];
__constant__ float float_c_kern_IL[MAX_FILTER_WIDTH];
__constant__ float float_c_kern_IH[MAX_FILTER_WIDTH];

__constant__ double double_c_kern_L[MAX_FILTER_WIDTH];
__constant__ double double_c_kern_H[MAX_FILTER_WIDTH];
__constant__ double double_c_kern_IL[MAX_FILTER_WIDTH];
__constant__ double double_c_kern_IH[MAX_FILTER_WIDTH];
