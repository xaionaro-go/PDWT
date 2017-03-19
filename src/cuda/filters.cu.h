#ifndef FILTERS_CU_H
#define FILTERS_CU_H

#define MAX_FILTER_WIDTH 40

// Filter buffers in cuda constant memory
extern __constant__ float float_c_kern_L[MAX_FILTER_WIDTH];
extern __constant__ float float_c_kern_H[MAX_FILTER_WIDTH];
extern __constant__ float float_c_kern_IL[MAX_FILTER_WIDTH];
extern __constant__ float float_c_kern_IH[MAX_FILTER_WIDTH];

extern __constant__ double double_c_kern_L[MAX_FILTER_WIDTH];
extern __constant__ double double_c_kern_H[MAX_FILTER_WIDTH];
extern __constant__ double double_c_kern_IL[MAX_FILTER_WIDTH];
extern __constant__ double double_c_kern_IH[MAX_FILTER_WIDTH];

#endif //FILTERS_CU_H

