#ifndef COMMON_H
#define COMMON_H

#include "utils.h"


__global__ void w_kern_soft_thresh(DTYPE* c_h, DTYPE* c_v, DTYPE* c_d,
    DTYPE beta, int Nr, int Nc);
__global__ void w_kern_soft_thresh_1d(DTYPE* c_d, DTYPE beta, int Nr, int Nc);
__global__ void w_kern_soft_thresh_appcoeffs(DTYPE* c_a, DTYPE beta, int Nr,
    int Nc);
__global__ void w_kern_proj_linf(DTYPE* c_h, DTYPE* c_v, DTYPE* c_d,
    DTYPE beta, int Nr, int Nc);
__global__ void w_kern_proj_linf_1d(DTYPE* c_d, DTYPE beta, int Nr, int Nc);
__global__ void w_kern_proj_linf_appcoeffs(DTYPE* c_a, DTYPE beta, int Nr,
    int Nc);


__global__ void w_kern_hard_thresh(DTYPE* c_h, DTYPE* c_v, DTYPE* c_d,
    DTYPE beta, int Nr, int Nc);
__global__ void w_kern_hard_thresh_appcoeffs(DTYPE* c_a, DTYPE beta, int Nr,
    int Nc);

__global__ void w_kern_circshift(DTYPE* d_image, DTYPE* d_out, int Nr, int Nc,
    int sr, int sc);

void w_call_soft_thresh(DTYPE** d_coeffs, DTYPE beta, w_info winfos,
    int do_thresh_appcoeffs, int normalize);
void w_call_hard_thresh(DTYPE** d_coeffs, DTYPE beta, w_info winfos,
    int do_thresh_appcoeffs, int normalize);
void w_call_proj_linf(DTYPE** d_coeffs, DTYPE beta, w_info winfos,
    int do_thresh_appcoeffs);
void w_shrink(DTYPE** d_coeffs, DTYPE beta, w_info winfos, int do_thresh_appcoeffs);
void w_call_circshift(DTYPE* d_image, DTYPE* d_image2, w_info winfos, int sr,
    int sc, int inplace = 1);

DTYPE** w_create_coeffs_buffer(w_info winfos);
void w_free_coeffs_buffer(DTYPE** coeffs, int nlevels);
void w_copy_coeffs_buffer(DTYPE** dst, DTYPE** src, w_info winfos);

DTYPE** w_create_coeffs_buffer_1d(w_info winfos);
void w_free_coeffs_buffer_1d(DTYPE** coeffs, int nlevels);
void w_copy_coeffs_buffer_1d(DTYPE** dst, DTYPE** src, w_info winfos);

__global__ void w_kern_hard_thresh_1d(DTYPE* c_d, DTYPE beta, int Nr, int Nc);
__global__ void w_kern_soft_thresh_1d(DTYPE* c_d, DTYPE beta, int Nr, int Nc);

void w_add_coeffs(DTYPE** dst, DTYPE** src, w_info winfos, DTYPE alpha=1.0f);
void w_add_coeffs_1d(DTYPE** dst, DTYPE** src, w_info winfos, DTYPE alpha=1.0f);

#endif
