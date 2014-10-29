/************************************************************************************************
* Implementing Singular Value Decomposition on GPU using CUDA using algorithm 			*
* given in IPDPS '09 paper "Singular Value Decomposition on GPU using CUDA"			*
*												*
* Copyright (c) 2009 International Institute of Information Technology, Hyderabad.		*
* All rights reserved.										*
*												*
* Permission to use, copy, modify and distribute this software and its documentation for 	*
* educational purpose is hereby granted without fee, provided that the above copyright		*
* notice and this permission notice appear in all copies of this software and that you do 	*
* not sell the software.									*
*												*
* THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR	*
* OTHERWISE.											*
* 												*
* Created by Sheetal Lahabar.									*
************************************************************************************************/

#ifndef _CUSVD_H_
#define _CUSVD_H_

/* Header files included */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>
#include <time.h>
#include <unistd.h>

#include <cutil.h>
#include "cublas.h"

cublasStatus status;

/* *****************************************************************************
 * This function call computes the Singular Value Decomposition of a real
 * dense matrix d_A with leading dimension M. It decomposes d_A as
 *
 *	d_A = d_MU * Diag(Sigma) * d_MVT
 *
 * Diag(Sigma) is a diagonal matrix with the elements of Sigma on the diagonal.
 *
 * It calls the subroutines cubidiagonal() and cudiagonal().
 *
 * It computes the left and right orthogonal matrices and stores in d_MU(MxM)
 * and d_MVT(NxN) respectively.
 * The singular values of d_A(MxN) is stored in Sigma vector on the CPU.
 *
 * Input parameters: 
 * 	  d_A - Matrix A on the device
 * 	 M, N - Dimensions of A
 *
 * Output parameters:
 *	 d_MU - Left orthogonal matrix on the device
 *	d_MVT - Right orthogonal matrix on the device
 *      Sigma - Singular Values of A on the CPU
 * *****************************************************************************/

bool cusvd(int M, int N, float *d_A, float *d_MU, float *d_MVT, double *Sigma);

#endif
