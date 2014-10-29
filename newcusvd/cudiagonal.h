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

#ifndef _CUDIAGONAL_H_
#define _CUDIAGONAL_H_

/* Header files included */

#include <stdio.h>
#include <stdlib.h>

#include "cutil.h"
#include "device_functions.h"
#include "cublas.h"

/*****************************************************************************************
 * This function diagonalizes the bidiagonal matrix B, such that 
 * B = dU * Diag(Sigma) * dV(T). 
 * 
 * The elements of Sigma are the singular values of matrix A. 
 *
 * Input parameters: 
 * 	M, N       - Dimensions of original matrix A.
 * 	diagonal  - Diagonal elements of B on the CPU.
 * 	superdiag - Superdiagonal elements of B on the CPU.
 *
 * Output parameters:
 *	dU     - Left householder matrix on the device obtained after diagonalizing B.
 *	dV     - Right householder matrix on the device obtained after diagonalizing B.
 *	diagonal - Singular values of A obtained after diagonalizing B.
 * ***************************************************************************************/

bool cudiagonal(int M, int N, double *diagonal, double *superdiag, float *dU, float *dV);

#endif
