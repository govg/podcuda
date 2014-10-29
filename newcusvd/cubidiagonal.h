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

#ifndef _CUBIDIAGONAL_H_
#define _CUBIDIAGONAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/timeb.h>
#include <time.h>
#include <unistd.h>

#include <cutil.h>
#include "cublas.h"

/*********************************************************************
 * This function returns 1 if eta>=0 and -1 for eta <0.
 * ******************************************************************/

float sign(float eta);

/*********************************************************************
 * This function initializes the mxm matrix with identity.
 * ******************************************************************/

float* iden(int m);

/*********************************************************************
 * This function initializes the mxm matrix with zeroes.
 * ******************************************************************/

float * zeros(int m);

/**********************************************************************************
 * This function computes the bidiagonal matrix B from a dense matrix d_A. 
 *	It decomposes d_A such that d_A = d_Q * B * d_P
 *
 * Assumes that M>=N.
 *	
 * The diagonal and superdiagonal elements of B are copied in diagonal and 
 * superdiag vectors on the CPU. The householder vectors reside in d_Q 
 * and d_P on the device.
 *
 * Input parameters: 
 *	M, N - Matrix dimensions of d_A
 *	d_A  - Matrix A on the device
 *
 * Output parameters:
 * 	d_Q       - Left Householder matrix on the device after bidiagonalizing d_A.
 * 	d_P       - Right Householder matrix on the device after bidiagonalizing d_A.
 * 	diagonal  - Diagonal elements of B on the CPU.
 *      superdiag - Superdiagonal elements of B on the CPU.
 *
 * *********************************************************************************/

bool cubidiagonal(int M, int N, float *d_A, float *d_Q, float *d_P, double *diagonal, double *superdiag);

#endif
