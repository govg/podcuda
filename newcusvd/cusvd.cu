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

#ifndef _CUDASVD_CU_
#define _CUDASVD_CU_

#include "cusvd.h"
#include "cubidiagonal.cu"
#include "cudiagonal.cu"

bool cusvd(int M, int N, float *d_A, float *d_MU, float *d_MVT, double *Sigma)
{
	int i = 0;
	bool result1, result2;

	float *dU, *dV;
	float *tempinitM, *tempinitN;
	double *diagonal, *superdiag;

	printf("Allocating resources, initializing variables\n");	

 	 diagonal =     (double*)malloc(sizeof(double)*N);
	superdiag = (double*)malloc(sizeof(double)*(N-1));

	tempinitM = (float*)malloc(sizeof(float)*M*M);
	tempinitN = (float*)malloc(sizeof(float)*N*N);

	for(i=0; i<M*M; i++)
		tempinitM[i]=0;

	for(i=0; i<N*N; i++)
		tempinitN[i]=0;

	for(i=0; i<M; i++)
		tempinitM[i*M+i]=1;

	for(i=0; i<N; i++)
		tempinitN[i*N+i]=1;

	CUDA_SAFE_CALL(cudaMalloc((void**)&d_Q, M*M*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_P, N*N*sizeof(float)));
     
        CUDA_SAFE_CALL(cudaMalloc((void**)&dU, M*M*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dV, N*N*sizeof(float)));

	CUDA_SAFE_CALL(cudaMemcpy(d_Q, tempinitM, M*M*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_P, tempinitN, N*N*sizeof(float), cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpy(dU, tempinitM, M*M*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dV, tempinitN, N*N*sizeof(float), cudaMemcpyHostToDevice));

 	result1 = cubidiagonal(M, N, d_A, d_Q, d_P, diagonal, superdiag);
  	result2 =   cudiagonal(M, N, diagonal, superdiag, dU, dV);

	cublasSgemm('n', 'n', M, M, M, 1, d_Q, M, dU, M, 0, d_MU, M);
        cublasSgemm('t', 'n', N, N, N, 1, dV, N, d_P, N, 0, d_MVT, N);

/*	
        float *fsigma = (float*)malloc(sizeof(float)*N);
	for(int i=0; i<N; i++)	
	{
		Sigma[i] = diagonal[i];
		fsigma[i] = (float)Sigma[i];
	}

	float *d_middle, *zero, *check2;	
	for(int i=0; i<M*M; i++)
		zero[i] = 0;

	zero = (float*)malloc(sizeof(float)*M*M);
	check2 = (float*)malloc(sizeof(float)*M*M);

        CUDA_SAFE_CALL(cudaMalloc((void**)&d_middle, M * M * sizeof(float)));
        CUDA_SAFE_CALL(cudaMemcpy(d_middle, zero, M * M * sizeof(float), cudaMemcpyHostToDevice));

        for(int i=0; i < M ;i++)
                CUDA_SAFE_CALL(cudaMemcpy(&d_middle[i*M+i], &fsigma[i], sizeof(float), cudaMemcpyHostToDevice));

        cublasSgemm('n', 'n', N, N, N, 1,     d_MU, N, d_middle, N, 0, d_middle, N);
        cublasSgemm('n', 'n', N, N, N, 1, d_middle, N,    d_MVT, N, 0, d_middle, N);

        CUDA_SAFE_CALL(cudaMemcpy(check2, d_middle, sizeof(float)*N*N, cudaMemcpyDeviceToHost));

        for(int i=0; i < M; i++)
                for(int j=0; j < N; j++)
                        printf("%f\n", check2[i*N+j]);
*/
	return result1 * result2;
}
#endif
