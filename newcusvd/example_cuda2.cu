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
* Tested on CUDA 2.0                                                                            *
************************************************************************************************/

#ifndef _EXAMPLE_CU_
#define _EXAMPLE_CU_

#include "example.h"

//Include the below file in your main program
#include "cusvd.cu"

float *initialize(int ind)
{
	int i = 0, j = 0, l = 0;
        float *temp = (float*)malloc(sizeof(float) * ind * ind);

        for(i=0 ; i < ind ; i++)
        {
                for(j=0 ; j < ind ; j++)
                {
                        if(i==j)
                                temp[l++] = 1;
                        else
                                temp[l++] = 0;
                }
        }
        return temp;
}

int main(int argc, char** argv)
{
	bool   result;
	double *Sigma; 

	//M>=N and M and N are a multiple of 32
	int M = 512, N = 512;
	float *A, *U, *VT, *d_A, *d_U, *d_VT;

	//Step 1 - Read A in column major order 
	A = (float*)malloc(sizeof(float) * M * N);

	FILE *fp = fopen("data", "r");
	for(i=0 ; i < M * N ; i++)
	{
		fscanf(fp,"%f", &A[i]);
	}
	fclose(fp);

	//Step 2
	Sigma = (double*)malloc(sizeof(double)*N);

	//Step 3
	CUT_DEVICE_INIT(argc, argv);
	status = cublasInit();

	if(status != CUBLAS_STATUS_SUCCESS)
	{
		fprintf(stderr, "Error in initialization");
		return EXIT_FAILURE;		
	}

	//Step 4
	
	status = CUDA_SAFE_CALL(cublasAlloc(M*N*sizeof(float), sizeof(float), (void**)&d_A));
	status = CUDA_SAFE_CALL(cublasAlloc(M*M*sizeof(float), sizeof(float), (void**)&d_U));	
	status = CUDA_SAFE_CALL(cublasAlloc(N*N*sizeof(float), sizeof(float), (void**)&d_VT));
	

	//Step 5
	U = initialize(M);
	VT = initialize(N);

	
	status = CUDA_SAFE_CALL(cublasSetMatrix(M, N, sizeof(float),  A, M,  d_A, M));
	status = CUDA_SAFE_CALL(cublasSetMatrix(M, N, sizeof(float),  U, M,  d_U, M));
	status = CUDA_SAFE_CALL(cublasSetMatrix(M, N, sizeof(float), VT, M, d_VT, M));
	

	//Step 6
	timer = 0;
	CUT_SAFE_CALL(cutCreateTimer(&timer));
	CUT_SAFE_CALL(cutStartTimer(timer));
	
	result = cusvd(M, N, d_A, d_U, d_VT, Sigma);

	CUT_SAFE_CALL(cutStopTimer(timer));
	printf("SVD processing time: %f (ms)\n", cutGetTimerValue(timer));	
	CUT_SAFE_CALL(cutDeleteTimer(timer));

/*	
 	printf("Copy and print VT matrix\n");	
	CUDA_SAFE_CALL(cudaMemcpy(VT, d_VT, sizeof(float)*N*N, cudaMemcpyDeviceToHost));
	for(int i=0; i < N; i++)
		for(int j=0; j < N; j++)
			printf("%f\n", check2[i*N+j]);
*/

	//Step 7
	free(A);
	CUDA_SAFE_CALL(cudaFree(d_A));
	CUDA_SAFE_CALL(cudaFree(d_U));
	CUDA_SAFE_CALL(cudaFree(d_VT));
	CUT_EXIT(argc, argv);
	return 0;
}
#endif
