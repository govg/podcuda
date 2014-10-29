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
#ifndef _CUBIDIAGONAL_CU_
#define _CUBIDIAGONAL_CU_

#include "cubidiagonal.h"

int l = 0, m, n, count = 0;
int j = 0, k = 0, Nb, size = 0, i = 0, kMax = 0, var = 0;
int time1, time2, time3, time4;

struct timeb tr;

//CPU variables
float *d_Q, *d_P;
float *Wmat, *Zmat;
float *e1, *idenU, *idenV, *In;

//Device variables
float *d_Zmat, *d_Wmat, *d_Umat, *d_Vmat;
float *d_w1, *d_z1, *d_v1, *d_x1;

float *d_U, *d_V;
float *d_qw, *d_pw;

float *d_Aor, *d_W;
float *d_e1, *d_zinput;
float *d_iden, *iden1, *d_temp, *d_zero;

float dot, alpha = 0 , beta = 0;
float one = 1, s, t, d_eta, alpha1 = -1;
float norm = 0, eta, v = 0, g = -1, sigma = 0, value = 0;

cublasStatus status1;

float sign(float eta)
{
	if(eta >= 0)
		return 1;
	else
		return -1;
}

float *iden(int m)
{
	int i = 0, j = 0, l = 0;
	float *a = (float*)malloc(sizeof(float) * m * m);

	for(i=0 ; i < m ; i++)
	{
		for(j=0 ; j < m ; j++)
		{
			if(i==j)
				a[l++]=1;
			else
				a[l++]=0;
		}
	}	
	return a;
}

float *zeros(int m)
{
	int i = 0;
	float *a = (float*)malloc(sizeof(float) * m);

	for(i=0 ; i < m ; i++)
		a[i] = 0;
	return a;
}

bool cubidiagonal(int M, int N, float *d_A, float *d_Q1, float *d_P1, double *diagonal, double *superdiag)
{
	bool result = 1;
 
	m = M;
	n = N;

//	if(n <= 128)
//		Nb = 1;
//	else
		Nb = 16;

	kMax = n/Nb;	

	idenU  = (float*)malloc(sizeof(float) * m * m);
	idenV  = (float*)malloc(sizeof(float) * n * n);

	Wmat   = (float*)malloc(sizeof(float) * m * m);
	Zmat   = (float*)malloc(sizeof(float) * n * Nb);

        In     = (float*)malloc(sizeof(float) * n * n);

	idenU  = iden(m);
	idenV  = iden(n);

	Wmat   = zeros(m * m);
	Zmat   = zeros(n * Nb);
	In     = zeros(n * n);
	
	if(m > n)
		size = m;
	else 
		size = n;

	iden1  = (float*)malloc(sizeof(float) * size);
	for(i=0; i < size ; i++)
		iden1[i] = 1;

	e1 = zeros(size);

	status1 = cublasAlloc(m*Nb, sizeof(d_Wmat[0]),   (void**)&d_Wmat);
	status1 = cublasAlloc(n*Nb, sizeof(d_Zmat[0]),   (void**)&d_Zmat);

	status1 = cublasAlloc(m*Nb,   sizeof(d_qw[0]),   (void**)&d_qw);
	status1 = cublasAlloc(n*Nb,   sizeof(d_pw[0]),   (void**)&d_pw);

	status1 = cublasAlloc( m*m,    sizeof(d_U[0]),    (void**)&d_U);
	status1 = cublasAlloc( n*n,    sizeof(d_V[0]),    (void**)&d_V);

	status1 = cublasAlloc( m*m, sizeof(d_Umat[0]),    (void**)&d_Umat);
	status1 = cublasAlloc( n*n, sizeof(d_Vmat[0]),    (void**)&d_Vmat);

	status1 = cublasAlloc( m*m,    sizeof(d_Q[0]),    (void**)&d_Q);
	status1 = cublasAlloc( n*n,    sizeof(d_P[0]),    (void**)&d_P);
	
	status1 = cublasAlloc(size,   sizeof(d_e1[0]),    (void**)&d_e1);
	status1 = cublasAlloc(size,   sizeof(d_w1[0]),    (void**)&d_w1);
	status1 = cublasAlloc(size,   sizeof(d_v1[0]),    (void**)&d_v1);
	status1 = cublasAlloc(size,   sizeof(d_x1[0]),    (void**)&d_x1);
	status1 = cublasAlloc(size,   sizeof(d_z1[0]),    (void**)&d_z1);
	status1 = cublasAlloc(size, sizeof(d_iden[0]),    (void**)&d_iden);
	status1 = cublasAlloc(size, sizeof(d_zero[0]),    (void**)&d_zero);
	status1 = cublasAlloc(size, sizeof(d_temp[0]),    (void**)&d_temp);
	status1 = cublasAlloc(size, sizeof(d_temp[0]),    (void**)&d_zinput);
	status1 = cublasAlloc(size, sizeof(d_temp[0]),    (void**)&d_Aor);
	status1 = cublasAlloc(size,    sizeof(d_W[0]),    (void**)&d_W);
	
	status1 = cublasSetMatrix(m, m, sizeof(d_U[0]), Wmat, m, d_U, m);
	status1 = cublasSetMatrix(n, n, sizeof(d_V[0]), In, n, d_V, n);

	status1 = cublasSetMatrix(m, m, sizeof(d_Q[0]), idenU, m, d_Q, m);
	status1 = cublasSetMatrix(n, n, sizeof(d_P[0]), idenV, n, d_P, n);

	status1 = cublasSetMatrix(m, Nb, sizeof(d_qw[0]), Wmat , m, d_qw, m);
	status1 = cublasSetMatrix(Nb, n, sizeof(d_pw[0]), Zmat, Nb, d_pw, Nb);

	status1 = cublasSetMatrix(m, Nb, sizeof(d_Wmat[0]), Wmat, m, d_Wmat, m);
	status1 = cublasSetMatrix(Nb, n, sizeof(d_Zmat[0]), Zmat, Nb, d_Zmat, Nb);

	status1 = cublasSetMatrix(m, m, sizeof(d_Umat[0]), Wmat, m, d_Umat, m);
	status1 = cublasSetMatrix(n, n, sizeof(d_Vmat[0]), In, n, d_Vmat, n);

	status1 = cublasSetVector(size,   sizeof(d_e1[0]),    e1, 1,   d_e1, 1);
	status1 = cublasSetVector(size, sizeof(d_iden[0]), iden1, 1, d_iden, 1);
	status1 = cublasSetVector(size, sizeof(d_zero[0]),  Wmat, 1, d_zero, 1);
	status1 = cublasSetVector(size, sizeof(d_zero[0]),  Wmat, 1,  d_Aor, 1);
	status1 = cublasSetVector(size, sizeof(d_zero[0]),  Wmat, 1, d_temp, 1);
	
	printf("Intialization Complete\n");
	printf("Computing bidiagonal matrix B\n");

	ftime(&tr);
	time1 = tr.time;
	time2 = tr.millitm;

	for(i=0; i < kMax ; i++)
	{
		status1 = cublasSetMatrix( m, Nb, sizeof(d_Wmat[0]), Wmat, m, d_Wmat, m );
		status1 = cublasSetMatrix( Nb, n, sizeof(d_Zmat[0]), Zmat, Nb, d_Zmat, Nb );

		   norm = cublasSnrm2( m - (Nb * i), &d_A[ m * Nb * i + Nb * i ] , 1 );
		status1 = cublasGetVector( 1, sizeof(float), &d_A[ m * Nb * i + Nb * i ], 1, &eta,1 );

		v = sign(eta) * norm;
		var = -1 * sign(eta);

		if(v != 0)
			sigma = (eta + v)/v;
		else 
			sigma = 0;

		value = v / (eta + v);
		status1 = cublasSetVector( 1, sizeof(float), &value, 1, &d_Umat[ m*Nb*i + Nb*i ], 1 );
		cublasSaxpy( m - (Nb*i), 1/(eta+v), &d_A[ m * Nb *i + Nb * i ], 1, &d_Umat[ m * Nb * i + Nb * i ], 1 );

		//Q---
		cublasSgemm( 'n', 'n', m, 1,  m - (Nb * i), sigma, &d_Q[ m * Nb * i ], m, &d_Umat[ m * Nb * i + Nb * i ], m, beta, &d_qw[0], m );
                cublasSgemm( 'n', 'n', m, 1, 1, -1, &d_qw[0], m, &d_Umat[ m * Nb * i + Nb * i ], 1, 1, &d_Q[ m * Nb * i ] , m );
		//Q---

		//value = norm * g;
		value = norm * var;
		status1 = cublasSetVector( 1 , sizeof(float), &value, 1, &d_A[ m * Nb * i + Nb * i ], 1 );

		cublasSgemm( 't','n', n - (Nb*i+1), 1, m - (Nb*i), -1 * sigma, &d_A[ m * Nb * i + Nb * i + m ], m, &d_Umat[ m * Nb * i + Nb * i ], m, beta, &d_v1[0], size );
		cublasSaxpy( n - (Nb*i+1), 1, &d_A[ m * Nb * i + Nb*i + m ], m, &d_v1[0], 1 );
		cublasScopy( n - (Nb*i+1), &d_v1[0], 1, &d_temp[0], 1 );	
		
		norm = cublasSnrm2( n - (Nb*i+1), &d_v1[0], 1 );
		status1 = cublasGetVector( 1, sizeof(float), &d_v1[0], 1, &e1[0], 1 );
		eta = e1[0];

		cublasSgemm( 't','t', n - (Nb*i+1), 1, 1 ,1, &d_A[ m * Nb * i + Nb*i + m ], m, &d_iden[0], n, -1, &d_v1[0], size );	
		cublasScopy( n - (Nb*i+1), &d_v1[0], 1, &d_Zmat[ Nb * i * Nb + Nb ], Nb );

		v = sign(eta) * norm;
		var = -1 * sign(eta);	

		if(v !=0)
			sigma = (eta + v)/v;
		else
			sigma = 0;
	
		value = v / (eta + v);
		status1 = cublasSetVector( 1 , sizeof(float), &value, 1, &d_Vmat[ n * Nb * i + Nb * i + n ], n );
		cublasSaxpy( n - (Nb*i+1), 1/(eta+v), &d_temp[0], 1, &d_Vmat[ n * Nb * i + Nb * i + n ], n );

		status1 = cublasSetVector( 1 , sizeof(float), &value, 1, &d_e1[0], 1 );
		cublasSaxpy( n-(Nb*i+1), 1/(eta+v), &d_temp[0], 1, &d_e1[0], 1 );

		//P---
		cublasSgemm( 't', 't', n, 1, n - (Nb*i+1), sigma, &d_P[ Nb * i + 1 ], n, &d_Vmat[ n * Nb * i + Nb * i + n ], n, beta, &d_W[0], size );
                cublasScopy( n, &d_W[0], 1, &d_pw[0], Nb );
                cublasSgemm( 'n', 't', 1, n, 1, -1, &d_Vmat[ n * Nb * i + Nb * i + n ], n, &d_W[0], size, 1, &d_P[ Nb * i +1 ], n ); 				//P---
		
		//value = norm * g;
		value = norm * var;

		cublasSgemm( 'n', 't',  m - (Nb*i), 1, n - (Nb*i+1), sigma, &d_A[ m * Nb * i + Nb * i + m ], m , &d_Vmat[ n * Nb * i + Nb * i + n ], n, beta, &d_Wmat[ Nb * i ], m );

		if(m * Nb * i + Nb * i + m < (m * n))
			status1 = cublasSetVector( 1, sizeof(float), &value, 1, &d_A[ m * Nb * i + Nb * i + m ], 1 );

		dot = cublasSdot ( n - (Nb*i+1), &d_Zmat[ Nb * i * Nb + Nb ], Nb, &d_Vmat[ n * Nb * i + Nb * i + n ], n );
		dot = dot * sigma;
		cublasSaxpy( n - (Nb*i+1), -1 * dot, &d_Vmat[ n * Nb * i + Nb * i + n ], n, &d_Zmat[ Nb * i * Nb + Nb ], Nb );

		alpha1 = -1;

		for(k=1 ; k < Nb ; k++)
		{
			cublasScopy( size - (Nb * i) , &d_zero[0], 1, &d_e1[0], 1 );			

			cublasSgemm( 'n', 'n', m - (Nb*i+k), 1, k, -1, &d_Umat[ Nb * m * i + Nb * i + k ], m, &d_Zmat[ Nb * i * Nb  + k * Nb ], Nb, 1, &d_A[ m * Nb * i + Nb * i + k * m + k ], m );
			cublasSgemm( 'n', 'n', m - (Nb*i+k), 1, k, -1, &d_Wmat[ i * Nb + k ], m, &d_Vmat[ n * Nb * i + Nb * i  + k * n ], n, 1, &d_A[ m * Nb * i + Nb * i + k * m + k ], m );
		
			norm = cublasSnrm2( m - (Nb*i+k) , &d_A[ m * Nb * i + Nb * i + k * m + k ] , 1 );
			status1 = cublasGetVector( 1, sizeof(float), &d_A[ m * Nb * i + Nb * i + k * m + k ], 1, &eta, 1 );

			v = sign(eta) * norm;
			var = -1 * sign(eta);

			if(v!=0)
				sigma = (eta + v)/v;
			else
				sigma = 0;

			value = v / (eta + v);
			status1 = cublasSetVector( 1 , sizeof(float), &value, 1, &d_e1[0], 1 );
			cublasSaxpy( m - (Nb*i+k), 1/(eta+v), &d_A[ m * Nb * i + Nb * i + k * m + k ], 1, &d_e1[0], 1 );

			//value = norm * g;
			value = norm * var;

			status1 = cublasSetVector( 1, sizeof(float), &value, 1, &d_A[ m * Nb * i + Nb * i + k * m + k ], 1 );
			cublasScopy( m - (Nb*i+k), &d_e1[0], 1, &d_Umat[ m * Nb * i + Nb * i + k * m + k ], 1 );
			
			//Q---
			cublasSgemm( 'n', 'n', m, 1, m - (Nb*i+k), sigma, &d_Q[ m * Nb * i + k * m ], m, &d_e1[0], size, beta, &d_W[0], size );
                        cublasSgemm( 'n', 't', m, m - (Nb*i+k), k, 1, &d_qw[0], m, &d_Umat[ Nb * i * m + Nb * i + k ], m, beta, &d_U[0], m );
                        cublasSgemm( 'n', 'n', m, 1, m - (Nb*i+k), -sigma, &d_U[0], m, &d_e1[0], size, 1, &d_W[0], size );

                        cublasScopy( m, &d_W[0], 1, &d_qw[ k * m ], 1 );

                        cublasSgemm( 'n', 't', m, 1, k, -1, &d_qw[0], m, &d_Umat[ m * Nb * i + Nb * i + k ], m, 1, &d_Q[ m * Nb * i + m * k ], m );
                        cublasSgemm( 'n', 't', m, 1, 1, -1, &d_W[0], size, &d_e1[0], size, 1, &d_Q[ m * Nb * i + m * k ], m );	
			//Q---
			
			cublasScopy( n - (Nb*i+k+1), &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m, &d_temp[0], 1 );

			cublasSgemm( 'n', 'n', 1, n - (Nb*i+k+1), k, -1, &d_Umat[ Nb * m * i + Nb * i + k ], m, &d_Zmat[ Nb * Nb  * i  + (k+1) * Nb ], Nb, 1, &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m ); 
			cublasSgemm( 'n', 'n', 1, n - (Nb*i+k+1), k ,-1, &d_Wmat[ Nb * i + k ], m, &d_Vmat[ n * Nb * i + Nb * i  + (k+1) * n ], n, 1, &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m );

			float s1 = 1 - sigma; 

			cublasScopy( n - (Nb*i+k+1), &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m, &d_Aor[0], 1 );

			cublasSgemm( 't', 'n', k, 1, m - (Nb*i+k+1), sigma, &d_Umat[ Nb * m * i + Nb * i + k + 1 ], m, &d_e1[1], size, 0, &d_U[0], m );
			cublasSgemm( 't', 'n', n - (Nb*i+k+1), 1, k, 1, &d_Zmat[ Nb * Nb * i + (k+1) * Nb ], Nb, &d_U[0], m, 0, &d_V[0], n );	
			cublasSgemm( 't', 'n', k, 1, m - (Nb*i+k+1), sigma, &d_Wmat[ Nb * i + k + 1 ], m, &d_e1[1], size, 0, &d_U[0], m );	
			cublasSgemm( 't', 'n', n - (Nb*i+k+1), 1, k, 1, &d_Vmat[ n * Nb * i + Nb * i + (k+1)* n ], n, &d_U[0], m, 1, &d_V[0], n );

			cublasSaxpy ( n - (Nb*i+k+1), s1, &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m, &d_V[0], 1 );

			cublasScopy( n - (Nb*i+k+1), &d_temp[0], 1, &d_A[ m * Nb * i + Nb * i + k * m + k + m ], m );
			cublasSgemm( 't', 't', n - (Nb*i+k+1), 1, m - (Nb*i+k+1), -1 * sigma, &d_A[ m * Nb * i + Nb * i + (k + 1) * m + k + 1 ], m, &d_e1[1], 1, beta, &d_v1[0], size );

			cublasSaxpy( n - (Nb*i+k+1), 1, &d_v1[0], 1, &d_V[0], 1 );
			cublasScopy( n - (Nb*i+k+1), &d_V[0], 1, &d_temp[0], 1 );

			norm = cublasSnrm2( n-(Nb*i+k+1), &d_temp[0], 1 );                           
			status1 = cublasGetVector( 1, sizeof(float), &d_temp[0], 1, &e1[0], 1 );     
			eta = e1[0];                                                               
               											   
			v = sign(eta) * norm;
			var = -1 * sign(eta);

			if(v!=0)
				sigma = (eta + v)/v;
			else
				sigma = 0;

			value = v / (eta + v);

			cublasScopy( size - (Nb * i), &d_zero[0], 1, &d_e1[0],  1 );			
			status1 = cublasSetVector( 1, sizeof(float), &value, 1, &d_e1[0], 1 );
			cublasSaxpy( n-(Nb*i+k+1), 1/(eta+v), &d_temp[0], 1, &d_e1[0], 1 );

			//P---
			if(Nb + 1 + Nb * (i-1) + k + 1 <= n)
			{
				cublasSgemm( 't', 'n', n, 1, n - (Nb*i+k+1), sigma, &d_P[ Nb * i + k + 1], n, &d_e1[0], size, beta, &d_W[0], size );
                                cublasSgemm( 't', 'n', n, n - (Nb*i+k+1), k, 1, &d_pw[0], Nb, &d_Vmat[ Nb * i * n + Nb * i + k * n + n ], n, beta, &d_V[0], n );
                                cublasSgemm( 'n', 'n', n, 1, n - (Nb*i+k+1), -sigma, &d_V[0], n, &d_e1[0], size, 1, &d_W[0], size );

                                cublasScopy( n, &d_W[0], 1, &d_pw[k], Nb );

                                cublasSgemm( 't', 'n', 1, n, k, -1, &d_Vmat[ n * Nb * i + Nb * i + k * n + n ], n, &d_pw[0], Nb, 1, &d_P[ Nb * i + k + 1 ], n );
                                cublasSgemm( 't', 't', 1, n, 1, -1, &d_e1[0], size, &d_W[0], size, 1, &d_P[ Nb * i + k + 1 ], n );
			}
			//P---

			
			cublasSgemm( 'n', 't', m - (Nb*i+k+1), 1, n - (Nb*i+k+1), 1, &d_A[ m * Nb * i + Nb*i + (k+1) * m + k + 1 ], m, &d_e1[0] , 1, beta, &d_Wmat[ Nb * i + k * m + k + 1 ], m );
			cublasSgemm( 'n', 't', k, 1, n - (Nb*i+k+1), 1, &d_Zmat[ Nb * Nb * i  + (k + 1)* Nb ], Nb, &d_e1[0], 1, 0, &d_U[0], m );
			cublasSgemm( 'n', 'n', m - (Nb*i+k+1), 1, k, -sigma, &d_Umat[ m * Nb * i + Nb * i + k + 1 ], m, &d_U[0], m, sigma, &d_Wmat[ Nb * i + k * m + k + 1 ], m );

			dot  =  sigma * cublasSdot ( n-(Nb*i+k+1), &d_Aor[0], 1, &d_e1[0], 1 );

			cublasSgemm( 'n', 't', k, 1, n - (Nb*i+k+1), 1, &d_Vmat[ n * Nb * i  + Nb * i + (k + 1) * n ], n, &d_e1[0], 1, 0, &d_U[0], m );
			cublasSgemm( 'n', 'n', m - (Nb*i+k+1), 1, k, -sigma, &d_Wmat[ Nb * i + k + 1 ], m, &d_U[0], m, 1, &d_Wmat[ Nb * i + k * m + k + 1 ], m );


			status1 = cublasSetVector( 1, sizeof(float), &dot, 1, &d_Wmat[ Nb * i + k * m + k ], 1 );
			cublasSgemm( 'n', 'n', n - (Nb*i+k+1), 1, 1, 1, &d_Aor[0], size, &d_iden[0], 1, -1, &d_temp[0], size ) ;


			cublasScopy( n - (Nb*i+k+1), &d_temp[0], 1, &d_Zmat[ Nb * i * Nb +  Nb * (k+1)  + k ], Nb );		 			     		cublasScopy( n - (Nb*i+k+1), &d_e1[0], 1, &d_Vmat[ n * Nb * i + Nb * i + (k + 1) * n + k ], n );
                        
			//value = norm * g;
			value = norm * var;
						
			if( m * Nb * i + Nb * i + (k + 1) * m + k < (m * n) ) 
				status1 = cublasSetVector( 1 , sizeof(float), &value, 1, &d_A[ m * Nb * i + Nb * i + (k + 1)* m + k ], 1 );

			dot = cublasSdot ( n - (Nb*i+k+1), &d_Zmat[ Nb * i * Nb +  Nb * (k+1)  + k ], Nb, &d_e1[0], 1 );
			dot = dot * -sigma;
			e1[0] = dot;

			cublasSaxpy( n - (Nb*i+k+1), e1[0], &d_e1[0], 1, &d_Zmat[ Nb * i * Nb +  Nb * (k+1)  + k ], Nb );
		}

		//Update A
		cublasSgemm( 'n', 'n', m - (Nb * (i+1)), n - (Nb * (i+1)), Nb, -1, &d_Umat[ Nb * i * m + Nb * i + Nb ], m, &d_Zmat[ Nb * Nb * (i+1) ], Nb, 1, &d_A[ m * Nb * (i+1)  + Nb * (i+1) ], m ) ;	
		cublasSgemm( 'n', 'n', m - (Nb * (i+1)), n - (Nb * (i+1)), Nb, -1, &d_Wmat[ Nb * (i + 1) ], m, &d_Vmat[ Nb * n * i + n * Nb  + Nb * i ], n, 1, &d_A[ m * Nb * (i+1)  + Nb * (i+1) ], m );

		//Update Q
		cublasSgemm( 'n', 't', m, m - (Nb * (i+1)), Nb, -1, &d_qw[0], m, &d_Umat[ m * Nb * i + Nb * i + Nb ], m, 1, &d_Q[ m * Nb * (i+1) ], m );
		//Update P
                if( n - (Nb+1) - Nb * (i-1) -1  >= Nb )                 
                        cublasSgemm( 't', 'n', n - (Nb * (i+1) + 1), n, Nb, -1, &d_Vmat[ n * Nb * (i+1) + Nb * i + n ], n, &d_pw[0], Nb, 1, &d_P[ Nb * (i+1) + 1], n );
      	}
	
	cudaThreadSynchronize();
	ftime(&tr);
	time3 = tr.time;
	time4 = tr.millitm;

	printf("Bidiagonalization complete: Time required %d %d \n", time3 - time1, time4 - time2);

	float *tempA;
	tempA = (float*)malloc(sizeof(float)*m*n);

	status1 = cublasGetMatrix(m, n, sizeof(float), &d_A[0], m, tempA, m);

	//Copy diagonal elements to CPU
        for(j=0; j < n ; j++)
                diagonal[j] = (double)tempA[(m+1)*j];

	//Copy superdiagonal elements to CPU
	for(j=0; j < n-1; j++)
		superdiag[j] = (double)tempA[m + (m+1)*j];

	CUDA_SAFE_CALL(cudaFree(d_w1));
	CUDA_SAFE_CALL(cudaFree(d_z1));
	CUDA_SAFE_CALL(cudaFree(d_v1));
	CUDA_SAFE_CALL(cudaFree(d_x1));
	CUDA_SAFE_CALL(cudaFree(d_e1));

	CUDA_SAFE_CALL(cudaFree(d_Umat));
	CUDA_SAFE_CALL(cudaFree(d_Vmat));
	CUDA_SAFE_CALL(cudaFree(d_Zmat));
	CUDA_SAFE_CALL(cudaFree(d_Wmat));

	CUDA_SAFE_CALL(cudaFree(d_U));
	CUDA_SAFE_CALL(cudaFree(d_V));

	CUDA_SAFE_CALL(cudaFree(d_qw));
	CUDA_SAFE_CALL(cudaFree(d_pw));

	CUDA_SAFE_CALL(cudaFree(d_Aor));
	CUDA_SAFE_CALL(cudaFree(d_W));
	CUDA_SAFE_CALL(cudaFree(d_iden));
	CUDA_SAFE_CALL(cudaFree(d_zero));
	CUDA_SAFE_CALL(cudaFree(d_temp));
	CUDA_SAFE_CALL(cudaFree(d_zinput));

	return result;
}
#endif 
