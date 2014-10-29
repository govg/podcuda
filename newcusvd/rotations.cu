/************************************************************************************************
* Implementing Singular Value Decomposition on GPU using CUDA using algorithm                   *
* given in IPDPS '09 paper "Singular Value Decomposition on GPU using CUDA"                     *
*                                                                                               *
* Copyright (c) 2009 International Institute of Information Technology, Hyderabad.              *
* All rights reserved.                                                                          *
*                                                                                               *
* Permission to use, copy, modify and distribute this software and its documentation for        *
* educational purpose is hereby granted without fee, provided that the above copyright          *
* notice and this permission notice appear in all copies of this software and that you do       *
* not sell the software.                                                                        *
*                                                                                               *
* THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR        *
* OTHERWISE.                                                                                    *
*                                                                                               *
* Created by Sheetal Lahabar.                                                                   *
************************************************************************************************/

/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer listed
in this license in the documentation and/or other materials
provided with the distribution.

- Neither the name of the copyright holders nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#include "cublas.h"
#include <cutil.h>
#include "rotations.h"

#ifndef _ROTATIONS_CU_
#define _ROTATIONS_CU_

#define xBS 64
#define BS 64
#define CS 192

/*************************************************************************
 * Kernel applies a series of elementary rotations to a matrix in the 
 * forward(top to bottom) direction.  
 * 
 * The algorithm multiplies the matrix by a sequence of rotation
 * transformations which are given by arrays d_c and d_s. The vectors d_c 
 * and d_s of length (M2-M1+1 or N2-N1+1) fit in the shared memory. 
 * 
 * Not the whole matrix but only a part of it is transformed (rows from M1
 * to M2 of dU or rows from N1 to N2 of dV). Only the elements of this submatrix are 
 * changed.
 * 
 * Input parameters:
 *     g_data      - Matrix to be transformed (dU or dV)
 *     d_c, d_s    - Coefficient vectors on the device
 * 
 * Output parameters:
 *     g_data      - Transformed matrix
 * 
 * Utility subroutine.
 *************************************************************************/

__global__ void forwardonce(float *g_data, float *d_c, float *d_s, int width, int height)
{
	__shared__ float a[BS];
	__shared__ float b[BS];

	__shared__ float c[CS];
	__shared__ float s[CS];

	int i=0;

	int tx   = threadIdx.x;
	int bx   = blockIdx.x;
	int Dimx = blockDim.x;

	float  ctemp=1, stemp=1;

	int    datacopy = CS / Dimx;

	int    ind      = __mul24(bx , Dimx) + tx;

	for(i=0; i < datacopy; i++)
	{
		c[ i * BS + tx ] = d_c[ i * BS + tx ];
		s[ i * BS + tx ] = d_s[ i * BS + tx ];
		__syncthreads();
	}

	a[tx]    = g_data[ind];
	__syncthreads();
	float a1, a2;
	float a1temp, a2temp;

	for(i=0 ; i < height-1; i++)
	{
		ctemp = c[i];
		stemp = s[i];

		b[tx] = g_data[ind + (i+1) * width];
		__syncthreads();

		a1temp = a[tx];
		a2temp = b[tx];

		a1     = ctemp * a1temp + stemp * a2temp;
		a2     = ctemp * a2temp - stemp * a1temp;
		__syncthreads();

		g_data[ind + (i) * width] = a1;
		__syncthreads();

		a[tx] = a2;
		__syncthreads();
	}

	g_data[ind + (i) * width] = a2;
	__syncthreads();
}

/****************************************************************************
 * Kernel applies a series of elementary rotations to a matrix in the 
 * backward(bottom to top) direction.  
 *
 * The algorithm multiplies the matrix by a sequence of rotation
 * transformations which are given by arrays d_c and d_s. The vectors d_c 
 * and d_s of length (M2-M1+1 or N2-N1+1) fit in the shared memory. 
 * 
 * Not the whole matrix but only a part of it is transformed (rows from M1
 * to M2 of dU and rows from N1 to N2 of dV). Only the elements of 
 * this submatrix are changed.
 * 
 * Input parameters:
 *    g_data      - Matrix to be transformed dU or dV
 *    d_c, d_s    - Coefficient vectors
 * 
 * Output parameters:
 *    g_data      - Transformed matrix
 *  
 * Utility subroutine.
 ****************************************************************************/

__global__ void backwardonce(float *g_data, float* d_c, float *d_s, int width, int height)
{
	__shared__ float a[BS];
	__shared__ float b[BS];

	__shared__ float c[CS];
	__shared__ float s[CS];

	int i=0;

	int tx   = threadIdx.x;
	int bx   = blockIdx.x;
	int Dimx = blockDim.x;

	float  ctemp=1, stemp=1;

	int    datacopy = CS / Dimx;

	int    ind      = __mul24(bx , Dimx) + tx;

	for(i=0; i < datacopy; i++)
	{
		c[i*BS + tx] = d_c[height - 2 - i*BS - tx];
		s[i*BS + tx] = d_s[height - 2 - i*BS - tx];
		__syncthreads();
	}

	a[tx]    = g_data[ind];
	__syncthreads();

	float a1, a2;
	float a1temp, a2temp;

	for(i=0 ; i < height-1; i++)
	{
		ctemp = c[i];
		stemp = s[i];

		b[tx] = g_data[ind - (i+1) * width];
		__syncthreads();

		a1temp = b[tx];
		a2temp = a[tx];

		a1     = ctemp * a1temp + stemp * a2temp;
		a2     = ctemp * a2temp - stemp * a1temp;
		__syncthreads();

		g_data[ind - i * width] = a2;
		__syncthreads();

		a[tx] = a1;
		__syncthreads();
	}

	g_data[ind - i * width] = a1;
	__syncthreads();

}

/*************************************************************************
 * Kernel applies a series of elementary rotations to a matrix in the 
 * backward(bottom to top) direction.  
 * 
 * The algorithm multiplies the matrix by a sequence of rotation
 * transformations which are given by arrays d_c and d_s. The vectors d_c 
 * and d_s of length (M2-M1+1 or N2-N1+1) do not fit in the shared memory.
 * They are loaded in batches. 
 * 
 * Not the whole matrix but only a part of it is transformed (rows from M1
 * to M2, rows from N1 to N2). Only the elements of this submatrix are 
 * changed.
 * 
 * Input parameters:
 * g_data      - Matrix to be transformed
 * d_c, d_s    - Coefficient vectors
 * 
 * Output parameters:
 * g_data is transformed.
 *  
 * Utility subroutine.
 ****************************************************************************/

__global__ void backward(float *g_data, float *d_c, float *d_s, int width, int height)
{
	__shared__ float a[BS];
	__shared__ float b[BS];

	__shared__ float c[CS];
	__shared__ float s[CS];

	int i=0, j=0, k=0;

	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int Dimx = blockDim.x;

	int elements = (512 - (Dimx * 2))/2;
	float ctemp, stemp;
	
	int ind = __mul24(bx, Dimx) + tx;
	float a1, a2, flag = 0;
	int calls;
	float a1temp, a2temp;

        int datacopy = CS / Dimx;

        if((height - 1) % elements == 0)
        {
                flag = 0;
        }
        else
        {
                flag = 1;
        }

        calls = (height - 1) / elements;

        a[tx] = g_data[ind];
        __syncthreads();

	for(i=0 ; i < calls ; i++)
        {
                for(k=0 ; k < datacopy ; k++)
                {
                        c[k * BS + tx] = d_c[ height - 2 - i * CS - k * BS - tx];
                        __syncthreads();
                        s[k * BS + tx] = d_s[ height - 2 - i * CS - k * BS - tx];
                        __syncthreads();
                }
                for(j=0; j < elements ; j++)
                {
                        ctemp = c[j];
                        stemp = s[j];

                        b[tx] = g_data[ind - i * CS * width - (j+1) * width];
                        __syncthreads();

                        a1temp = b[tx];
                        a2temp = a[tx];
                        __syncthreads();

                        a1 = ctemp * a1temp + stemp * a2temp;
                        a2 = ctemp * a2temp - stemp * a1temp;
                        __syncthreads();

                        g_data[ind - i * CS * width - j * width] = a2;
                        __syncthreads();

                        a[tx] = a1;
                        __syncthreads();
                }
        }
	if(flag == 0)
        {
                g_data[ind - (i-1) * CS * width - j * width] = a1;
                __syncthreads();
        }
        else
        {
                for(k=0; k < datacopy; k++)
                {
                        c[ k * BS + tx ] = d_c[ height - 2 - i * CS - k * BS - tx ];
                        __syncthreads();
                        s[ k * BS + tx ] = d_s[ height - 2 - i * CS - k * BS - tx ];
                        __syncthreads();
                }
                for( j = 0; j < height - 1 - (calls * elements); j++)
                {
                        ctemp = c[j];
                        stemp = s[j];

                        b[tx] = g_data[ind - i * CS * width - (j + 1) * width];
                        __syncthreads();

                        a1temp = b[tx];
                        a2temp = a[tx];
			__syncthreads();

                        a1     = ctemp * a1temp + stemp * a2temp;
                        a2     = ctemp * a2temp - stemp * a1temp;
                        __syncthreads();

                        g_data[ ind - i * CS * width - j * width ] = a2;
                        __syncthreads();

                        a[tx] = a1;
                        __syncthreads();
                }

                g_data[ind - i * CS * width - j * width] = a1;
                __syncthreads();
        }
}

/*************************************************************************
 * Kernel applies a series of elementary rotations to a matrix in the 
 * forward(top to bottom) direction.  
 * 
 * The algorithm multiplies the matrix by a sequence of rotation
 * transformations which are given by arrays d_c and d_s. The vectors d_c 
 * and d_s of length (M2-M1+1 or N2-N1+1) do not fit in the shared memory. 
 * They are loaded in batches to the shared memory.
 * 
 * Not the whole matrix but only a part of it is transformed (rows from M1
 * to M2, rows from N1 to N2). Only the elements of this submatrix are 
 * changed.
 * 
 * Input parameters:
 *   g_data      - Matrix to be transformed
 *   d_c, d_s    - Coefficient vectors
 *  
 * Output parameters:
 *   g_data is transformed.
 *   
 *   Utility subroutine.
 **************************************************************************/


__global__ void forward(float *g_data, float *d_c, float *d_s, int width, int height)
{
	__shared__ float a[BS];
	__shared__ float b[BS];

	__shared__ float c[CS];
	__shared__ float s[CS];

	int i=0, j=0, k=0;

	int tx    =  threadIdx.x;
	int bx    =  blockIdx.x;
	int Dimx  =  blockDim.x;

	int elements = (512 - (Dimx * 2))/2;

	float ctemp, stemp;

	int ind = __mul24(bx, Dimx) + tx;

	float a1, a2, flag=0;
	int calls = 0;
	float a1temp, a2temp;

	int datacopy = CS / Dimx;

	if((height - 1) % elements == 0)
	{
		flag = 0;
	}
	else
	{
		flag = 1;
	}

	calls = (height - 1) / elements;

	a[tx] = g_data[ind];
	__syncthreads();

	for(i=0 ; i < calls ; i++)
	{
		for(k=0 ; k < datacopy ; k++)
		{
			c[k * BS + tx] = d_c[i * CS + k * BS + tx];
			__syncthreads();
			s[k * BS + tx] = d_s[i * CS + k * BS + tx];
			__syncthreads();
		}

		for(j=0; j < elements ; j++)
		{
			ctemp = c[j];
			stemp = s[j];

			b[tx] = g_data[ind + i * CS * width + (j+1) * width];
			__syncthreads();
			a1temp = a[tx];
			a2temp = b[tx];
			__syncthreads();

			a1 = ctemp * a1temp + stemp * a2temp;
			a2 = ctemp * a2temp - stemp * a1temp;
			__syncthreads();

			g_data[ind + i * CS * width + j * width] = a1;
			__syncthreads();

			a[tx] = a2;
			__syncthreads();
		}
	}
	if(flag == 0)
	{
		g_data[ind + (i-1) * CS * width + j * width] = a2;
		__syncthreads();
	}
	else
	{
		for(k=0; k < datacopy; k++)
		{
			c[ k * BS + tx ] = d_c[ i * CS + k * BS + tx ];
			__syncthreads();
			s[ k * BS + tx ] = d_s[ i * CS + k * BS + tx ];
			__syncthreads();
		}
		for( j = 0; j < height - 1 - (calls * elements); j++)
		{
			ctemp = c[j];
			stemp = s[j];

			b[tx] = g_data[ind + i * CS * width + (j + 1) * width];
			__syncthreads();

			a1temp = a[tx];
			a2temp = b[tx];
			a1     = ctemp * a1temp + stemp * a2temp;
			a2     = ctemp * a2temp - stemp * a1temp;
			__syncthreads();

			g_data[ ind + i * CS * width + j * width ] = a1;
			__syncthreads();

			a[tx] = a2;
			__syncthreads();
		}

		g_data[ind + i * CS * width + j * width] = a2;
		__syncthreads();
	}
}

/***********************************************************************************
 * Application of a sequence of  elementary rotations to a matrix
 *
 * The algorithm post-multiplies the matrix by a sequence of rotation
 * transformations which is given by arrays C and S. Depending on the value
 * of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
 * rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.
 *
 * Not the whole matrix but only a part of it is transformed (rows from M1
 * to M2, columns from N1 to N2). Only the elements of this submatrix are changed.
 *
 * Input parameters:
 *     IsForward   -   the sequence of the rotation application.
 *         M1,M2   -   the range of rows to be transformed.
 *     N1, N2      -   the range of columns to be transformed.
 *     C,S         -   transformation coefficients.
 *                     Array whose index ranges within [1..N2-N1].
 *     WORK        -   working array whose index ranges within [M1..M2].
 *     mat         -   matrix to be transformed on the device    
 *     dC          -   temporary memory array on the device
 *     dd          -   temporaty memory array on the device
 *
 * Output parameters:
 *     mat         -   transformed matrix on the device
 *                                                                 
 * Utility subroutine.
 ************************************************************************************/

float* applyrotationsfromtheleft(bool isforward,
		int m1,
		int m2,
		int n1,
		int n2,
		const ap::real_1d_array& c,
		const ap::real_1d_array& s,
		ap::real_1d_array& work,
		float* mat, float* dC, float* dd)
{
	double ctemp;
	double stemp;

	if( m1>m2||n1>n2 )
	{
		return mat;
	}

	if( isforward )
	{
		if( n1!=n2 )
		{
			int i=0;
			int actlength = m2-m1+1;
			float *cs, *ss;

			dim3 threads;
			dim3 grid;

			cs = (float*)malloc(sizeof(float) * (actlength+xBS));
			ss = (float*)malloc(sizeof(float) * (actlength+xBS));

			threads.x = xBS;
			threads.y = 1;
		
			grid.x = (n2-n1+1)/xBS;
			grid.y = 1;

			for(i=1;i<m2-m1+1;i++)
                        {
                                ctemp = c(i);
                                stemp = s(i);

                                if(!(ctemp!=1 || stemp!=0))
                                {
                                        cs[i-1] = 1;
                                        ss[i-1] = 0;
                                }
                                else 
                                {
                                        cs[i-1]=ctemp;
                                        ss[i-1]=stemp;  
                                }
                        }       

			CUDA_SAFE_CALL(cudaMemcpy(dd, &ss[0], (m2-m1+xBS)*sizeof(float), cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL(cudaMemcpy(dC, &cs[0], (m2-m1+xBS)*sizeof(float), cudaMemcpyHostToDevice));  

			if(m2-m1 <= CS)
			{
				forwardonce<<<grid, threads>>>(&mat[m1*(n2-n1+1)], dC, dd, n2-n1+1, m2-m1+1);
			}
			else
			{ 
				forward<<<grid, threads>>>(&mat[m1*(n2-n1+1)], dC, dd, n2-n1+1, m2-m1+1);
			}

			/*	for(j = m1; j <= m2-1; j++)
				{
					ctemp = c(j-m1+1);
					stemp = s(j-m1+1);
					if( ctemp!=1 || stemp!=0 )
					{
						jp1 = j+1;
						ap::vmove(&work(n1), &a(jp1, n1), ap::vlen(n1,n2), ctemp);
						ap::vsub(&work(n1), &a(j, n1), ap::vlen(n1,n2), stemp);
						ap::vmul(&a(j, n1), ap::vlen(n1,n2), ctemp);
						ap::vadd(&a(j, n1), &a(jp1, n1), ap::vlen(n1,n2), stemp);
						ap::vmove(&a(jp1, n1), &work(n1), ap::vlen(n1,n2));
					}
				}
			*/
			return mat;
		}
		else
		{
			//Special Case
			printf("Special case :(\n");

			/*
  			  for(j = m1; j <= m2-1; j++)
			  {
				  ctemp = c(j-m1+1);
				  stemp = s(j-m1+1);
				  if( ctemp!=1||stemp!=0 )
				  {
					  temp = a(j+1,n1);
					  a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
					  a(j,n1) = stemp*temp+ctemp*a(j,n1);
			  	}
			  }
			*/
		}
		return mat;
	}
	else
	{
		if( n1!=n2 )
		{
			int i=0;
			int actlength = m2-m1+1;

			float *cs, *ss;
			cs = (float*)malloc(sizeof(float) * (actlength+xBS));
			ss = (float*)malloc(sizeof(float) * (actlength+xBS));

			dim3 threads;
			dim3 grid;

			grid.x = (n2-n1+1)/xBS;	
			grid.y = 1;

			threads.x = xBS;
			threads.y = 1;

			for(i=1;i<m2-m1+1;i++)
			{
				ctemp = c(i);
				stemp = s(i);

				if(!(ctemp!=1 || stemp!=0))
				{
					cs[i-1] = 1;
					ss[i-1] = 0;
				}
				else
				{
					cs[i-1]=ctemp;
					ss[i-1]=stemp;
				}

			}

			CUDA_SAFE_CALL(cudaMemcpy(dd, &ss[0], (m2-m1+xBS)*sizeof(float), cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL(cudaMemcpy(dC, &cs[0], (m2-m1+xBS)*sizeof(float), cudaMemcpyHostToDevice));

			if(m2-m1 <= CS)	
			{
				backwardonce<<<grid, threads>>>(&mat[ m2 * (n2-n1+1) ], dC, dd, n2-n1+1, m2-m1+1);
			}
			else
			{
				backward<<<grid, threads>>>(&mat[ m2 * (n2-n1+1) ], dC, dd, n2-n1+1, m2-m1+1);
			}

			/*for(j = m2-1; j >= m1; j--)
			  {
			  ctemp = c(j-m1+1);
			  stemp = s(j-m1+1);

			  if( ctemp!=1||stemp!=0 )
			  {
			  jp1 = j+1;
			  ap::vmove(&work(n1), &a(jp1, n1), ap::vlen(n1,n2), ctemp);
			  ap::vsub(&work(n1), &a(j, n1), ap::vlen(n1,n2), stemp);
			  ap::vmul(&a(j, n1), ap::vlen(n1,n2), ctemp);
			  ap::vadd(&a(j, n1), &a(jp1, n1), ap::vlen(n1,n2), stemp);
			  ap::vmove(&a(jp1, n1), &work(n1), ap::vlen(n1,n2));
			  }

			  }*/

			return mat;
		}
		else
		{
			printf("Special case :(\n");
			//Special Case
			/*
			   for(j = m2-1; j >= m1; j--)
			   {
			   	ctemp = c(j-m1+1);
	  		   	stemp = s(j-m1+1);
	    		        if( ctemp!=1||stemp!=0 )
			        {
				   temp = a(j+1,n1);
				   a(j+1,n1) = ctemp*temp-stemp*a(j,n1);
				   a(j,n1) = stemp*temp+ctemp*a(j,n1);
			        }
			   }
			 */

		}

	}
	return mat;
}


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    WORK        -   working array whose index ranges within [M1..M2].
    mat         -   matrix to be transformed on the device
    dC          -   temporary memory array on the device
    dd          -   temporary memory array on the device
Output parameters:
    mat         -   transformed matrix on the device.

Utility subroutine.
*************************************************************************/

float* applyrotationsfromtheright(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_1d_array& work, float* mat, float* dC, float* dd)
{
    return mat;

/*  int j;
    int jp1;
    double ctemp;
    double stemp;
    double temp;

    if( isforward )
    {
        if( m1!=m2 )
        {
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(work.getvector(m1, m2), a.getcolumn(jp1, m1, m2), ctemp);
                    ap::vsub(work.getvector(m1, m2), a.getcolumn(j, m1, m2), stemp);
                    ap::vmul(a.getcolumn(j, m1, m2), ctemp);
                    ap::vadd(a.getcolumn(j, m1, m2), a.getcolumn(jp1, m1, m2), stemp);
                    ap::vmove(a.getcolumn(jp1, m1, m2), work.getvector(m1, m2));
                }
            }
        }
        else
        {
	    //Special case
            for(j = n1; j <= n2-1; j++)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
    else
    {
        if( m1!=m2 )
        {
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    jp1 = j+1;
                    ap::vmove(work.getvector(m1, m2), a.getcolumn(jp1, m1, m2), ctemp);
                    ap::vsub(work.getvector(m1, m2), a.getcolumn(j, m1, m2), stemp);
                    ap::vmul(a.getcolumn(j, m1, m2), ctemp);
                    ap::vadd(a.getcolumn(j, m1, m2), a.getcolumn(jp1, m1, m2), stemp);
                    ap::vmove(a.getcolumn(jp1, m1, m2), work.getvector(m1, m2));
                }
            }
        }
        else
        {
	    //Special Case
            for(j = n2-1; j >= n1; j--)
            {
                ctemp = c(j-n1+1);
                stemp = s(j-n1+1);
                if( ctemp!=1||stemp!=0 )
                {
                    temp = a(m1,j+1);
                    a(m1,j+1) = ctemp*temp-stemp*a(m1,j);
                    a(m1,j) = stemp*temp+ctemp*a(m1,j);
                }
            }
        }
    }
    */
}


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/
void generaterotation(double f, double g, double& cs, double& sn, double& r)
{
    double f1;
    double g1;

    if( g==0 )
    {
        cs = 1;
        sn = 0;
        r = f;
    }
    else
    {
        if( f==0 )
        {
            cs = 0;
            sn = 1;
            r = g;
        }
        else
        {
            f1 = f;
            g1 = g;
            r = sqrt(ap::sqr(f1)+ap::sqr(g1));
            cs = f1/r;
            sn = g1/r;
            if( fabs(f)>fabs(g)&&cs<0 )
            {
                cs = -cs;
                sn = -sn;
                r = -r;
            }
        }
    }
}
#endif


