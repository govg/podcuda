#########################################################################################
#											#
#	                                                                                #
#	Software implements Singular Value Decomposition on GPU using CUDA              #
#	     										#
#	Copyright (c) 2009 International Institute of Information Technology. 		#
#	All rights reserved.								#
#											#
#	Permission to use, copy, modify and distribute this software and its 		#
#	documentation for research purpose is hereby granted without fee, 		#
#	provided that the above copyright notice and this permission notice appear	#
#	in all copies of this software and that you do not sell the software.		#
#											#
#	THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND, 		#
#	EXPRESS, IMPLIED OR OTHERWISE.							#
#											#
#	Please report any issues to Sheetal Lahabar (sheetal@research.iiit.ac.in)	#
#											#
#	Please cite following paper, if you use this software for research purpose	#
#											#
#	"Singular Value Decomposition on GPU using CUDA"				#
#	Sheetal Lahabar and P. J. Narayanan						#
#	Proceedings of 23rd IEEE International Parallel and Distributed Processing 	#
#	Symposium, May 2009								#
#											#
#	The software uses the diagonalization routine for a bidiagonal matrix    	#
#	downloaded from									#	
#	http://www.alglib.net/matrixops/other/bdsvd.php. The copyright rights		#
#	are included in all the source files used.					#	
#											#
#											#
#											#
#########################################################################################

=========================================================================================

1. Using the Code in the Linux Environment 
	1. Download and untar cuSvd.tar from the website.
	   It contains cusvd directory.

2. Make sure you have CUDA installed on your system.
3. Include cusvd.cu in your main program.
   For more information look into example code.
4. make to compile the code
5. It has been tested on NVIDIA GTX 280, GTX 8800 and Tesla S1070

=========================================================================================

2. Example Code
	Example code is given as a help to understand how to use the code.
	1. example.cu is the main program.
	2. make to compile the example code.
	3. Executable file "cusvd" is created.

=========================================================================================

3. Brief description on how to use the code.

   Our code will solve the Singular Value Decomposition problem for large real dense 
   matrix A (MxN) with leading dimension M and M>=N and M>=64 and N>=64.
   It decomposes the given matrix A (MxN) as A = U * Sigma * V(T). 
   
   1) It assumes that the matrix dimensions are a multiple of 32. It the matrix dimensions 
      are not a multiple of 32 then pad the matrix with zeroes to the next multiple of 32.
   2) For a given matrix, M>=N.
   3) For M < N, Compute the SVD of A transpose.
   4) The individual timings excludes the cost of transfer from CPU to the GPU and 
      initializing the matrices. However, the SVD processing time includes the time for 
      transferring and initializing the matrix.

   Steps executed in example.cu
   Step 1. Allocate memory for matrix A (MxN)(Assumptions: M>=N, Dimensions are a multiple 
           of 32) on the CPU. Read the elements in column major format in A. 
   Step 2. Allocate memory for Sigma on the CPU.
   Step 3. Initialize the device. 
   Step 4. Allocate Device memory: d_A for A(MxN), d_U for U(MxM) and d_VT for VT(NxN).
   Step 5. Intialize d_A on the device and unitary matrices d_U and d_VT on device.
   Step 6. Function call cusvd() computes the SVD of A.    
   Step 7. Free the allocated device memory.
 
=========================================================================================

4. Overall Singular Value Decomposition algorithm
   
   1. Include cusvd.cu to your main program.  
  
   Invoke the following series of functions calls:
   2. Sigma = (double*)malloc(sizeof(double)*N);
   3. status = cublasAlloc(M*N*sizeof(float), sizeof(float), (void**)&d_A));
   4. status = cublasAlloc(M*M*sizeof(float), sizeof(float), (void**)&d_U));
   5. status = cublasAlloc(N*N*sizeof(float), sizeof(float), (void**)&d_VT));
   6. status = cublasSetMatrix(M, M, sizeof(float), identityM, M, d_U, M));
   7. status = cublasSetMatrix(N, N, sizeof(float), identityN, N, d_VT, N));
   8. status = cublasSetMatrix(M, N, sizeof(float), A, M, d_A, M));
   9. result = cusvd(M, N, d_A, d_U, d_VT, Sigma);
  10. cudaFree(d_A);
  11. cudaFree(d_U);
  10. cudaFree(d_VT);

=========================================================================================

5. Using the Code in Windows Environment 

   The step process required for using the code on windows remains same as in linux 
   environment. 

========================================================================================= 

