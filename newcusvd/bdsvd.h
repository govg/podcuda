/************************************************************************************************
* Implementing Singular Value Decomposition on GPU using CUDA using algorithm                   *
* given in IPDPS '09 paper "Singular Value Decomposition on GPU using CUDA"                     *
*                                                                                               *
* Copyright (c) 2009 International Institute of Information Technology, Hyderabad.              *
* All rights reserved.                                                                          *
*                                                                                               *
* Permission to use, copy, modify and distribute this software and its documentation for        *
* educational purpose is hereby granted without fee, provided that the above copyright         *
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

#ifndef _BDSVD_H_
#define _BDSVD_H_

#include "ap.h"
#include <sys/timeb.h>
#include <time.h>

#include "rotations.h"

/*************************************************************************
Singular value decomposition of a bidiagonal matrix (extended algorithm)

The algorithm performs the singular value decomposition  of  a  bidiagonal
matrix B (upper or lower) representing it as B = Q*S*P^T, where Q and  P -
orthogonal matrices, S - diagonal matrix with non-negative elements on the
main diagonal, in descending order.

The  algorithm  finds  singular  values.  In  addition,  the algorithm can
calculate  matrices  Q  and P.

The feature of the algorithm is its ability to find  all  singular  values
including those which are arbitrarily close to 0  with  relative  accuracy
close to  machine precision. If the parameter IsFractionalAccuracyRequired
is set to True, all singular values will have high relative accuracy close
to machine precision. If the parameter is set to False, only  the  biggest
singular value will have relative accuracy  close  to  machine  precision.
The absolute error of other singular values is equal to the absolute error
of the biggest singular value.

Input parameters:
    D       -   main diagonal of matrix B.
                Array whose index ranges within [0..N-1].
    E       -   superdiagonal (or subdiagonal) of matrix B.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix B.
    IsUpper -   True, if the matrix is upper bidiagonal.
    IsFractionalAccuracyRequired -
                accuracy to search singular values with.
    NRU     -   number of rows in matrix U.
    NCVT    -   number of columns in matrix VT.
    dU      -   Matrix Q on the GPU.
    dV      -   Matrix P on the GPU.
    dC      -   Temporary vector on the device.
    dd      -   Temporary vector on the device.

Output parameters:
    D       -   singular values of matrix B in descending order.
    dU      -   if NRU>0, contains matrix Q, i.e. Left householder matrix 
                after diagonalization.
    dV      -   if NCVT>0, contains matrix PT, i.e. Right householder matrix
                after diagonalization.
Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged (rare case).

Additional information:
    The type of convergence is controlled by the internal  parameter  TOL.
    If the parameter is greater than 0, the singular values will have
    relative accuracy TOL. If TOL<0, the singular values will have
    absolute accuracy ABS(TOL)*norm(B).
    By default, |TOL| falls within the range of 10*Epsilon and 100*Epsilon,
    where Epsilon is the machine precision. It is not  recommended  to  use
    TOL less than 10*Epsilon since this will  considerably  slow  down  the
    algorithm and may not lead to error decreasing.
History:
    * 31 March, 2007.
        changed MAXITR from 6 to 12.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1999.
*************************************************************************/
/*bool rmatrixbdsvd(ap::real_1d_array& d,
     ap::real_1d_array e,
     int n,
     bool isupper,
     bool isfractionalaccuracyrequired,
     ap::real_2d_array& u,
     int nru,
     ap::real_2d_array& c,
     int ncc,
     ap::real_2d_array& vt,
     int ncvt, float* dU,float* dV, float* CC, float** values, float* C, float* dC, float* dd);
*/

bool rmatrixbdsvd(ap::real_1d_array& d,
     ap::real_1d_array e,
     int n,
     bool isupper,
     bool isfractionalaccuracyrequired,
     int nru,
     int ncvt, 
     float *dU, 
     float *dV, 
     float *dC, 
     float *dd);


/*************************************************************************
Obsolete 1-based subroutine. See RMatrixBDSVD for 0-based replacement.

History:
    * 31 March, 2007.
        changed MAXITR from 6 to 12.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1999.
*************************************************************************/

bool bidiagonalsvddecomposition(ap::real_1d_array& d,
     ap::real_1d_array e,
     int n,
     bool isupper,
     bool isfractionalaccuracyrequired,
     int ustart,
     int nru,
     int cstart,
     int ncc, 
     int vstart,
     int ncvt, 
     float *dU, 
     float *dV, 
     float *dC, 
     float *dd);

/*************************************************************************
 * Kernel returns -1*a for the input vector a.
 * **********************************************************************/

__global__ void positive(float *a);

/*************************************************************************
 * Kernel swaps vectors a and b.
 * **********************************************************************/

__global__ void swapfunc(float *a, float *b);

/*************************************************************************
 * Kernel returns 
 * mm0 = mm0 * cosr - sinr * mm1;
 * mm1 = mm1 * cosr + sinr * mm0;
 *
 * for the input vectors mm0 and mm1 and scalars cosr and sinr. 
 * **********************************************************************/

__global__ void vttwobytwo(float *mm1, float *mm0, float cosr, float sinr, int n);

#endif
