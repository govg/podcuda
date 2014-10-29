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

#ifndef _ROTATIONS_H_
#define _ROTATIONS_H_

#include "ap.h"
#include "cublas.h"
#include <cutil.h>

/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    WORK        -   working array whose index ranges within [N1..N2].
    mat         -   Matrix on the device.
    dC          -   Temporary vector on the device
    dd          -   Temporary vector on the device 

Output parameters:
    mat         -   Transformed matrix

Utility subroutine.
*************************************************************************/

float * applyrotationsfromtheleft(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s, ap::real_1d_array& work, float *mat, float *dC, float *dd);

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
    mat         -   Matrix on the device. 
    dC          -   Temporary vector on the device.
    dd          -   Temporary vector on the device.


Output parameters:
    mat         -   Transformed matrix.

Utility subroutine.
*************************************************************************/

float* applyrotationsfromtheright(bool isforward,
     int m1,
     int m2,
     int n1,
     int n2,
     const ap::real_1d_array& c,
     const ap::real_1d_array& s,
     ap::real_1d_array& work, float *mat, float *dC, float *dd);


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/

void generaterotation(double f, double g, double& cs, double& sn, double& r);

#endif
