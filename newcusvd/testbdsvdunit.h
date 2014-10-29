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

/************************************************************************************************
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

************************************************************************************************/

#ifndef _TESTBDSVDUNIT_H_
#define _TESTBDSVDUNIT_H_

#include "ap.h"

/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
**********************************************************************/

bool testbdsvd(int M, int N, bool silent, float *dU, float *dV, float *dC, float *dd, double *diagonal, double *superdiag);

/*************************************************************************
Silent unit test
*************************************************************************/

bool testbdsvdunit_test_silent();

/*************************************************************************
Unit test
*************************************************************************/

bool testbdsvdunit_test();

#endif
