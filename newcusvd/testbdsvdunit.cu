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
* **********************************************************************************************/

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

#ifndef _TESTBDSVDUNIT_CU_
#define _TESTBDSVDUNIT_CU_

#include <stdio.h>
#include <stdlib.h>
#include "cublas.h"
#include <cutil.h>
#include <sys/timeb.h>
#include <time.h>

#include "testbdsvdunit.h"
#include "bdsvd.h"

struct timeb tt1;
void fillidentity(ap::real_2d_array& a, int n);
void fillsparsede(ap::real_1d_array& d,
     ap::real_1d_array& e,
     int n,
     double sparcity);
void getbdsvderror(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     bool isupper,
     const ap::real_2d_array& u,
     const ap::real_2d_array& c,
     const ap::real_1d_array& w,
     const ap::real_2d_array& vt,
     double& materr,
     double& orterr,
     bool& wsorted);
void testbdsvdproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     double& materr,
     double& orterr,
     bool& wsorted,
     bool& wfailed);

/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/

bool testbdsvd(int M, int N, bool silent, float *dU, float *dV, float *dC, float *dd, double *diagonal, double *superdiag)
{
    ap::real_1d_array d;
    ap::real_1d_array e;
    ap::real_2d_array mempty;
    
    int n;
    int maxn;
    int i;
    bool failcase;
    
    maxn = N;
    int maxm = M;

    d.setbounds(0, maxn-1);
    e.setbounds(0, maxn-2);

    n = N;
    int m = maxm;
    
    for(i=0; i < n; i++)
	    d(i) = diagonal[i];

    for(i=0; i < n-1; i++)
	    e(i) = superdiag[i];	

    int n1 = 0, n2 = 0, n3 = 0, n4 = 0;

    printf("Computing Diagonal matrix Sigma\n");
    ftime(&tt1);
    n1 = tt1.time;
    n2 = tt1.millitm;

    failcase = rmatrixbdsvd(d, e, n, true, true, m, n, dU, dV, dC, dd);
    cudaThreadSynchronize();

    ftime(&tt1);
    n3 = tt1.time;
    n4 = tt1.millitm;

    printf("Time required for diagonalization %d %d \n", (n3-n1), (n4-n2));
    for(i=0; i < n; i++)
	    diagonal[i] = d(i);

    return failcase;

    //
    // special case: zero divide matrix
    // unfixed LAPACK routine should fail on this problem
    //

    /*    
	  n = 7;
	  d(0) = -6.96462904751731892700e-01;
	  d(1) = 0.00000000000000000000e+00;
	  d(2) = -5.73827770385971991400e-01;
	  d(3) = -6.62562624399371191700e-01;
	  d(4) = 5.82737148001782223600e-01;
	  d(5) = 3.84825263580925003300e-01;
	  d(6) = 9.84087420830525472200e-01;
	  e(0) = -7.30307931760612871800e-02;
	  e(1) = -2.30079042939542843800e-01;
	  e(2) = -6.87824621739351216300e-01;
	  e(3) = -1.77306437707837570600e-02;
	  e(4) = 1.78285126526551632000e-15;
	  e(5) = -4.89434737751289969400e-02;
	  rmatrixbdsvd(d, e, n, true, false, mempty, 0, mempty, 0, mempty, 0);
    */

    //
    // zero matrix, several cases
    //

   /*
    for(i = 0; i <= maxn-1; i++)
    {
	    d(i) = 0;
    }
    for(i = 0; i <= maxn-2; i++)
    {
	    e(i) = 0;
    }
    for(n = 1; n <= maxn; n++)
    {
	    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
    }

    //
    // Dense matrix
    //
    for(n = 1; n <= maxn; n++)
    {
	    for(pass = 1; pass <= 10; pass++)
	    {
		    for(i = 0; i <= maxn-1; i++)
		    {
			    d(i) = 2*ap::randomreal()-1;
		    }
		    for(i = 0; i <= maxn-2; i++)
		    {
			    e(i) = 2*ap::randomreal()-1;
		    }
		    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
	    }
    }

    //
    // Sparse matrices, very sparse matrices, incredible sparse matrices
    //
    for(n = 1; n <= maxn; n++)
    {
	    for(pass = 1; pass <= 10; pass++)
	    {
		    fillsparsede(d, e, n, 0.5);
		    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
		    fillsparsede(d, e, n, 0.8);
		    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
		    fillsparsede(d, e, n, 0.9);
		    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
		    fillsparsede(d, e, n, 0.95);
		    testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
	    }
    }

    //
    // report
    //

    failr = double(failcount)/double(succcount+failcount);
    waserrors = materr>threshold||orterr>threshold||!wsorted||failr>failthreshold;

    if( !silent )
    {
	    printf("TESTING BIDIAGONAL SVD DECOMPOSITION\n");
	    printf("SVD decomposition error:                 %5.3le\n",
			    double(materr));
	    printf("SVD orthogonality error:                 %5.3le\n",
			    double(orterr));
	    printf("Singular values order:                   ");
	    if( wsorted )
	    {
		    printf("OK\n");
	    }
	    else
	    {
		    printf("FAILED\n");
	    }
	    printf("Always converged:                        ");
	    if( !wfailed )
	    {
		    printf("YES\n");
	    }
	    else
	    {
		    printf("NO\n");
		    printf("Fail ratio:                              %5.3lf\n",
				    double(failr));
	    }
	    printf("Fail matrix test:                        ");
	    if( !failcase )
	    {
		    printf("AS EXPECTED\n");
	    }
	    else
	    {
		    printf("CONVERGED (UNEXPECTED)\n");
	    }
	    printf("Threshold:                               %5.3le\n",
			    double(threshold));
	    if( waserrors )
	    {
		    printf("TEST FAILED\n");
	    }
	    else
	    {
		    printf("TEST PASSED\n");
	    }
	    printf("\n\n");
    }
    result = !waserrors;
    return result;
    */
}

void fillidentity(ap::real_2d_array& a, int n)
{
/*	int i;
	int j;
	a.setbounds(0, n-1, 0, n-1);
	for(i = 0; i <= n-1; i++)
	{
		for(j = 0; j <= n-1; j++)
		{
			if( i==j )
			{
				a(i,j) = 1;
			}
			else
			{
				a(i,j) = 0;
			}
		}
	}
*/
}

void fillsparsede(ap::real_1d_array& d,
		ap::real_1d_array& e,
		int n,
		double sparcity)
{
/*	int i;
	int j;

	d.setbounds(0, n-1);
	e.setbounds(0, ap::maxint(0, n-2));
	for(i = 0; i <= n-1; i++)
	{
		if( ap::randomreal()>=sparcity )
		{
			d(i) = 2*ap::randomreal()-1;
		}
		else
		{
			d(i) = 0;
		}
	}
	for(i = 0; i <= n-2; i++)
	{
		if( ap::randomreal()>=sparcity )
		{
			e(i) = 2*ap::randomreal()-1;
		}
		else
		{
			e(i) = 0;
		}
	}
*/
}

void getbdsvderror(const ap::real_1d_array& d,
		const ap::real_1d_array& e,
		int n,
		bool isupper,
		const ap::real_2d_array& u,
		const ap::real_2d_array& c,
		const ap::real_1d_array& w,
		const ap::real_2d_array& vt,
		double& materr,
		double& orterr,
		bool& wsorted)
{
/*      int i;
        int j;
        int k;
        double locerr;
        double sm;

	//
	// decomposition error
	//

	locerr = 0;
	for(i = 0; i <= n-1; i++)
	{
	for(j = 0; j <= n-1; j++)
	{
	sm = 0;
	for(k = 0; k <= n-1; k++)
	{
	sm = sm+w(k)*u(i,k)*vt(k,j);
	}
	if( isupper )
	{
	if( i==j )
	{
	locerr = ap::maxreal(locerr, fabs(d(i)-sm));
	}
	else
	{
	if( i==j-1 )
	{
	locerr = ap::maxreal(locerr, fabs(e(i)-sm));
	}
	else
	{
	locerr = ap::maxreal(locerr, fabs(sm));
	}
	}
	}
	else
	{
	if( i==j )
	{
	locerr = ap::maxreal(locerr, fabs(d(i)-sm));
	}
	else
	{
	if( i-1==j )
	{
	locerr = ap::maxreal(locerr, fabs(e(j)-sm));
	}
	else
	{
	locerr = ap::maxreal(locerr, fabs(sm));
	}
	}
	}
	}
	}
	materr = ap::maxreal(materr, locerr);

	//
	// check for C = U'
	// we consider it as decomposition error
	//

	locerr = 0;
	for(i = 0; i <= n-1; i++)
	{
		for(j = 0; j <= n-1; j++)
		{
		locerr = ap::maxreal(locerr, fabs(u(i,j)-c(j,i)));
		}
	}	
	materr = ap::maxreal(materr, locerr);

	//
	// orthogonality error
	//

	locerr = 0;
	for(i = 0; i <= n-1; i++)
	{
		for(j = i; j <= n-1; j++)
		{
			sm = ap::vdotproduct(u.getcolumn(i, 0, n-1), u.getcolumn(j, 0, n-1));
			if( i!=j )
			{
				locerr = ap::maxreal(locerr, fabs(sm));
			}
			else
			{
				locerr = ap::maxreal(locerr, fabs(sm-1));
			}
			sm = ap::vdotproduct(&vt(i, 0), &vt(j, 0), ap::vlen(0,n-1));
			if( i!=j )
			{
				locerr = ap::maxreal(locerr, fabs(sm));
			}
			else
			{
				locerr = ap::maxreal(locerr, fabs(sm-1));
			}
		}
	}
	orterr = ap::maxreal(orterr, locerr);

	//
	// values order error
	//

	for(i = 1; i <= n-1; i++)
	{
		if( w(i)>w(i-1) )
		{
			wsorted = false;
		}
	}
*/
}


void testbdsvdproblem(const ap::real_1d_array& d,
		const ap::real_1d_array& e,
		int n,
		double& materr,
		double& orterr,
		bool& wsorted,
		bool& wfailed)
{
 /*   
		   ap::real_2d_array u;
		   ap::real_2d_array vt;
		   ap::real_2d_array c;
		   ap::real_1d_array w;
		   int i;
		   int j;
		   int k;
		   double v;
		   double mx;

		   mx = 0;
		   for(i = 0; i <= n-1; i++)
		   {
			   if( fabs(d(i))>mx )
		       	   {
 				   mx = fabs(d(i));
			   }
		   }
		   
		   for(i = 0; i <= n-2; i++)
		   {
		   	if( fabs(e(i))>mx )
		   	{
		   		mx = fabs(e(i));
		   	}
		   }
		   if( mx==0 )
		   {
		  	 mx = 1;
		   }

	//
	// Upper BDSVD tests
	//

	w.setbounds(0, n-1);
	fillidentity(u, n);
	fillidentity(vt, n);
	fillidentity(c, n);
	for(i = 0; i <= n-1; i++)
	{
		w(i) = d(i);
	}
	if( !rmatrixbdsvd(w, e, n, true, false, u, n, c, n, vt, n) )
	{
		failcount = failcount+1;
		wfailed = true;
		return;
	}
	getbdsvderror(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);
	fillidentity(u, n);
	fillidentity(vt, n);
	fillidentity(c, n);
	for(i = 0; i <= n-1; i++)
	{
		w(i) = d(i);
	}
	if( !rmatrixbdsvd(w, e, n, true, true, u, n, c, n, vt, n) )
	{
		failcount = failcount+1;
		wfailed = true;
		return;
	}
	getbdsvderror(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);

	//
	// Lower BDSVD tests
	//
	w.setbounds(0, n-1);
	fillidentity(u, n);
	fillidentity(vt, n);
	fillidentity(c, n);
	for(i = 0; i <= n-1; i++)
	{
		w(i) = d(i);
	}
	if( !rmatrixbdsvd(w, e, n, false, false, u, n, c, n, vt, n) )
	{
		failcount = failcount+1;
		wfailed = true;
		return;
	}
	getbdsvderror(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);

	fillidentity(u, n);
	fillidentity(vt, n);
	fillidentity(c, n);

	for(i = 0; i <= n-1; i++)
	{
		w(i) = d(i);
	}
	if( !rmatrixbdsvd(w, e, n, false, true, u, n, c, n, vt, n) )
	{
		failcount = failcount+1;
		wfailed = true;
		return;
	}
	getbdsvderror(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);

	//
	// update counter
	//
	succcount = succcount+1;
*/
}

/*************************************************************************
  Silent unit test
 *************************************************************************/
bool testbdsvdunit_test_silent()
{
	bool result = 1;
	//result = testbdsvd(true);
	return result;
}


/*************************************************************************
  Unit test
 *************************************************************************/
bool testbdsvdunit_test()
{
	bool result = 1;
	//result = testbdsvd(false);
	return result;
}
#endif

