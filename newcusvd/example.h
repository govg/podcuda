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

#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

/* Header files included */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/timeb.h>

#include "cutil.h"
#include "device_functions.h"
#include "cublas.h"

unsigned int timer = 0;
using namespace std;

/********************************************
 * This function initializes the matrix with
 * identity.
********************************************/
float *initialize(int M);

#endif
