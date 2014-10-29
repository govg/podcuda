/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from magma_z_init.cpp normal z -> d, Tue Sep  2 12:38:35 2014
       @author Hartwig Anzt
*/

#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ostream>
#include <assert.h>
#include <stdio.h>
#include "magmasparse_d.h"
#include "magma.h"
#include "mmio.h"



using namespace std;








/**
    Purpose
    -------

    Initialize a magma_d_vector.


    Arguments
    ---------

    @param
    x           magma_d_vector
                vector to initialize   

    @param
    mem_loc     magma_location_t
                memory for vector 

    @param
    num_rows    magma_int_t
                desired length of vector      

    @param
    values      double
                entries in vector


    @ingroup magmasparse_daux
    ********************************************************************/

magma_int_t 
magma_d_vinit(    magma_d_vector *x, 
                  magma_location_t mem_loc,
                  magma_int_t num_rows, 
                  double values ){

    x->memory_location = Magma_CPU;
    x->num_rows = num_rows;
    x->nnz = num_rows;
    if( mem_loc == Magma_CPU ){
        x->memory_location = Magma_CPU;

        magma_dmalloc_cpu( &x->val, num_rows );
        if ( x->val == NULL )
            return MAGMA_ERR_HOST_ALLOC;
        for( magma_int_t i=0; i<num_rows; i++)
             x->val[i] = values; 
        return MAGMA_SUCCESS;  
    }
    else if( mem_loc == Magma_DEV ){
        x->memory_location = Magma_DEV;

        if (MAGMA_SUCCESS != magma_dmalloc( &x->val, x->num_rows)) 
            return MAGMA_ERR_DEVICE_ALLOC;

        magmablas_dlaset(MagmaFull, num_rows, 1, values, values, x->val, num_rows);

        return MAGMA_SUCCESS; 
    }
    return MAGMA_SUCCESS; 
}



   


