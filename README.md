##PODCUDA

This library aims to be a CUDA based implementation of the POD.

###POD (Proper Orthogonal Decomposition)

[POD](http://en.wikipedia.org/wiki/Principal_components_analysis) or Proper Orthogonal Decomposition is a tool used in data analysis.

 

In simple terms, it works by generating a lower dimensional approximate of
noisy high dimensional data. It is closely related to PCA, SVD and a lot of
other such methods

###CUDA

[CUDA]( http://www.nvidia.com/object/cuda_home_new.html) or Compute Unified Device Architecture, is a parallel computing platform
created by NVIDIA corp. It allows using general GPUs present in mid-high end
computers for computation. The advantage of using this over traditional CPU
based computation is that there are generally around hundreds of compute cores
on a normal GPU. While current generation of GPUs do not allow better than
single (double in case of high end) precision computation, it still is good
enough for a lot of trivially parallelizable computation.

<hr>
The primary focus of this library will be ease of use, and speed. CUDA does
have its limitations, and I will try to work around it as much as it allows.


The target is to be finished with single precision computation in the coming
month, so that work on double precision versions can be started by the start
of October.


<hr>

###Related work and references : 

Related decompositions : 

* [Eigen Decomposition](http://en.wikipedia.org/wiki/Matrix_decomposition#Eigendecomposition)
* [K-L Theorem](http://en.wikipedia.org/wiki/Karhunen%E2%80%93Lo%C3%A8ve_theorem)

Helpful tutorials and introductions : 

* http://web.stanford.edu/~amsallem/CA-CME345-Ch4.pdf
* http://www.iisc.ernet.in/currsci/apr102000/tutorial2.pdf


For queries, contact me at govind.93@gmail.com.

I also lurk on Quakenet and Freenode, govg in both.
