      Program convert_to_binary !run with mpirun -np 1 ./a.out
 
      implicit none
      include 'mpif.h'

      integer,parameter :: M=1001, N=401, Nf=1984

      integer          :: i,j,k,nn
      double precision :: vort(M,N),d
      character*100    :: filename

      integer :: ierr, rank, nproc

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)


      open(1,file='fileread')
     
      do k=1,Nf
         read(1,*)filename
	 nn=len_trim(filename)

         write(*,*)filename(1:nn)

	 open(10,file=filename(1:nn))
	 do j=1,N
	 do i=1,M
	    read(10,*)vort(i,j)
	 enddo
	 enddo

	 open(20,file=filename(1:5)//'snapshot.bin',form='UNFORMATTED')
         write(20)vort 
	 close(20)
      enddo

      call MPI_FINALIZE(ierr)

      End program
