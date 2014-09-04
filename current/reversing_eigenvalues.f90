      program reversing_eigenvelue
      implicit none

       integer, parameter:: Nf=1984
       integer :: i,j
       double precision :: x1(Nf),x2(Nf),a(Nf),b(Nf)
       double precision :: x1v(Nf,Nf),x2v(Nf,Nf),av(Nf,Nf),bv(Nf,Nf)
	   
	
!---------------------------------------------------------------------------   
      ! open(1,file="eigenvalue1.dat")
       open(2,file="eigenvalue2.dat")
	   
       do i=1,Nf
       !read(1,*), x1(i)
       read(2,*), x2(i)
       end do
       !close(1)
       close(2)
	   
       do i=0,Nf-1
      ! a(i+1)= x1(Nf-i)
       b(i+1)= x2(Nf-i)
       end do
	   
       ! OPEN(3,FILE="eigenvalue_AAT.dat")
        OPEN(4,FILE="eigenvalue_BBT.dat")	

        do i=1,Nf
        !   write(3,*),a(i)
           write(4,*), b(i)
        end do	
       ! close(3)
        close(4)

        write(*,*) "end of eigenvalue reversing"

!------------------------------------------------------------------------------


      ! open(10,file="eigenvector1.dat")
       open(20,file="eigenvector2.dat")
	   
       do j=1,Nf
          do i=1,Nf
       !   read(10,*), x1v(i,j)
          read(20,*), x2v(i,j)
          end do
       end do
      ! close(10)
       close(20)

      do j=0,Nf-1
           do i=1,Nf
	 !  av(i,j+1)= x1v(i,Nf-j)
	   bv(i,j+1)= x2v(i,Nf-j)
           end do
       end do
	   
       ! OPEN(11,FILE="eigenvector_AAT.dat")
        OPEN(12,FILE="eigenvector_BBT.dat")	

       do j=1,Nf
        do i=1,Nf
       ! write(11,*), av(i,j)
        write(12,*), bv(i,j)
        end do	
      end do
       ! close(11)
        close(12)
		
		
      end 	   
