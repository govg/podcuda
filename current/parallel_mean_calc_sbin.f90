
        module data_all

        implicit none
        include 'mpif.h'

        integer,parameter :: M=1001 , N=401  , Nf=1984

        double precision,dimension(1:M,1:N) :: psi,vort
        double precision,dimension(1:M,1:N) :: mpsi,mvort,mlvort,mlpsi,mlvort_temp, mlpsi_temp
        double precision,dimension(1:Nf,1:Nf) :: AAt,BBt
        double precision,dimension(1:M,1:N) :: x,h1
        double precision,dimension(1:M,1:N) :: y,h2

        double precision :: dzi,deta,idzi,ideta

        character*20 :: fn(1:Nf),fnn(1:Nf),nn1
        character*20 :: st(1:Nf),stt(1:Nf)

        end module
 !-------------------------------------------------------------------------------------------
 !------------------------------------------------------------------------------------------- 


        program read_mean_file

        use data_all
        implicit none

        integer :: i,j,k,nn(Nf),id,local_Nf,local_strt,local_end,rm,local_rm
        integer ::rank,nproc,ierr,status(MPI_STATUS_SiZE)
       
        call MPI_INIT(ierr)
        
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)        
        
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        
        local_Nf=Nf/nproc
        rm = mod(Nf,nproc)

        if (rank.lt.(nproc-rm)) then

           local_strt=1+(rank*local_Nf)
           local_end=local_strt+local_Nf-1
           
        else

           local_strt=1+(rank*local_Nf) + (rank-(nproc-rm))
           local_end =local_strt+local_Nf

        endif

        print*,local_strt,local_end, rank

        mvort = 0.0d0
!        mpsi  = 0.0d0
        
        if(rank.eq.0) then
          open(11,file="fileread")
          DO i = 1, Nf 
            ! read(11,*)
             read(11,*) fn(i)
             nn(i)=len_trim(fn(i))
             nn1=fn(i)
             st(i) = nn1(1:nn(i)) 
             stt(i)= st(i)
             write(*,*)stt(i)
          ENDDO 
          close(11)
       endif
          
         
       DO i=1,Nf

          call MPI_BCAST(nn(i), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(fn(i), nn(i), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(st(i), nn(i), MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)
          call MPI_BCAST(stt(i),nn(i), MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)

       ENDDO

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       mlvort=0.0d0
       mlpsi=0.0d0
 
        DO k = local_strt,local_end
                
                open(k,file=st(k),form='unformatted')
                
                read(k)vort

		print*,vort(1,2)

         	DO j = 1, N
          		DO i = 1, M
                                mlvort(i,j)=mlvort(i,j)+vort(i,j)
          		ENDDO
         	ENDDO
         	close(k)

        	print *,"MEAN ",stt(k)," IS OVER" 

       ENDDO
       
       DO i=1,M
 
          call MPI_REDUCE(mlvort(i,1:N), mvort(i,1:N), N, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(mlpsi(i,1:N), mpsi(i,1:N), N, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

       ENDDO  
        
       if(rank.eq.0)then

        do j=1,N
           do i=1,M

           mvort(i,j) = mvort(i,j)/ Nf
           print*,mvort(i,j)
    !       mpsi(i,j)  = mpsi(i,j)/ Nf

	   end do
	end do

        open(1,file="P_MEAN_psi_vor.test")
        write(1,*)"VARIABLES = x, Y, PSI, VORT" 
        WRITE(1,*)"ZONE I = ",M,"J = ",N,"F = POINT" 
        DO j = 1, N
         DO i = 1, M
          write(1,10)  mvort(i,j)
          10 format (4(E1  7.10,1X))
         ENDDO
        ENDDO
        close(1)

        !call system ("gzip MEAN_psi_vor")

        open(1,file="P_MEAN_psi_vor.bin",form="UNFORMATTED")
        write(1)  mvort
        close(1)
        
      endif

     call MPI_FINALIZE(ierr)

     end
