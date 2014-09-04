        !-------------------------------------------------------------------------------------------
        !-------------------------------------------------------------------------------------------
        module data_global
         implicit none 

         include 'mpif.h'

         integer,parameter :: M=1001 , N=401 , Nf=1984

         double precision,dimension(1:M,1:N) :: psi1,psi2,vort1,vort2
         double precision,dimension(1:M,1:N) :: mpsi,mvort
         double precision,dimension(1:Nf,1:Nf) :: AAt,BBt
         double precision,dimension(1:M,1:N) :: x,h1
         double precision,dimension(1:M,1:N) :: y,h2

         double precision :: dzi,deta,idzi,ideta

         character*20 :: fn(1:Nf),fnn(1:Nf),nn1
         character*20 :: st(1:Nf),stt(1:Nf),sr
         
         integer      :: rm(1:100),cm(1:100)          !No of elements must be greater than or equal to (nproc-1)!
         integer      :: istart=0,n_runs
        end module

        !-------------------------------------------------------------------------------------------
        !-------------------------------------------------------------------------------------------

 
        program Matrix_finding
       
         use data_global
         implicit none
  
         double precision,dimension(1:M,1:N) :: xx,yy
         double precision,dimension(0:N)     :: sum_psi,sum_vor

         double precision                    :: ax(1:4),bx(1:4),Intg_psi,Intg_vor 
         
         integer                             :: rank, nproc, ierr, dest, source, status(MPI_STATUS_SIZE)
         integer                             :: i,j,k,r,c,nn(Nf),p
         integer                             :: local_strt_t, local_end_t, local_strt_b, local_end_b
         integer                             :: local_strt, local_end, s, e,rem,local_n,rrem,frem
         integer                             :: k_local,r_local,c_local,c_strt

         integer, dimension(1:32)            :: local_strt0,local_end0,local_strt_t0,local_end_t0
         integer, dimension(1:32)            :: local_strt_b0,local_end_b0

         
         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

         !----------------------Defining Local Start and Local End---------------------! 
         
         rem=mod(Nf,2*nproc)
         local_n=Nf/(2*nproc)    
         rrem=mod(rem,nproc)
         frem=mod(Nf,nproc)

         !Done so as to evenly distribute the work among the processors; handles remainder (Nf%nproc.ne.0)
         
         if(rank .lt. frem)then

           s=rank*(Nf/nproc)+1+rank
           e=s+(Nf/nproc)-1+1

         else

           s=rank*(Nf/nproc)+1+frem
           e=s+(Nf/nproc)-1

         endif


         if(rem.lt.nproc)then

         if(rank .lt. rem)then

           local_strt_t= (rank*local_n + 1) + rank               ! in the top half 
           local_end_t = (local_strt_t + local_n -1) +1

           local_end_b= Nf-(rank*local_n)                        ! in the bottom half    
           local_strt_b = local_end_b-local_n+1

         else

           local_strt_t= (rank*local_n + 1) + rem
           local_end_t = (local_strt_t + local_n -1) 

           local_end_b= Nf-(rank*local_n)         
           local_strt_b = local_end_b-local_n+1
   
         endif
         
         else

           if(rank .lt.rrem)then

           local_strt_t= (rank*local_n + 1) + 2*rank               ! in the top half 
           local_end_t = (local_strt_t + local_n -1) +2

           local_end_b= Nf-(rank*local_n)                        ! in the bottom half    
           local_strt_b = local_end_b-local_n+1

           else

           local_strt_t= (rank*local_n + 1) + rrem+rank
           local_end_t = (local_strt_t + local_n -1)+1

           local_end_b= Nf-(rank*local_n)
           local_strt_b = local_end_b-local_n+1
          
           endif
           
         endif

         print*,rank,local_strt_t, local_end_t, 'from the top' ,s,e

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
         print*,rank,local_strt_b, local_end_b, 'from the bottom'

         !-------------------------------------------------------------------------------!
         
         !----------Getting and broadcasting grid data, file details and mean------------!

         if(rank.eq.0)then
         
           call grid_properties
          !call grid_generation_cylinder
          
           open(11,file="fileread")       
        
           DO i = 1, Nf 
         !     read(11,*)                  !To skip a line!
              read(11,*) fn(i)
              nn(i)=len_trim(fn(i))
              nn1=fn(i)
              st(i) = nn1(1:nn(i)) 
              stt(i)= st(i)
              write(*,*)st(i)

           END DO
  
         close(11)
        
         call read_mean_exist_file
        
        endif
       
        call MPI_BCAST(dzi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(deta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        DO i=1,Nf

           call MPI_BCAST(nn(i), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(fn(i), nn(i), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(st(i), nn(i), MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)
           call MPI_BCAST(stt(i),nn(i), MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)

        END DO
        
        DO j=1,N

           call MPI_BCAST(h1(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(h2(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

           call MPI_BCAST(mpsi(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(mvort(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        END DO
       
        dest=0

        if (rank.ne.0)then

          call MPI_SEND(local_strt_t, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)
          call MPI_SEND(local_end_t, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)

          call MPI_SEND(local_strt_b, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)
          call MPI_SEND(local_end_b, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)

        endif

        if (rank.eq.0)then

           DO i=1,nproc-1
              
              source =i
              
              call MPI_RECV(local_strt_t0(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)
              call MPI_RECV(local_end_t0(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)

              call MPI_RECV(local_strt_b0(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)
              call MPI_RECV(local_end_b0(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)

           END DO
      
        endif

       !----------------------------------------------------------------------------------!

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        AAt = 0.0d0
        BBt = 0.0d0
        
        DO r=s,e

          ! call system ("gunzip "//fn(r))

        END DO

        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        n_runs=1
        write(sr(3:3),'(I1)')mod(rank,10)
        write(sr(2:2),'(I1)')(mod(rank,100))/10
        write(sr(1:1),'(I1)')rank/100
        nn=len_trim(sr)
        print*,sr
        
        if (istart.eq.1) then

           open(rank+10,file='running_status_'//sr(1:3))
           read(rank+10,*)k_local,r_local,c_local
           close(rank+10)

       else

           k_local=1

       endif

       open(rank+200,file='data'//sr(1:3)//'.bin',position='append')


        DO k=k_local,2,1

           if(k.eq.1)then                                   !top half!

             local_strt= local_strt_t
             local_end = local_end_t                    

           else                                             !bottom half!

             local_strt= local_strt_b                  
             local_end = local_end_b
             
           endif
           
           if (istart.eq.1.and.n_runs.eq.1) local_strt=r_local

           DO r = local_strt, local_end, 1

                open(rank,file=st(r),form="UNFORMATTED")
   
                !read(rank,*)
                !read(rank,*)
                !read(rank,*)
                !read(rank,*)
                !read(rank,*)
                !read(rank,*)
      

                read(rank)vort1
              !  print*,vort1(M,N)
                close(rank)
                
                if(istart.eq.1.and.n_runs.eq.1)then
                  c_strt=c_local
                else
                  c_strt=1
                endif
         
                DO c = c_strt, r, 1

                   open(rank,file=stt(c),form="UNFORMATTED")
                   
                  ! read(rank,*)
                  ! read(rank,*)
                  ! read(rank,*)
                  ! read(rank,*)
                  ! read(rank,*)
                  ! read(rank,*)

                   read(rank)vort2
                   close(rank)
 
                   !-----------------------Finding the elements of the correlation matrix-----------------------!

         !               Intg_psi = 0.0d0
                        Intg_vor = 0.0d0

          
                     DO j = 2,N	
           
                        DO i = 1, M-1
 
         !                  ax(1) = (psi(i  ,j  ,1)-mpsi(i  ,j  ))*(psi(i  ,j  ,2)- mpsi(i  ,j  ))*h1(i  ,j  )*h2(i  ,j  )
         !                  ax(2) = (psi(i  ,j-1,1)-mpsi(i  ,j-1))*(psi(i  ,j-1,2)- mpsi(i  ,j-1))*h1(i  ,j-1)*h2(i  ,j-1)
         !                  ax(3) = (psi(i+1,j  ,1)-mpsi(i+1,j  ))*(psi(i+1,j  ,2)- mpsi(i+1,j  ))*h1(i+1,j  )*h2(i+1,j  )
         !                  ax(4) = (psi(i+1,j-1,1)-mpsi(i+1,j-1))*(psi(i+1,j-1,2)- mpsi(i+1,j-1))*h1(i+1,j-1)*h2(i+1,j-1)

                           bx(1) = (vort1(i  ,j  )-mvort(i  ,j  ))*(vort2(i  ,j  )-mvort(i  ,j  ))*h1(i  ,j  )*h2(i  ,j  )
                           bx(2) = (vort1(i  ,j-1)-mvort(i  ,j-1))*(vort2(i  ,j-1)-mvort(i  ,j-1))*h1(i  ,j-1)*h2(i  ,j-1)
                           bx(3) = (vort1(i+1,j  )-mvort(i+1,j  ))*(vort2(i+1,j  )-mvort(i+1,j  ))*h1(i+1,j  )*h2(i+1,j  )
                           bx(4) = (vort1(i+1,j-1)-mvort(i+1,j-1))*(vort2(i+1,j-1)-mvort(i+1,j-1))*h1(i+1,j-1)*h2(i+1,j-1)

        !                   Intg_psi = Intg_psi + 0.25d0 * (ax(1)+ax(2)+ax(3)+ax(4)) * deta * dzi
                           Intg_vor = Intg_vor + 0.25d0 * (bx(1)+bx(2)+bx(3)+bx(4)) * deta * dzi

                        END DO

                    END DO   

        !            Intg_psi = Intg_psi / Nf
                    Intg_vor = Intg_vor / Nf                  
                  
          
         !         AAt(r,c) = Intg_psi 
                  BBt(r,c) = Intg_vor 
                  
                                    
        !          open(rank+10,file='running_status3_'//sr(1:3))
        !          write(rank+10,*)k,r,c+1
        !          close(rank+10)

                  open(rank+100,file='running_status_'//sr(1:3))
                  write(rank+100,*)k,r,c+1
                  close(rank+100)

          !        write(rank+200,*)AAt(r,c),BBt(r,c),r,c
                   write(rank+200,*)BBt(r,c),k,r,c

                  !----------------------------------------------------------------------------------------------------------!

                  write(*,*)r,c,"IS DONE", rank

              END DO
              
              n_runs=n_runs+1

            END DO

        END DO
        
        close(rank+200)
       
       ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

        DO r=s,e

         !  call system ("gzip "//st(r))

        END DO

        !--------------------------------------------Message Passing-----------------------------------------------! 
       ! open(rank+20,file='data'//sr(1:3))
        
       ! DO k=1,2

       !     if(k.eq.1)then                                   !top half!

       !      local_strt= local_strt_t
       !      local_end = local_end_t

       !      if (rank.eq.0)then

       !         DO i=1,nproc-1

       !            local_strt0(i)= local_strt_t0(i)
       !            local_end0(i) = local_end_t0(i)

       !         END DO

       !      endif

       !    else                                             !bottom half!

       !      local_strt= local_strt_b
       !      local_end = local_end_b

       !       if (rank.eq.0)then

       !         DO i=1,nproc-1

       !            local_strt0(i)= local_strt_b0(i)
       !            local_end0(i) = local_end_b0(i)

       !         END DO

      !       endif

      !     endif

      !        DO r = local_strt, local_end

      !           DO c=1,r 
                 
      !              dest   =0
                    
        !            read(20+rank,*)AAt(r,c),BBt(r,c)
      !               read(20+rank,*)BBt(r,c)

      !              if (rank.ne.0) then

      !               call MPI_SEND(r, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)               !Position!
      !               call MPI_SEND(c, 1, MPI_INTEGER, dest, rank, MPI_COMM_WORLD, ierr)
                   !  call MPI_SEND(AAt(r,c), 1, MPI_DOUBLE_PRECISION, dest, rank, MPI_COMM_WORLD, ierr)   !Data!
      !               call MPI_SEND(BBt(r,c), 1, MPI_DOUBLE_PRECISION, dest, rank, MPI_COMM_WORLD, ierr)

      !              endif

      !         END DO

      !       END DO
             
      !     if(rank.eq.0)then 
             
      !       DO i=1,nproc-1
    
      !          DO r = local_strt0(i), local_end0(i)

      !             DO c=1,r

      !                 source=i

      !                 call MPI_RECV(rm(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)
      !                 call MPI_RECV(cm(i), 1, MPI_INTEGER,source, i, MPI_COMM_WORLD, status, ierr)
                    !   call MPI_RECV(AAt(rm(i),cm(i)), 1, MPI_DOUBLE_PRECISION, source, i,  MPI_COMM_WORLD, status, ierr)
      !                 call MPI_RECV(BBt(rm(i),cm(i)), 1, MPI_DOUBLE_PRECISION, source, i,  MPI_COMM_WORLD, status, ierr)    

      !              END DO   
                  
      !            END DO
              
      !         END DO

      !      endif

      !   END DO

         !--------------------------------------------------------------------------------------------------------------!

      !  if(rank.eq.0) then
       
      !    DO r = 1, Nf

      !       DO c = 1, Nf

      !          if (c .gt. r) then

                   !AAt(r,c) = AAt(c,r)              
      !             BBt(r,c) = BBt(c,r)

      !          endif

      !       ENDDO

      !    ENDDO
        
         ! open(11,file="P_MATRIX_AAT1.dat")
      !    open(12,file="P_MATRIX_BBT1.dat")

      !    DO c = 1, Nf

      !       DO r = 1, Nf

                !write(11,*) AAt(r,c)
      !          write(12,*) BBt(r,c)

      !       ENDDO

      !    ENDDO

         ! CLOSE(11)
      !    CLOSE(12)
 
      !  end if
        
        call MPI_FINALIZE(ierr)

        end program Matrix_finding
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
		
      subroutine grid_generation_cylinder
       
       use data_global
       IMPLICIT NONE

       INTEGER :: i,j,ip1,im1,jp1,jm1
       DOUBLE PRECISION :: rj,chn,pi
       DOUBLE PRECISION :: x_eta(M,N),x_zi(M,N)
       DOUBLE PRECISION :: y_eta(M,N),y_zi(M,N)

      
       dzi=1.0d0/(M-1); deta=1.0d0/(N-1)
       pi=4*datan(1.0d0)

       do j=1,N
        chn = ((j-1)*1.0d0)/(N-1)
        rj = 0.50d0+20.d0*(1-(dtanh(3.20d0*(1-chn)))/(dtanh(3.20d0)))
        do i=1,M
          x(i,j) = rj*cos(2*pi*(i-1)/(M-1))
          y(i,j) = rj*sin(2*pi*(i-1)/(M-1))
        enddo
       enddo
       
       DO j=2,N-1
         DO i=2,M
           im1=i-1
           ip1=i+1
           IF(i.EQ.M) ip1=2
           jm1=j-1
           jp1=j+1
           x_zi(i,j)=(x(ip1,j)-x(im1,j))*0.50d0/dzi
           y_zi(i,j)=(y(ip1,j)-y(im1,j))*0.50d0/dzi
           x_eta(i,j)=(x(i,jp1)-x(i,jm1))*0.50d0/deta
           y_eta(i,j)=(y(i,jp1)-y(i,jm1))*0.50d0/deta
         END DO
       END DO

       DO i=2,M
            im1=i-1
            ip1=i+1
            IF(i.EQ.M) ip1=2
            x_zi(i,1)=(x(ip1,1)-x(im1,1))*0.50d0/dzi
            y_zi(i,1)=(y(ip1,1)-y(im1,1))*0.50d0/dzi

            x_eta(i,1)=(x(i,2)-x(i,1))/deta
            y_eta(i,1)=(y(i,2)-y(i,1))/deta

            x_zi(i,N)=(x(ip1,N)-x(im1,N))*0.50d0/dzi
            y_zi(i,N)=(y(ip1,N)-y(im1,N))*0.50d0/dzi
            x_eta(i,N)=(x(i,N)-x(i,N-1))/deta
            y_eta(i,N)=(y(i,N)-y(i,N-1))/deta
       END DO

       DO j=1,N
         x_zi(1,j)=x_zi(M,j)
         y_zi(1,j)=y_zi(M,j)
         x_eta(1,j)=x_eta(M,j)
         y_eta(1,j)=y_eta(M,j)
       END DO

       DO j=1,N
         DO i=1,M
           h1(i,j)=dsqrt(x_zi(i,j)*x_zi(i,j)+y_zi(i,j)*y_zi(i,j))
           h2(i,j)=dsqrt(x_eta(i,j)*x_eta(i,j)+y_eta(i,j)*y_eta(i,j))
         END DO
       END DO 
       write(*,*) "end of grid generation"

      end subroutine grid_generation_cylinder


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------  

        subroutine read_mean_exist_file
        
         use data_global
         implicit none

         integer :: i,j
        
         !call system ("gunzip MEAN_psi_vor.gz") 

         open(1,file="P_MEAN_psi_vor.bin",form="UNFORMATTED")
       !  open(1,file="P_MEAN_psi_vor.test")

         !   read(1,*)
         !   read(1,*) 
        ! do j=1,N
	! do i=1,M
        !    read(1,*) mvort(i,j)
	     read(1)mvort
!	     print*, mvort(1,1)
        ! enddo
	! enddo

         close(1)

         !call system ("gzip MEAN_psi_vor")

         write(*,*) "end of reading mean existing file"

        end subroutine read_mean_exist_file
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

       subroutine grid_properties

         use data_global
         IMPLICIT NONE

         INTEGER :: i,j,ip1,im1,jp1,jm1
         DOUBLE PRECISION :: rj,chn,pi
         DOUBLE PRECISION :: x_eta(M,N),x_zi(M,N)
         DOUBLE PRECISION :: y_eta(M,N),y_zi(M,N)

         double precision :: d
     
         open(1,file='GRID_POD.dat',status='old')

         read(1,*)
         read(1,*)
       
         DO j=1,N
            DO i=1,M
               read(1,*)x(i,j),y(i,j),d,d,d
            END DO
         END DO
        
        dzi=1.0d0/(M-1)
        deta=1.0d0/(N-1)      

       DO j=2,N-1
         DO i=2,M
           im1=i-1
           ip1=i+1
           IF(i.EQ.M) ip1=2
           jm1=j-1
           jp1=j+1
           x_zi(i,j)=(x(ip1,j)-x(im1,j))*0.50d0/dzi
           y_zi(i,j)=(y(ip1,j)-y(im1,j))*0.50d0/dzi
           x_eta(i,j)=(x(i,jp1)-x(i,jm1))*0.50d0/deta
           y_eta(i,j)=(y(i,jp1)-y(i,jm1))*0.50d0/deta
         END DO
       END DO

       DO i=2,M
            im1=i-1
            ip1=i+1
            IF(i.EQ.M) ip1=2
            x_zi(i,1)=(x(ip1,1)-x(im1,1))*0.50d0/dzi
            y_zi(i,1)=(y(ip1,1)-y(im1,1))*0.50d0/dzi

            x_eta(i,1)=(x(i,2)-x(i,1))/deta
            y_eta(i,1)=(y(i,2)-y(i,1))/deta

            x_zi(i,N)=(x(ip1,N)-x(im1,N))*0.50d0/dzi
            y_zi(i,N)=(y(ip1,N)-y(im1,N))*0.50d0/dzi
            x_eta(i,N)=(x(i,N)-x(i,N-1))/deta
            y_eta(i,N)=(y(i,N)-y(i,N-1))/deta
       END DO

       DO j=1,N
         x_zi(1,j)=x_zi(M,j)
         y_zi(1,j)=y_zi(M,j)
         x_eta(1,j)=x_eta(M,j)
         y_eta(1,j)=y_eta(M,j)
       END DO

       DO j=1,N
         DO i=1,M
           h1(i,j)=dsqrt(x_zi(i,j)*x_zi(i,j)+y_zi(i,j)*y_zi(i,j))
           h2(i,j)=dsqrt(x_eta(i,j)*x_eta(i,j)+y_eta(i,j)*y_eta(i,j))
         END DO
       END DO
 

       end subroutine grid_properties
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
