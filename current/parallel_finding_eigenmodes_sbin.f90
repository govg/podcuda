!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
        module data_global

        implicit none
        include 'mpif.h'

        integer,parameter :: M=1001 , N=401 ,Nr=401, Nf=1984 , Nmode=32  
        integer,parameter :: NI=10 , Md=M/NI                            !not used!

        double precision,dimension(1:M,1:N,1:Nmode) :: spodmode,vpodmode,nspodmode,nvpodmode
        double precision,dimension(1:M,1:N)         :: psi,vor,psio 
        double precision,dimension(1:M,1:N)         :: mpsi,mvor
        double precision,dimension(1:Nf,1:Nf)       :: U1,U2
        double precision,dimension(1:Nf,1:Nmode)    :: sco,vco
        double precision,dimension(1:Nmode)         :: snorm,vnorm
        double precision,dimension(1:Nf)            :: eigv,eigs         !not used! 
        double precision,dimension(1:M,1:N)         :: x,h1,y,h2
      
        double precision :: dzi,deta,idzi,ideta,time

        double precision, parameter :: dtime=0.2d0,t0=0.2d0

        integer :: istart,iend,igstart,igend,nn
        integer :: ierr, rank, nproc
        integer :: local_strt, local_end, local_Nf
        integer :: strt_m,end_m,local_Nm,rm,rmm


        character*20 :: fn(1:Nf),fn1(1:Nf)
        character*20 :: st(1:Nf),nn1,stt(1:Nf),sttt(1:Nf)

        end module

!-------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------- 

        program POS_cal_mode

        use data_global
        implicit none

        double precision,dimension(1:N) :: sum_psi,sum_vor
        double precision :: ax(1:4),bx(1:4),Intg_psi,Intg_vor  

        integer :: i,j,r,c,k,ig

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

        local_Nf= Nf/nproc
        rm=mod(Nf,nproc)

        local_Nm= Nmode/nproc
        rmm=mod(Nmode,nproc)

        if (rank.lt.(nproc-rm)) then          !Modes divided among processors; handles remainder(Nf%nproc ne 0)!
           
           local_strt=1+ rank*local_Nf
           local_end =(rank+1)*local_Nf

        else
           
           local_strt= 1+rank*local_Nf+(rank-(nproc-rm))
           local_end =local_strt+local_Nf

        endif

        if(rank.lt.(nproc-rmm))then
          
          strt_m=1+ rank*local_Nm
          end_m =(rank+1)*local_Nm
          
        else

          strt_m=1+ rank*local_Nm+(rank-(nproc-rmm))
          end_m =strt_m+local_Nm

        endif

        print*,local_strt,local_end,rank
        print*,strt_m,end_m,rank


        call grid_properties           !Airfoil grid data file reqd!
        call read_mean_exist_file

        write(*,*) "mean read"

           if(rank.eq.0) then

           open(11,file="fileread")

           DO i = 1, Nf 

              !read(11,*)

              read(11,*) fn(i)
              nn=len_trim(fn(i))
              nn1=fn(i)
              fn1(i)=fn(i)
              st(i) = nn1(1:nn)
              stt(i)= st(i)  

           !   write(*,*) st(i)

           END DO 

           close(11)

           !open(3,file="eigenvalue_AAT.dat")
           !open(4,file="eigenvalue_BBT.dat")

           !DO r = 1, Nf

             ! read(3,*) eigs(r)
	     ! read(4,*) eigv(r)

          !ENDDO

          !close(3)
          !close(4)

          write(*,*) "eigen values read"
        
         ! open(1,file="eigenvector_AAT.dat")
          open(2,file="eigenvector_BBT.dat")

          DO c = 1, Nf
             DO r = 1, Nf

          !      read(1,*) U1(r,c) 
                read(2,*) U2(r,c)

             END DO
          END DO

          close(1)
          close(2)

          write(*,*) "end of reading vectors"
        
       end if

        call MPI_BCAST(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dzi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(deta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        DO i=1,Nf

           call MPI_BCAST(fn(i), nn, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(st(i), nn, MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)
           call MPI_BCAST(stt(i),nn, MPI_CHARACTER, 0,  MPI_COMM_WORLD, ierr)

           call MPI_BCAST(U1(1:Nf,i), Nf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(U2(1:Nf,i), Nf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        END DO

        DO j=1,N

           call MPI_BCAST(h1(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(h2(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

           call MPI_BCAST(mpsi(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
           call MPI_BCAST(mvor(1:M,j), M, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        END DO

        
!----------------------------------------------------------------------------------!
!     MATRIX MULTIPLICATION IN TERMS OF TENSORIAL OUTER PRODUCT AND CONTRACTION	   !
!            METHOD :: U' * A = Sigma * V' : V POD MODES MATRIX	                   !
!----------------------------------------------------------------------------------!

        spodmode = 0.0d0
        vpodmode = 0.0d0

        print*,'done'
        
        !DO k = local_strt, local_end

           !call system ("gunzip "//fn(k))

       ! END DO

        DO k = 1,Nf

           open(10,file=stt(k),form='unformatted')

           !write(*,*) stt(k)

!           read(10,*)
!           read(10,*)
!           read(10,*)
!           read(10,*)
!           read(10,*)
!           read(10,*)


            read(10) vor

            close(10)

           ! call system ("gzip "//stt(k)) 
         
            DO r = strt_m,end_m

               DO j = 1, Nr

                  DO i = 1, M

                   !  spodmode(i,j,r) = spodmode(i,j,r) + U1(k,r) * (psi(i,j) - mpsi(i,j))
                     vpodmode(i,j,r) = vpodmode(i,j,r) + U2(k,r) * (vor(i,j) - mvor(i,j))

                  END DO

                END DO
            END DO

            print*,k,'is done'

        END DO

        write(*,*) "end of pod mode calculations"
         
!-----------------------------------------------------------------------------------------------! 
!                                NORMALISING W.R.T ITS NORM					!
!-----------------------------------------------------------------------------------------------!


        snorm = 0.0d0
        vnorm = 0.0d0
          
        write(*,*) "entering normalising norm do loop", rank

        DO r = strt_m, end_m

            !Intg_psi = 0.0d0
            Intg_vor = 0.0d0
          
            DO j = 2, Nr
               
               DO i = 1, M-1

                !  ax(1) = spodmode(i  ,j  ,r) * spodmode(i  ,j  ,r) * h1(i  ,j  ) * h2(i  ,j  )
                !  ax(2) = spodmode(i+1,j  ,r) * spodmode(i+1,j  ,r) * h1(i+1,j  ) * h2(i+1,j  )
                !  ax(3) = spodmode(i  ,j-1,r) * spodmode(i  ,j-1,r) * h1(i  ,j-1) * h2(i  ,j-1)
                !  ax(4) = spodmode(i+1,j-1,r) * spodmode(i+1,j-1,r) * h1(i+1,j-1) * h2(i+1,j-1)

                  bx(1) = vpodmode(i  ,j  ,r) * vpodmode(i  ,j  ,r) * h1(i  ,j  ) * h2(i  ,j  )
                  bx(2) = vpodmode(i+1,j  ,r) * vpodmode(i+1,j  ,r) * h1(i+1,j  ) * h2(i+1,j  )
                  bx(3) = vpodmode(i  ,j-1,r) * vpodmode(i  ,j-1,r) * h1(i  ,j-1) * h2(i  ,j-1)
                  bx(4) = vpodmode(i+1,j-1,r) * vpodmode(i+1,j-1,r) * h1(i+1,j-1) * h2(i+1,j-1)

                  Intg_psi = Intg_psi + 0.25d0 * (ax(1)+ax(2)+ax(3)+ax(4)) * dzi * deta
                  Intg_vor = Intg_vor + 0.25d0 * (bx(1)+bx(2)+bx(3)+bx(4)) * dzi * deta
                 
                  !print*,Intg_vor

                END DO
         
            END DO
         
           ! snorm(r) = Intg_psi
            vnorm(r) = Intg_vor        !the number of norm is equal to no. of the modes choosen
         
        ENDDO

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        DO r = strt_m,end_m

            DO j = 1, Nr

               DO i = 1, M

             !     nspodmode(i,j,r) = spodmode(i,j,r) 
                  nvpodmode(i,j,r) = vpodmode(i,j,r)!/vnorm(r)

               END DO

            END DO

        END DO 

!--------------------------------------------------------------------------------------------!
!                         FINDING THE POD CONSTANTS a_k(t)                                   !
!                   METHOD 1 : a_k = (phi_k , V) ; (.) INNER PRODUCT    		     !
!--------------------------------------------------------------------------------------------!
       
        write(*,*) "entering finding ak loop" , rank
       
        DO k = 1,Nf

              !call system ("gunzip "//fn1(k))

               !write(*,*), stt(k)

               open(16,file=stt(k),form='unformatted')

         !      write(*,*)stt(k)

             !  read(16,*)
             !  read(16,*)
             !  read(16,*)
             !  read(16,*)
             !  read(16,*)
             !  read(16,*)

             read(16) vor

             close(16)

               !call system ("gzip "//sttt(k))
           
            DO r = strt_m, end_m

               ! Intg_psi = 0.0d0
                Intg_vor = 0.0d0
          
                DO j = 2, N
          
                   DO i = 1, M-1

                !      ax(1)= (psi(i  ,j  )-mpsi(i  ,j  ))*spodmode(i  ,j  ,r)*h1(i  ,j  )*h2(i  ,j  )
                !      ax(2)= (psi(i+1,j  )-mpsi(i+1,j  ))*spodmode(i+1,j  ,r)*h1(i+1,j  )*h2(i+1,j  )
                !      ax(3)= (psi(i  ,j-1)-mpsi(i  ,j-1))*spodmode(i  ,j-1,r)*h1(i  ,j-1)*h2(i  ,j-1)
                !      ax(4)= (psi(i+1,j-1)-mpsi(i+1,j-1))*spodmode(i+1,j-1,r)*h1(i+1,j-1)*h2(i+1,j-1)                              
          
                      bx(1)= (vor(i  ,j  )-mvor(i  ,j  ))*vpodmode(i  ,j  ,r)*h1(i  ,j  )*h2(i  ,j  ) 
                      bx(2)= (vor(i+1,j  )-mvor(i+1,j  ))*vpodmode(i+1,j  ,r)*h1(i+1,j  )*h2(i+1,j  ) 
                      bx(3)= (vor(i  ,j-1)-mvor(i  ,j-1))*vpodmode(i  ,j-1,r)*h1(i  ,j-1)*h2(i  ,j-1) 
                      bx(4)= (vor(i+1,j-1)-mvor(i+1,j-1))*vpodmode(i+1,j-1,r)*h1(i+1,j-1)*h2(i+1,j-1) 
          
                !      Intg_psi = Intg_psi + 0.25d0 * (ax(1)+ax(2)+ax(3)+ax(4)) * dzi * deta
                      Intg_vor = Intg_vor + 0.25d0 * (bx(1)+bx(2)+bx(3)+bx(4)) * dzi * deta

                   END DO
          
                END DO

              ! sco(k,r) = Intg_psi/snorm(r)
               vco(k,r) = Intg_vor/vnorm(r)

               write(*,*) k,r,"IS DONE" , rank

            END DO

        END DO

        DO k = local_strt, local_end

           !call system ("gzip "//st(k))

        END DO 

        
        call write_pod_data
        
        call MPI_FINALIZE(ierr)
        
        end program POS_cal_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

      subroutine grid_generation

       USE data_global
       IMPLICIT NONE

       INTEGER :: i,j,ip1,im1,jp1,jm1

       DOUBLE PRECISION :: rj,chn,pi
       DOUBLE PRECISION :: x_eta(M,N),x_zi(M,N)
       DOUBLE PRECISION :: y_eta(M,N),y_zi(M,N)

       dzi=1.0d0/(M-1); deta=1.0d0/(N-1)
       pi=4*datan(1.0d0)

       DO j=1,N

          chn = ((j-1)*1.0d0)/(N-1)
          rj = 0.50d0+20.d0*(1-(dtanh(3.20d0*(1-chn)))/(dtanh(3.20d0)))
          DO i=1,M

             x(i,j) = rj*cos(2*pi*(i-1)/(M-1))
             y(i,j) = rj*sin(2*pi*(i-1)/(M-1))

          END DO

       END DO
       
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

 end subroutine grid_generation

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------  

        subroutine read_mean_exist_file

          use data_global
          implicit none
          integer :: i,j
        
        
           open(15,file="P_MEAN_psi_vor.bin",form='unformatted')

           read(15) mvor

           close(15)
       
           write(*,*) "end of reading the mean value file"

        end subroutine read_mean_exist_file

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

 
        subroutine write_pod_data
        
          use data_global
          implicit none

          double precision :: dt
          integer :: i,j,k,r
          character*14 :: str
        
           dt = 0.2d0 
           write(*,*) "printing subroutine"
         
           DO r = strt_m, end_m
         
              str(01:01) = ACHAR(48 + int(r/100)-10*int(r/1000))
              str(02:02) = ACHAR(48 + int(r/10)-10*int(r/100))
              str(03:03) = ACHAR(48 + int(r/1)-10*int(r/10))

              open(1,file="POD_MODE_phi_"//str(1:3))

              write(1,*)"VARIABLES= X, Y, VMODE"
              write(1,*)"ZONE I = ", M, "J = ", N, "F = POINT"

              DO j = 1, N

                 DO i = 1, M

                     write(1,*)x(i,j), y(i,j), nvpodmode(i,j,r)

                 END DO

              END DO

            close(1) 
         
            open(2,file="P_POD_CONSTANTS_"//str(1:3))

            write(2,*)"VARIABLES= T, SVO"

                time=t0

                DO k = 1, Nf    
                   write(2,*) time , vco(k,r)     !change according to the time range selected!
                   time=time+dtime
                END DO

            close(2) 

        END DO 
        
        end subroutine write_pod_data

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
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


!-------------------------------------------------------------------------------------------
!                                     THE END                                              
!-------------------------------------------------------------------------------------------
