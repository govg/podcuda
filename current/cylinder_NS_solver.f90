!------------------------------------------------------------------------------!
!|                           *********************                            |!
!|                            !  PROGRAM MAGNUS !                             |!
!|                           *********************                            |!
!------------------------------------------------------------------------------!

       MODULE VARIABLES
       IMPLICIT NONE
       INCLUDE "mpif.h"   

       INTEGER :: rank,error,ierr
       INTEGER :: ivrt,pass,irs,istart
       INTEGER *8 :: NEQ,NEQP                                      

       INTEGER, PARAMETER :: MAX_FFT=512,M1=9,mlimit=10
       INTEGER :: COM_FFT_CR,FILT_FREQ,itrcount
       DOUBLE PRECISION, PARAMETER :: FFT_MIN_AMP=100000.0d0,TF=0.95d0
 
       INTEGER, PARAMETER :: imax=1001,jmax=401,nproc=16,plimit=150 
       INTEGER, PARAMETER :: n=imax,nm1=n-1,m=(jmax-1)/nproc+1
       INTEGER :: nin=392,nout=608,ibg,ibg1,ibg2
       INTEGER, PARAMETER :: arr_lim=(n+4)/2
       INTEGER, PARAMETER :: arr_lim1=arr_lim,arr_lim2=arr_lim
       REAL, PARAMETER :: alpha_fil=0.2,Fp=0.00,Afk=0.000d0
!       REAL, PARAMETER :: alpha_fil1d=0.49,expos=0.17

       REAL *8 :: victim1(-11:m+12),var(-11:m+12),fun(-11:m+12)        
       REAL *8 :: victim2(1:arr_lim),var_np(1:arr_lim),fun_np(1:arr_lim)
       REAL *8 :: fa(n,m),str,pi,deta,dzi,para,eps,res_sum
!       REAL *8 :: xx1,xx2,xst,s_zi(0:n),ds_zi(0:n),betaL
 
       DOUBLE PRECISION :: vrtn_rhs(0:n,0:m),vrtnf(0:n+1,0:m),resxF(0:n,0:mlimit) 
       DOUBLE PRECISION :: Xco(-2:n+2,-5:m+6),Yco(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: Xgo(-2:n+2,jmax),Ygo(-2:n+2,jmax)
       DOUBLE PRECISION :: h11(-2:n+2,-5:m+6),h22(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: h11g(-2:n+2,-5:jmax),h22g(-2:n+2,-5:jmax)

       DOUBLE PRECISION :: h22bh11_im1(-2:n+2,-5:m+6),h22bh11_ip1(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: h11bh22_jm1(-2:n+2,-5:m+6),h11bh22_jp1(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: h22bh11(-2:n+2,-5:m+6),h11bh22(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: dzi_sq_inv,deta_sq_inv,re_inv,h11h22_inv(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: two_dzi_inv,two_deta_inv,dh11bdn(-2:n+2,1),dh22bdn(-2:n+2,1)  
       DOUBLE PRECISION :: dh11bdz(-2:n+2,1),dh22bdz(-2:n+2,1)    

       DOUBLE PRECISION :: psn(-2:n+2,-5:m+6),psnp(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: psng(-2:n+2,-5:jmax),vrtgl(-2:n+2,-5:jmax)
       DOUBLE PRECISION :: ometag(-2:n+2,1:jmax),psnpg(-2:n+2,-5:jmax)
       DOUBLE PRECISION :: vrto(-2:n+2,-11:m+12),vrtn(-2:n+2,-11:m+12)
       DOUBLE PRECISION :: vrt_base(-2:n+2,-5:m+6),vrt_rk(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: f_q(1:4),fac(1:4),temp(-2:n+2,-5:m+6),tem(-2:n+2,-11:m+12)

       REAL *8 :: omega0,delt,omega1,sf,omega,gee  
       REAL *8 :: rlx,jrel,tolerance,tolerance2,step1,step2,dt1
       REAL *8 :: sample_dt,tmax,dr1,dt,re,ro,theta_mean,t,time,prestime

       DOUBLE PRECISION :: prn(-2:n+2,-5:jmax),pro(-2:n+2,-5:m+6)
       DOUBLE PRECISION :: u(-2:n+2,-5:jmax),v(-2:n+2,-5:jmax)
       DOUBLE PRECISION :: vhold(-2:n+2,-5:jmax),vhnew(-2:n+2,-5:jmax)
       DOUBLE PRECISION :: bvaluesp(-2:n+2),bvalueop(-2:n+2)

       DOUBLE PRECISION :: rhsp(-2:n+2,-5:jmax),respx(0:n,0:plimit)
       DOUBLE PRECISION :: omxi(-2:n+2,-5:m+6),ometa(-2:n+2,-11:m+12)
       DOUBLE PRECISION :: a(-2:n+2),b(-2:n+2),c(-2:n+2),q1c(-2:n+2)
       DOUBLE PRECISION :: w1c(-2:n+2),r1c(-2:n+2),p1c(-2:n+2),u1c(-2:n+2)
       DOUBLE PRECISION :: f(-2:n+2),v1c(-2:n+2),vrtn_fil(-2:n+2,1:m+12)
       DOUBLE PRECISION :: p(-11:m+12),q(-11:m+12),r(-11:m+12)
       DOUBLE PRECISION :: phi(-11:m+12),d(-11:m+12),w(-11:m+12)
       
       REAL *8, PARAMETER :: xl = 0.3793894912d0               ! #
       REAL *8, PARAMETER :: xr = 0.3793894912d0               ! #
       REAL *8, PARAMETER :: am1 = -0.787786895d0              ! #
       REAL *8, PARAMETER :: ap1 = 0.787786895d0               ! #
       REAL *8, PARAMETER :: am2 = -0.045801298d0              ! #
       REAL *8, PARAMETER :: ap2 = 0.045801298d0               ! #
       REAL *8, PARAMETER :: beta1 = 0.06d0, beta2 = 0.11d0     ! #

       REAL *8, PARAMETER :: one_by_six = 1.0d0/6.0d0
       REAL *8, PARAMETER :: one_by_three = 1.0d0/3.0d0

       REAL *8 :: AP(6,n,plimit),APM(6,n,plimit)
       DOUBLE PRECISION :: cl,cd,cm,ke,ze1
       DOUBLE PRECISION :: Abig(6,-2:n+2,-5:m+6),resx(-2:n+2,-5:m+6)

       END MODULE VARIABLES  

       MODULE VARIABLESf
       USE VARIABLES
       IMPLICIT NONE
       SAVE

       INTEGER *8, PARAMETER :: l=(n-1)*(m)
       DOUBLE PRECISION :: u4(l),w4(l)
       DOUBLE PRECISION :: a4(l),b4(l),c4(l),d4(l),e4(l)
       DOUBLE PRECISION :: p4(l),q4(l),r4(l)

       END MODULE VARIABLESf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
       
       PROGRAM MAGNUS

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,k,rhsnum,status,process,ip1,im1,jp1,jm1
       INTEGER :: source,jmin,root,temp1,npr,kco,jstart,jend,ibm
       REAL *8 :: dpsi,psid(0:n),zid,fx,bt
       DOUBLE PRECISION :: t1

       CHARACTER*14 :: filerhs, filestr, filevrt
       CHARACTER :: dummyc,g

       dimension status(MPI_STATUS_SIZE)

       CALL MPI_INIT(error)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npr,error)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)

       open(1,file='DATA')
            read(1,*) re,tmax,dt,sample_dt
            read(1,*) ro,dr1,str
            read(1,*) para,rlx,jrel,theta_mean
            read(1,*) irs,istart
            read(1,*) omega0,delt,omega1,sf,gee
            read(1,*) tolerance, tolerance2,step1, step2, dt1
       close(1)
 
       filerhs(1:3) = "RHS" ; filestr(1:3) = "SVP" 

!       betaL=Fp*Re 
       itrcount = 1
       CALL INITIA

!        if(rank.eq.0)write(*,*)"initia"

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       DO WHILE(t.LE.tmax) !Main Time Loop Starts

        IF (rank.eq.0) THEN ! unsteady blowing & suction at
        CALL BOPSI
        write(*,*) t,"TIME"
        END IF

        ivrt= ivrt+1
        CALL RUNGE

!        IF(((MOD(ivrt,100)).EQ.0).and.(rank.EQ.0)) THEN
!        OPEN(1,FILE='vort.data',POSITION='append')
!           WRITE(1,*) t,vrtn(700,1),vrtn(800,1),vrtn(900,1),vrtn(1000,1),vrtn(1100,1),vrtn(1200,1),&
!                vrtn(1300,1),vrtn(1400,1),vrtn(1500,1),vrtn(1600,1),vrtn(1700,1),vrtn(1800,1),&
!                vrtn(1900,1),vrtn(2000,1),vrtn(2100,1),vrtn(2200,1),vrtn(2300,1),vrtn(2400,1),&
!                vrtn(2500,1),vrtn(2600,1),vrtn(2700,1),vrtn(2800,1),vrtn(2900,1),vrtn(3000,1),&
!                vrtn(3100,1),vrtn(3200,1),vrtn(3300,1),vrtn(3400,1),vrtn(3500,1),vrtn(3600,1),&
!                vrtn(3700,1),vrtn(3800,1)
!        CLOSE(1)
!       
!        OPEN(1,FILE='vort_near_exc.data',POSITION='append')
!           WRITE(1,*) t,vrtn(470,1),vrtn(500,1),vrtn(525,1),vrtn(550,1),vrtn(575,1),vrtn(600,1),&
!           vrtn(625,1),vrtn(650,1),&
!           vrtn(670,1),vrtn(695,1),vrtn(725,1),vrtn(750,1),vrtn(775,1),vrtn(825,1),&
!           vrtn(850,1),vrtn(875,1)
!        CLOSE(1)
!        ENDIF

        IF (dabs(t-prestime).LT.dt/2) THEN  !! PRESSURE FILE WRITE SECTION !!
          ! write(*,*) t,"TIME"
           time = t
        
           CALL MPI_BARRIER(mpi_comm_world,ierr)

           jstart=1
           jend=m-1
           root=0 !send the data from each process to root and create datafile

           IF (rank.EQ.nproc-1) jend = m

           IF (rank.NE.0) THEN
              temp1=rank               
              DO j=jstart,jend
                call mpi_send(psn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
              END DO
           END IF

           IF (rank.EQ.0) THEN
              DO j=jstart,jend
                 DO i=1,n
                   psng(i,j)=psn(i,j)
                 END DO
              END DO

              DO temp1=1,nproc-1
                 jend = m-1
                 IF (temp1.EQ.nproc-1) jend = m
                 source=temp1
                 DO j=jstart,jend
                    call mpi_recv(psng(1:n,j+source*(m-1)),n,mpi_double_precision, &
                    source,temp1,mpi_comm_world,status,ierr)
                 END DO
              END DO
           END IF          

           ! VORTICITY DATA TRANSFER !
           jstart=1
           jend=m-1
           IF (rank.EQ.nproc-1) jend = m
           root = 0

           IF (rank.NE.0) THEN
              temp1=rank
              DO j=jstart,jend
               call mpi_send(vrtn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
              END DO
           END IF

           IF (rank.EQ.0) THEN
              DO j=jstart,jend
                 DO i=1,n
                    vrtgl(i,j)=vrtn(i,j)
                 END DO
              END DO

              DO temp1=1,nproc-1
                 jend = m-1
                 IF (temp1.EQ.nproc-1) jend = m
                 source=temp1
                 DO j=jstart,jend
                    call mpi_recv(vrtgl(1:n,j+source*(m-1)),n,mpi_double_precision, &
                    source,temp1,mpi_comm_world,status,ierr)
                 END DO
              END DO
           END IF

          jstart=1
          jend=m-1
          root=0 !send the data from each process to root and create datafile

          IF(rank.EQ.nproc-1) jend = m

          IF(rank.NE.0) THEN
            temp1=rank               
            DO j=jstart,jend
             call mpi_send(psnp(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
            END DO
          END IF

          IF(rank.EQ.0) THEN
            DO j=jstart,jend
               DO i=1,n
                 psnpg(i,j)=psnp(i,j)
               END DO
            END DO

            DO temp1=1,nproc-1
               jend = m-1
               IF(temp1.EQ.nproc-1) jend = m
               source=temp1
               DO j=jstart,jend
                  call mpi_recv(psnpg(1:n,j+source*(m-1)),n,mpi_double_precision, &
                  source,temp1,mpi_comm_world,status,ierr)
               END DO
            END DO
          END IF          

!          IF (rank.EQ.0) THEN
!              CALL PRESSURE
!              prestime = prestime+sample_dt
!          END IF

          CALL MPI_BARRIER(mpi_comm_world,ierr) 
        END IF


        IF((MOD(ivrt,irs)).EQ.0) THEN !! TEMP. FILE WRITE SECTION !!  

            CALL MPI_BARRIER(mpi_comm_world,ierr)

            jstart=1
            jend=m-1
            root=0 !send the data from each process to root and create datafile

            IF(rank.EQ.nproc-1) jend = m

            IF(rank.NE.0) THEN
              temp1=rank               
              DO j=jstart,jend
               call mpi_send(psn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
              END DO
            END IF

            IF(rank.EQ.0) THEN
              DO j=jstart,jend
                 DO i=1,n
                   psng(i,j)=psn(i,j)
                 END DO
              END DO

              DO temp1=1,nproc-1
                 jend = m-1
                 IF(temp1.EQ.nproc-1) jend = m
                 source=temp1
                 DO j=jstart,jend
                    call mpi_recv(psng(1:n,j+source*(m-1)),n,mpi_double_precision, &
                    source,temp1,mpi_comm_world,status,ierr)
                 END DO
              END DO
            END IF          

!!!!!!!!!!!!!!!!!! VORTICITY DATA TRANSFER !!!!!!!!!!!!!!!!!!!!!!
          jstart=1
          jend=m-1
          IF(rank.EQ.nproc-1) jend = m
          root = 0

          IF(rank.NE.0) THEN
             temp1=rank
             DO j=jstart,jend
               call mpi_send(vrtn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
             END DO
          END IF

          IF(rank.EQ.0) THEN
             DO j=jstart,jend
                DO i=1,n
                   vrtgl(i,j)=vrtn(i,j)
                END DO
             END DO

          DO temp1=1,nproc-1
               jend = m-1
               IF(temp1.EQ.nproc-1) jend = m
             source=temp1
             DO j=jstart,jend
                call mpi_recv(vrtgl(1:n,j+source*(m-1)),n,mpi_double_precision, &
                source,temp1,mpi_comm_world,status,ierr)
             END DO
          END DO
          END IF

          IF(rank.EQ.0) THEN
          WRITE(*,*) 'WRITING IN TEMP FILE'
          OPEN(2,file='FTEMPSVP',access='sequential')

          WRITE(2,*) "#",t ,ivrt
          WRITE(2,*) "#",re,omega,gee,para,prestime
          WRITE(2,*) "TITLE = STREAM VORTICITY  CONTOURS"
          WRITE(2,*) "VARIABLES = X,Y,Psi,Vorticity "
          WRITE(2,*) "ZONE T=S I=",imax," J=",jmax," F = POINT"

          jstart =1 
          DO k=0,nproc-1
             jend = m-1
             IF(k.EQ.nproc-1) jend = m
             DO j=jstart,jend
                DO i =1,n
                   WRITE(2,*)xgo(i,j+k*(m-1)),ygo(i,j+k*(m-1)),psng(i,j+k*(m-1)),vrtgl(i,j+k*(m-1)),0.0d0
                END DO
             END DO
          END DO
          CLOSE(2)
          write(*,*) "WRITTEN"
          
          END IF

        END IF

        IF((MOD(ivrt,irs+1)).EQ.0) THEN !! TEMP. FILE WRITE SECTION !!  

            CALL MPI_BARRIER(mpi_comm_world,ierr)

            jstart=1
            jend=m-1
            root=0 !send the data from each process to root and create datafile

            IF(rank.EQ.nproc-1) jend = m

            IF(rank.NE.0) THEN
              temp1=rank               
              DO j=jstart,jend
               call mpi_send(psn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
              END DO
            END IF

            IF(rank.EQ.0) THEN
              DO j=jstart,jend
                 DO i=1,n
                   psng(i,j)=psn(i,j)
                 END DO
              END DO

              DO temp1=1,nproc-1
                 jend = m-1
                 IF(temp1.EQ.nproc-1) jend = m
                 source=temp1
                 DO j=jstart,jend
                    call mpi_recv(psng(1:n,j+source*(m-1)),n,mpi_double_precision, &
                    source,temp1,mpi_comm_world,status,ierr)
                 END DO
              END DO
            END IF          

!!!!!!!!!!!!!!!!!! VORTICITY DATA TRANSFER !!!!!!!!!!!!!!!!!!!!!!
          jstart=1
          jend=m-1
          IF(rank.EQ.nproc-1) jend = m
          root = 0

          IF(rank.NE.0) THEN
             temp1=rank
             DO j=jstart,jend
               call mpi_send(vrtn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
             END DO
          END IF

          IF(rank.EQ.0) THEN
             DO j=jstart,jend
                DO i=1,n
                   vrtgl(i,j)=vrtn(i,j)
                END DO
             END DO

          DO temp1=1,nproc-1
               jend = m-1
               IF(temp1.EQ.nproc-1) jend = m
             source=temp1
             DO j=jstart,jend
                call mpi_recv(vrtgl(1:n,j+source*(m-1)),n,mpi_double_precision, &
                source,temp1,mpi_comm_world,status,ierr)
             END DO
          END DO
          END IF

          IF(rank.EQ.0) THEN
          WRITE(*,*) 'WRITING IN TEMP FILE'
          OPEN(2,file='FTEMPSVP1',access='sequential')

          WRITE(2,*) "#",t ,ivrt
          WRITE(2,*) "#",re,omega,gee,para,prestime
          WRITE(2,*) "TITLE = STREAM VORTICITY  CONTOURS"
          WRITE(2,*) "VARIABLES = X,Y,Psi,Vorticity "
          WRITE(2,*) "ZONE T=S I=",imax," J=",jmax," F = POINT"

          jstart =1 
          DO k=0,nproc-1
             jend = m-1
             IF(k.EQ.nproc-1) jend = m
             DO j=jstart,jend
                DO i =1,n
                   WRITE(2,*)xgo(i,j+k*(m-1)),ygo(i,j+k*(m-1)),psng(i,j+k*(m-1)),vrtgl(i,j+k*(m-1)),0.0d0
                END DO
             END DO
          END DO
          CLOSE(2)
          END IF

        END IF

        t1 = step1
        DO WHILE(t1.LE.step2)
          IF( dabs(t1-t).LE.dt/2.0d0) THEN

            CALL MPI_BARRIER(mpi_comm_world,ierr)

            jstart=1
            jend=m-1
            root=0 !send the data from each process to root and create datafile

            IF(rank.EQ.nproc-1) jend = m

            IF(rank.NE.0) THEN
              temp1=rank               
              DO j=jstart,jend
               call mpi_send(psn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
              END DO
            END IF

            IF(rank.EQ.0) THEN
              DO j=jstart,jend
                 DO i=1,n
                   psng(i,j)=psn(i,j)
                 END DO
              END DO

              DO temp1=1,nproc-1
                 jend = m-1
                 IF(temp1.EQ.nproc-1) jend = m
                 source=temp1
                 DO j=jstart,jend
                    call mpi_recv(psng(1:n,j+source*(m-1)),n,mpi_double_precision, &
                    source,temp1,mpi_comm_world,status,ierr)
                 END DO
              END DO
            END IF          

!!!!!!!!!!!!!!!!!! VORTICITY DATA TRANSFER !!!!!!!!!!!!!!!!!!!!!!
          jstart=1
          jend=m-1
          IF(rank.EQ.nproc-1) jend = m
          root = 0

          IF(rank.NE.0) THEN
             temp1=rank
             DO j=jstart,jend
               call mpi_send(vrtn(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
             END DO
          END IF

          IF(rank.EQ.0) THEN
             DO j=jstart,jend
                DO i=1,n
                   vrtgl(i,j)=vrtn(i,j)
                END DO
             END DO

          DO temp1=1,nproc-1
               jend = m-1
               IF(temp1.EQ.nproc-1) jend = m
             source=temp1
             DO j=jstart,jend
                call mpi_recv(vrtgl(1:n,j+source*(m-1)),n,mpi_double_precision, &
                source,temp1,mpi_comm_world,status,ierr)
             END DO
          END DO
          END IF
!!!!!!!!!!!!!!!!!! OLD psn data transfer !!!!!!!!!!!!!!!
          jstart=1
          jend=m-1
          root=0 !send the data from each process to root and create datafile

          IF(rank.EQ.nproc-1) jend = m

          IF(rank.NE.0) THEN
            temp1=rank               
            DO j=jstart,jend
             call mpi_send(psnp(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
            END DO
          END IF

          IF(rank.EQ.0) THEN
            DO j=jstart,jend
               DO i=1,n
                 psnpg(i,j)=psnp(i,j)
               END DO
            END DO

            DO temp1=1,nproc-1
               jend = m-1
               IF(temp1.EQ.nproc-1) jend = m
               source=temp1
               DO j=jstart,jend
                  call mpi_recv(psnpg(1:n,j+source*(m-1)),n,mpi_double_precision, &
                  source,temp1,mpi_comm_world,status,ierr)
               END DO
            END DO
          END IF          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          IF(rank.EQ.0) THEN
!             rhsnum = int(t*100 + 0.01d0)
!             filerhs(4:4) = ACHAR(48 + int(rhsnum/10000)); rhsnum = mod(rhsnum,10000)
!             filerhs(5:5) = ACHAR(48 + int(rhsnum/1000)); rhsnum = mod(rhsnum,1000)
!             filerhs(6:6) = ACHAR(48 + int(rhsnum/100)); rhsnum = mod(rhsnum,100)
!             filerhs(7:7) = "."
!             filerhs(8:8) = ACHAR(48+int(rhsnum/10)); filerhs(9:9) = ACHAR(48 + mod(rhsnum,10))
!             filerhs(10:13) = ".dat"; filestr(4:13) = filerhs(4:13)

        filerhs(01:01) = ACHAR(48 + int(t / 100) - 10 * int(t / 1000))
        filerhs(02:02) = ACHAR(48 + int(t / 10) - 10 * int(t / 100))
        filerhs(03:03) = ACHAR(48 + int(t / 1) - 10 * int(t / 10))
        filerhs(04:04) = "."
        filerhs(05:05) = ACHAR(48 + int(t * 10) - 10 * int(t * 1))
        filerhs(06:06) = ACHAR(48 + int(t * 100) - 10 * int(t * 10))
        filerhs(07:07) = ACHAR(48 + int(t * 1000) - 10 * int(t * 100))
        filerhs(08:08) = ACHAR(48 + int(t * 10000) - 10 * int(t * 1000))
        filerhs(09:09) = ACHAR(48 + int(t * 100000) - 10 * int(t * 10000))
        filerhs(10:10) = ACHAR(48 + int(t * 1000000) - 10 * int(t * 100000))
        filerhs(11:14) = ".dat"; filestr(4:14) = filerhs(4:14)
 
             OPEN(3,file=filerhs)
             WRITE(3,3)g,t,ivrt
             WRITE(3,4)g,re,omega,gee,para,prestime
             WRITE(3,*) "TITLE = STREAM VORTICITY CONTOURS"
             WRITE(3,*) "VARIABLES = X,Y,Psi,Vorticity,OLD_psi "
             WRITE(3,*) "ZONE T=S I=",n," J=",m," F = POINT"
             WRITE(3,10) t
 
             jstart =1 
             DO k=0,nproc-1
                jend = m-1
                IF(k.EQ.nproc-1) jend = m
                DO j=jstart,jend
                   DO i =1,n
                      WRITE(3,*)xgo(i,j+k*(m-1)),ygo(i,j+k*(m-1)),psng(i,j+k*(m-1)),vrtgl(i,j+k*(m-1)),psnpg(i,j+k*(m-1))
                   END DO
                END DO
             END DO
             CLOSE(3)
          END IF
          call system("gzip "//filerhs)
          END IF
           t1 = t1 + dt1
        END DO

        t = t + dt
        itrcount = itrcount + 1
6         FORMAT(1x,' SOLN AT TIME=',e12.5,3x,'ITR # ',i9)
3         FORMAT(T1,A1,1X,E13.6,1X,I10)
4         FORMAT(T1,A1,1X,F7.1,1X,F4.1,1X,F3.1,1X,F4.1,1X,E13.6)
10        FORMAT(2X,'AUXDATA time =  "',F6.2,' "')    
       ENDDO 
    
       CALL MPI_FINALIZE(error)
   
       END PROGRAM MAGNUS
 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE INITIA

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,pfios,ip1,im1,jp1,jm1,jstart,jend,jj,stride,kount,root,k,status
       integer :: source,temp1,jstartN,jendN
       DOUBLE PRECISION :: rj, theta, asq, sum, psp,yy,chn,dtheta,radio,mssg
       DOUBLE PRECISION :: x_eta(n,-5:m+6),x_zi(n,-5:m+6),dum
       DOUBLE PRECISION :: y_eta(n,-5:m+6),y_zi(n,-5:m+6),psidum
       CHARACTER :: g, dummyc
       DIMENSION status(MPI_STATUS_SIZE)

       dzi=1.0d0/dble(n-1); deta=1.0d0/dble(jmax-1)
       prestime = dt
       pi=4*datan(1.0d0); g="#"

       WRITE(*,1) re,tmax,dt
1      FORMAT(1x,'Re=',E12.5,1x,'tmax=',E12.5,1x,'dt=',E12.5)

       OPEN(1, file = 'FINALGRID.dat') ! Every processor will read this file!
       READ(1,*)
       READ(1,*)
       DO j=1,jmax
        DO i=1,n
         READ(1,*) Xgo(i,j),Ygo(i,j)
        END DO
       END DO

       jstart = 1
       jend = m
       DO i=1,n
          DO j=jstart,jend
             jj=(m-1)*rank+j
             xco(i,j)=xgo(i,jj)
             yco(i,j)=ygo(i,jj) !Corresponding XY coordinates are assigned!
          END DO
       END DO

       if(rank.eq.nproc-1)  write(*,*)Xco(nin,m),Yco(nin,m),nin
       if(rank.eq.nproc-1)  write(*,*)Xco(nout,m),Yco(nout,m),nout

       DO j=jstart,jend
         xco(n,j)=xco(1,j);xco(n+1,j)=xco(2,j);xco(n+2,j)=xco(3,j)
         xco(0,j)=xco(n-1,j); xco(-1,j)=xco(n-2,j); xco(-2,j)=xco(n-3,j)

         yco(n,j)=yco(1,j);yco(n+1,j)=yco(2,j);yco(n+2,j)=yco(3,j)
         yco(0,j)=yco(n-1,j); yco(-1,j)=yco(n-2,j);yco(-2,j)=yco(n-3,j)
       END DO

!!!!!!!!!!!!!!!!!!!!!!!!!DATA TRANSFER BLOCK!!!!!!!!!!!!!!!!!!!!!!
       temp = xco
       call msg_pass
       do i=-2,n+2
          xco(i,0) = temp(i,0); xco(i,m+1) = temp(i,m+1)
          xco(i,-1) = temp(i,-1);xco(i,m+2) = temp(i,m+2)
          xco(i,-2) = temp(i,-2);xco(i,m+3) = temp(i,m+3)
          xco(i,-3) = temp(i,-3);xco(i,m+4) = temp(i,m+4)
          xco(i,-4) = temp(i,-4);xco(i,m+5) = temp(i,m+5)
          xco(i,-5) = temp(i,-5);xco(i,m+6) = temp(i,m+6)
       end do

       temp = yco
       call msg_pass
       do i=-2,n+2
          yco(i,0) = temp(i,0); yco(i,m+1) = temp(i,m+1)
          yco(i,-1) = temp(i,-1);yco(i,m+2) = temp(i,m+2)
          yco(i,-2) = temp(i,-2);yco(i,m+3) = temp(i,m+3)
          yco(i,-3) = temp(i,-3);yco(i,m+4) = temp(i,m+4)
          yco(i,-4) = temp(i,-4);yco(i,m+5) = temp(i,m+5)
          yco(i,-5) = temp(i,-5);yco(i,m+6) = temp(i,m+6)
       end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       IF(rank.eq.0) THEN
!         s_zi(1)=0.0d0 
!         DO i=2,n
!            j=1
!            ds_zi(i)=dsqrt((xco(i,j)-xco(i-1,j))**2+(yco(i,j)-yco(i-1,j))**2)
!            s_zi(i)=s_zi(i-1)+ds_zi(i)
!         END DO
!
!         DO i=2,N-1
!            IF (xco(i,1).GE.expos) THEN
!               ibg=i
!               EXIT
!            ENDIF
!         ENDDO
!
!         ibg1= ibg-11
!         ibg2= ibg+11
!         xx1= s_zi(ibg1)
!         xx2= s_zi(ibg2)
!         xst= 0.50d0*(xx1+xx2)
!       END IF


       jstart =1 ; jend = m
       IF (rank.EQ.0) jstart = 2
       IF (rank.EQ.nproc-1) jend = m-1
 
       DO j=jstart,jend
         DO i=2,n
           im1=i-1
           ip1=i+1
           IF(i.EQ.n) ip1=2
           jm1=j-1
           jp1=j+1
           x_zi(i,j)=(Xco(ip1,j)-Xco(im1,j))*0.50d0/dzi
           y_zi(i,j)=(Yco(ip1,j)-Yco(im1,j))*0.50d0/dzi
           x_eta(i,j)=(Xco(i,jp1)-Xco(i,jm1))*0.50d0/deta
           y_eta(i,j)=(Yco(i,jp1)-Yco(i,jm1))*0.50d0/deta
         END DO
       END DO

       DO i=2,n
         IF(rank.EQ.0) THEN
            im1=i-1
            ip1=i+1
            IF(i.EQ.n) ip1=2
            x_zi(i,1)=(Xco(ip1,1)-Xco(im1,1))*0.50d0/dzi
            y_zi(i,1)=(Yco(ip1,1)-Yco(im1,1))*0.50d0/dzi
            x_eta(i,1)=(Xco(i,2)-Xco(i,1))/deta
            y_eta(i,1)=(Yco(i,2)-Yco(i,1))/deta
         END IF

         IF(rank.EQ.nproc-1) THEN
            im1=i-1
            ip1=i+1
            IF(i.EQ.n) ip1=2
            x_zi(i,m)=(Xco(ip1,m)-Xco(im1,m))*0.50d0/dzi
            y_zi(i,m)=(Yco(ip1,m)-Yco(im1,m))*0.50d0/dzi
            x_eta(i,m)=(Xco(i,m)-Xco(i,m-1))/deta
            y_eta(i,m)=(Yco(i,m)-Yco(i,m-1))/deta
         END IF
       END DO

       DO j=1,m
         x_zi(1,j)=x_zi(n,j)
         y_zi(1,j)=y_zi(n,j)
         x_eta(1,j)=x_eta(n,j)
         y_eta(1,j)=y_eta(n,j)
       END DO

       DO j=1,m
         DO i=1,n
           h11(i,j)=dsqrt(x_zi(i,j)*x_zi(i,j)+y_zi(i,j)*y_zi(i,j))
           h22(i,j)=dsqrt(x_eta(i,j)*x_eta(i,j)+y_eta(i,j)*y_eta(i,j))
         END DO
       END DO 

       IF(RANK.EQ.0) THEN
         DO i=1,n 
            dh11bdn(i,1)=(h11(i,2)-h11(i,1))/deta
            dh22bdn(i,1)=(h22(i,2)-h22(i,1))/deta
         END DO

!         DO i=2,n-1
!            dh11bdz(i,1)=(h11(i+1,1)-h11(i-1,1))/(2.0d0*dzi)
!            dh22bdz(i,1)=(h22(i+1,1)-h22(i-1,1))/(2.0d0*dzi)
!         END DO

!         dh11bdz(1,1)=(h11(2,1)-h11(n-1,1))/(2.0d0*dzi)
!         dh22bdz(1,1)=(h22(2,1)-h22(n-1,1))/(2.0d0*dzi)

!         dh11bdz(n,1)=dh11bdz(1,1)
!         dh22bdz(n,1)=dh22bdz(1,1) 

!         DO i=1,n
!            dh11bdn(i,1)=(4.0d0*h11(i,2)-3.0d0*h11(i,1)-h11(i,3))/(2.0d0*deta)
!            dh22bdn(i,1)=(4.0d0*h22(i,2)-3.0d0*h22(i,1)-h22(i,3))/(2.0d0*deta)
!         END DO

         DO i=2,n-1
            dh11bdz(i,1)=(h11(i+1,1)-h11(i-1,1))/(2.0d0*dzi)
            dh22bdz(i,1)=(h22(i+1,1)-h22(i-1,1))/(2.0d0*dzi)
         END DO

         dh11bdz(1,1)=(h11(2,1)-h11(n-1,1))/(2.0d0*dzi)
         dh22bdz(1,1)=(h22(2,1)-h22(n-1,1))/(2.0d0*dzi)

         dh11bdz(n,1)=dh11bdz(1,1)
         dh22bdz(n,1)=dh22bdz(1,1) 
       END IF
 
       DO j=1,m
         h11(0,j)  = h11(n-1,j)
         h11(-1,j) = h11(n-2,j)
         h11(-2,j) = h11(n-3,j)
         h11(n,j)  = h11(1,j)
         h11(n+1,j)= h11(2,j)
         h11(n+2,j)= h11(3,j)

         h22(0,j)  = h22(n-1,j)
         h22(-1,j) = h22(n-2,j)
         h22(-2,j) = h22(n-3,j)
         h22(n,j)  = h22(1,j)
         h22(n+1,j)= h22(2,j)
         h22(n+2,j)= h22(3,j)
       END DO

!!!!!!!!!!!!!!!!!!!!!!!!!DATA TRANSFER BLOCK!!!!!!!!!!!!!!!!!!!!!!
       temp = h11
       CALL msg_pass
       DO i=-2,n+2
          h11(i,0) = temp(i,0); h11(i,m+1) = temp(i,m+1)
          h11(i,-1) = temp(i,-1);h11(i,m+2) = temp(i,m+2)
          h11(i,-2) = temp(i,-2);h11(i,m+3) = temp(i,m+3)
          h11(i,-3) = temp(i,-3);h11(i,m+4) = temp(i,m+4)
          h11(i,-4) = temp(i,-4);h11(i,m+5) = temp(i,m+5)
          h11(i,-5) = temp(i,-5);h11(i,m+6) = temp(i,m+6)
       END DO

       temp = h22
       CALL msg_pass
       DO i=-2,n+2
          h22(i,0) = temp(i,0); h22(i,m+1) = temp(i,m+1)
          h22(i,-1) = temp(i,-1);h22(i,m+2) = temp(i,m+2)
          h22(i,-2) = temp(i,-2);h22(i,m+3) = temp(i,m+3)
          h22(i,-3) = temp(i,-3);h22(i,m+4) = temp(i,m+4)
          h22(i,-4) = temp(i,-4);h22(i,m+5) = temp(i,m+5)
          h22(i,-5) = temp(i,-5);h22(i,m+6) = temp(i,m+6)
       END DO
!!!!!!!!!!!!!!!!!! Transfer of h11 and h22 to h11g and h22g !!!!!!!!!!!!!!!!

       CALL MPI_BARRIER(mpi_comm_world,ierr)

       jstart=1
       jend=m-1
       root=0 !send the data from each process to root and create datafile

       IF(rank.EQ.nproc-1) jend=m

       IF (rank.NE.0) THEN
          temp1=rank               
          DO j=jstart,jend
             call mpi_send(h11(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
          END DO
       END IF

       IF (rank.EQ.0) THEN
          DO j=jstart,jend
             DO i=1,n
                h11g(i,j)=h11(i,j)
             END DO
          END DO

          DO temp1=1,nproc-1
             jend = m-1
             IF (temp1.EQ.nproc-1) jend=m
             source=temp1
             DO j=jstart,jend
                call mpi_recv(h11g(1:n,j+source*(m-1)),n,mpi_double_precision, &
                source,temp1,mpi_comm_world,status,ierr)
             END DO
          END DO
       END IF          
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       jstart=1
       jend=m-1
       root=0 !send the data from each process to root and create datafile

       IF(rank.EQ.nproc-1) jend = m

       IF (rank.NE.0) THEN
          temp1=rank               
          DO j=jstart,jend
             call mpi_send(h22(1:n,j),n,mpi_double_precision,root,temp1,mpi_comm_world,ierr)
          END DO
       END IF

       IF (rank.EQ.0) THEN
          DO j=jstart,jend
             DO i=1,n
                h22g(i,j)=h22(i,j)
             END DO
          END DO

          DO temp1=1,nproc-1
             jend = m-1
             IF (temp1.EQ.nproc-1) jend = m
             source=temp1
             DO j=jstart,jend
                call mpi_recv(h22g(1:n,j+source*(m-1)),n,mpi_double_precision, &
                source,temp1,mpi_comm_world,status,ierr)
             END DO
          END DO
       END IF

       IF (rank.EQ.0) THEN
          DO j = 1,jmax
             h11g(n,j) = h11g(1,j)
             h11g(0,j) = h11g(n-1,j)
             h11g(-1,j)= h11g(n-2,j)
             h11g(-2,j)= h11g(n-3,j)
             h11g(n+1,j)= h11g(2,j)
             h11g(n+2,j)= h11g(3,j)

             h22g(n,j) = h22g(1,j)
             h22g(0,j) = h22g(n-1,j)
             h22g(-1,j)= h22g(n-2,j)
             h22g(-2,j)= h22g(n-3,j)
             h22g(n+1,j)= h22g(2,j)
             h22g(n+2,j)= h22g(3,j)
          END DO     
       END IF   
 
!!!!!!!!!!!!!!!!ADDITIONAL BLOCK TO REDUCE COMPUTATIONS!!!!!!!!!!!

       dzi_sq_inv = 1.0d0/(dzi*dzi)
       deta_sq_inv = 1.0d0/(deta*deta)
       two_dzi_inv = 1.0d0/(2.0d0*dzi)
       two_deta_inv = 1.0d0/(2.0d0*deta)
       re_inv = 1.0d0/re

       DO i=2,n
          ip1=i+1;im1=i-1
          IF (i==n) ip1=2
          jstartN = 1
          jendN = m
          IF (rank.EQ.0) jstartN = 2
          IF (rank.EQ.nproc-1) jendN = m-1

          DO j=jstartN,jendN
             h11h22_inv(i,j) = 1.0d0/(h11(i,j)*h22(i,j))

             h22bh11_im1(i,j) = 0.5*(h22(im1,j)/h11(im1,j)+h22(i,j)/h11(i,j))
             h22bh11_ip1(i,j) = 0.5*(h22(ip1,j)/h11(ip1,j)+h22(i,j)/h11(i,j))
             h22bh11(i,j) = -(h22bh11_im1(i,j)+h22bh11_ip1(i,j))

             h11bh22_jm1(i,j) = 0.5*(h11(i,j)/h22(i,j)+h11(i,j-1)/h22(i,j-1))
             h11bh22_jp1(i,j) = 0.5*(h11(i,j)/h22(i,j)+h11(i,j+1)/h22(i,j+1))
             h11bh22(i,j) = -(h11bh22_jm1(i,j)+h11bh22_jp1(i,j))
          END DO
       END DO

!!!!!!!!!!BLOCK COMPLETE!!!!!!!!!!!!!!!!!!!!!!

       CALL POISSONCOEFF

       IF(istart.EQ.0) THEN

       jstart=1;jend=m
       IF (rank.EQ.0) jstart = 2

       IF(rank.EQ.0) THEN
          DO i=1,n
             psn(i,1) = 15.0d0
             psnp(i,1) = psn(i,1)
          END DO
       END IF       

       DO i=1,n
          DO j=jstart,jend
             psn(i,j)=15.0d0+Yco(i,j)*dcos(theta_mean)+Xco(i,j)*dsin(theta_mean)
             psnp(i,j) = psn(i,j)
          END DO
       END DO


       DO i=1,n
         DO j=jstart,jend
          vrtn(i,j) = 0.0d0; vrto(i,j) = vrtn(i,j)
         END DO
       END DO

      DO j=jstart,jend
         vrtn(0,j) = vrtn(n-1,j);vrtn(-1,j) = vrtn(n-2,j)
         vrtn(-2,j) = vrtn(n-3,j);vrtn(n,j) = vrtn(1,j)
         vrtn(n+1,j) = vrtn(2,j);vrtn(n+2,j) = vrtn(3,j)

         vrto(0,j) = vrto(n-1,j);vrto(-1,j) = vrto(n-2,j)
         vrto(-2,j) = vrto(n-3,j);vrto(n,j) = vrto(1,j)
         vrto(n+1,j) = vrto(2,j);vrto(n+2,j) = vrto(3,j)

         psn(n,j)=psn(1,j);psn(n+1,j)=psn(2,j);psn(n+2,j)=psn(3,j)
         psn(0,j)=psn(n-1,j); psn(-1,j)=psn(n-2,j); psn(-2,j)=psn(n-3,j)
      END DO

      temp = psn
      call msg_pass
      do i=-2,n+2
      psn(i,0) = temp(i,0); psn(i,m+1) = temp(i,m+1)
      psn(i,-1) = temp(i,-1);psn(i,m+2) = temp(i,m+2)
      psn(i,-2) = temp(i,-2);psn(i,m+3) = temp(i,m+3)
      psn(i,-3) = temp(i,-3);psn(i,m+4) = temp(i,m+4)
      psn(i,-4) = temp(i,-4);psn(i,m+5) = temp(i,m+5)
      psn(i,-5) = temp(i,-5);psn(i,m+6) = temp(i,m+6)
      end do

      CALL VORTICITYPASS
      CALL VORTICITY_OLD_PASS

      CALL STREAMSOLVER 

      DO j=1,m
         psn(n,j)=psn(1,j);psn(n+1,j)=psn(2,j);psn(n+2,j)=psn(3,j)
         psn(0,j)=psn(n-1,j); psn(-1,j)=psn(n-2,j); psn(-2,j)=psn(n-3,j)
      END DO
                 
      temp = psn
      CALL msg_pass
      DO i=-2,n+2
      psn(i,0) = temp(i,0); psn(i,m+1) = temp(i,m+1)
      psn(i,-1) = temp(i,-1);psn(i,m+2) = temp(i,m+2)
      psn(i,-2) = temp(i,-2);psn(i,m+3) = temp(i,m+3)
      psn(i,-3) = temp(i,-3);psn(i,m+4) = temp(i,m+4)
      psn(i,-4) = temp(i,-4);psn(i,m+5) = temp(i,m+5)
      psn(i,-5) = temp(i,-5);psn(i,m+6) = temp(i,m+6)
      END DO
 
      IF(rank.eq.0) THEN
         CALL BOVRT
      END IF
      END IF
     
      IF(istart.NE.0) THEN

        open(2,file='FTEMPSVP',form='formatted',access='sequential')
           print*,'Reading from temp'
           read(2,*)dummyc,t ,ivrt
3          format(T1,A1,1X,E13.6,1X,I10)
           read(2,*)dummyc,re,omega,gee,para,prestime
4          format(T1,A1,1X,F7.1,1X,F4.1,1X,F3.1,1X,F4.1,1X,E13.6)
5          format(//)
           read(2,*)
           read(2,*)
           read(2,*)
           DO j =1,jmax
              DO i=1,n 
                 read(2,*) Xgo(i,j),Ygo(i,j),psng(i,j),vrtgl(i,j),dum
              END DO
           END DO
        close(2)

        jstart = 1
        jend = m

        DO i=1,n
           DO j=jstart,jend
              jj=(m-1)*rank+j

              xco(i,j)=xgo(i,jj)
              yco(i,j)=ygo(i,jj)
              psn(i,j)=psng(i,jj)
              vrtn(i,j)=vrtgl(i,jj)

              vrto(i,j) = vrtn(i,j) 
           END DO
        END DO

        print *,'prestime =',prestime,'Re=',re
        t = t + dt
      END IF

      IF (rank.EQ.0) THEN
         DO i=0,n
            DO j=0,jmax
               prn(i,j)=1.0d0
            END DO
         END DO
      END IF
 
      temp = psn
      call msg_pass
      do i=-2,n+2
         psn(i,0) = temp(i,0); psn(i,m+1) = temp(i,m+1)
         psn(i,-1) = temp(i,-1);psn(i,m+2) = temp(i,m+2)
         psn(i,-2) = temp(i,-2);psn(i,m+3) = temp(i,m+3)
         psn(i,-3) = temp(i,-3);psn(i,m+4) = temp(i,m+4)
         psn(i,-4) = temp(i,-4);psn(i,m+5) = temp(i,m+5)
         psn(i,-5) = temp(i,-5);psn(i,m+6) = temp(i,m+6)
      end do

       IF (rank.EQ.0) THEN
          CALL BOVRT !EXTRA ADDITION TO CALCULATE BOVRT AT IMPULSIVE START!
       END IF

       CALL VORTICITYPASS
       CALL VORTICITY_OLD_PASS

      RETURN

      END SUBROUTINE INITIA

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE BOVRT

       USE VARIABLES  
       IMPLICIT NONE

       INTEGER :: i,j,ip1,im1
       DOUBLE PRECISION :: term1,term2,term3,term4,psn_eta
       double precision :: ftsi,stsi,ttsi,zid,aa,bb,cc

        aa = 5.0d0*15.1875d0
        bb = 4.0d0*35.4375d0
        cc = 3.0d0*20.25d0

        do i = 1,N-1
           j = 1
!           if ((s_zi(i).ge.xx1).and.(s_zi(i).le.xx2)) then
!              if (s_zi(i) < xst) then
!                 zid= (s_zi(i)-xx1)/(xst-xx1)
!                 vrtn(i,j)=Afk*(aa*zid**4-bb*zid**3+cc*zid**2)*dsin(betaL*t)/(xst-xx1)
!                 vrtn(i,j)=vrtn(i,j)+(-1.0d0/(deta*h22(i,j))**2)*(2.0d0*(psn(i,j+1)-psn(i,j)))
!              else
!                 zid= (xx2-s_zi(i))/(xx2-xst)
!                 vrtn(i,j)=Afk*(aa*zid**4-bb*zid**3+cc*zid**2)*dsin(betaL*t)/(xx2-xst)
!                 vrtn(i,j)=vrtn(i,j)+(-1.0d0/(deta*h22(i,j))**2)*(2.0d0*(psn(i,j+1)-psn(i,j)))
!              endif
!           else
              vrtn(i,j)= (-1.0d0/(deta*h22(i,j))**2)*(2.0d0*(psn(i,j+1)-psn(i,j)))
!           endif
        enddo 

       j=1
       DO i=1,n-1
          vrto(i,1)=vrtn(i,1)
       END DO

       vrtn(n,1)=vrtn(1,1); vrtn(0,1)=vrtn(n-1,1)
       vrto(n,1)=vrto(1,1); vrto(0,1)=vrto(n-1,1)

       RETURN
       END SUBROUTINE BOVRT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE BOPSI

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j
       DOUBLE PRECISION :: zid,aa,bb,cc,funct

       aa= 15.1875d0/6.0d0
       bb= 35.4375d0/5.0d0
       cc= 20.25d0/4.0d0

       do i= 1,n
         j= 1
!         if ((s_zi(i).ge.xx1).and.(s_zi(i).le.xx2)) then
!          if(s_zi(i) < xst) then
!           zid= (s_zi(i)-xx1)/(xst-xx1)
!           psn(i,j)= 15.0d0-(xst-xx1)*Afk*(aa*zid**6-bb*zid**5+cc*zid**4)*dsin(betaL*t)
!          else
!           zid= (xx2-s_zi(i))/(xx2-xst)
!           psn(i,j)= 15.0d0-(xx2-xst)*Afk*(aa*zid**6-bb*zid**5+cc*zid**4)*dsin(betaL*t)
!          endif
!         else
          psn(i,j)= 15.0d0
!         endif
       enddo

       
       END SUBROUTINE BOPSI
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE RUNGE

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,ip1,im1,vn,jstart,jend,k,pass1,cnt,freq
       REAL *8 :: victim(-2:n+2),sum1
       REAL :: a0_fil2,a1_fil2,a0_fil4,a1_fil4,a2_fil4

!       a0_fil2 = 0.50d0+2.0d0*alpha_fil
!       a1_fil2 = 0.250d0+alpha_fil
!
!       a0_fil4 = (5.0d0+12.0d0*alpha_fil)/8.0d0
!       a1_fil4 = 0.250d0+alpha_fil
!       a2_fil4 = (-1.0d0+4*alpha_fil)/16.0d0

       f_q(1) = 0.259674099d0
       f_q(2) = 0.499969d0
       fac(1) = 0.0d0
       fac(2) = 0.0d0
       fac(3) = 1.0d0

!       f_q(1)=0.5d0;f_q(2)=0.5d0;f_q(3)=1.0d0;f_q(4)=1.0d0
!       fac(1)=one_by_six;fac(2)=one_by_three;fac(3)=fac(2);fac(4)=fac(1)
        
       IF(rank.eq.0)THEN
         jstart = 1
         jend = m+6
       ELSEIF(rank.eq.nproc-1)THEN
         jstart = -5
         jend =  m
       ELSE
         jstart = -5
         jend = m+6
       END IF

       DO j = jstart,jend
          DO i = -2,n+2
             psnp(i,j) = psn(i,j)
             vrt_base(i,j) = vrtn(i,j)
             vrt_rk(i,j) = 0.0d0
          END DO
          psnp(0,j)=psnp(n-1,j) 
       END DO

       DO pass = 1,3
        CALL PERIODIC_BOUNDARY
        CALL VTE  
        IF(pass.EQ.3) THEN
          IF(rank.EQ.nproc-1) CALL OBC
        END IF
       END DO

!       IF((rank.NE.nproc-1)) THEN 
!       DO vn = 1,m
!         victim2(1:266)=vrtn(1:266,vn)
!         CALL npfilter
!         vrtn(1:266,vn)=victim2(1:266)
!         DO i=528,263,-1
!           victim2(528-i+1) = vrtn(i,vn)
!         END DO
!         CALL npfilter
!         DO i=1,266
!           vrtn(528+1-i,vn)=victim2(i)
!         END DO
!         vrtn(n,vn)=vrtn(1,vn);vrtn(0,vn)=vrtn(n-1,vn)
!       END DO
!      END IF
! 
!      IF(((MOD(ivrt,2)).EQ.0).and.(rank.EQ.0)) THEN
!          do j=1,mlimit
!             do i=1,n
!                vrtn_rhs(i,j) = vrtn(i,j)
!                vrtnf(i,j) = vrtn(i,j)
!             end do
!          vrtnf(0,j) = vrtn(0,j)
!          vrtn_rhs(0,j) = vrtn(0,j)
!          end do
!
!          DO i=1,n-1
!             j=2
!             vrtn_rhs(i,j)=a0_fil2*vrtn(i,j)+0.50d0*a1_fil2*(vrtn(i+1,j)+vrtn(i-1,j)+ &
!             vrtn(i,j-1)+vrtn(i,j+1))
!             j=mlimit-1
!             vrtn_rhs(i,j)=a0_fil2*vrtn(i,j)+0.50d0*a1_fil2*(vrtn(i+1,j)+vrtn(i-1,j)+ &
!             vrtn(i,j-1)+vrtn(i,j+1))
!          END DO
!
!          DO j=2,mlimit-1
!             i=1
!             vrtn_rhs(i,j)=a0_fil2*vrtn(i,j)+0.50d0*a1_fil2*(vrtn(i+1,j)+vrtn(i-1,j)+ &
!             vrtn(i,j-1)+vrtn(i,j+1))
!
!             i=n-1
!             vrtn_rhs(i,j)=a0_fil2*vrtn(i,j)+0.50d0*a1_fil2*(vrtn(i+1,j)+vrtn(i-1,j)+ &
!             vrtn(i,j-1)+vrtn(i,j+1))
!          END DO
!
!          DO i=2,n-2
!             DO j=3,mlimit-2
!                vrtn_rhs(i,j)=a0_fil4*vrtn(i,j)+0.50d0*a1_fil4*(vrtn(i+1,j)+vrtn(i-1,j)+ &
!                vrtn(i,j-1)+vrtn(i,j+1))+0.50d0*a2_fil4*(vrtn(i+2,j)+vrtn(i-2,j)+ &
!                vrtn(i,j-2)+vrtn(i,j+2))
!             END DO
!          END DO
!       

!          CALL FILTER_XY
!
!          do j=2,mlimit-1
!             do i=1,n-1
!                vrtn(i,j) = vrtnf(i,j)
!             end do
!             vrtn(n+1,j) = vrtn(2,j); vrtn(-1,j) = vrtn(n-2,j)
!             vrtn(0,j) = vrtn(n-1,j); vrtn(n,j)=vrtn(1,j)
!          end do
!
!      END IF 
         
        CALL PERIODIC_BOUNDARY

        CALL STREAMSOLVER

        IF (rank.eq.0) THEN
           CALL BOVRT
        END IF

        CALL PERIODIC_BOUNDARY

        RETURN        
      END SUBROUTINE RUNGE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE STREAMSOLVER

       USE VARIABLES
       IMPLICIT NONE
        
       INTEGER :: i,j,k,iconverge,ip1,im1,jstart,jend,jp1,jm1,ii,jj,counter
       INTEGER :: kcount,tag,dest,source,kount,stride,status,kco,process
       REAL *8 :: rho0,rho1,w0,w1cg,betacg,alphacg,transpro,rho1gg,transpro1,transpro2
       REAL *8 :: var123,var1234,var12345,var123g,var1234g,var12345g,sum
       REAL *8,DIMENSION(-2:n+2,-5:m+6) :: pcg0=0.0d0,pcg1=0.0d0,vcg0=0.0d0,vcg1=0.0d0
       REAL *8,DIMENSION(-2:n+2,-5:m+6) :: scg=0.0d0,tcg=0.0d0,y=0.0d0
       REAL *8,DIMENSION(-2:n+2,-5:m+6) :: kinvt=0.0d0,kinvs=0.0d0
       REAL *8,DIMENSION(-2:n+2,-5:m+6) :: z=0.0d0,res0=0.0d0 
       DOUBLE PRECISION :: maxerr,maxerrg

       rho0=1.0d0;alphacg=1.0d0;w0=1.0d0;kcount=0
       counter=1

!       open(20,file='stream_itr.dat')
!       IF(rank.eq.0) THEN
!       sum=0.0d0 
!       DO i=1,n
!          psn(i,1)=psn(i,2)+((h22(i,1)*deta)**2)*vrtn(i,1)/2.0d0
!          sum=sum+psn(i,1) 
!       END DO
!       DO i=1,n
!          psn(i,1)=dble(sum/n)
!       END DO 
!       write(20,*)counter, psn(1,1)
!       psn(n,1)= psn(1,1)
!       END IF
       
       temp = psn 
       CALL PASS_temp          
       DO i=1,n
          psn(i,0)=temp(i,0); psn(i,m)=temp(i,m)
       END DO

       jstart = 1;jend = m
       DO j=jstart,jend
         psn(n,j)=psn(1,j);psn(n+1,j)=psn(2,j);psn(n+2,j)=psn(3,j)
         psn(0,j)=psn(n-1,j); psn(-1,j)=psn(n-2,j); psn(-2,j)=psn(n-3,j)
       END DO

       CALL CONVERGE(iconverge,kcount)

       res0 = resx     
 
       DO WHILE( iconverge.EQ.0 )
           counter=counter+1
!          IF(rank.EQ.0) write(*,*) 'RES_Total = ',res_sum,iconverge
          jstart = 1
          jend = m-1

          IF (rank.EQ.0) jstart = 2

          rho1 = 0.0d0
          DO i = 1,n-1
             DO j = jstart,jend
                rho1 = rho1 + res0(i,j)*resx(i,j)
             END DO
          END DO

          CALL MPI_ALLREDUCE(rho1,rho1gg,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          MPI_COMM_WORLD,error)
       
          betacg = (rho1gg/rho0)*(alphacg/w0)
          DO j=jstart,jend
             DO i=1,n-1
                pcg1(i,j)=resx(i,j)+betacg*(pcg0(i,j)-w0*vcg0(i,j))
             END DO
          END DO
            
          DO j = jstart,jend
             pcg1(0,j) = pcg1(n-1,j)
             pcg1(-1,j) = pcg1(n-2,j)
             pcg1(-2,j) = pcg1(n-3,j)
             pcg1(n,j) = pcg1(1,j)
             pcg1(n+1,j) = pcg1(2,j)
             pcg1(n+2,j) = pcg1(3,j)
          END DO

          temp = pcg1 
          CALL PASS_temp
          DO i=1,n
             pcg1(i,0) = temp(i,0); pcg1(i,m) = temp(i,m)
          END DO

          CALL AMULTIPLY(vcg1,pcg1)

          var123= TRANSPRO(rank,n,m,res0,vcg1)

          CALL MPI_ALLREDUCE(var123,var123g,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          MPI_COMM_WORLD,error)

          alphacg = rho1gg/var123g

          DO j = jstart,jend
             DO i = 1, n-1
                scg(i,j)=resx(i,j)-alphacg*vcg1(i,j)
             END DO
          END DO

          DO j = jstart,jend
             scg(0,j) = scg(n-1,j); scg(n,j) = scg(1,j)
          END DO

          temp = scg
          CALL PASS_temp
          DO i=1,n
             scg(i,0) = temp(i,0); scg(i,m) = temp(i,m)
          END DO

          CALL AMULTIPLY(tcg,scg)
        
          var1234= TRANSPRO(rank,n,m,tcg,tcg)

          CALL MPI_ALLREDUCE(var1234,var1234g,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          MPI_COMM_WORLD,error)          
       
          var12345= TRANSPRO(rank,n,m,tcg,scg)

          CALL MPI_ALLREDUCE(var12345,var12345g,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          MPI_COMM_WORLD,error)

          w1cg = var12345g/var1234g
          DO j = jstart,jend
             DO i = 1, n-1
                psn(i,j) = psn(i,j) + alphacg*pcg1(i,j) + w1cg*scg(i,j)
             END DO
          END DO

!          IF(rank.eq.0) THEN
!          sum=0.0d0
!          DO i=1,n
!          psn(i,1)=psn(i,2)+((h22(i,1)*deta)**2)*vrtn(i,1)/2.0d0
!          sum=sum+psn(i,1)
!          END DO
!          DO i=1,n
!          psn(i,1)=dble(sum/n)
!          END DO
!          write(20,*) counter,psn(1,1)
!          psn(n,1)= psn(1,1)
!          END IF

          temp=psn
          CALL PASS_temp
          DO i=1,n
             psn(i,0)=temp(i,0); psn(i,m)=temp(i,m)
          END DO

          DO j = 0,m
             psn(n,j)=psn(1,j);psn(0,j)=psn(n-1,j);psn(n+1,j)=psn(2,j) 
          END DO

          CALL CONVERGE(iconverge,kcount)       
          pcg0 = pcg1; vcg0 = vcg1; 
          rho0 = rho1gg; w0 = w1cg
       END DO

       RETURN

       END SUBROUTINE STREAMSOLVER

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE CONVERGE(iconverge,kcount)

       USE VARIABLES
       IMPLICIT NONE
        
       INTEGER,INTENT(INOUT) :: iconverge,kcount
       INTEGER :: i,j,ip1,im1,jp1,jm1,jstart,jend,ii,jj
       REAL *8 :: a1,b1,c1,d1,e1,res_max,trash(-2:N+2,-5:m+6),sum

!       IF(rank.eq.0) THEN
!          sum=0.0d0
!          DO i=1,n
!          psn(i,1)=psn(i,2)+((h22(i,1)*deta)**2)*vrtn(i,1)/2.0d0
!          sum=sum+psn(i,1)
!          END DO
!          DO i=1,n
!          psn(i,1)=dble(sum/n)
!          END DO
!          psn(n,1)= psn(1,1)
!       END IF

       kcount = kcount+1
       res_max = 0.0d0
       jstart = 1
       jend = m-1
       IF(rank.EQ.0) jstart = 2
       DO j = jstart,jend
          jp1 = j+1; jm1 = j-1
          DO i = 1,n-1
            ip1 = i+1; im1 = i-1
            a1 = Abig(1,i,j)*psn(i,jm1)
            b1 = Abig(2,i,j)*psn(im1,j)
            c1 = Abig(3,i,j)*psn(i,j)
            d1 = Abig(4,i,j)*psn(ip1,j)
            e1 = Abig(5,i,j)*psn(i,jp1)
            resx(i,j) = -vrtn(i,j)/Abig(6,i,j)-(a1+b1+c1+d1+e1)
            IF(DABS(resx(i,j))>res_max) THEN
              res_max = DABS(resx(i,j))
              ii = i ; jj= j
            END IF
          END DO
       END DO

!       IF(rank.EQ.0) THEN
!       OPEN(1,FILE='rank0_residue.dat',POSITION='append')
!         write(1,*)res_max,ii,jj
!       CLOSE(1)
!       END IF

!       IF(rank.EQ.1) THEN
!       OPEN(1,FILE='rank1_residue.dat',POSITION='append')
!         write(1,*)res_max,ii,jj
!       CLOSE(1)
!       END IF

!       IF(rank.EQ.2) THEN
!       OPEN(1,FILE='rank2_residue.dat',POSITION='append')
!         write(1,*)res_max,ii,jj
!       CLOSE(1)
!       END IF
 
!       IF(rank.EQ.3) THEN
!       OPEN(1,FILE='rank3_residue.dat',POSITION='append')
!         write(1,*)res_max,ii,jj
!       CLOSE(1)
!       END IF

       CALL MPI_ALLREDUCE(res_max,res_sum,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
            MPI_COMM_WORLD,error)
            iconverge = 1
            IF (res_sum.GT.tolerance) iconverge = 0
            !IF(rank.EQ.0) write(*,*) 'RES_Total = ',res_sum,iconverge,kcount
       RETURN
       END SUBROUTINE CONVERGE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       FUNCTION TRANSPRO(rankN,n2,m2,tres1,tres2)
       
       IMPLICIT NONE
 
       INTEGER,INTENT(IN) :: rankN,n2,m2
       REAL *8,INTENT(IN),DIMENSION(-2:n2+2,-5:m2+6) :: tres1, tres2
       INTEGER :: i,j,jstart,jend
       REAL *8 :: transpro 
       
       transpro = 0.0d0

       jstart = 1
       jend = m2-1
       IF(rankN.EQ.0) jstart = 2

       DO i = 1,n2-1
          DO j = jstart,jend
             transpro = transpro + tres1(i,j)*tres2(i,j)
          END DO
       END DO
       
       END FUNCTION TRANSPRO

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  
       SUBROUTINE TOCSRFROMS(x, xcsr, tog)
       
       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(INOUT),DIMENSION(0:n,0:m) :: x
       REAL *8,INTENT(INOUT),DIMENSION(NEQ) :: xcsr
       INTEGER,INTENT(IN) :: tog
       INTEGER :: i,j,eqn

       if(tog.EQ.0) then
        do i = 1,n-1
         do j = 2, m-1
          eqn = (j-2)*(n-1) + i
          xcsr(eqn) = x(i,j)
         enddo
        enddo
       endif

       if(tog.EQ.1) then 
        do i = 1,n-1
         do j = 2, m-1
          eqn = (j-2)*(n-1) + i
          x(i,j) = xcsr(eqn)
         enddo
        enddo
       endif
       
       RETURN
         
       END SUBROUTINE TOCSRFROMS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE FORWARDS(y1,x1)

       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(IN),DIMENSION(NEQ) :: y1
       REAL *8,INTENT(OUT),DIMENSION(NEQ) :: x1
       REAL *8 :: dotproduct
       INTEGER :: i,k1,k2

       x1(1) = y1(1)
       do i = 2, NEQ
!        k1 = IAL(i); k2 = IAL(i+1) - 1
!        x1(i) = y1(i) - DOTPRODUCT(AL(k1:k2),x1(JAL(k1:k2)),k2-k1+1)
       enddo
       
       END SUBROUTINE FORWARDS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE BACKWARDS(y1)

       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(INOUT),DIMENSION(NEQ) :: y1
       REAL *8 :: dotproduct
       INTEGER :: i,k1,k2

!       y1(NEQ) = y1(NEQ)/AU(IAU(NEQ))
       do i = NEQ-1,1,-1
!        k1 = IAU(i); k2 = IAU(i+1) - 1
!        y1(i) = ( y1(i) - DOTPRODUCT(AU(k1+1:k2),y1(JAU(k1+1:k2)),k2-k1) )/AU(k1)
       enddo
       
       END SUBROUTINE BACKWARDS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


       SUBROUTINE AMULTIPLY(ty,tx)
       
       USE VARIABLES
       IMPLICIT NONE
       
       REAL *8,INTENT(IN),DIMENSION(-2:n+2,-5:m+6) :: tx
       REAL *8,INTENT(OUT),DIMENSION(-2:n+2,-5:m+6) :: ty
       INTEGER :: i,j, im1, ip1, jm1, jp1,jstart,jend
       REAL *8 :: a1,b1,c1,d1,e1

       jstart = 1
       jend = m-1
       IF(rank.EQ.0) jstart = 2

       DO i = 1,n-1
          im1 = i-1; ip1 = i+1; 
          DO j = jstart, jend
             jm1 = j-1; jp1 = j+1
             a1 = Abig(1,i,j)*tx(i,jm1)
             b1 = Abig(2,i,j)*tx(im1,j)
             c1 = Abig(3,i,j)*tx(i,j)
             d1 = Abig(4,i,j)*tx(ip1,j)
             e1 = Abig(5,i,j)*tx(i,jp1)
             ty(i,j) = a1+b1+c1+d1+e1
          END DO
       END DO

       END SUBROUTINE AMULTIPLY

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE VTE

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,ip1,im1,jstart,jend
       DOUBLE PRECISION :: theq,rhsf,rhs12,rhs3,rhs4,rhs1,rhs2
       DOUBLE PRECISION :: f1d,f2d,u1,v1
       DOUBLE PRECISION :: beta
       DOUBLE PRECISION :: a1,b1,c1
       DOUBLE PRECISION :: term1,term2,term3,term4

       CALL COMPACT

       DO i=2,n

          ip1=i+1;im1=i-1
          IF (i==n) ip1=2

          jstart = 1
          jend = m

          IF (rank .EQ. 0) jstart = 2
          IF (rank.EQ.nproc-1) jend = m-1

        DO j=jstart,jend

          u1 = (psn(i,j+1)-psn(i,j-1))*two_deta_inv
          v1 = -(psn(ip1,j)-psn(im1,j))*two_dzi_inv

          rhs1 = (h22bh11_im1(i,j)*vrto(im1,j)+h22bh11(i,j)*vrto(i,j)+ &
                 h22bh11_ip1(i,j)*vrto(ip1,j))*(dzi_sq_inv)

          rhs2 = (h11bh22_jm1(i,j)*vrto(i,j-1)+h11bh22(i,j)*vrto(i,j)+ &
                 h11bh22_jp1(i,j)*vrto(i,j+1))*(deta_sq_inv)
         
          rhs3 = -u1*omxi(i,j)

          rhs4 = -v1*ometa(i,j)

          rhsf=rhs4+rhs3+(rhs1+rhs2)*re_inv
       
! ***********RK4 CONVERSION STARTS************** !
          theq = dt*rhsf*h11h22_inv(i,j)

          vrt_rk(i,j) = vrt_rk(i,j)+fac(pass)*theq
         
          IF (pass.LT.3) THEN
             vrtn(i,j) = vrt_base(i,j)+f_q(pass)*theq 
          ELSE
             vrtn(i,j) = vrt_base(i,j)+vrt_rk(i,j) 
          END IF      
! ***********RK4 CONVERSION ENDS**************** !
 
          vrtn(1,j) = vrtn(n,j)
        END DO
       END DO

       !!!!!!!!!!!!!!!!!!!!!!
       tem = vrtn
       CALL PASS_tem
       DO i=1,n
          vrtn(i,m)=tem(i,m)
       END DO
       !!!!!!!!!!!!!!!!!!!!!!

       DO j=jstart,jend
          vrtn(n,j)=vrtn(1,j); vrtn(0,j)=vrtn(n-1,j)
       END DO

       DO i=1,n
        DO j=jstart,jend
         vrto(i,j)=vrtn(i,j)
        END DO
       END DO

       CALL VORTICITYPASS
       CALL VORTICITY_OLD_PASS

       RETURN
       END SUBROUTINE VTE     

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE COMPACT

       USE VARIABLES
       IMPLICIT NONE
       
       INTEGER *4 :: i,j,jj,jstart,jend,jstart2,jend2,jstart1,jend1
       INTEGER *4 :: ip1,ip2,im1,im2,temp1,source,status,root
       DOUBLE PRECISION :: s,alfa,u1,v1,xi(-11:m+12),vrtc(-2:n+2,-11:m+12)

       !!!!!!!!!!!!!!!!!!!!!!
       tem = vrtn
       CALL PASS_tem
       DO i=1,n
          vrtn(i,m)=tem(i,m)
       END DO
       !!!!!!!!!!!!!!!!!!!!!!
      
       DO i=1,n
          DO j=1,m
             omxi(i,j)=0.0d0
             ometa(i,j)=0.0d0
          END DO
       END DO
  
       CALL PERIODIC_BOUNDARY 

       jstart1 = 1
       jend1 = m
       IF(rank.EQ.nproc-1) jend1 = m-1
       IF(rank.EQ.0) jstart1 = 2

       !****  XI   DERIVATIVE
       DO j=jstart1,jend1
          vrtn(0,j)=vrtn(n-1,j)
          vrtn(-1,j)=vrtn(n-2,j)
          vrtn(n+1,j)=vrtn(2,j)
          DO i=1,n-1
             f(i)=(vrtn(i-2,j)*am2+vrtn(i-1,j)*am1+vrtn(i+1,j)*ap1+&
                  vrtn(i+2,j)*ap2)*(n-1)
          END DO

          CALL TRIPERIODIC 

          DO i=1,n-1
             omxi(i,j) = f(i)
          END DO
          omxi(n,j) = omxi(1,j) 
       END DO

      !* ADDING DISSIPATION *!

       jstart1 = 1
       jend1 = m
       IF(rank.EQ.nproc-1) jend1 = m-1
       IF(rank.EQ.0) jstart1 = 2

       DO i = 1,n-1
          ip2 = i+2
          ip1 = i+1
          im2 = i-2
          im1 = i-1

          IF(i.EQ.1) THEN
            im2 = n-2
            im1 = n-1
          END IF

          IF(i.EQ.2) THEN
            im2 = n-1
            im1 = 1
          END IF

          IF(i.EQ.n) THEN
            ip2 = 3
            ip1 = 2
          END IF
          IF(i.EQ.n-1) THEN
            ip1 = 1
            ip2 = 2
          END IF

          DO j = jstart1,jend1
             u1=(psn(i,j+1)-psn(i,j-1))*0.5d0*(jmax-1)
             IF (u1 .ge. 0) THEN
                 alfa=.06d0
             ELSE
                 alfa=-.06d0
             END IF
             omxi(i,j)=omxi(i,j)+alfa*(vrtn(ip2,j)-4.d0*vrtn(ip1,j)+ &
                      6.d0*vrtn(i,j)-4.d0*vrtn(im1,j)+vrtn(im2,j))*(n-1)
          END DO
       END DO

       DO j = jstart1,jend1
          omxi(n,j) = omxi(1,j)
       END DO

       jstart1=1
       jend1=m
       IF(rank.EQ.nproc-1) jend1=m-1
       IF(rank.EQ.0) jstart1=2

       DO j=jstart1,jend1
          omxi(n,j) = omxi(1,j)
       END DO

!* ETA  DERIVATIVE *!

       CALL VORTICITYPASS
 
       IF(rank.EQ.0)THEN
          jstart = 1
          jend = m+12
       ELSEIF(rank.EQ.nproc-1)THEN
          jstart = -11
          jend =  m
       ELSE
         jstart = -11
         jend = m+12
       END IF

      DO i=1,n-1
         DO j=jstart+2,jend-2
            d(j)= (vrtn(i,j-2)*am2+vrtn(i,j-1)*am1+vrtn(i,j+1)*ap1+&
                   vrtn(i,j+2)*ap2)*(jmax-1)
         END DO

         d(jstart)=(jmax-1)*(-1.5d0*vrtn(i,jstart)+2.0d0*vrtn(i,jstart+1)- &
                     0.5d0*vrtn(i,jstart+2)) 
         d(jstart+1)=(jmax-1)*((2.0d0*beta1-1)*vrtn(i,jstart)*one_by_three-(8*beta1*one_by_three+0.5d0)*&
                     vrtn(i,jstart+1)+(4.0d0*beta1+1)*vrtn(i,jstart+2)-&
                     (8.0d0*beta1*one_by_three+one_by_six)*vrtn(i,jstart+3)+&
                     2.0d0*beta1*vrtn(i,jstart+4)*one_by_three) 
         d(jend-1)=(jmax-1)*((1.d0-2.d0*beta2)*vrtn(i,jend)*one_by_three+(8*beta2*one_by_three+0.5d0)* &
                    vrtn(i,jend-1)-(4*beta2+1)*vrtn(i,jend-2)+(8*beta2*one_by_three+one_by_six)* &
                    vrtn(i,jend-3)-2*beta2*vrtn(i,jend-4)*one_by_three)       
         d(jend)=(jmax-1)*(1.5d0*vrtn(i,jend)-2.0d0*vrtn(i,jend-1)+0.5d0*vrtn(i,jend-2))

         !* SOLVING NON-PERIODIC TRIDIAGONAL SYSTEM *!
         xi(jstart) = d(jstart)/phi(jstart)
         DO j = jstart+1,jend
            xi(j) = (d(j)-p(j)*xi(j-1))/phi(j)
         END DO
         w(jend) = xi(jend)

         DO j = jend-1,jstart,-1
            w(j) = xi(j)-r(j)/phi(j)*w(j+1)
         END DO
         DO j = jstart,jend
            ometa(i,j) = w(j)
         END DO
      END DO

      DO i = 1,n-1
         DO j = jstart,jend
            jj = (jstart+jend)-j
            vrtc(i,jj) = vrtn(i,j)
         END DO
      END DO

      DO i=1,n-1
         DO j=jstart+2,jend-2
            d(j)= (vrtc(i,j-2)*am2+vrtc(i,j-1)*am1+vrtc(i,j+1)*ap1+&
                   vrtc(i,j+2)*ap2)*(jmax-1)
         END DO

         d(jstart)=(jmax-1)*(-1.5d0*vrtc(i,jstart)+2.0d0*vrtc(i,jstart+1)- &
                    0.5d0*vrtc(i,jstart+2))
         d(jstart+1)=(jmax-1)*((2.0d0*beta1-1)*vrtc(i,jstart)*one_by_three-(8*beta1*one_by_three+0.5d0)* &
                      vrtc(i,jstart+1)+(4.0d0*beta1+1)*vrtc(i,jstart+2)- &
                     (8.0d0*beta1*one_by_three+one_by_six)*vrtc(i,jstart+3)+   &
                      2.0d0*beta1*vrtc(i,jstart+4)*one_by_three)

         d(jend-1)=(jmax-1)*((1.d0-2.d0*beta2)*vrtc(i,jend)*one_by_three+(8*beta2*one_by_three+0.5d0)* &
                    vrtc(i,jend-1)-(4*beta2+1)*vrtc(i,jend-2)+(8*beta2*one_by_three+one_by_six)*  &
                    vrtc(i,jend-3)-2*beta2*vrtc(i,jend-4)*one_by_three)
         d(jend)=(jmax-1)*(1.5d0*vrtc(i,jend)-2.0d0*vrtc(i,jend-1)+0.5d0*vrtc(i,jend-2))
       
         xi(jstart) = d(jstart)/phi(jstart)
         DO j = jstart+1,jend
            xi(j) = (d(j)-p(j)*xi(j-1))/phi(j)
         END DO
         w(jend) = xi(jend)
         DO j = jend-1,jstart,-1
            w(j) = xi(j)-r(j)/phi(j)*w(j+1)
         END DO

         DO j = jstart,jend
            jj = (jstart+jend)-j
            ometa(i,j)=(ometa(i,j)-w(jj))*0.5d0
         END DO
      END DO

      temp = psn
      call msg_pass
      do i=-2,n+2
         psn(i,0) = temp(i,0); psn(i,m+1) = temp(i,m+1)
         psn(i,-1) = temp(i,-1);psn(i,m+2) = temp(i,m+2)
         psn(i,-2) = temp(i,-2);psn(i,m+3) = temp(i,m+3)
         psn(i,-3) = temp(i,-3);psn(i,m+4) = temp(i,m+4)
         psn(i,-4) = temp(i,-4);psn(i,m+5) = temp(i,m+5)
         psn(i,-5) = temp(i,-5);psn(i,m+6) = temp(i,m+6)
      end do

       IF(rank.EQ.0)THEN
          jstart = 1
          jend = m+6
       ELSEIF(rank.EQ.nproc-1)THEN
          jstart = -5
          jend =  m
       ELSE
         jstart = -5
         jend = m+6
       END IF

      im1=i-1 
      DO i=1,n-1
         DO j=jstart+2,jend-2
            IF (i.EQ.1) psn(0,j)=psn(n-1,j)  
            v1=-(psn(i+1,j)-psn(i-1,j))*0.5d0*(n-1)
            IF (v1.GE.0.d0) THEN
               alfa=.06d0
            ELSE
               alfa=-.06d0
            END IF
            ometa(i,j)=ometa(i,j)+alfa*(vrtn(i,j+2)-4.d0*vrtn(i,j+1)  &
                    +6.d0*vrtn(i,j)-4.d0*vrtn(i,j-1)+vrtn(i,j-2))*(jmax-1)
         END DO
      END DO

      DO i=1,n-1
         j=jstart+1
         IF (i.EQ.1) psn(0,j)=psn(n-1,j)
         v1=-(psn(i+1,j)-psn(i-1,j))*0.5d0*(n-1)
         IF (v1.GE.0.d0) THEN
            ometa(i,j)=(vrtn(i,j)-vrtn(i,j-1))*(jmax-1)
         ELSE
            ometa(i,j)=(vrtn(i,j+1)-vrtn(i,j))*(jmax-1)
         END IF

         j=jend-1
         v1= -(psn(i+1,j)-psn(i-1,j))*0.5d0*(n-1)
         IF (v1.GE.0.d0) THEN
            ometa(i,j)=(vrtn(i,j)-vrtn(i,j-1))*(jmax-1)
         ELSE
            ometa(i,j)=(vrtn(i,j+1)-vrtn(i,j))*(jmax-1)
         END IF
      END DO

      DO j=jstart,jend
         ometa(n,j)=ometa(1,j)
         omxi(n,j)=omxi(1,j)
      END DO

!!!!!!!!!!!!!!!!!!!!!!
       tem =  ometa
       CALL PASS_tem
       DO i=1,n
          ometa(i,0)=tem(i,0);ometa(i,m)=tem(i,m)
       END DO
!!!!!!!!!!!!!!!!!!!!!!

      RETURN
         
      END SUBROUTINE COMPACT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE TRIPERIODIC

       USE VARIABLES 
       IMPLICIT NONE

       INTEGER *4 :: j
       REAL *8 :: sum

       v1c(1) = f(1)/q1c(1)
       DO j = 2,nm1-1
          v1c(j) = (f(j)-p1c(j)*v1c(j-1))/q1c(j)
       END DO 

       sum = 0.0d0
       DO j =1,nm1-2
          sum=sum+r1c(j)*v1c(j)
       END DO

       v1c(nm1) = (f(nm1)-sum-p1c(nm1)*v1c(nm1-1))/q1c(nm1)
        
       f(nm1) = v1c(nm1)
       f(nm1-1) = v1c(nm1-1)-u1c(nm1-1)*f(nm1)
       DO j =2,nm1-1
          f(nm1-j)=v1c(nm1-j)-u1c(nm1-j)*f(nm1+1-j)-w1c(nm1-j)*f(nm1)
       END DO

       RETURN

       END SUBROUTINE TRIPERIODIC

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE OBC

       USE VARIABLES
       IMPLICIT NONE

       INTEGER ::i, im1, ip1,j
       REAL *8 ::rhs1,rhs2,asq,theta,rj,psp,yy,psn0(n),psn01(n),psn02(n)
       REAL *8 ::as,bs,cs,ds,t1,t2,t3,psn1(n)

!       DO i=1,nin-1
!          psn(i,m) =15.0d0+Yco(i,m)*dcos(theta_mean)+Xco(i,m)*dsin(theta_mean)
!       END DO

!       DO i=nout+1,n
!          psn(i,m) = 15.0d0+Yco(i,m)*dcos(theta_mean)+Xco(i,m)*dsin(theta_mean)
!       END DO

       psn(0,m) = psn(n-1,m)

       !*  SOLVING THE WAVE EQUATION AT THE BOUNDARY *!
       DO i = nin-1,nout+1
          ip1 = i+1; im1 = i-1
          v(i,m-1)=-(psnp(ip1,m-1)-psnp(im1,m-1))/(2.0d0*dzi*h11(i,m-1))
          v(i,m)=-(psnp(ip1,m)-psnp(im1,m))/(2.0d0*dzi*h11(i,m))
       END DO

       DO i = nin-1,nout+1
          v(i,m) = v(i,m)-v(i,m)*dt*(v(i,m)-v(i,m-1))/(h22(i,m)*deta)
       END DO

       psn0(nout+2)=psn(nout+2,m) ; psn0(nout+1)=psn(nout+1,m)
       DO i = nout,nin,-1
          psn0(i) = psn0(i+2) + v(i+1,m)*2.0d0*dzi*h11(i+1,m)
       END DO

       psn1(nin-2)=psn(nin-2,m) ; psn1(nin-1)=psn(nin-1,m)
       DO i = nin,nout
          psn1(i) = psn1(i-2) - v(i-1,m)*2.0d0*dzi*h11(i-1,m)
       END DO

       DO i = nin,nout
          psn(i,m) = 0.5d0*(psn1(i)+psn0(i))
       END DO

       DO i = nin,nout
          rhs1 = h22(i+1,m)/(2.0d0*h11(i+1,m))*(psn(i+2,m)-psn(i,m)) &
                 -h22(i-1,m)/(2.0d0*h11(i-1,m))*(psn(i,m)-psn(i-1,m))

          rhs2 = h11(i,m)/h22(i,m)*(psn(i,m)-psn(i,m-1)) &
                -h11(i,m-1)/(2.0d0*h22(i,m-1))*(psn(i,m)-psn(i,m-2))

          vrtn(i,m) = -(1.0d0/((h11(i,m)*h22(i,m)*(dzi**2)))*rhs1 + &
                        1.0d0/(h11(i,m)*h22(i,m)*(deta**2))*rhs2)
       END DO

       DO i=1,n
          vrto(i,m) = vrtn(i,m)
       END DO

       END SUBROUTINE OBC

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE PRESSURE

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,jp1,jm1,ip1,im1,k,iconverge,tempnum,kcount
       REAL *8 :: ze,zw,zn,zs,dpdeta(m),as,bs,cs
       REAL *8 :: rho0,rho1,w0,w1cg,betacg,alphacg,transpro
       REAL *8 :: term1,term2,term3,prntem,RHM
       REAL *8,DIMENSION(0:n,0:plimit) :: pcg0=0.0d0,pcg1=0.0d0,vcg0=0.0d0,vcg1=0.0d0
       REAL *8,DIMENSION(0:n,0:plimit) :: scg=0.0d0,tcg=0.0d0,y=0.0d0,z=0.0d0
       REAL *8,DIMENSION(0:n,0:plimit) :: kinvt=0.0d0,kinvs=0.0d0
       REAL *8,DIMENSION(0:n,0:plimit) :: res0=0.0d0 
       REAL *8,DIMENSION(NEQP) :: pcg1csr,ycsr,zcsr,tcgcsr,scgcsr
       REAL *8,DIMENSION(NEQP) :: kinvtcsr,kinvscsr
       pcg1csr=0.0d0;ycsr=0.0d0;zcsr=0.0d0;tcgcsr=0.0d0
       scgcsr=0.0d0;kinvtcsr=0.0d0;kinvscsr=0.0d0
        
       !* CALCULATIONS OF VELOCITIES *!
       DO i =1,n-1
          ip1=i+1; im1=i-1
          IF (i.EQ.1) im1=n-1; IF (i.EQ.n-1) ip1=1

          DO j =2,plimit-1
             jp1 =j+1; jm1 =j-1
             u(i,j)=(psng(i,jp1)-psng(i,jm1))/(2.0d0*deta*h22g(i,j))
             v(i,j)=-(psng(ip1,j)-psng(im1,j))/(2.0d0*dzi*h11g(i,j))
          END DO

          u(i,1)=ro*omega
          v(i,1)=-(psng(ip1,1)-psng(im1,1))/(2.0d0*dzi*h11g(i,1))
          u(i,plimit)=(psng(i,plimit)-psng(i,plimit-1))/(deta*h22g(i,plimit))
          v(i,plimit)=-(psng(ip1,plimit)-psng(im1,plimit))/(2.0d0*dzi*h11g(i,plimit))
       END DO

       DO j =1,plimit
          u(0,j)=u(n-1,j);u(n,j)=u(1,j)
          v(0,j)=v(n-1,j);v(n,j)=v(1,j)
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !* CALCULATIONS OF RHS *!
       DO j = 2,plimit-1
          jp1 = j+1; jm1 = j-1
          DO i = 1,n-1
             ip1 = i+1; im1 = i-1
             IF(i.eq.1) im1 = n-1 ; IF(i.EQ.n-1) ip1 = 1
             IF (i==1) THEN
                ze = (h22g(i,j)+h22g(ip1,j))*(v(ip1,j)+v(i,j))*(vrtgl(ip1,j)+vrtgl(i,j))
                zw = (h22g(i,j)+h22g(n-1,j))*(v(n-1,j)+v(i,j))*(vrtgl(n-1,j)+vrtgl(i,j))
                zn = (h11g(i,j)+h11g(i,jp1))*(u(i,jp1)+u(i,j))*(vrtgl(i,jp1)+vrtgl(i,j))
                zs = (h11g(i,j)+h11g(i,jm1))*(u(i,jm1)+u(i,j))*(vrtgl(i,jm1)+vrtgl(i,j))
             ELSE
                ze = (h22g(i,j)+h22g(ip1,j))*(v(ip1,j)+v(i,j))*(vrtgl(ip1,j)+vrtgl(i,j))
                zw = (h22g(i,j)+h22g(im1,j))*(v(im1,j)+v(i,j))*(vrtgl(im1,j)+vrtgl(i,j))
                zn = (h11g(i,j)+h11g(i,jp1))*(u(i,jp1)+u(i,j))*(vrtgl(i,jp1)+vrtgl(i,j))
                zs = (h11g(i,j)+h11g(i,jm1))*(u(i,jm1)+u(i,j))*(vrtgl(i,jm1)+vrtgl(i,j))
             END IF
             rhsp(i,j)=((ze-zw)/(8.0d0*dzi)-(zn-zs)/(8.0d0*deta))
          END DO
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !* CALCULATIONS OF NEUMANN BOUNDARY CONDITIONS *!

       DO j =1,plimit
          DO i =1,n-1
             ip1=i+1; im1=i-1;
             IF(i.EQ.1) im1 = n-1 ; IF(i.EQ.n-1) ip1 = 1
             vhnew(i,j)=-(psng(ip1,j)-psng(im1,j))/(2*dzi)
             vhold(i,j)=-(psnpg(ip1,j)-psnpg(im1,j))/(2*dzi)
          END DO
       END DO

       DO i = 1,n-1
          ip1=i+1; im1=i-1;
          IF(i.EQ.1) im1=n-1 ; IF(i.EQ.n-1) ip1=1

          j=plimit-1
          term1 = 1.0d0/(4.0d0*re*dzi)*(vrtgl(ip1,j)-vrtgl(im1,j) + &
                  vrtgl(ip1,j+1)-vrtgl(im1,j+1))

          term2 = -0.125d0*(h11g(i,j)+h11g(i,j+1))*(u(i,j) + &
                  u(i,j+1))*(vrtgl(i,j)+vrtgl(i,j+1))

          term3 = -(vhnew(i,j)-vhold(i,j)+vhnew(i,j+1)-vhold(i,j+1))/(2.0d0*dt)

          cs = 0.50d0*(h11g(i,j)/h22g(i,j)+h11g(i,j+1)/h22g(i,j+1))

          bvalueop(i)=(term1+term2+term3)*deta/cs

          j=1
          term1= 1.0d0/(4.0d0*re*dzi)*(vrtgl(ip1,j)-vrtgl(im1,j) + &
                 vrtgl(ip1,j+1)-vrtgl(im1,j+1))

          term2= - 0.125d0*(h11g(i,j)+h11g(i,j+1))*(u(i,j) + &
                 u(i,j+1))*(vrtgl(i,j)+vrtgl(i,j+1))

          term3= -(vhnew(i,j)-vhold(i,j)+vhnew(i,j+1)-vhold(i,j+1))/(2.0d0*dt)

          cs = 0.50d0*(h11g(i,j)/h22g(i,j) + h11g(i,j+1)/h22g(i,j+1))

          bvaluesp(i)=(term1+term2+term3)*deta/cs
       END DO

       !* ACTUAL BICGSTAB FOR PRESSURE *!

       rho0=1.0d0; alphacg=1.0d0; w0=1.0d0
       kcount=0

       !* BOUNDARY CONDITIONS EMPLOYED *!
       DO i = 1,n-1
          prn(i,1)=prn(i,2)-bvaluesp(i)
          prn(i,plimit)=prn(i,plimit-1)+bvalueop(i)
          prn(1,plimit)=1.5d0
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL CONVERGEP(iconverge,kcount)
       res0 = respx
       DO WHILE(iconverge.EQ.0)

       rho1=TRANSPRO(rank,n,plimit,res0,respx)

!       rho1 = 0.0d0
!       DO i = 1,n-1
!          DO j = 2,plimit-1
!             rho1=rho1+res0(i,j)*respx(i,j)
!          END DO
!       END DO

       betacg=(rho1/rho0)*(alphacg/w0)

          DO j =2,plimit-1
             DO i =1,n-1
                pcg1(i,j)=respx(i,j)+betacg*(pcg0(i,j)-w0*vcg0(i,j))
             END DO
          END DO

          DO j = 2,plimit-1
             pcg1(0,j) = pcg1(n-1,j); pcg1(n,j) = pcg1(1,j)
          END DO

          CALL AMULTIPLYP(vcg1,pcg1)

          alphacg = rho1/TRANSPRO(rank,n,plimit,res0,vcg1)

          DO j = 2,plimit-1
             DO i = 1,n-1
                scg(i,j) = respx(i,j)-alphacg*vcg1(i,j)
             END DO
          END DO

          DO j = 2,plimit-1
             scg(0,j)=scg(n-1,j); scg(n,j)=scg(1,j)
          END DO

          CALL AMULTIPLYP(tcg,scg)

          w1cg=TRANSPRO(rank,n,plimit,tcg,scg)/TRANSPRO(rank,n,plimit,tcg,tcg)

          DO j = 2,plimit-1
             DO i = 1,n-1
                prn(i,j)=prn(i,j)+alphacg*pcg1(i,j)+w1cg*scg(i,j)
             END DO
          END DO

          !* BOUNDARY CONDITIONS EMPLOYED *!
          DO i = 1,n-1
             prn(i,1) = prn(i,2) - bvaluesp(i)
             prn(i,plimit) = prn(i,plimit-1) + bvalueop(i)
             prn(1,plimit) = 1.5d0
          END DO

          CALL CONVERGEP(iconverge,kcount)

          pcg0 = pcg1; vcg0 = vcg1
          rho0 = rho1; w0 = w1cg
       END DO

       DO i = 1,n
          DO j = 1,plimit
             prn(i,j) = prn(i,j) - 0.5d0*(u(i,j)**2 + v(i,j)**2)
          END DO
       END DO

       CALL CLCD

       RETURN

       END SUBROUTINE PRESSURE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE CONVERGEP(iconverge,kcount)

       USE VARIABLES
       IMPLICIT NONE
        
       INTEGER,INTENT(INOUT) :: iconverge,kcount
       INTEGER :: i,j,ip1,im1,jp1,jm1
       REAL *8 :: a1,b1,c1,d1,e1

       kcount = kcount + 1

       DO j=1,plimit
          prn(n,j) = prn(1,j) 
          prn(0,j) = prn(n-1,j) 
       END DO

       DO j =2,plimit-1
          jp1=j+1; jm1=j-1
          DO i =1,n-1
             ip1 = i+1; im1 = i-1
             a1 = AP(1,i,j)*prn(i,jm1) 
             b1 = AP(2,i,j)*prn(im1,j)
             c1 = AP(3,i,j)*prn(i,j)
             d1 = AP(4,i,j)*prn(ip1,j)
             e1 = AP(5,i,j)*prn(i,jp1)
             respx(i,j) = rhsp(i,j)/AP(6,i,j)-(a1+b1+c1+d1+e1)
          END DO
       END DO

       iconverge = 1
       DO j =2,plimit-1
          DO i =1,n-1
             IF(dabs(respx(i,j)).GT.tolerance2) THEN
               iconverge = 0; 
               print 50, dabs(respx(i,j)),i,j,kcount
50             format('RESPRE= ',E12.5,2X,'I=',I3,2X,'J=',I3,2X,'COUNT=',I7)
               EXIT
             END IF
          END DO
          IF(iconverge.EQ.0) THEN
!            write(*,*) iconverge,respx(i,j),i,j    
           EXIT
          END IF
       END DO
   
       IF (kcount.GT.5000) THEN 
          iconverge = 1
          open (1,file='convp',position='append')
               write(1,*) t
          close(1)
       END IF

       RETURN

       END SUBROUTINE CONVERGEP

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE TOCSRFROM(x,xcsr,tog)
       
       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(INOUT),DIMENSION(0:n,0:plimit) :: x
       REAL *8,INTENT(INOUT),DIMENSION(NEQP) :: xcsr
       INTEGER,INTENT(IN) :: tog
       INTEGER :: i,j,eqn

       IF (tog.EQ.0) THEN
          DO i = 1,n-1
             DO j = 2,plimit-1
                eqn=(j-2)*(n-1)+i
                xcsr(eqn)=x(i,j)
             END DO
          END DO
       END IF

       IF (tog.EQ.1) THEN 
          DO i = 1,n-1
             DO j = 2,plimit-1
                eqn=(j-2)*(n-1)+i
                x(i,j) = xcsr(eqn)
             END DO
          END DO
       END IF
       
       RETURN
         
       END SUBROUTINE TOCSRFROM

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE FORWARD(y1,x1)

       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(IN),DIMENSION(NEQP) :: y1
       REAL *8,INTENT(OUT),DIMENSION(NEQP) :: x1
       REAL *8 :: dotproduct
       INTEGER :: i,k1,k2

       x1(1) = y1(1)
       DO i = 2,NEQP
!          k1=IALP(i); k2=IALP(i+1)-1
!          x1(i)=y1(i)-DOTPRODUCT(ALP(k1:k2),x1(JALP(k1:k2)),k2-k1+1)
       END DO
       
       END SUBROUTINE FORWARD

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 
       SUBROUTINE BACKWARD(y1)

       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(INOUT),DIMENSION(NEQP) :: y1
       REAL *8 :: dotproduct
       INTEGER :: i,k1,k2

!       y1(NEQP) = y1(NEQP)/AUP(IAUP(NEQP))
       do i = NEQP-1,1,-1
!          k1=IAUP(i); k2=IAUP(i+1)-1
!          y1(i)=(y1(i)-DOTPRODUCT(AUP(k1+1:k2),y1(JAUP(k1+1:k2)),k2-k1))/AUP(k1)
       enddo
       
       END SUBROUTINE BACKWARD

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       FUNCTION DOTPRODUCT(xx,yy,dimen)

       IMPLICIT NONE
        
       INTEGER,INTENT(IN) :: dimen
       REAL *8,INTENT(IN),DIMENSION(dimen) :: xx,yy 
       REAL *8 :: dotproduct
       INTEGER :: i

       dotproduct = 0.0d0
       do i = 1,dimen
!        dotproduct = dotproduct + xx(i)*yy(i)
       enddo

       END FUNCTION DOTPRODUCT
       
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE AMULTIPLYP(ty,tx)

       USE VARIABLES
       IMPLICIT NONE
       
       REAL *8,INTENT(IN),DIMENSION(0:n,0:plimit) :: tx
       REAL *8,INTENT(OUT),DIMENSION(0:n,0:plimit) :: ty
       INTEGER :: i,j, im1, ip1, jm1, jp1
       REAL *8 :: a1,b1,c1,d1,e1

       DO i = 1,n-1
          im1 = i-1; ip1 = i+1; 
          DO j = 2, plimit-1
             jm1 = j-1; jp1 = j+1
             a1 = APM(1,i,j)*tx(i,jm1) 
             b1 = APM(2,i,j)*tx(im1,j)
             c1 = APM(3,i,j)*tx(i,j)
             d1 = APM(4,i,j)*tx(ip1,j)
             e1 = APM(5,i,j)*tx(i,jp1)
             ty(i,j) = a1+b1+c1+d1+e1
          END DO
       END DO

       END SUBROUTINE AMULTIPLYP

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      subroutine filter(victim)  !6th order filter
      USE VARIABLES       
      implicit none
      integer :: i,istartf,iendf
      real*8 :: q0,q1,q2,q3,alpha1,a12,a22,a32,a42,a52,a62
      real*8 :: a72,a13,a23,a33,a43,a53,a63,a73
      real*8 :: c1,c2,c3,c4,c5
      real*8 :: victim(-2:n+2)
      REAL *8 :: w1c_f,r1c_f,p1c_f,q1c_f,u1c_f  ! #                P
      REAL *8 :: v1c_f,alpha_f,f_f              ! #                P
      common/filt/w1c_f(n-2),r1c_f(n-2),p1c_f(n),q1c_f(n),u1c_f(n),v1c_f(n),alpha_f
      common/filt1/f_f(n)
   
      q0 = 11.0d0 / 16.0d0 + 5.0d0 * alpha_f / 8.0d0  !6th order filter
      q1 = 15.0d0 / 32.0d0 + 17.0d0 * alpha_f / 16.0d0
      q2 = -3.0d0 / 16.0d0 + 3.0d0 * alpha_f / 8.0d0
      q3 = 1.0d0 / 32.0d0  - alpha_f / 16.0d0

      a12 = 1.0d0 / 64.0d0 + 31.0d0 * alpha_f / 32.0d0
      a22 = 29.0d0 / 32.0d0 + 3.0d0 * alpha_f / 16.0d0
      a32 = 15.0d0 / 64.0d0 + 17.0d0 * alpha_f / 32.0d0
      a42 = - 5.0d0 / 16.0d0 + 5.0d0 * alpha_f / 8.0d0
      a52 = 15.0d0 / 64.0d0 - 15.0d0 * alpha_f / 32.0d0
      a62 = - 3.0d0 / 32.0d0 + 3.0d0 * alpha_f / 16.0d0
      a72 = 1.0d0 / 64.0d0 - alpha_f / 32.0d0

      a13 = -1.0d0 / 64.0d0 + 1.0d0 * alpha_f / 32.0d0
      a23 = 3.0d0 / 32.0d0 + 13.0d0 * alpha_f / 16.0d0
      a33 = 49.0d0 / 64.0d0 + 15.0d0 * alpha_f / 32.0d0
      a43 = 5.0d0 / 16.0d0 + 3.0d0 * alpha_f / 8.0d0
      a53 = -15.0d0 / 64.0d0 + 15.0d0 * alpha_f / 32.0d0
      a63 = 3.0d0 / 32.0d0 - 3.0d0 * alpha_f / 16.0d0
      a73 = -1.0d0 / 64.0d0 + alpha_f / 32.0d0

      c1 = 15.0d0 / 16.0d0
      c2 = 4.0d0 / 16.0d0
      c3 = -6.0d0 / 16.0d0
      c4 = c2
      c5 = -1.0d0 / 16.0d0

      istartf = 1
      iendf = n

      victim(0) = victim(iendf-1)
      victim(-1) = victim(iendf-2)
      victim(-2) = victim(iendf-3)

      victim(iendf) = victim(1)
      victim(iendf+1) = victim(2)
      victim(iendf+2) = victim(3)

      do i = istartf , iendf-1 

      f_f(i) = q3 * 0.5d0 * (victim(i-3) + victim(i+3)) &
      + q2 * 0.5d0 * ( victim(i-2) + victim(i+2) ) &
      + q1 * 0.5d0 * ( victim(i+1) + victim(i-1) ) &
      + q0 * victim(i)

      end do

      call TRIPERIODIC_FIL

      victim(1:n-1) = f_f(1:n-1)

      end subroutine filter

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

       SUBROUTINE TRIPERIODIC_FIL

       USE VARIABLES 
       IMPLICIT NONE

       INTEGER *4 :: j
       REAL *8 :: sum
       REAL *8 :: w1c_f,r1c_f,p1c_f,q1c_f,u1c_f  ! #                P
       REAL *8 :: v1c_f,alpha_f                  ! #                P
       REAL *8 :: f_f
       common /filt/ w1c_f(n-2),r1c_f(n-2),p1c_f(n),q1c_f(n),u1c_f(n),v1c_f(n),alpha_f
       common /filt1/ f_f(n)

       v1c_f(1) = f_f(1)/q1c_f(1)
              do j = 2 , nm1 - 1
        v1c_f(j) = (f_f(j) - p1c_f(j) * v1c_f(j-1)) / q1c_f(j)
       end do 

       sum = 0.0d0
       do j = 1 , nm1 - 2
        sum = sum + r1c_f(j) * v1c_f(j)
       end do

       v1c_f(nm1) = (f_f(nm1) - sum - p1c_f(nm1) * v1c_f(nm1-1)) / q1c_f(nm1)
        
       f_f(nm1) = v1c_f(nm1)
       f_f(nm1-1) = v1c_f(nm1-1) - u1c_f(nm1-1) * f_f(nm1)
       do j = 2 , nm1 - 1
        f_f(nm1-j) = v1c_f(nm1-j) - u1c_f(nm1-j) * f_f(nm1+1-j) - w1c_f(nm1-j) * f_f(nm1)
       end do

       RETURN

       END SUBROUTINE TRIPERIODIC_FIL

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
       SUBROUTINE POISSONCOEFF 

       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,k,j1,j2,jrow,eqn,meqn,uptr((jmax-2)*(n-1))
       INTEGER :: jstart,jend,check,enum,jw,jj,iw((jmax-2)*(n-1))
       DOUBLE PRECISION :: as,bs,cs,sum,t1,es,ds,fs
       DOUBLE PRECISION :: invdzis,invdetas
       DOUBLE PRECISION :: w1c_f,r1c_f,p1c_f,q1c_f,u1c_f
       DOUBLE PRECISION :: v1c_f,alpha_f
       common/filt/w1c_f(n-2),r1c_f(n-2),p1c_f(n),q1c_f(n),u1c_f(n),v1c_f(n),alpha_f

       IF (rank.EQ.0) THEN
          jstart = 2
          jend = m+5
       ELSEIF (rank.EQ.nproc-1) THEN
          jstart = -4
          jend =  m-1
       ELSE
          jstart = -4
          jend = m+5
       END IF

       invdzis = 1.0d0/(dzi*dzi); invdetas = 1.0d0/(deta*deta)

       DO j = jstart, jend
          DO i = 1,n-1
             as = 0.50d0*(h11(i,j-1)/h22(i,j-1) + h11(i,j)/h22(i,j))
             cs = 0.50d0*(h11(i,j+1)/h22(i,j+1) + h11(i,j)/h22(i,j))
             bs = -(as+cs)
        
             IF (i==1) THEN
                ds = 0.50d0*(h22(n-1,j)/h11(n-1,j) + h22(1,j)/h11(1,j))
                es = 0.50d0*(h22(2,j)/h11(2,j) + h22(1,j)/h11(1,j))
             ELSE
                ds = 0.50d0*(h22(i-1,j)/h11(i-1,j) + h22(i,j)/h11(i,j))
                es = 0.50d0*(h22(i+1,j)/h11(i+1,j) + h22(i,j)/h11(i,j))
             END IF

             fs = -(ds+es)
             Abig(1,i,j) = as*invdetas
             Abig(2,i,j) = ds*invdzis
             Abig(3,i,j) = fs*invdzis + bs*invdetas
             Abig(4,i,j) = es*invdzis
             Abig(5,i,j) = cs*invdetas
             Abig(1:5,i,j) = Abig(1:5,i,j)/(h11(i,j)*h22(i,j))
             Abig(6,i,j) = Abig(3,i,j)
             Abig(1:5,i,j) = Abig(1:5,i,j)/Abig(6,i,j)
          END DO
       END DO

       !* FOR PRESSURE POISSON SUBROUTINE *!
       IF (rank.EQ.0) THEN
          DO j = 2,plimit-1
             DO i = 1,n-1
                as=0.50d0*(h11g(i,j-1)/h22g(i,j-1)+h11g(i,j)/h22g(i,j))
                cs=0.50d0*(h11g(i,j+1)/h22g(i,j+1)+h11g(i,j)/h22g(i,j))
                bs=-(as+cs)
                ds=0.50d0*(h22g(i-1,j)/h11g(i-1,j)+h22g(i,j)/h11g(i,j))

                IF (i==1) THEN
                   ds = 0.50d0*(h22g(n-1,j)/h11g(n-1,j) + h22g(i,j)/h11g(i,j))
                   es = 0.50d0*(h22g(i+1,j)/h11g(i+1,j) + h22g(i,j)/h11g(i,j))
                ELSE
                   ds = 0.50d0*(h22g(i-1,j)/h11g(i-1,j) + h22g(i,j)/h11g(i,j))
                   es = 0.50d0*(h22g(i+1,j)/h11g(i+1,j) + h22g(i,j)/h11g(i,j))
                END IF

                fs = -(ds+es)
                AP(1,i,j) = as*invdetas
                AP(2,i,j) = ds*invdzis
                AP(3,i,j) = fs*invdzis+bs*invdetas
                AP(4,i,j) = es*invdzis
                AP(5,i,j) = cs*invdetas
                AP(6,i,j) = AP(3,i,j)

                IF (j.GT.1.AND.j.LT.(plimit)) AP(6,i,j)=AP(3,i,j)
                IF (j.EQ.2) AP(6,i,j)=AP(3,i,j)+AP(1,i,j)
                IF (j.EQ.(plimit-1).AND.i.NE.1) AP(6,i,j)=AP(3,i,j)+AP(5,i,j)

                AP(1:5,i,j)=AP(1:5,i,j)/AP(6,i,j)

                APM(1,i,j)=AP(1,i,j);APM(2,i,j)=AP(2,i,j);APM(3,i,j)=AP(3,i,j)
                APM(4,i,j)=AP(4,i,j);APM(5,i,j)=AP(5,i,j);APM(6,i,j)=AP(6,i,j)

             END DO
          END DO

          DO i =1,n-1
             APM(3,i,2)=APM(3,i,2)+APM(1,i,2)
             APM(1,i,2)=0.0d0
             IF(i.NE.1) APM(3,i,plimit-1)=APM(3,i,plimit-1)+APM(5,i,plimit-1)
             APM(5,i,plimit-1)=0.0d0
          END DO

       END IF

        
       IF (rank.EQ.0) THEN
          jstart = 1
          jend = m+6
       ELSEIF (rank.EQ.nproc-1) THEN
          jstart = -5
          jend =  m
       ELSE
          jstart = -5
          jend = m+6
       END IF
      
       !* FOR COMPACT SCHEME *!
       !* XI DERIVATIVE *!
           
       p1c(1) = 0.0d0
       u1c(nm1) = 0.0d0
       
       DO i=1,n-1
          a(i) = xl; b(i) = 1.0d0; c(i) = xr
       END DO

       q1c(1) = b(1);  u1c(1) = c(1)/b(1);  
       w1c(1) = a(1)/b(1); r1c(1) = c(nm1)

       DO j=2,nm1-2
          p1c(j) = a(j)
          q1c(j) = b(j)-p1c(j)*u1c(j-1)
          u1c(j) = c(j)/q1c(j)
          w1c(j) = - p1c(j)*w1c(j-1)/q1c(j)
          r1c(j) = - r1c(j-1)*u1c(j-1)
       END DO

       p1c(nm1-1) = a(nm1-1)
       q1c(nm1-1) = b(nm1-1)-p1c(nm1-1)*u1c(nm1-2)
       u1c(nm1-1) = (c(nm1-1)-p1c(nm1-1)*w1c(nm1-2))/q1c(nm1-1)
       p1c(nm1) = a(nm1)-r1c(nm1-2)*u1c(nm1-2)
       sum = 0.0d0

       DO j=1,nm1-2
          sum=sum+r1c(j)*w1c(j)
       END DO

       q1c(nm1)=b(nm1)-sum-p1c(nm1)*u1c(nm1-1)

      !* ETA DERIVATIVE *!

       IF(rank.EQ.0) THEN
         jstart = 1
         jend = m+12
       ELSEIF(rank.EQ.nproc-1) THEN
         jstart = -11
         jend =  m
       ELSE
         jstart = -11
         jend = m+12
       END IF
 
       DO j = jstart,jend
          p(j) = xl
          q(j) = 1.0d0
          r(j) = xr
       END DO

       r(jstart) = 0.0d0; p(jstart+1) = 0.0d0; r(jstart+1) = 0.0d0 
       p(jend-1) = 0.0d0; r(jend-1) = 0.0d0; p(jend) = 0.0d0 
       phi(jstart) = q(jstart)

       DO j = jstart+1,jend
          phi(j) = q(j)-p(j)*r(j-1)/phi(j-1)
       END DO

!!!filter...........................................

       alpha_f = 0.49
       p1c_f(1) = 0.0d0
       u1c_f(nm1) = 0.0d0

       DO i=1,n-1
          a(i) = alpha_f; b(i) = 1.0d0; c(i) = alpha_f
       END DO
       q1c_f(1) = b(1);  u1c_f(1) = c(1)/b(1);  
       w1c_f(1) = a(1)/b(1); r1c_f(1) = c(nm1)

       DO j = 2 , nm1 - 2
          p1c_f(j) = a(j)
          q1c_f(j) = b(j)-p1c_f(j)*u1c_f(j-1)
          u1c_f(j) = c(j)/q1c_f(j)
          w1c_f(j) = - p1c_f(j)*w1c_f(j-1)/q1c_f(j)
          r1c_f(j) = - r1c_f(j-1)*u1c_f(j-1)
       END DO

       p1c_f(nm1-1) = a(nm1-1)
       q1c_f(nm1-1) = b(nm1-1) - p1c_f(nm1-1) * u1c_f(nm1-2)
       u1c_f(nm1-1) = (c(nm1-1) - p1c_f(nm1-1) * w1c_f(nm1-2)) / q1c_f(nm1-1)
       p1c_f(nm1) = a(nm1) - r1c_f(nm1-2) * u1c_f(nm1-2)
       sum = 0.0d0

       DO j=1,nm1-2
          sum=sum+r1c_f(j)*w1c_f(j)
       END DO
       q1c_f(nm1)=b(nm1)-sum-p1c_f(nm1)*u1c_f(nm1-1)

       RETURN
       END SUBROUTINE POISSONCOEFF
        
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
       SUBROUTINE CLCD
       ! Calculates loads and moment coefficients

       USE VARIABLES
       IMPLICIT NONE

       INTEGER *4 :: i, ip1
       REAL *8 :: xp(n-1), yp(n-1)
       REAL *8 :: press(n-1), shear(n-1), dx, dy
       REAL *8 :: dnf_p(n-1), dnf_s(n-1), daf_p(n-1), daf_s(n-1)
       REAL *8 :: dnf(n-1), daf(n-1)
       REAL *8 :: nf_p, nf_s, af_p, af_s, nf, af
       REAL *8 :: cm_le, cm_rc

       ! Initializations
       nf_p=0.0d0
       nf_s=0.0d0
       af_p=0.0d0
       af_s=0.0d0
       nf=0.0d0
       af=0.0d0
       cl=0.0d0
       cd=0.0d0
       cm_le=0.0d0
       cm_rc=0.0d0

       do i=1,n
           Xco(i,1) = -Xco(i,1)
       end do

       ! Calculation of coordinates of midpoint of each surface panel
       DO i=1,n-1
         ip1=i+1

         xp(i)=0.5d0*(Xco(i,1)+Xco(ip1,1))
         yp(i)=0.5d0*(Yco(i,1)+Yco(ip1,1))
       END DO

       ! Determination of press and shear at the midpoint of each surface panel

       DO i=1,n-1
         ip1=i+1

         press(i)=0.5d0*(prn(i,1)+prn(ip1,1))
         shear(i)=0.5d0*(-vrtn(i,1)-vrtn(ip1,1))/re
       END DO

       ! Determination of axial and normal forces
       DO i=1,n-1
         ip1=i+1

         dx=Xco(i+1,1)-Xco(i,1)
         dy=Yco(i+1,1)-Yco(i,1)

         dnf_p(i)=-press(i)*dx
         daf_p(i)=press(i)*dy

         dnf_s(i)=shear(i)*dy
         daf_s(i)=shear(i)*dx

         dnf(i)=dnf_p(i)+dnf_s(i)
         daf(i)=daf_p(i)+daf_s(i)
       END DO

       ! Writing axial forces on each panel to file
       OPEN(1,FILE='AXIAL',STATUS='replace')
       DO i=1,n
         WRITE(1,*) i,daf_s(i),daf_p(i)
       END DO
       CLOSE(1)

       DO i=1,n-1
         nf_p=nf_p+dnf_p(i)
         af_p=af_p+daf_p(i)

         nf_s=nf_s+dnf_s(i)
          af_s=af_s+daf_s(i)

         nf=nf+dnf(i)
         af=af+daf(i)
       END DO

       ! Determination of cl, cd, cm_le and cm_rc
       cl=af*dsin(theta_mean)+nf*dcos(theta_mean)
       cd=af*dcos(theta_mean)-nf*dsin(theta_mean)

       DO i=1,n-1
         cm_le=cm_le+dnf(i)*(xp(i))-daf(i)*(yp(i))
         cm_rc=cm_rc+dnf(i)*(xp(i))-daf(i)*(yp(i))
       END DO

       do i=1,n
           Xco(i,1) = -Xco(i,1)
       end do

       ! Writing loads and moment coefficients to files
       OPEN(1,FILE='CP',POSITION='append',STATUS='unknown')
       DO i=1,n
         WRITE(1,*) Xco(i,1),2.0d0*prn(i,1)
       END DO
       CLOSE(1)

       OPEN(1,FILE='CP_TIME',POSITION='append',STATUS='unknown')
       DO i=1,n
         WRITE(1,*) i,time,2.0d0*prn(i,1)
       END DO
       CLOSE(1)

       OPEN(1,FILE='CL',POSITION='append',STATUS='unknown')
       WRITE(1,*) time,2.0d0*cl
       CLOSE(1)

       OPEN(1,FILE='CD',POSITION='append',STATUS='unknown')
       WRITE(1,*) time,2.0d0*cd
       CLOSE(1)
       
       OPEN(1,FILE='CM_LE',POSITION='append',STATUS='unknown')
       WRITE(1,*) time,2.0d0*cm_le
       CLOSE(1)

       OPEN(1,FILE='CM_ROT',POSITION='append',STATUS='unknown')
       WRITE(1,*) time,2.0d0*cm_rc
       CLOSE(1)

       END SUBROUTINE CLCD

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE PERIODIC_BOUNDARY
       USE VARIABLES
       IMPLICIT NONE

       INTEGER :: i,j,jstart,jend

       IF(rank.eq.0)THEN
       jstart = 1
       jend = m+6
       ELSEIF(rank.eq.nproc-1)THEN
       jstart = -5
       jend =  m
       ELSE
       jstart = -5
       jend = m+6
       END IF

       DO j = jstart,jend
         psn(0,j) = psn(n-1,j)
         psn(-1,j) = psn(n-2,j)
         psn(-2,j) = psn(n-3,j)
         psn(n,j) = psn(1,j)
         psn(n+1,j) = psn(2,j)
         psn(n+2,j) = psn(3,j)
       END DO
       
       IF(rank.eq.0)THEN
       jstart = 1
       jend = m+12
       ELSEIF(rank.eq.nproc-1)THEN
       jstart = -11
       jend =  m
       ELSE
       jstart = -11
       jend = m+12
       END IF

       DO j = jstart,jend 
         vrtn(0,j) = vrtn(n-1,j)
         vrtn(-1,j) = vrtn(n-2,j)
         vrtn(-2,j) = vrtn(n-3,j)
         vrtn(n,j) = vrtn(1,j)
         vrtn(n+1,j) = vrtn(2,j)
         vrtn(n+2,j) = vrtn(3,j)

         vrto(0,j) = vrto(n-1,j)
         vrto(-1,j) = vrto(n-2,j)
         vrto(-2,j) = vrto(n-3,j)
         vrto(n,j) = vrto(1,j)
         vrto(n+1,j) = vrto(2,j)
         vrto(n+2,j) = vrto(3,j)
       END DO

       END SUBROUTINE PERIODIC_BOUNDARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine msg_pass
      use VARIABLES
      implicit none
      integer :: i,j
      integer :: status(mpi_status_size),source,dest,tag,count,sendtag,recvtag

!Ist block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!IInd block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)

      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!
!IIIrd block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!IV th block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!Vth block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if

!!!!!!!!!!!
!VIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,-5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!
      end subroutine msg_pass

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine msg_pass1
      use VARIABLES
      implicit none
      integer :: i,j
      integer :: status(mpi_status_size),source,dest,tag,count,sendtag,recvtag

!Ist block
       IF(rank.eq.0) THEN
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       ELSE IF (rank.EQ.nproc-1) THEN
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr) 
       ELSE  
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

         sendtag=0   ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)  
         Call MPI_Recv(temp(1:n,m),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       END IF
!!!!!!!!


      end subroutine msg_pass1

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine msg_pass2
      use VARIABLES
      implicit none
      integer :: i,j
      integer :: status(mpi_status_size),source,dest,tag,count,sendtag,recvtag

!Ist block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(temp(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(temp(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(temp(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(temp(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if

      end subroutine msg_pass2      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine msg_pass3
      use VARIABLES
      implicit none
      integer :: i,j
      integer :: status(mpi_status_size),source,dest,tag,count,sendtag,recvtag

!Ist block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!IInd block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)

      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-1),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-2),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!
!IIIrd block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-2),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-3),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!IV th block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-3),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-4),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!Vth block
      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-4),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-5),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if

!!!!!!!!!!!
!VIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-5),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-6),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!
!VIIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+7),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,8),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-6),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,8),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-7),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+7),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!VIIIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-8),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+8),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-7),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,9),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-7),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,9),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-8),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+8),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!IXth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-9),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+9),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-8),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,10),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-8),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,10),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-9),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+9),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!Xth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-10),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+10),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-9),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,11),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-9),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,11),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-10),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+10),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!XIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-11),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+11),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-10),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,12),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-10),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,12),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-11),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+11),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

!XIIth block

      if(rank.eq.0) then
         sendtag=0 ; recvtag=1
         source=rank+1;dest=rank+1
         Call MPI_Send(tem(1:n,m-12),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+12),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      else if(rank.eq.nproc-1) then
         sendtag=1  ;  recvtag=0
         source=rank-1   ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-11),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,13),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
      else
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,-11),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,13),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         sendtag=0    ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-12),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m+12),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
      end if
!!!!!!!!!!!

      end subroutine msg_pass3

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine msg_pass4
      use VARIABLES
      implicit none
      integer :: i,j
      integer :: status(mpi_status_size),source,dest,tag,count,sendtag,recvtag

!Ist block
       IF(rank.eq.0) THEN
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(tem(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(tem(1:n,m),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       ELSE IF (rank.EQ.nproc-1) THEN
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr) 
       ELSE  
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(tem(1:n,0),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(tem(1:n,1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

         sendtag=0   ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(tem(1:n,m-1),n,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)  
         Call MPI_Recv(tem(1:n,m),n,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       END IF
!!!!!!!!
      end subroutine msg_pass4

!----------------------------------------------------------------------
!----------------------------------------------------------------------
 

      subroutine npfilter  !6th order filter
      USE VARIABLES
      implicit none
      integer :: i,istartf,iendf
      real*8 :: q0,q1,q2,q3,alpha1
      real*8 :: eta1f,q02,q12,q03,q13,q23
      real*8 :: q04,q14,q24,q34

      alpha1 = 0.450d0
      eta1f = 0.001d0

      q0 = 22.0d0/32.0d0 + 10.0d0 * alpha1 / 16.0d0 - 10.0d0*eta1f  !6th order filter with eta
      q1 = 15.0d0/32.0d0 + 17.0d0 * alpha1 / 16.0d0 + 15.0d0*eta1f
      q2 = -6.0d0/32.0d0 + 6.0d0 * alpha1 / 16.0d0 - 6.0d0*eta1f
      q3 = 1.0d0/32.0d0 - alpha1/16.0d0 + eta1f

      q02=0.50d0+alpha1 !2nd ordwer filter
      q12=0.50d0+alpha1

      q03=5.0d0/8.0d0 + 3.0d0*alpha1/4.0d0 !4th order filter
      q13=1.0d0/2.0d0 + alpha1
      q23=-1.0d0/8.0d0 + alpha1/4.0d0

      q04 = 11.0d0 / 16.0d0 + 5.0d0 * alpha1 / 8.0d0  !6th order filter
      q14 = 15.0d0 / 32.0d0 + 17.0d0 * alpha1 / 16.0d0
      q24 = -3.0d0 / 16.0d0 + 3.0d0 * alpha1 / 8.0d0
      q34 = 1.0d0 / 32.0d0  - alpha1 / 16.0d0

      istartf = 5
      iendf = arr_lim-4

      fun_np(1)=victim2(1) ! Unfiltered
      fun_np(arr_lim)=victim2(arr_lim) ! Unfiltered

      !!!!!!!! Composite filter block starts !!!!!!!!!!!!!!!!!!!

      fun_np(2)= q12*0.5d0*(victim2(3)+victim2(1))+q02*victim2(2)
      fun_np(arr_lim-1)= q12*0.5d0*(victim2(arr_lim-2)+victim2(arr_lim))+ &
                      q02*victim2(arr_lim-1)


      fun_np(3)= q23*0.5d0*(victim2(5)+victim2(1))+q13*0.5d0*(victim2(4)+ &
              victim2(2))+q03*victim2(3)

      fun_np(arr_lim-2)= q23*0.5d0*(victim2(arr_lim-4)+victim2(arr_lim))+ &
                      q13*0.5d0*(victim2(arr_lim-3)+victim2(arr_lim-1))+ &
                      q03*victim2(arr_lim-2)


      fun_np(4)= q34*0.5d0*(victim2(7)+victim2(1)) + q24*0.5d0*(victim2(6)+victim2(2))+ &
              q14*0.5d0*(victim2(5)+ victim2(3))+q04*victim2(4)
      fun_np(arr_lim-3)= q34*0.5d0*(victim2(arr_lim)+victim2(arr_lim-6)) + &
                      q24*0.5d0*(victim2(arr_lim-1)+victim2(arr_lim-5))+ &
                      q14*0.5d0*(victim2(arr_lim-2)+victim2(arr_lim-4))+ &
                      q04*victim2(arr_lim-3)


      !!!!!!!! Composite filter block ends !!!!!!!!!!!!!!!!!!!

      do i = istartf,iendf
         fun_np(i) = eta1f*victim2(i+4) + (q3*0.5d0-5.0d0*eta1f)*victim2(i+3)+ &
                (q2*0.5d0+10.0d0*eta1f)*victim2(i+2)+ &
                (q1*0.5d0-10.0d0*eta1f)*victim2(i+1)+(q0+5.0d0*eta1f)*victim2(i)+ &
                (q1*0.5d0-eta1f)*victim2(i-1)+q2*0.5d0*victim2(i-2)+ &
                 q3*0.5d0*victim2(i-3)
      end do

 
      call npTRIPERIODIC_FIL2

      victim2(1:arr_lim) = var_np(1:arr_lim)

      end subroutine npfilter

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
       SUBROUTINE npTRIPERIODIC_FIL2

       USE VARIABLES
       IMPLICIT NONE

       INTEGER *4 :: j
       REAL *8 :: sum
       REAL *8 :: a9(arr_lim),b9(arr_lim),c9(arr_lim)
       REAL *8 :: beta9(arr_lim),gamm(arr_lim),alphat

       a9(1) = 0.0d0
       b9(1) = 1.0d0
       c9(1) = 0.0d0
       a9(arr_lim) = 0.0d0
       b9(arr_lim) = 1.0d0
       c9(arr_lim) = 0.0d0
       alphat = 0.450d0

       do j = 2,arr_lim-1
        a9(j) = alphat
        b9(j) = 1.0d0
        c9(j) = alphat
       end do

       beta9(1) = b9(1)
       do j = 2,arr_lim
           beta9(j) = b9(j) - ((a9(j)*c9(j-1))/beta9(j-1))
       end do

       gamm(1) = fun_np(1)/b9(1)

        do j=2,arr_lim
           gamm(j) = (fun_np(j)-a9(j)*gamm(j-1))/beta9(j)
       end do

       var_np(arr_lim) = gamm(arr_lim)

       do j=arr_lim-1,1,-1
           var_np(j) = gamm(j)-((c9(j)*var_np(j+1))/beta9(j))
       end do

       RETURN

       END SUBROUTINE npTRIPERIODIC_FIL2

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      Subroutine filter2  !6th order filter
      USE VARIABLES
      implicit none
      integer :: i,istartf,iendf,jstart,jend
      real*8 :: q0,q1,q2,q3,alpha1,a12,a22,a32,a42,a52,a62
      real*8 :: a72,a13,a23,a33,a43,a53,a63,a73
      real*8 :: c1,c2,c3,c4,c5
      real*8 :: eta1f,q02,q12,q03,q13,q23
      real*8 :: q04,q14,q24,q34

      IF(rank.EQ.0)THEN
        jstart = 1
        jend = m+12
      ELSEIF(rank.EQ.nproc-1)THEN
        jstart = -11
        jend =  m
      ELSE
        jstart = -11
        jend = m+12
      END IF

      alpha1 = 0.4d0
      eta1f = 0.00

      q0 = 22.0d0/32.0d0 + 10.0d0 * alpha1 / 16.0d0 - 10.0d0*eta1f  !6th order filter with eta
      q1 = 15.0d0/32.0d0 + 17.0d0 * alpha1 / 16.0d0 + 15.0d0*eta1f
      q2 = -6.0d0/32.0d0 + 6.0d0 * alpha1 / 16.0d0 - 6.0d0*eta1f
      q3 = 1.0d0/32.0d0 - alpha1/16.0d0 + eta1f

      q02=0.50d0+alpha1 !2nd ordwer filter
      q12=0.50d0+alpha1

      q03=5.0d0/8.0d0 + 3.0d0*alpha1/4.0d0 !4th order filter
      q13=1.0d0/2.0d0 + alpha1
      q23=-1.0d0/8.0d0 + alpha1/4.0d0

      q04 = 11.0d0 / 16.0d0 + 5.0d0 * alpha1 / 8.0d0  !6th order filter
      q14 = 15.0d0 / 32.0d0 + 17.0d0 * alpha1 / 16.0d0
      q24 = -3.0d0 / 16.0d0 + 3.0d0 * alpha1 / 8.0d0
      q34 = 1.0d0 / 32.0d0  - alpha1 / 16.0d0

      istartf = jstart+4
      iendf = jend-4
      fun(jstart)=victim1(jstart) ! Unfiltered
      fun(jend)=victim1(jend) ! Unfiltered

      !!!!!!!! Composite filter block starts !!!!!!!!!!!!!!!!!!!
      fun(jstart+1)= q12*0.5d0*(victim1(jstart+2)+victim1(jstart))+q02*victim1(jstart+1)
      fun(jend-1)= q12*0.5d0*(victim1(jend-2)+victim1(jend))+q02*victim1(jend-1)


!      fun(jstart+2)= q23*0.5d0*(victim1(jstart+4)+victim1(jstart))+q13*0.5d0*(victim1(jstart+3)+ &
!               victim1(jstart+1))+q03*victim1(jstart+2)

!      fun(jend-2)= q23*0.5d0*(victim1(jend-4)+victim1(jend))+q13*0.5d0*(victim1(jend-3)+ &
!               victim1(jend-1))+q03*victim1(jend-2)


!      fun(jstart+3)= q34*0.5d0*(victim1(jstart+6)+victim1(jstart))+q24*0.5d0*(victim1(jstart+5)+ &
!                     victim1(jstart+1))+q14*0.5d0*(victim1(jstart+4)+ victim1(jstart+2))+q04*victim1(jstart+3)

!      fun(jend-3)= q34*0.5d0*(victim1(jend)+victim1(jend-6)) + q24*0.5d0*(victim1(jend-1)+victim1(jend-5))+ &
!                   q14*0.5d0*(victim1(jend-2)+ victim1(jend-4))+q04*victim1(jend-3)


      !!!!!!!! Composite filter block ends !!!!!!!!!!!!!!!!!!!

      do i =jstart+2,jend-2
         fun(i) = q23*0.5d0*(victim1(i+2)+victim1(i-2))+q13*0.5d0*(victim1(i+1)+ &
               victim1(i-1))+q03*victim1(i)
      end do

!         fun(i) = eta1f*victim1(i+4) + (q3*0.5d0-5.0d0*eta1f)*victim1(i+3)+ &
!                (q2*0.5d0+10.0d0*eta1f)*victim1(i+2)+ &
!                (q1*0.5d0-10.0d0*eta1f)*victim1(i+1)+(q0+5.0d0*eta1f)*victim1(i)+ &
!                (q1*0.5d0-eta1f)*victim1(i-1)+q2*0.5d0*victim1(i-2)+ &
!                 q3*0.5d0*victim1(i-3)
!      end do

      call TRIPERIODIC_FIL2

      victim1(jstart:jend) = var(jstart:jend)
      end subroutine filter2
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
       SUBROUTINE TRIPERIODIC_FIL2

       USE VARIABLES
       IMPLICIT NONE
       INTEGER *4 :: j,jstart,jend
       REAL *8 :: sum
       REAL *8 :: a9(-11:m+12),b9(-11:m+12),c9(-11:m+12),beta9(-11:m+12),gamm(-11:m+12),alphat

       IF(rank.EQ.0)THEN
        jstart = 1
        jend = m+12
      ELSEIF(rank.EQ.nproc-1)THEN
        jstart = -11
        jend =  m
      ELSE
        jstart = -11
        jend = m+12
      END IF

       a9(jstart) = 0.0d0
       b9(jstart) = 1.0d0
       c9(jstart) = 0.0d0
       a9(jend) = 0.0d0
       b9(jend) = 1.0d0
       c9(jend) = 0.0d0
       alphat = 0.4

       do j = jstart+1,jend-1
          a9(j) = alphat
          b9(j) = 1.0d0
          c9(j) = alphat
       end do

       beta9(jstart) = b9(jstart)
       do j = jstart+1,jend
           beta9(j) = b9(j) - ((a9(j)*c9(j-1))/beta9(j-1))
       end do

       gamm(jstart) = fun(jstart)/b9(jstart)

       do j=jstart+1,jend
           gamm(j) = (fun(j)-a9(j)*gamm(j-1))/beta9(j)
       end do

       var(jend) = gamm(jend)

       do j=jend-1,jstart,-1
           var(j) = gamm(j) - ((c9(j)*var(j+1))/beta9(j))
       end do

       RETURN

       END SUBROUTINE TRIPERIODIC_FIL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE VORTICITYPASS
      USE VARIABLES
      IMPLICIT NONE

      !DOUBLE PRECISION, DIMENSION(1:n,1:13) :: right_send_vor,left_send_vor
      !DOUBLE PRECISION, DIMENSION(1:n,1:13) :: right_recv_vor,left_recv_vor 
      DOUBLE PRECISION, DIMENSION(1:13) :: right_send_vor,left_send_vor
      DOUBLE PRECISION, DIMENSION(1:13) :: right_recv_vor,left_recv_vor
      INTEGER :: intype,right_tag=5,left_tag=6,jstart,jend
      INTEGER :: source,dest,tag,count,sendtag,recvtag,i,j
      INTEGER :: status(MPI_STATUS_SIZE)

      !CALL MPI_TYPE_CONTIGUOUS(n*13,MPI_double_precision,intype,ierr)
      CALL MPI_TYPE_CONTIGUOUS(13,MPI_double_precision,intype,ierr)
      CALL MPI_TYPE_COMMIT(intype,ierr)
      
      DO i = 1,n
      DO j = 0,12
          !DO i = 1,n
             IF (rank == 0) THEN
                !right_send_vor(i,1+j)=vrtn(i,m-12+j) 
                right_send_vor(1+j)=vrtn(i,m-12+j)
             ELSE IF(rank == (nproc-1)) THEN
                !left_send_vor(i,1+j)=vrtn(i,1+j)
                left_send_vor(1+j)=vrtn(i,1+j)
             ELSE
                !right_send_vor(i,1+j)=vrtn(i,m-12+j)
                !left_send_vor(i,1+j)=vrtn(i,1+j) 
                right_send_vor(1+j)=vrtn(i,m-12+j)
                left_send_vor(1+j)=vrtn(i,1+j)
             END IF
          !END DO
       END DO
       

       IF (rank==0) THEN
          CALL MPI_SEND(right_send_vor,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv_vor,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
       ELSE IF(rank==(nproc-1)) THEN
          CALL MPI_SEND(left_send_vor,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv_vor,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
       ELSE
          CALL MPI_SEND(right_send_vor,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv_vor,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
          CALL MPI_SEND(left_send_vor,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv_vor,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
       END IF
 
       DO j = 0,12
          !DO i = 1,n
             IF (rank == 0) THEN
                !vrtn(i,m+j) = right_recv_vor(i,1+j)
                vrtn(i,m+j) = right_recv_vor(1+j)
             ELSE IF(rank == (nproc-1)) THEN
                !vrtn(i,-11+j) = left_recv_vor(i,1+j)
                vrtn(i,-11+j) = left_recv_vor(1+j)
             ELSE
                !vrtn(i,m+j) = right_recv_vor(i,1+j)
                !vrtn(i,-11+j) = left_recv_vor(i,1+j) 
                vrtn(i,m+j) = right_recv_vor(1+j)
                vrtn(i,-11+j) = left_recv_vor(1+j)
             END IF
          !END DO
       END DO
      END DO
         
      CALL MPI_TYPE_FREE(intype,ierr)

      END SUBROUTINE VORTICITYPASS
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE VORTICITY_OLD_PASS
      USE VARIABLES
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(1:13) :: right_send_vor,left_send_vor
      DOUBLE PRECISION, DIMENSION(1:13) :: right_recv_vor,left_recv_vor 

      INTEGER :: intype,right_tag=5,left_tag=6,jstart,jend
      INTEGER :: source,dest,tag,count,sendtag,recvtag,i,j
      INTEGER :: status(MPI_STATUS_SIZE)

      CALL MPI_TYPE_CONTIGUOUS(13,MPI_double_precision,intype,ierr)
      CALL MPI_TYPE_COMMIT(intype,ierr)

       DO i = 1,n 
       DO j = 0,12
         ! DO i = 1,n
             IF (rank == 0) THEN
                right_send_vor(1+j)=vrto(i,m-12+j)
             ELSE IF(rank == (nproc-1)) THEN
                left_send_vor(1+j)=vrto(i,1+j)
             ELSE
                right_send_vor(1+j)=vrto(i,m-12+j)
                left_send_vor(1+j)=vrto(i,1+j)
             END IF
          !END DO
       END DO

       IF (rank==0) THEN
          CALL MPI_SEND(right_send_vor,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv_vor,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
       ELSE IF(rank==(nproc-1)) THEN
          CALL MPI_SEND(left_send_vor,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv_vor,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
       ELSE
          CALL MPI_SEND(right_send_vor,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv_vor,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
          CALL MPI_SEND(left_send_vor,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv_vor,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
       END IF
 
       DO j = 0,12
          !DO i = 1,n
             IF (rank == 0) THEN
                vrto(i,m+j) = right_recv_vor(1+j)
             ELSE IF(rank == (nproc-1)) THEN
                vrto(i,-11+j) = left_recv_vor(1+j)
             ELSE
                vrto(i,m+j) = right_recv_vor(1+j)
                vrto(i,-11+j) = left_recv_vor(1+j)
             END IF
          !END DO
       END DO

      END DO
      CALL MPI_TYPE_FREE(intype, ierr)

      END SUBROUTINE VORTICITY_OLD_PASS
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PASS_temp
      USE VARIABLES
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(1:n) :: right_send,left_send
      DOUBLE PRECISION, DIMENSION(1:n) :: right_recv,left_recv 

      INTEGER :: intype,right_tag=5,left_tag=6,jstart,jend
      INTEGER :: source,dest,tag,count,sendtag,recvtag,i,j
      INTEGER :: status(MPI_STATUS_SIZE)

      CALL MPI_TYPE_CONTIGUOUS(n,MPI_double_precision,intype,ierr)
      CALL MPI_TYPE_COMMIT(intype,ierr)

      DO i = 1,n
         IF (rank == 0) THEN
            right_send(i)=temp(i,m-1)
         ELSE IF(rank == (nproc-1)) THEN
            left_send(i)=temp(i,1)
         ELSE
            right_send(i)=temp(i,m-1)
            left_send(i)=temp(i,1)
         END IF
      END DO

      IF (rank==0) THEN
          CALL MPI_SEND(right_send,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
      ELSE IF(rank==(nproc-1)) THEN
          CALL MPI_SEND(left_send,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
      ELSE
          CALL MPI_SEND(right_send,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
          CALL MPI_SEND(left_send,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
      END IF
 
      DO i = 1,n
         IF (rank == 0) THEN
            temp(i,m) = right_recv(i)
         ELSE IF(rank == (nproc-1)) THEN
            temp(i,0) = left_recv(i)
         ELSE
            temp(i,m) = right_recv(i)
            temp(i,0) = left_recv(i)
         END IF
      END DO

      CALL MPI_TYPE_FREE(intype, ierr)

      END SUBROUTINE PASS_temp
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PASS_tem
      USE VARIABLES
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(1:n) :: right_send,left_send
      DOUBLE PRECISION, DIMENSION(1:n) :: right_recv,left_recv 

      INTEGER :: intype,right_tag=5,left_tag=6,jstart,jend
      INTEGER :: source,dest,tag,count,sendtag,recvtag,i,j
      INTEGER :: status(MPI_STATUS_SIZE)

      CALL MPI_TYPE_CONTIGUOUS(n,MPI_double_precision,intype,ierr)
      CALL MPI_TYPE_COMMIT(intype,ierr)

      DO i = 1,n
         IF (rank == 0) THEN
            right_send(i)=tem(i,m-1)
         ELSE IF(rank == (nproc-1)) THEN
            left_send(i)=tem(i,1)
         ELSE
            right_send(i)=tem(i,m-1)
            left_send(i)=tem(i,1)
         END IF
      END DO

      IF (rank==0) THEN
          CALL MPI_SEND(right_send,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
      ELSE IF(rank==(nproc-1)) THEN
          CALL MPI_SEND(left_send,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
      ELSE
          CALL MPI_SEND(right_send,1,intype,rank+1,left_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(right_recv,1,intype,rank+1,right_tag,MPI_COMM_WORLD,status,ierr)
          CALL MPI_SEND(left_send,1,intype,rank-1,right_tag,MPI_COMM_WORLD,ierr)
          CALL MPI_RECV(left_recv,1,intype,rank-1,left_tag,MPI_COMM_WORLD,status,ierr)
      END IF
 
      DO i = 1,n
         IF (rank == 0) THEN
            tem(i,m) = right_recv(i)
         ELSE IF(rank == (nproc-1)) THEN
            tem(i,0) = left_recv(i)
         ELSE
            tem(i,m) = right_recv(i)
            tem(i,0) = left_recv(i)
         END IF
      END DO
      CALL MPI_TYPE_FREE(intype, ierr)

      END SUBROUTINE PASS_tem
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE FFT

      USE VARIABLES
      IMPLICIT NONE

      INTEGER *4 :: i,j,N1,L,LE,LE2,IP,ND2,k,NM11
      REAL *8 :: AMP(MAX_FFT),FREQ(MAX_FFT)
      DOUBLE COMPLEX :: U1,S,T1
      DOUBLE COMPLEX :: sample(MAX_FFT)

      DO i = 1,MAX_FFT
         sample(i) = 0.0d0
      END DO

      DO i=1,n
         sample(i)=vrtn(i,20)
      END DO

      N1=2**M1
      pi=ATAN(1.0d0)*4.0d0

      DO 20 L= 1,M1
         LE = 2**(M1+1-L)
         LE2= LE/2
         U1 = (1.0,0.0)
         S = CMPLX(COS(pi/FLOAT(LE2)),-SIN(pi/FLOAT(LE2)))
         DO 20 j=1,LE2
            DO 10 i=j,N1,LE
               IP = I+LE2
               T1=sample(i)+sample(IP)
               sample(IP) =(sample(i)-sample(IP))*U1
10             sample(i) = T1
20             U1 = U1*S
               ND2 = N1/2
               NM11 = N1-1
               j =1
               DO 50 I=1,NM11
                  IF (I.GE.J) GOTO 30
                     T1 = sample(J)
                     sample(j)=sample(i)
                     sample(i)=T1
30                   K=ND2
40                   IF(K.GE.j) GOTO 50
                     J=J-k
                     K=k/2
                     GOTO 40
50                   J = J+k

      DO i = 1, N1/2
         FREQ(i) = 2*pi*(i-1)/N1
         Amp(i) = SQRT((REAL(sample(i)))**2+(AIMAG(sample(i)**2))**2)
      END DO

      IF (Amp(N1/2).gt.100) THEN
         FILT_FREQ = 2
      ELSE
         FILT_FREQ = 1
      END IF

!      write(*,*) FILT_FREQ, Amp(N1/2), "FILTE_FREQ,AMP"


      RETURN
      END SUBROUTINE FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE FILTER_XY

       USE VARIABLES
       IMPLICIT NONE
       
       INTEGER :: i,j,k,iconverge, ip1, im1
       INTEGER :: kcount
       REAL *8 :: rho0, rho1, w0, w1cg, betacg, alphacg, transpro_fil
       REAL *8,DIMENSION(0:n,0:mlimit)::pcg0=0.0d0,pcg1=0.0d0,vcg0=0.0d0,vcg1=0.0d0
       REAL *8,DIMENSION(0:n,0:mlimit)::scg=0.0d0,tcg=0.0d0,y=0.0d0,z=0.0d0
       REAL *8,DIMENSION(0:n,0:mlimit)::kinvt=0.0d0,kinvs=0.0d0,res01=0.0d0
       REAL *8,DIMENSION(1000) :: pcg1csr,ycsr,zcsr,tcgcsr,scgcsr
       REAL *8,DIMENSION(1000) :: kinvtcsr,kinvscsr
       pcg1csr = 0.0d0; ycsr = 0.0d0; zcsr = 0.0d0; tcgcsr = 0.0d0
       scgcsr = 0.0d0; kinvtcsr = 0.0d0; kinvscsr = 0.0d0


       rho0 = 1.0d0; alphacg  = 1.0d0; w0 = 1.0d0
       kcount = 0;
       do j = 1,mlimit
        vrtnf(n,j) = vrtnf(1,j); vrtnf(0,j) = vrtnf(n-1,j); vrtnf(n+1,j) = vrtnf(2,j)
       enddo
       CALL CONVERGE_FIL(iconverge,kcount)
       res01 = resxF

       do while( iconverge.EQ.0 )
        rho1 = TRANSPRO_FIL(n,mlimit,res01,resxF)
        betacg = (rho1/rho0)*(alphacg/w0)
        do j = 2,mlimit-1
         do i = 1, n-1
          pcg1(i,j) = resxF(i,j) + betacg*(pcg0(i,j) - w0*vcg0(i,j))
         enddo
        enddo

        do j = 2,mlimit-1
         pcg1(0,j) = pcg1(n-1,j); pcg1(n,j) = pcg1(1,j)
        enddo

        CALL AMULTIPLY_FIL(vcg1,pcg1)

        alphacg = rho1/TRANSPRO_FIL(n,mlimit,res01, vcg1)

        do j = 2,mlimit-1
         do i = 1, n-1
          scg(i,j) = resxF(i,j) - alphacg*vcg1(i,j)
         enddo
        enddo
        do j = 2,mlimit-1
         scg(0,j) = scg(n-1,j); scg(n,j) = scg(1,j)
        enddo

        CALL AMULTIPLY_FIL(tcg,scg)

        w1cg = TRANSPRO_FIL(n,mlimit, tcg, scg)/TRANSPRO_FIL(n,mlimit, tcg, tcg)

        do j = 2,mlimit-1
         do i = 1, n-1
          vrtnf(i,j) = vrtnf(i,j) + alphacg*pcg1(i,j) + w1cg*scg(i,j)
         enddo
        enddo

        do j = 1,mlimit
         vrtnf(n,j) = vrtnf(1,j); vrtnf(0,j) = vrtnf(n-1,j); vrtnf(n+1,j) = vrtnf(2,j)
        enddo

        CALL CONVERGE_FIL(iconverge,kcount)
        pcg0 = pcg1; vcg0 = vcg1;
        rho0 = rho1; w0 = w1cg
       enddo

       RETURN

       END SUBROUTINE FILTER_XY

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE CONVERGE_FIL(iconverge,kcount)

       USE VARIABLES
       IMPLICIT NONE

       INTEGER,INTENT(INOUT) :: iconverge,kcount
       INTEGER :: i,j,ip1,im1,jp1,jm1
       REAL *8 :: a1,b1,c1,d1,e1

       kcount = kcount + 1

       do j = 2,mlimit-1
        jp1 = j+1; jm1 = j-1
        do i = 1,n-1
         ip1 = i+1; im1 = i-1
         a1 = alpha_fil*vrtnf(i,jm1)
         b1 = alpha_fil*vrtnf(im1,j)
         c1 = vrtnf(i,j)
         d1 = alpha_fil*vrtnf(ip1,j)
         e1 = alpha_fil*vrtnf(i,jp1)
         resxF(i,j) = vrtn_rhs(i,j)-(a1+b1+c1+d1+e1)
        end do
       end do

       iconverge = 1
       do j = 2,mlimit-1
        do i = 1,n-1
         if(dabs(resxF(i,j)).GT.tolerance) then
          iconverge = 0
!          print 50, dabs(resxF(i,j)),i,j,kcount
 50       format('RES_FIL = ',E12.5,2X,'I =',I3,2X,'J =',I3,2X,'COUNT =',I7)
          exit
         endif
        end do
        if(iconverge.EQ.0) exit
       end do

       RETURN

       END SUBROUTINE CONVERGE_FIL

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
       FUNCTION TRANSPRO_FIL(n2,m2,tres1,tres2)

       IMPLICIT NONE

       INTEGER,INTENT(IN) :: n2,m2
       REAL *8,INTENT(IN),DIMENSION(0:n2,0:m2) :: tres1, tres2
       INTEGER :: i,j
       REAL *8 :: transpro_fil

       transpro_fil = 0.0d0
       do i = 1,n2-1
        do j = 2,m2-1
         transpro_fil = transpro_fil + tres1(i,j)*tres2(i,j)
        enddo
       enddo
       END FUNCTION TRANSPRO_FIL

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

       SUBROUTINE AMULTIPLY_FIL(ty,tx)

       USE VARIABLES
       IMPLICIT NONE

       REAL *8,INTENT(IN),DIMENSION(0:n,0:mlimit) :: tx
       REAL *8,INTENT(OUT),DIMENSION(0:n,0:mlimit) :: ty
       INTEGER :: i,j, im1, ip1, jm1, jp1
       REAL *8 :: a1,b1,c1,d1,e1

       do i = 1,n-1
        im1 = i-1; ip1 = i+1;
        do j = 2, mlimit-1
         jm1 = j-1; jp1 = j+1
         a1 = alpha_fil*tx(i,jm1)
         b1 = alpha_fil*tx(im1,j)
         c1 = tx(i,j)
         d1 = alpha_fil*tx(ip1,j)
         e1 = alpha_fil*tx(i,jp1)
         ty(i,j) = a1 + b1 + c1 + d1 + e1
        enddo
       enddo

       END SUBROUTINE AMULTIPLY_FIL

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

                                                                                           

