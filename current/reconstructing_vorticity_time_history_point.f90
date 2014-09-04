
      Module Var

       INTEGER,PARAMETER                        :: M=1001,N=401,Nf=1957,Nm=5,loc=4504 !171546
       DOUBLE PRECISION, DIMENSION(1:Nf)        :: w,wd,wr,w_eig,w_mode
       DOUBLE PRECISION, DIMENSION(1:Nf)        :: a,t
       DOUBLE PRECISION                         :: x,y,wm,wt
       CHARACTER*15                             :: filename
       integer                                  :: u        

      End Module 



      Program Reconstruction

        USE Var
        IMPLICIT NONE
        
        INTEGER             :: i,j,k,r,s_len
        DOUBLE PRECISION    :: dummy,dummy1
        CHARACTER*15         :: str,str1
   

        !--------------------------------Finding the Disturbance Vorticity------------------------------!

        OPEN(1,FILE='fileread')
        OPEN(12,FILE='P_MEAN_psi_vor1')

        READ(12,*)
        READ(12,*)

        DO k=1,Nf

       !  READ(1,*)         
         READ(1,*)filename 
       
         OPEN(11,FILE=filename)

         print*,'opened ',filename

         OPEN(13,FILE='90_disturbance_vorticity_0.504_0.dat')

       !  READ(11,*)
       !  READ(11,*)
       !  READ(11,*)
       !  READ(11,*)
       !  READ(11,*)
       !  READ(11,*)

         DO i=1,loc

            READ(11,*)w(k)
           
              
           if (k.eq.1)then
             
           
               READ(12,*)x,y, wm
            

            endif
              
           ENDDO
          
           print*, x,y, wm
 
           wd(k)=w(k)-wm

           s_len=len_trim(filename)
           read(filename(1:5),*)t(k)

           WRITE(13,*)t(k),wd(k)

     ENDDO
     
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)

        
        !------------------------------------------------------------------------------------------------------!
        !----------------Reconstruction of Vorticity at (.50052,6.24E-02)(line no:691)-------------------------!
        wr=0.0d0        
       DO r=1,Nm

           str1(01:01) = ACHAR(48 + int(r/100)-10*int(r/1000))
           str1(02:02) = ACHAR(48 + int(r/10)-10*int(r/100))
           str1(03:03) = ACHAR(48 + int(r/1)-10*int(r/10))

         !  OPEN(11,FILE='P_POD_CONSTANTS_'//str1(1:3))

           OPEN(11,FILE='POD_'//str1(3:3))
           OPEN(12,FILE='POD_MODE_phi_'//str1(1:3))

           OPEN(13,FILE='SLE5_reconstructed_vorticity_0.5_0.dat')
 
          ! READ(11,*)           

           READ(12,*)
           READ(12,*)
        
          DO i=1,loc

            READ(12,*)x,y,wt

          ENDDO
          
          print*,x,y

          DO k=1,Nf

              READ(11,*)t(k),a(k)
              w_eig(k)=a(k)*wt

              print*,a(k),wt,r
              
              wr(k)=wr(k)+w_eig(k)

          END DO


          CLOSE(11)
          CLOSE(12)

        ENDDO

        DO k=1,Nf

           WRITE(13,*)t(k),wr(k)          

        END DO
 
                   
        !--------------------------------------------------------------------!

      End Program Reconstruction
