      Program jigzaw
       implicit none

       integer          :: i,j,k,l,r,c,Nf=1984,Nl=246140,Np=8
       double precision :: BBT(1:1984,1:1984),d
       character*3        :: sr

       do i=0,Np-1

          write(sr(3:3),'(I1)')mod(i,10)
          write(sr(2:2),'(I1)')(mod(i,100))/10
          write(sr(1:1),'(I1)')i/100

          print*,sr

          open(2, file='data'//sr//'.bin')
          
          do j=1,Nl
             read(2,*)d,k,r,c
             BBT(r,c)=d
          enddo

          close(2)

       enddo

         DO r = 1, Nf

             DO c = 1, Nf

                if (c .gt. r) then
            
                   BBt(r,c) = BBt(c,r)

                endif

             ENDDO

          ENDDO

         ! open(11,file="P_MATRIX_AAT1.dat")
          open(12,file="P_MATRIX_BBT1.dat")

          DO c = 1, Nf

             DO r = 1, Nf

                !write(11,*) AAt(r,c)
                write(12,*) BBt(r,c)

             ENDDO

          ENDDO

          CLOSE(11)
          CLOSE(12)

       
      End Program jigzaw
