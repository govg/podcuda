      Program file_conditioning
      
       implicit none
       
       double precision:: d1,d2,psi,vort
       
       integer         :: i,j,Nf=12
       character*20    :: str
       
       open(1,file='fileread')

       do i=1,Nf
          
          read(1,*)
          read(1,*)str

          open(2,file=str(1:17))
        
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)
          read(2,*)

          do j=1,1001*401

             read(2,*)d1,d2,psi,vort
          
             open(3,file=str(1:5)//'snapshot')
             write(3,*)vort
             

          enddo

          print*,str
          close(2)
          close(3)
          
       enddo

      End Program file_conditioning
