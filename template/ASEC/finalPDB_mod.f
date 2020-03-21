      Program Update_final
      implicit real*8 (a-h,o-z)
      character line25*25,line47*47,line49*49
c      dimension indice(i)

      open(1,file='Initial.xyz',status='old')
      open(2,file='Final.xyz',status='old')
      open(3,file='Insert.xyz',status='unknown')
c      open(3,file='template_inp2tk',status='old')
      
CCCCCCCCC Number os atoms in the protein (numato)
      numato=numero
CCCCCCCCC
      read(1,*)
      read(2,*)
      write(3,'(i6)')numato
      do i=1,numato
         read(1,'(A,i4)')line47,indice
         read(2,'(A,i4,A)')line49,ind,line25
         write(3,'(A,i4,A)')line49,indice,line25
      enddo
      end
      


      
