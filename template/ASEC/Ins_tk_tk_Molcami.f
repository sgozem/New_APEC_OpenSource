      Program Update_final
      implicit real*8 (a-h,o-z)
      character coord_last*34,line80*80,line13*13,line34*34,
     &line35*35
      dimension ichr(finall),coord_last(numero)

      open(1,file='last_tk.xyz',status='old')
      open(2,file='new_tk.xyz',status='old')
      open(3,file='Insert-tk.xyz',status='unknown')
      open(4,file='qmmmatoms',status='old')
CCCCCCCCC Number os atoms in the protein (numato)
      numatoms=numero
      ifin=finall

      do i=1,ifin
         read(4,*)ichr(i)
      enddo

      read(1,*)
      do i=1,numatoms
         read(1,'(A,A)')line13,coord_last(i)
      enddo
      
      read(2,*)
      write(3,'(i6)')numatoms
      do i=1,numatoms
         k=0
         do j=1,ifin
            if (i.eq.ichr(j)) then
               k=1
            endif
         enddo
         if (k.eq.0) then
            read(2,'(A)')line80
            write(3,'(A)')line80
         else
            read(2,'(A,A,A)')line13,line34,line35
            write(3,'(A,A,A)')line13,coord_last(i),line35
         endif
      enddo
      end
      


      
