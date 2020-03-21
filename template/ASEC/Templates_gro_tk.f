      Program Gro2tk
      implicit real*8 (a-h,o-z)
      dimension coordgro(numero,3),coordtk(numero,3),
     &          isel(numero),igro2tk(numero),iseltk(numero)
      character line13*13,line22*22
      open(1,file='final_Config.gro',status='old')
      open(2,file='coordinates_tk.xyz',status='old')
      open(3,file='template_tk2gro',status='unknown')
      open(4,file='template_gro2tk',status='unknown')
c      open(5,file='list_gro',status='old')
c      open(6,file='list_tk',status='unknown')

CCCCCCCCC Number of atoms in the protein (numatoms)

      numatoms=numero
      
CCCCCCCCC
      read(1,*)
      read(1,*)
      read(2,*)
      do i=1,numatoms
         read(1,'(A,f6.3,2x,f6.3,2x,f6.3)')line22,coordgro(i,1),
     &   coordgro(i,2),coordgro(i,3)
         read(2,'(A,f6.2,6x,f6.2,6x,f6.2,6x,A31)')line13,coordtk(i,1),
     &   coordtk(i,2),coordtk(i,3)
      enddo
      
      write(3,'(A30)')' Tinker atoms vs Gromacs atoms'           
      do i=1,numatoms
         do j=1,numatoms
          if ((nint(coordtk(i,1)*100).eq.nint(coordgro(j,1)*1000)).and.
     &        (nint(coordtk(i,2)*100).eq.nint(coordgro(j,2)*1000)).and.
     &        (nint(coordtk(i,3)*100).eq.nint(coordgro(j,3)*1000))) then
                write(3,'(i5,3x,i5)')i,j
            endif
         enddo
      enddo
      
      write(4,'(A30)')'Gromacs atoms vs Tinker atoms'
      do i=1,numatoms
         do j=1,numatoms
          if ((nint(coordgro(i,1)*1000).eq.nint(coordtk(j,1)*100)).and.
     &        (nint(coordgro(i,2)*1000).eq.nint(coordtk(j,2)*100)).and.
     &        (nint(coordgro(i,3)*1000).eq.nint(coordtk(j,3)*100))) then
              write(4,'(i5,3x,i5)')i,j
              igro2tk(i)=j
            endif            
         enddo
      enddo

c      read(5,'(i5)')numsel
c      read(5,*)
c      do i=1,numsel
c         read(5,'(i5)')ii
c         isel(i)=ii
c      enddo 
c      k=1
c      kk=0
c      do i=1,numatoms
c         if (isel(k).ne.i) then
c            kk=kk+1
c            iseltk(kk)=igro2tk(i)
c         else
c            k=k+1
c         endif
c      enddo
c
c      do i=1,kk-1
c         do j=i+1,kk
c           if (iseltk(j).lt.iseltk(i)) then
c              mayor=iseltk(i)
c              iseltk(i)=iseltk(j)
c              iseltk(j)=mayor
c           endif
c         enddo
c      enddo
c      do i=1,kk
c         write(6,'(A6,1x,i5)')'ACTIVE',iseltk(i)
c      enddo
      end
 
