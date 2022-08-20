      Program makeforce
      implicit real*8 (a-h,o-z)
      character line52*52,line49*49,line11*11,vdw*57,line13*13,
     &          atomtype*52,line6*6 
      dimension atomtype(2999),charges(numero),indx2(numero+1), 
     &          indx3(3000),iused3(numero),sigma(100),epsil(100),
     &          iused2(numero),ifxxall(tailall),itype(tailall),
     &          chargefxx0(numero)
      open(1,file='coordinates_tk.xyz',status='old')
      open(2,file='atom.dat',status='old')
      open(3,file='vdw.dat',status='old')
      open(4,file='list_tk.dat',status='old')
      open(5,file='new_atom.dat',status='unknown')
      open(6,file='new_vdw.dat',status='unknown')
      open(7,file='charges.dat',status='old')
      open(8,file='new_charges.dat',status='unknown')

CCCCCCCCCCCCCCCC
C     This program gives you the additional charges (new_charges) and atoms types
C     (new_vdw) to be added to the new_melacu51.prm
C
C     icont: number of atoms which are moving along the molecular dynamics
C     numindx2: number of atom types on the melacu51.prm file
C     numcharges: number of charges on the melacu51.prm force field
C     numatoms: number of atoms on the xyz tinker file (includes the link atom)
C     chnorm: for charge normalization
C   
CCCCCCCCCCCCCCCC

      numindx2=atomos
      numcharges=cargas
      numatoms=numero
      numvdw=vander
      chnorm=100.0d0

      ifxx=tailall

      read(1,*)
      do i=1,numatoms
         read(1,'(A,i4)')line49,indx2(i)
      enddo

      do i=1,numindx2
         read(2,'(A,i4,3x,i2,A)')line11,ind2,ind3,line52
         atomtype(ind2)=line52
         indx3(ind2)=ind3
      enddo

      do i=1,numvdw
         read(3,'(A13,i2,14x,f9.6,2x,f9.6)')line13,ind3,sigma(i),
     &       epsil(i)
      enddo
      
      do i=1,numcharges
         read(7,'(A,i4,14x,f7.4)')line11,ind2,charge
         charges(ind2)=charge
      enddo

      if (ifxx.gt.0) then
         open(9,file='chargefxx0',status='old')
         do i=1,ifxx
            read(9,*)line6,iatom,tempcharge,itype(i)
            ifxxall(i)=iatom
            chargefxx0(iatom)=tempcharge
         enddo
      endif
c This loop check if the new atomtype indx2 of atom icont was
c already used before. It is syncronized with the atomtypes
c added by ASEC.f in the xyz file.
c This does not apply to the ASEC points of the tail,
c which are added sequantially brecause each of them has a 
c different charge
      read(4,'(17x,i6)')icont
      newindx2=2101
      do i=1,icont
         k=0
         kk=0
         read(4,'(i6)')ind1
         kkk=0
         if (ifxx.gt.0) then
            do j=1,ifxx
               if (ind1.eq.ifxxall(j)) then
                  kkk=1
               endif
            enddo
         endif

         if (kkk.eq.0) then
            iused3(i)=indx3(indx2(ind1))
            if (i.ge.2) then
               do j=1,i-1
                  if (iused3(i).eq.iused3(j)) then
                     k=1
                  endif
               enddo
            endif
            iused2(i)=indx2(ind1)
            if (i.ge.2) then
               do j=1,i-1
                  if (iused2(i).eq.iused2(j)) then
                     kk=1
                  endif
               enddo
            endif
            if (kk.eq.0) then
              write(5,'(A11,i4,2x,i3,A,4x,i1)')'atom       ',newindx2,
     &              indx3(indx2(ind1))+100,atomtype(indx2(ind1)),0
              write(8,'(A11,i4,14x,f9.6)')'charge     ',
     &              newindx2,charges(indx2(ind1))/chnorm
              newindx2=newindx2+1
            endif
            if (k.eq.0) then
               write(6,'(A12,i3,14x,f9.6,2x,f13.10)')'vdw         ',
     &               indx3(indx2(ind1))+100,sigma(indx3(indx2(ind1))),
     &               epsil(indx3(indx2(ind1)))/chnorm**2
            endif
         else
            write(5,'(A11,i4,2x,i3,A,4x,i1)')'atom       ',newindx2,
     &           indx3(indx2(ind1))+100,atomtype(indx2(ind1)),0
            write(8,'(A11,i4,14x,f9.6)')'charge     ',
     &           newindx2,chargefxx0(ind1)/chnorm
            write(6,'(A12,i3,14x,f9.6,2x,f13.10)')'vdw         ',
     &           indx3(indx2(ind1))+100,sigma(indx3(indx2(ind1))),
     &           epsil(indx3(indx2(ind1)))/chnorm**2
            newindx2=newindx2+1
         endif
      enddo

      ifixnewindx2=4001

      if (ifxx.gt.0) then
         do i=1,ifxx
            if (itype(i).eq.0) then
               write(5,'(A11,i4,2x,i3,A,4x,i1)')'atom       ',
     &              ifixnewindx2,indx3(indx2(ifxxall(i))),
     &                                atomtype(indx2(ifxxall(i))),0
               write(8,'(A11,i4,14x,f9.6)')'charge     ',
     &              ifixnewindx2,chargefxx0(ifxxall(i))
               ifixnewindx2=ifixnewindx2+1
            endif
         enddo
      endif
      end
 
