      Program makeforce
      implicit real*8 (a-h,o-z)
      character line12*12,line*80,line30*30,line39*39,intername*7
      dimension bondav(2),angleav(2),properav(2),ainproperav(2),
     &         alj14av(2),coul14av(2),aljsrav(2),coulsrav(2),
     &         apotentav(2),enerterm(terminos),intername(terminos)
      open(1,file='CAV_energies',status='old')
      open(2,file='molcas.out',status='old')
      open(3,file='Energy_CAV_TINKER',status='unknown')
c      open(7,file='molcas.out',status='old')
      
      iconf=configu
      ilines=lineas
      iterms=terminos
      iterms2=terminos2

c      energy=0.0d0

      igrnum=ilines/iterms
      
      do i=1,iterms
         enerterm(i)=0.0d0
      enddo
      do i=1,igrnum
         do j=1,iterms
            read(1,'(A,5x,f12.8)')intername(j),ener
            enerterm(j)=enerterm(j)+ener
         enddo
      enddo
      do i=1,iterms
         enerterm(i)=enerterm(i)/iconf
      enddo      
      energymm=sum(enerterm)

      write(3,'(A)')' '
      write(3,'(A)')' '
      write(3,'(A)')' Average MM energy terms  '
      do i=1,iterms
         write(3,'(A,1x,f12.8,A,3x,f20.8,A)')intername(i),enerterm(i),
     & '  Eh',enerterm(i)*627.5091809d0,' Kcal/mol'
      enddo
      write(3,'(A)')' '
      write(3,'(8x,f12.8,A)')energymm,'  Eh '
      write(3,'(f20.8,A)')energymm*627.5091809d0,'  kcal/mol '
      write(3,'(A)')' '
      write(3,'(A)')' '

c      do i=1,ilines
c         read(1,'(A,f12.8)')line12,ener
c         energy=energy+ener
c      enddo
c      energymm=energy/iconf

c      write(3,'(A)')' Average MM energy from MOLCAS-TINKER '
c      write(3,'(A)')' '
c      write(3,'(8x,f12.8,A)')energymm,'  Eh '
c      write(3,'(f20.8,A)')energymm*627.5091809d0,'  kcal/mol '

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Reading the first MM energies from molcas output

      l=0
      energyrest=0.0d0
      do while (l.eq.0)
         read(2,'(A80)')line
         k=index(line,'MM energy components passed:')
            if (k.eq.3) then
               do i=1,iterms2
                  read(2,'(A,f12.8)')line12,ener
                  energyrest=energyrest+ener
               enddo
               l=1
            endif         
      enddo
      close(2)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Reading the CASPT2 energy

      open(2,file='molcas.out',status='old')
      l=0
      do while (l.eq.0)
         read(2,'(A80)')line
         k= index(line,'  Total CASPT2 energies:')
            if (k.eq.1) then
croot_2               read(2,*)
croot_3               read(2,*)
               read(2,'(A,f15.8)')line39,energyqm
               l=1
            endif
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      write(3,'(A)')' Average MM energy from MOLCAS-TINKER '
c      write(3,'(A)')' '
c      write(3,'(8x,f12.8,A)')energymm,'  Eh '
c      write(3,'(f20.8,A)')energymm*627.5091809d0,'  kcal/mol '
  
      Total=energyqm-energyrest+energymm
      Totalkcal=Total*627.5091809d0
      TotalkJ=Total*4.184d0
      write(3,'(A)')' '
      write(3,'(A)')' " QM ENERGY (QM+ele) " '
      write(3,'(A)')' '
      write(3,'(5x,f15.8,A)')energyqm-energyrest,' Eh '
      write(3,'(f20.8,A)')(energyqm-energyrest)*627.5091809d0,
     &' kcal/mol '
      write(3,'(A)')' '
      write(3,'(A)')' '
      write(3,'(A)')' MOLCAS-TINKER Total Energy '
      write(3,'(A)')' '
      write(3,'(5x,f15.8,A)')Total,' Eh '
      write(3,'(f20.8,A)')Totalkcal,' kcal/mol '
      write(3,'(A)')' '
      write(3,'(A)')' '
c      write(6,'(A)')' ""TOTAL ENERGY OF THE CAVITY"" '
c      write(6,'(A)')' '
c      write(6,'(3x,f20.10,A)')Total,' kcal/mol '
c      write(6,'(3x,f20.10,A)')TotalEh,' Eh '
c      write(6,'(3x,f20.10,A)')TotalkJ,' kJ/mol '
    
 
      end

      
      
