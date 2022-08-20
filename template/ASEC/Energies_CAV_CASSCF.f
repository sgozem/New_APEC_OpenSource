      Program makeforce
      implicit real*8 (a-h,o-z)
      character line*80,line1*35,line91*91,line22*22,
     &          line126*126,line39*39,line26*26,line30*30,line47*47 
      dimension bondav(2),angleav(2),properav(2),ainproperav(2),
     &         alj14av(2),coul14av(2),aljsrav(2),coulsrav(2),
     &         apotentav(2)
C      open(1,file='md.log',status='old')
      open(1,file='md_1.log',status='old')
      open(2,file='md_2.log',status='old')
C      open(4,file='md_3.log',status='old')
C      open(5,file='Values_full.dat',status='unknown')
      open(6,file='Energy_CAV',status='unknown')
      open(7,file='molcas.out',status='old')
c      open(7,file='test',status='unknown')

C      write(5,'(A)')'     Bond            Angle            Proper Dih     
C     &    Inproper       LJ-14          Coulomb-14       LJ (SR)      
C     &Coulomb (SR)    Potential'
      do l=1,2
        bondsum=0.0d0
        anglesum=0.0d0
        propersum=0.0d0
        ainpropersum=0.0d0
        alj14sum=0.0d0
        coul14sum=0.0d0
        aljsrsum=0.0d0
        coulsrsum=0.0d0
        apotentsum=0.0d0     
        k=0
        do j=1,1300
           do while (k.eq.0)
           read(l,'(a)')line
           k= index(line,'Energies (kJ/mol)')
              if (k.eq.4) then
                read(l,'(a)')line
                read(l,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &                bond,angle,proper,ainproper,alj14
                read(l,'(a)')line                 
                read(l,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &            coul14,aljsr,coulsr,apotent
              if (j.gt.300) then 
c                 if (l.eq.1) then
c                 write(5,'(2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,
c     &            2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7)')bond*
c     &            0.2390057D0,
c     &            angle*0.2390057D0,proper*0.2390057D0,ainproper*
c     &            0.2390057D0,alj14*0.2390057D0,coul14*0.2390057D0,
c     &            aljsr*0.2390057D0,coulsr*0.2390057D0,apotent*
c     &            0.2390057D0
c                 endif
                 bondsum=bondsum+bond*0.2390057D0
                 anglesum=anglesum+angle*0.2390057D0
                 propersum=propersum+proper*0.2390057D0
                 ainpropersum=ainpropersum+ainproper*0.2390057D0
                 alj14sum=alj14sum+alj14*0.2390057D0
                 coul14sum=coul14sum+coul14*0.2390057D0
                 aljsrsum=aljsrsum+aljsr*0.2390057D0
                 coulsrsum=coulsrsum+coulsr*0.2390057D0
                 apotentsum=apotentsum+apotent*0.2390057D0
              endif
              endif
           enddo
           k=0
        enddo
        bondav(l)=bondsum/1000.0d0
        angleav(l)=anglesum/1000.0d0
        properav(l)=propersum/1000.0d0
        ainproperav(l)=ainpropersum/1000.0d0
        alj14av(l)=alj14sum/1000.0d0
        coul14av(l)=coul14sum/1000.0d0
        aljsrav(l)=aljsrsum/1000.0d0
        coulsrav(l)=coulsrsum/1000.0d0
        apotentav(l)=apotentsum/1000.0d0
      enddo
c      write(5,'(A)')' Average after termalization ' 
c      write(5,'(2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,
c     &      2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7)')bondav(1),
c     &      angleav(1),properav(1),ainproperav(1),alj14av(1),
c     &      coul14av(1),aljsrav(1),coulsrav(1),apotentav(1)

      aljsrCAV=aljsrav(1)-aljsrav(2)
      alj14CAV=alj14av(1)-alj14av(2)
      coulsrCAV=coulsrav(1)-coulsrav(2)
      coul14CAV=coul14av(1)-coul14av(2)
      bondCAV=bondav(1)-bondav(2)
      angleCAV=angleav(1)-angleav(2)
      properCAV=properav(1)-properav(2)
      ainproperCAV=ainproperav(1)-ainproperav(2)
      TotalCAV=aljsrCAV+alj14CAV+coulsrCAV+coul14CAV+bondCAV+
     &         angleCAV+properCAV+ainproperCAV
      TotalCAVEh=TotalCAV*0.0015936d0
      TotalCAVkJ=TotalCAV*4.184d0
      write(6,'(A)')' MM energy terms of the cavity '
      write(6,'(A)')' '
      write(6,'(8(3x,f20.10),A)')bondCAV,angleCAV,properCAV,
     & ainproperCAV,aljsrCAV,alj14CAV,coulsrCAV,coul14CAV,' kcal/mol '
      write(6,'(A)')' '
      write(6,'(A)')' Total MM energy of the cavity '
      write(6,'(3x,f20.10,A)')TotalCAV,' kcal/mol '
      write(6,'(3x,f20.10,A)')TotalCAVEh,' Eh '
      write(6,'(3x,f20.10,A)')TotalCAVkJ,' kJ/mol '

C Reading the last MM energies from molcas output

c      do
c         read(7,'(A80)',IOSTAT=io)line
c         if (io.gt.0) then
cc            write(*,*) 'Check charges.out. Something was wrong'
c            EXIT
c         else if (io.eq.0) then
c            k= index(line,'MM energy components passed:')
c            if (k.eq.3) then
c               energymm=0.0d0
c               do i=1,5
c                  read(7,'(A,f15.8)')line30,ener
c                  energymm=energymm+ener
c               enddo
c            endif
c         else if (io.lt.0) then
c             goto 1
c         endif
c 1    enddo
c      write(8,*)energymm
c      close(7)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Reading the first MM energies from molcas output

c      l=0
c      energymm=0.0d0
c      do while (l.eq.0)
c         read(7,'(A80)')line
c         k= index(line,'MM energy components passed:')
c            if (k.eq.3) then
c               do i=1,5
c                  read(7,'(A,f15.8)')line30,ener
c                  energymm=energymm+ener
c               enddo
c               l=1
c            endif         
c      enddo
c      close(7)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Reading the first CASSCF energy 

c      open(7,file='molcas.out',status='old')      
      l=0
      do while (l.eq.0)
         read(7,'(A80)')line
         k= index(line,'Final state energy')
            if (k.eq.7) then
               read(7,*)
               read(7,*)
               read(7,'(A,f15.8)')line47,energ
               energyqm=energ
               l=1
            endif
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Total=energyqm*627.5091809d0+TotalCAV
      TotalEh=Total*0.0015936d0
      TotalkJ=Total*4.184d0
      write(6,'(A)')' '
      write(6,'(A)')' MOLCAS-TINKER Total Energy '
      write(6,'(A)')' Chromophore-Protein interaction energy '
      write(6,'(A)')' '
      write(6,'(3x,f20.10,A)')energyqm*627.5091809d0,' kcal/mol '
      write(6,'(A)')' '
      write(6,'(A)')' '
      write(6,'(A)')' ""TOTAL ENERGY OF THE CAVITY"" '
      write(6,'(A)')' '
      write(6,'(3x,f20.10,A)')Total,' kcal/mol '
      write(6,'(3x,f20.10,A)')TotalEh,' Eh '
      write(6,'(3x,f20.10,A)')TotalkJ,' kJ/mol '
    
 
      end

      
      
