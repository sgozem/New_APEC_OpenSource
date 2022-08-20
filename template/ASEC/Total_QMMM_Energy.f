      Program makeforce
      implicit real*8 (a-h,o-z)
      character line*80,line1*35,line91*91,line22*22,
     &          line126*126,line39*39,line26*26,line30*30,line47*47,
     &           seed*1 
      dimension bondav(2),angleav(2),properav(2),ainproperav(2),
     &         alj14av(2),coul14av(2),aljsrav(2),coulsrav(2),
     &         apotentav(2),potential(5),collection(0:250010),
     &         values(250000)

      iseeds=2
      itotal=25000
      ithermal=1000
      open(1,file='Total_QMMM_Energy.dat',status='unknown')
      open(2,file='CASPT2.out',status='old')

      potentialav=0.0d0
      icont=0
      do l=1,iseeds
         write(seed,'(i1)')l
         open(l+2,file='md_rerun_'//seed//'.log',status='old')
         potentialsum=0.0d0
         k=0
         do j=1,itotal+1
            do while (k.eq.0)
            read(l+2,'(a)')line
            k= index(line,'Energies (kJ/mol)')
               if (k.eq.4) then
                 read(l+2,'(a)')line
                 read(l+2,'(a)')line
c                 read(l+2,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,
c     &                e13.6)')bond,angle,proper,ainproper,alj14
                 read(l+2,'(a)')line
                 read(l+2,'(a)')line               
c                 read(l+2,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,
c     &                 e13.6)')coul14,aljsr,coulsr,apotent
                 read(l+2,'(a)')line
                 read(l+2,'(2x,e13.6)')potent
               if (j.gt.ithermal+1) then
                  collection(icont)=potent*0.2390057d0
                  icont=icont+1 
                  potentialsum=potentialsum+potent
               endif
               endif
            enddo
            k=0
         enddo
         potentialav=potentialav+potentialsum/(itotal-ithermal)
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCC
      open(5,file='average',status='unknown')
c      do l=0,199999
c         write(5,*)collection(l)
c      enddo
      iconta2=0
      do i=0,1000,5
         suma=0.0d0
         iconta=0
         do kk=0,250000
            values(kk)=0.0d0
         enddo
         if (i.eq.0) then
            do j=0,199999
               suma=suma+collection(j)
               iconta=iconta+1
               values(iconta)=collection(j)
            enddo
         else
            do j=0,199999,i
               suma=suma+collection(j)
               iconta=iconta+1
               values(iconta)=collection(j)
            enddo
         endif
         average=suma/iconta

c Standard Deviation and mean error
         square=0.0d0
         do k=1,iconta
            square=square+(values(k)-average)**2
         enddo
         std=sqrt(square/iconta)
         error=std/sqrt(iconta*1.0d0)
         write(5,'(i6,3x,f18.10,3x,f18.10)')iconta,average,error
         
      enddo

      potentialav=potentialav/iseeds
      TotalMM=potentialav
      TotalMMEh=TotalMM*0.0003808d0
      TotalMMkcal=TotalMM*0.2390057d0
      write(1,'(A)')' '
      write(1,'(A)')' Total MM energy '
      write(1,'(3x,f20.10,A)')TotalMMkcal,' kcal/mol '
      write(1,'(3x,f20.10,A)')TotalMMEh,' Eh '
      write(1,'(3x,f20.10,A)')TotalMM,' kJ/mol '

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C Reading the CASPT2 energy

      open(2,file='CASPT2.out',status='old')
      l=0
      do while (l.eq.0)
         read(2,'(A80)')line
         k= index(line,'  Total CASPT2 energies:')
            if (k.eq.1) then
croot_2               read(2,*)
croot_3               read(2,*)
croot_4               read(2,*)
croot_5               read(2,*)
croot_6               read(2,*)
croot_7               read(2,*)
croot_8               read(2,*)
               read(2,'(A,f15.8)')line39,energ
               energyqm=energ
               l=1
            endif
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      TotalEh=energyqm+TotalMMEh
      Total=energyqm*627.5091809d0+TotalMMkcal
      TotalkJ=energyqm*2625.498413d0+TotalMM
      write(1,'(A)')' '
      write(1,'(A)')' MOLCAS-TINKER Total Energy '
      write(1,'(A)')' Chromophore-Protein interaction energy '
      write(1,'(A)')' '
      write(1,'(3x,f20.10,A)')energyqm*627.5091809d0,' kcal/mol '
      write(1,'(3x,f20.10,A)')energyqm, ' Eh '
      write(1,'(3x,f20.10,A)')energyqm*2625.498413d0,' kJ/mol '
      write(1,'(A)')' '
      write(1,'(A)')' '
      write(1,'(A)')' ""TOTAL QMMM ENERGY"" '
      write(1,'(A)')' '
      write(1,'(3x,f20.10,A)')Total,' kcal/mol '
      write(1,'(3x,f20.10,A)')TotalEh,' Eh '
      write(1,'(3x,f20.10,A)')TotalkJ,' kJ/mol '
    
 
      end

      
      
