      Program Autocorrelation
      implicit real*8 (a-h,o-z)
      character line*80,line1*35,line91*91,line22*22,
     &          line126*126,line47*47,line26*26 
      dimension aInterTot(3000000),aInterZero(3000000),
     &          aInter(3000000),c(3000)
c               bond(3000),angle(3000),proper(3000),ainproper(3000),
c     &         alj14(3000),coul14(3000),aljsr(3000),coulsr(3000),
c     &         apotent(3000),InterP(3000)
      open(1,file='md.log',status='old')
      open(2,file='md_zero.log',status='old')
      open(3,file='Correlation.dat',status='unknown')
      open(4,file='values.dat',status='unknown')
CCCCCCCC   
CC    long, number of fs i the MD
CCCCCCCC
      long=1000
      icorr=10
      
      k=0
      kk=10
      kkk=1
c      do j=1,long*10
      do j=1,13000
         do while (k.eq.0)
           read(1,'(a)')line
           k= index(line,'Energies (kJ/mol)')
           if (k.eq.4) then
               read(1,'(a)')line
               read(1,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &            bond,angle,proper,ainproper,alj14
               read(1,'(a)')line
               read(1,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &          coul14,aljsr,coulsr,apotent      
             if ((j.gt.3000).and.(mod(kk,10).eq.0)) then
c             if (j.gt.3000) then
               aInterTot(kkk)=aljsr+coulsr
c               aInterTot(kkk)=apotent
               kkk=kkk+1
c               write(4,*)apotent
             endif
             kk=kk+1
           endif
         enddo
         k=0
      enddo
      k=0
      kkk=1
      do j=1,1300
         do while (k.eq.0)
           read(2,'(a)')line
           k= index(line,'Energies (kJ/mol)')
           if (k.eq.4) then
               read(2,'(a)')line
               read(2,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &            bond,angle,proper,ainproper,alj14
               read(2,'(a)')line
               read(2,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &          coul14,aljsr,coulsr,apotent      
             if (j.gt.300) then
               aInterZero(kkk)=aljsr+coulsr
               kkk=kkk+1
             endif
           endif
         enddo
         k=0
      enddo
      do i=1,long
        aInter(i)=aInterTot(i)-aInterZero(i)
c        aInter(i)=aInterTot(i)
        write(4,*)aInter(i)
      enddo
      
CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC
      sum0=0.0d0
      sum01=0.0d0
      do i=1,long
         sum0=sum0+aInter(i)*aInter(i)
         sum01=sum01+aInter(i)
      enddo
      aver0=sum0/long
      aver01=sum01/long
c      write(4,*)aver0
      
      do it=1,icorr
         sum1=0.0d0
         sum2=0.0d0
         sum3=0.0d0
         do i=1,long-it
            sum1=sum1+aInter(i)*aInter(i+it)
            sum2=sum2+aInter(i)
            sum3=sum3+aInter(i+it)
         enddo
         aver1=sum1/(long-it)
         aver2=sum2/(long-it)
         aver3=sum3/(long-it)
         c(it)=(aver1-aver2*aver3)/(aver0-aver01*aver01)
         write(3,'(2x,i5,4x,f6.3)')it,c(it)
      enddo
      end

