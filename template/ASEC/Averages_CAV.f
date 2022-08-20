      Program makeforce
      implicit real*8 (a-h,o-z)
      character line*80,line1*35,line91*91,line22*22,
     &          line126*126,line47*47,line26*26 
      dimension bond(3000),angle(3000),proper(3000),ainproper(3000),
     &         alj14(3000),coul14(3000),aljsr(3000),coulsr(3000),
     &         apotent(3000)
      open(1,file='md.log',status='old')
      open(2,file='Values.dat',status='unknown')
      open(3,file='Averages.dat',status='unknown')      
      write(2,'(A)')'     Bond            Angle            Proper Dih     
     &    Inproper       LJ-14          Coulomb-14       LJ (SR)      
     &Coulomb (SR)    Potential'
      write(3,'(A)')'     Bond            Angle            Proper Dih     
     &    Inproper       LJ-14          Coulomb-14       LJ (SR)      
     &Coulomb (SR)    Potential'

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
      do j=1,2300
         do while (k.eq.0)
         read(1,'(a)')line
         k= index(line,'Energies (kJ/mol)')
            if (k.eq.4) then
               read(1,'(a)')line
               read(1,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &              bond(j),angle(j),proper(j),ainproper(j),alj14(j)
               read(1,'(a)')line
c               write(2,'(a)')line
               
               read(1,'(2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6,2x,e13.6)')
     &          coul14(j),aljsr(j),coulsr(j),apotent(j)
            if (j.gt.300) then   
               write(2,'(2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,
     &          2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7)')bond(j)*
     &          0.2390057D0,
     &          angle(j)*0.2390057D0,proper(j)*0.2390057D0,ainproper(j)*
     &          0.2390057D0,alj14(j)*0.2390057D0,coul14(j)*0.2390057D0,
     &          aljsr(j)*0.2390057D0,coulsr(j)*0.2390057D0,apotent(j)*
     &          0.2390057D0
               bondsum=bondsum+bond(j)*0.2390057D0
               anglesum=anglesum+angle(j)*0.2390057D0
               propersum=propersum+proper(j)*0.2390057D0
               ainpropersum=ainpropersum+ainproper(j)*0.2390057D0
               alj14sum=alj14sum+alj14(j)*0.2390057D0
               coul14sum=coul14sum+coul14(j)*0.2390057D0
               aljsrsum=aljsrsum+aljsr(j)*0.2390057D0
               coulsrsum=coulsrsum+coulsr(j)*0.2390057D0
               apotentsum=apotentsum+apotent(j)*0.2390057D0
            endif
            endif
         enddo
         k=0
      enddo
      bondav=bondsum/2000.0d0
      angleav=anglesum/2000.0d0
      properav=propersum/2000.0d0
      ainproperav=ainpropersum/2000.0d0
      alj14av=alj14sum/2000.0d0
      coul14av=coul14sum/2000.0d0
      aljsrav=aljsrsum/2000.0d0
      coulsrav=coulsrsum/2000.0d0
      apotentav=apotentsum/2000.0d0      
      write(3,'(2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7,
     & 2x,f14.7,2x,f14.7,2x,f14.7,2x,f14.7)')bondav,angleav,
     &properav,ainproperav,alj14av,coul14av,aljsrav,coulsrav,
     &apotentav
      end
