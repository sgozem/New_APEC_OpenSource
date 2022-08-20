      Program Update_final
      implicit real*8 (a-h,o-z)
      character line12*12,atom*4,res*3,SOL*4,ACI*4
      dimension iatnum(total),atom(total),coord(total,3),
     &          SOL(3),ACI(7)

      open(1,file='HETATM.pdb',status='old')
      open(2,file='HETATM.gro',status='unknown')
      
CCCCCCCCC Number os atoms in the protein (numato)
      numatom=total
      initres=residuo+1
      initnum=numero+1
CCCCCCCCC
c      write(3,'(i6)')numato
      do i=1,numatom
         read(1,'(12x,A,14x,3(f8.3))')atom(i),coord(i,1),
     &                              coord(i,2),coord(i,3)
      enddo
      SOL(1)='  OW'
      SOL(2)=' HW1'
      SOL(3)=' HW2'

      ACI(1)='HH31'
      ACI(2)=' CH3'
      ACI(3)='HH32'
      ACI(4)='HH33'
      ACI(5)='   C'
      ACI(6)=' OC1'
      ACI(7)=' OC2'

      icont=1
      do i=1,numatom
         if (atom(icont) == ' OT ') then
            do j=1,3
               write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'SOL',SOL(j),
     &              initnum,coord(icont,1)/10.0d0,coord(icont,2)/10.0d0,
     &              coord(icont,3)/10.0d0
               initnum=initnum+1
               icont=icont+1
            enddo 
            initres=initres+1
         endif
         if (atom(icont) == 'NA  ') then
            write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'NA ','  NA',
     &           initnum,coord(icont,1)/10.0d0,coord(icont,2)/10.0d0,
     &           coord(icont,3)/10.0d0
            initnum=initnum+1
            icont=icont+1
            initres=initres+1
         endif
         if (atom(icont) == 'CL  ') then
            write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'CL ','  CL',
     &           initnum,coord(icont,1)/10.0d0,coord(icont,2)/10.0d0,
     &           coord(icont,3)/10.0d0
            initnum=initnum+1
            icont=icont+1
            initres=initres+1
         endif
         if (atom(icont) == ' CT ') then

ccc The first HH31 atoms is the first one on the ACI

            write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'ACI',ACI(1),
     &          initnum,coord(icont+1,1)/10.0d0,coord(icont+1,2)/10.0d0,
     &          coord(icont+1,3)/10.0d0
               initnum=initnum+1
               icont=icont+1
            write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'ACI',ACI(2),
     &          initnum,coord(icont-1,1)/10.0d0,coord(icont-1,2)/10.0d0,
     &          coord(icont-1,3)/10.0d0
               initnum=initnum+1
               icont=icont+1
            do j=3,7
              write(2,'(I5,A,3x,A,I5,3(f8.3))')initres,'ACI',ACI(j),
     &          initnum,coord(icont,1)/10.0d0,coord(icont,2)/10.0d0,
     &          coord(icont,3)/10.0d0
               initnum=initnum+1
               icont=icont+1               
            enddo
            initres=initres+1
         endif
      if (icont.eq.numatom+1) then
         goto 20
      endif
 20   enddo
      end

