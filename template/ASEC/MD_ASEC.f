      Program Mean_conf
      implicit real*8 (a-h,o-z)
CCCCCCCCC
c  num: number of configurations in the ASEC
c  numero: is the maximum number of atoms which are moving in the dynamic
CCCCCCCCC

      dimension atomo(confs,numero),amino(confs,numero),
     &       iaminum(confs,numero),coord(confs,numero,3),
     &       idiff(numero),xsuma(numero),ysuma(numero),
     &       zsuma(numero),rms(confs),xmean(numero),
     &       ymean(numero),zmean(numero),square(confs),
     &       igro2tk(numero),idifftk(numero),
     &       itk2gro(numero)


      open(1,file='Selected_100.gro',status='old')
      open(2,file='list_gro.dat',status='unknown')
      open(3,file='ASEC_tk.xyz',status='unknown')
      open(4,file='Charges_gro.xyz',status='unknown')
      open(5,file='Best_Dynamic.xyz',status='unknown')
      open(6,file='rms',status='unknown')
      open(7,file='Best_Config.gro',status='unknown')      
      open(8,file='Dynamic.xyz',status='unknown')
      open(9,file='template_gro2tk',status='old')
      open(10,file='template_tk2gro',status='old')
      open(11,file='list_tk.dat',status='unknown')

      num=confs
      numatoms=numero
      iproatoms=protatoms

CCCCCCCCC 
c  Read the information from the dynamic
CCCCCCCCC

      do i=1,num
        if (i.eq.1) then
	  read(1,*)
	  read(1,*)
	else
	  read(1,*)a,b,c
	  read(1,*)
	  read(1,*)  
        endif        
	do ii=1,numatoms 
	      read(1,'(i5,2a5,i5,3f8.3,3f8.4)')iaminum(i,ii),
     &        amino(i,ii),atomo(i,ii),j,coord(i,ii,1),
     &        coord(i,ii,2),coord(i,ii,3)
	enddo
      enddo

CCCCCCCCC
c         Selection of the atoms participating in the dynamic. "icont" is the
c         number of these atoms and "idiff (ii)" is the number of each atom
CCCCCCCCC
      j=0
      icont=0
      do ii=1,numatoms
	do k=1,num-1
	  if ((coord(k,ii,1).ne.coord(k+1,ii,1)).or.
     &        (coord(k,ii,2).ne.coord(k+1,ii,2)).or.
     &        (coord(k,ii,3).ne.coord(k+1,ii,3))) then
             j=j+1
          endif
        enddo
        if (j.ne.0) then
	  write(2,'(i5,2x,2a5,i5)')ii,amino(1,ii),atomo(1,ii),
     &          iaminum(1,ii)
          icont=icont+1
          idiff(icont)=ii
          j=0
        endif
      enddo
      
CCCCCCCC 
c       RMS calculation of each configuration in relation to the 
c       mean configuration and selection of the closest one  
CCCCCCCC

      do l=1,numatoms
         xsuma(l)=0.0d0
         ysuma(l)=0.0d0
         zsuma(l)=0.0d0 
      enddo
      write(4,'(I6)')icont*num
      write(4,*)
      do l=1,icont
        if (idiff(l).le.iproatoms) then
          do m=1,num
             write(4,'(a2,3x,f8.3,2x,f8.3,2x,f8.3)')'Xx',coord(m,
     &             idiff(l),1)*10.0d0,coord(m,idiff(l),2)*10.0d0,
     &             coord(m,idiff(l),3)*10.0d0
             xsuma(l)=xsuma(l)+coord(m,idiff(l),1)*10.0d0
             ysuma(l)=ysuma(l)+coord(m,idiff(l),2)*10.0d0
             zsuma(l)=zsuma(l)+coord(m,idiff(l),3)*10.0d0
          enddo
        xmean(l)=xsuma(l)/num
        ymean(l)=ysuma(l)/num
        zmean(l)=zsuma(l)/num
        endif
      enddo
      do i=1,num
         square(i)=0.0d0
      enddo
      
      do i=1,num
	do l=1,icont
          square(i)=square(i)+
     &    (coord(i,idiff(l),1)*10.0d0-xmean(l))**2.0d0+
     &    (coord(i,idiff(l),2)*10.0d0-ymean(l))**2.0d0+
     &    (coord(i,idiff(l),3)*10.0d0-zmean(l))**2.0d0
	enddo
	rms(i)=dsqrt(square(i)/(icont*3.0d0))
	write(6,'(f9.5)')rms(i)
	if (i.eq.1) then
	   ibest=1
	   bestrms=rms(1)
	else   
	if (rms(i).le.bestrms) then
	   ibest=i
	   bestrms=rms(i)
	endif
	endif
      enddo
      write(6,*)
      write(6,'(a10)')'RMS MIN'
      write(6,'(I5,f9.5)')ibest,bestrms

CCCCCCCC
c     Writing on:
c     ASEC_tk.xyz: grid of charges on tinker format
c     list_tk: these are the moving atoms in tinker format, but excluding the MM atoms
c     Best_Dynamic.xyz: Best xyz configuration, showing the moving atoms only 
c     Best_Config.pdb: Best .pdb configuration
c     Best_Config.gro: Best .gro configuration
c     Dynamic.xyz: Dynamic, showing the moving atoms only
CCCCCCCC
      
      read(9,*)
      do i=1,numatoms
         read(9,'(i5,3x,i5)')j,k
         igro2tk(j)=k
      enddo
      do i=1,icont
         idifftk(i)=igro2tk(idiff(i))
      enddo
cccc
c     Ordering the list of moving atoms in the tk format
cccc
      minimun=1000000
      init=1
      do i=1,icont
         do j=init,icont
            if (idifftk(j).lt.minimun) then
               minimun=idifftk(j)
               k=j
            endif
         enddo
         temp=idifftk(i)
         idifftk(i)=idifftk(k)
         idifftk(k)=temp
         init=init+1
         minimun=1000000
      enddo
      
      write(11,'(A15,2x,i6)')'Number of atoms',icont
      do i=1,icont
         write(11,'(i6)')idifftk(i)
      enddo

      read(10,*)
      do i=1,numatoms
         read(10,'(i5,3x,i5)')j,k
         itk2gro(j)=k
      enddo

      write(3,'(I6)')icont*(num-1)
      write(3,*)
      do m=1,icont
         do l=1,num
           if (l.ne.ibest) then
              write(3,'(a2,2x,f8.3,4x,f8.3,4x,f8.3)')'Xx',
     &              coord(l,itk2gro(idifftk(m)),1)*10.0d0,
     &              coord(l,itk2gro(idifftk(m)),2)*10.0d0,
     &              coord(l,itk2gro(idifftk(m)),3)*10.0d0
           endif
         enddo
      enddo
 
CCCCCCCC

       write(5,'(I4)')icont
       write(5,*)      
      do l=1,icont
	write(5,'(a5,3x,f8.3,2x,f8.3,2x,f8.3)')atomo(ibest,idiff(l)),
     &        coord(ibest,idiff(l),1)*10.0d0,
     &        coord(ibest,idiff(l),2)*10.0d0,
     &        coord(ibest,idiff(l),3)*10.0d0
      enddo
            
CCCCCCCC

      write(7,'(a)')'Generated by Mean_conf'
      write(7,'(I5)')numatoms
      do i=1,numatoms
       write(7,'(i5,2a5,i5,3f8.3,3f8.4)')
     &      iaminum(ibest,i),amino(ibest,i),atomo(ibest,i),i,
     &      coord(ibest,i,1),coord(ibest,i,2),
     &      coord(ibest,i,3)
      enddo
      write(7,'(1x,f9.5,1x,f9.5,1x,f9.5)')a,b,c
      
CCCCCCCC

      do i=1,num
      write(8,'(I4)')icont
      write(8,*)
         do ii=1,icont
           write(8,'(a5,3x,f8.3,2x,f8.3,2x,f8.3)')atomo(i,idiff(ii)),
     &           coord(i,idiff(ii),1)*10.0d0,
     &           coord(i,idiff(ii),2)*10.0d0,
     &           coord(i,idiff(ii),3)*10.0d0            
         enddo
      enddo
      
      end

      
