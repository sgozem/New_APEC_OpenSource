      subroutine grid(irid,itid,natom,ntar,q,r,chlab,ngp,rgrd)
*
      implicit none
*
      include 'param.h'
*
      integer irid(matom),itid(mtar),natom,ngp,ntar
      double precision q(matom),r(matom,3),rgrd(mgp,3)
      character*2 chlab
*
      integer ia,iatom,ib,ic,igp,ipart,itar,ixyz
      double precision a,a2,b,b2,c,c2,d(3),qtot,rcc(mtar,3)
      double precision x,x2,y,y2,z,z2
*
      write(*,'(a)') 'Computing center of charge of each target on ' //
     &               chlab // ' side...'
      write(*,*)
      write(*,'(a)') 'Target residue       Q     CCX     CCY     CCZ'
*
      do ixyz = 1,3
         d(ixyz) = 0.0d0
      end do
      qtot = 0.0d0
      itar = 1
      do iatom = 1,natom
         if (irid(iatom).eq.itid(itar)) then
            do ixyz = 1,3
               d(ixyz) = d(ixyz)+q(iatom)*r(iatom,ixyz)
            end do
            qtot = qtot+q(iatom)
            if (iatom.eq.natom.or.irid(iatom).ne.irid(iatom+1)) then
               if (abs(qtot).lt.small) then
                  write(*,1000) irid(iatom),qtot
                  write(*,'(a)') 'Error: Charge too small, exiting'
                  stop
               end if
               do ixyz = 1,3
                  rcc(itar,ixyz) = d(ixyz)/qtot
               end do
               write(*,1100) itid(itar),qtot,(rcc(itar,ixyz),ixyz=1,3)
               do ixyz = 1,3
                  d(ixyz) = 0.0d0
               end do
               qtot = 0.0d0
               itar = itar+1
            end if
         end if
      end do
*
      write(*,*)
*
      write(*,'(a)') 'Constructing spatial grid for ' // chlab //
     &               ' side...'
      write(*,*)
      write(*,'(a)') 'Target residue   Initial grid size     ' //
     &               'Final grid size'
*
      igp = 0
      do itar = 1,ntar
         ipart = 0
         do ia = -np,np
            a = dble(ia)*rho/dble(np)
            a2 = a**2
            do ib = -np,np
               b = dble(ib)*rho/dble(np)
               b2 = b**2
               do ic = -np,np
                  c = dble(ic)*rho/dble(np)
                  c2 = c**2
                  if (a2+b2+c2.le.rho2) then
                     x = rcc(itar,1)+a
                     y = rcc(itar,2)+b
                     z = rcc(itar,3)+c
                     do iatom = 1,natom
                        x2 = (x-r(iatom,1))**2
                        y2 = (y-r(iatom,2))**2
                        z2 = (z-r(iatom,3))**2
                        if (x2+y2+z2.lt.sigma2) go to 100
                     end do
                     igp = igp+1
                     ipart = ipart+1
                     rgrd(igp,1) = x
                     rgrd(igp,2) = y
                     rgrd(igp,3) = z
 100                 continue
                  end if
               end do
            end do
         end do
         write(*,1200) itid(itar),(2*np+1)**3,ipart
      end do
      ngp = igp
*
      write(*,*)
*
      write(*,1300) 'Total grid size:',ngp
      write(*,*)
*
 1000 format(i14,f8.3)
 1100 format(i14,4f8.3)
 1200 format(i14,2i20)
 1300 format(a,i8)
*
      end
