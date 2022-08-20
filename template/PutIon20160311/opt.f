      subroutine opt(kcurr,natom,ngp,qion,qkeep,q,r,rgrd,rkeep,igplow,
     &               elow)
*
      implicit none
*
      include 'param.h'
*
      integer igplow,kcurr,natom,ngp,qion,qkeep(mion)
      double precision elow,q(matom),r(matom,3),rgrd(mgp,3)
      double precision rkeep(mion,3)
*
      integer iatom,igp,iion
      double precision ec,ewrk,kc,xyz,x2,y2,z2
*
* Coulomb's constant chosen for length in units of angstrom and
* energy in units of eV
*
      kc = qe/(4.0d0*pi*eps*ang)
*
* For each grid point, compute ion-protein and ion-ion contributions
* to electrostatic energy. Check for minimum by comparing current
* energy with the previously lowest one.
*
      elow = big
      igplow = 0
      do igp = 1,ngp
         ec = 0.0d0
         do iatom = 1,natom
            x2 = (rgrd(igp,1)-r(iatom,1))**2
            y2 = (rgrd(igp,2)-r(iatom,2))**2
            z2 = (rgrd(igp,3)-r(iatom,3))**2
            xyz = sqrt(x2+y2+z2)
            ewrk = kc*dble(qion)*q(iatom)/xyz
            ec = ec+ewrk
         end do
         do iion = 1,kcurr-1
            x2 = (rgrd(igp,1)-rkeep(iion,1))**2
            y2 = (rgrd(igp,2)-rkeep(iion,2))**2
            z2 = (rgrd(igp,3)-rkeep(iion,3))**2
            xyz = sqrt(x2+y2+z2)
            if (xyz.ge.sigma) then
               ewrk = kc*dble(qion)*qkeep(iion)/xyz
            else
               ewrk = big
            end if
            ec = ec+ewrk
         end do
         if (ec.lt.elow) then
            elow = ec
            igplow = igp
         end if
      end do
*
      end
