      subroutine place(natom,ngpec,ngpic,qec,qic,q,r,rgrdec,rgrdic)
*
      implicit none
*
      include 'param.h'
*
      integer natom,ngpec,ngpic,qec,qic
      double precision q(matom),r(matom,3),rgrdec(mgp,3),rgrdic(mgp,3)
*
      integer igplow,iion,iwrk,ixyz,kionec,kionic,kcurr,nionec,nionic
      integer qionec,qionic,qkeep(mion)
      double precision elow,rkeep(mion,3)
      logical zskip
      character*2 chtype
*
* Compute numbers of ions and their charges
*
      nionic = abs(qic)
      nionec = abs(qec)
      if (nionic.gt.0) then
         qionic = -qic/nionic
      end if
      if (nionec.gt.0) then
         qionec = -qec/nionec
      end if
*
* Place counterions one by one...
*
      write(*,'(a)') 'Figuring out optimal ion placement...'
      write(*,*)
*
      write(*,'(a)') 'Ion    Type  Energy      RX      RY      RZ'
*
      kionic = 0
      kionec = 0
*
 100  if (kionic.lt.nionic.or.kionec.lt.nionec) then
*
* ...on intracellular side
*
         zskip = nionec.gt.nionic.and.kionec.eq.0
         if (kionic.lt.nionic.and.(.not.zskip)) then
            kionic = kionic+1
            kcurr = kionic+kionec
            call opt(kcurr,natom,ngpic,qionic,qkeep,q,r,rgrdic,rkeep,
     &               igplow,elow)
            if (qionic.eq.1) then
               chtype = 'NA'
            else
               chtype = 'CL'
            end if
            write(*,1000) 'IC',kionic,chtype,elow,(rgrdic(igplow,ixyz),
     &                    ixyz=1,3)
            qkeep(kcurr) = qionic
            do ixyz = 1,3
               rkeep(kcurr,ixyz) = rgrdic(igplow,ixyz)
            end do
         end if
*
* ...and on extracellular side
*
         if (kionec.lt.nionec) then
            kionec = kionec+1
            kcurr = kionic+kionec
            call opt(kcurr,natom,ngpec,qionec,qkeep,q,r,rgrdec,rkeep,
     &               igplow,elow)
            if (qionec.eq.1) then
               chtype = 'NA'
            else
               chtype = 'CL'
            end if
            write(*,1000) 'EC',kionec,chtype,elow,(rgrdec(igplow,ixyz),
     &                    ixyz=1,3)
            qkeep(kcurr) = qionec
            do ixyz = 1,3
               rkeep(kcurr,ixyz) = rgrdec(igplow,ixyz)
            end do
         end if
*
      go to 100
      end if
*
      write(*,*)
*
* Write counterion positions in pdb format
*
      iwrk = 0
*
      do iion = 1,kcurr
         if (qkeep(iion).eq.-1) then
            iwrk = iwrk+1
            write(luout,1100) iwrk,iwrk,(rkeep(iion,ixyz),ixyz=1,3)
         end if
      end do
*
      do iion = 1,kcurr
         if (qkeep(iion).eq.1) then
            iwrk = iwrk+1
            write(luout,1200) iwrk,iwrk,(rkeep(iion,ixyz),ixyz=1,3)
         end if
      end do
*
 1000 format(a2,i3.3,4x,a2,4f8.3)
 1100 format('HETATM',i5,1x,' CL ',1x,'CL ',1x,'A',i4,4x,3f8.3)
 1200 format('HETATM',i5,1x,' NA ',1x,'NA ',1x,'A',i4,4x,3f8.3)
*
      end
