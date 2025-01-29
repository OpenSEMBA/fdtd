
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module PMLbodies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!17/08/15 update!!!!!!!!!!!
!!!Elimino el tratamiento de los campos magneticos de PMLbody para programar un multiPMLbody 
!!!solo teniendo en cuenta los parametros efectivos y sin actualizar los magneticos.
!!!Mangento en el fichero PMLbody_pre170815_noupdateababienH.F90 la version antigua
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module PMLbodies



   use Report

   use fdetypes
   implicit none
   private
   !structures needed by the PMLbody

   TYPE  ::  BerPML__t
      REAL (KIND=RKIND), pointer ::  field,Plus,Minu
      REAL (KIND=RKIND)          ::  gx2,P_be,P_ce,Psi,transversalDelta,del,posi
      integer (kind=4)  ::  minTotal,maxTotal
   END TYPE BerPML__t


   TYPE  ::  berpml_t
      integer (kind=4)   ::   NumNodes,orient
      type (BerPML__t), allocatable, dimension (:) :: nodes
   end type berpml_t

   !!!!!valiables locales
   type (berpml_t), save, target   ::  berpmlE,berpmlH

   integer (kind=4), parameter :: PMLorden = 2
   REAL (KIND=RKIND), parameter :: CoeffReflPML=1e-4

!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
   public AdvancePMLbodyE,AdvancePMLbodyH,InitPMLbodies,DestroyPMLbodies,StorefieldsPMLbodies,calc_pmlbodypar

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitPMLbodies(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,layoutnumber,size,g2,Gm2,ThereArePMLbodies,resume,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00

      type (SGGFDTDINFO), intent(IN)     ::  sgg
      REAL (KIND=RKIND)     , pointer, dimension ( : )   ::   g2,gm2
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   ::  &
      sggMiEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
      sggMiEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
      sggMiEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
      sggMiHx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
      sggMiHy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
      sggMiHz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      REAL (KIND=RKIND)   , intent(in) , target     :: &
      Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE),&
      Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE),&
      Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE),&
      Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE),&
      Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE),&
      Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      REAL (KIND=RKIND) , dimension (:)   , intent(in)   :: Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
                                                         &  Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
                                                         &  Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE), &
                                                            Idxe(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE), &
                                                            Idye(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE), &
                                                            Idze(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)

      REAL (KIND=RKIND)  :: sigma
      integer (kind=4), intent(in) :: layoutnumber,size
      logical, INTENT(IN)  :: resume
      logical, INTENT(OUT)  ::  ThereArePMLbodies
      integer (kind=4)  ::  jmed,j1,conta,k1,i1,orient
      integer (kind=4), dimension(0:sgg%nummedia) ::maxx,minx,maxy,miny,maxz,minz
      character(len=BUFSIZE) :: buff
      character (LEN=BUFSIZE)  ::  whoami
      type (BerPML__t), pointer :: PML_
      logical :: unstable

!
        eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
!         
!!!
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
      unstable=.false.
!
      ThereArePMLbodies=.FALSE.

      !precontaje


      minx=2**20 ; miny=minx; minz=minx; maxx=-minx; maxy=-miny; maxz=-minz;
      conta=0
      Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
         Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
            Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
               jmed=sggMiEx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEx) then
                       if (i1 < minx(jmed)) minx(jmed)=i1
                       if (j1 < miny(jmed)) miny(jmed)=j1
                       if (k1 < minz(jmed)) minz(jmed)=k1
                       if (i1 > maxx(jmed)) maxx(jmed)=i1
                       if (j1 > maxy(jmed)) maxy(jmed)=j1
                       if (k1 > maxz(jmed)) maxz(jmed)=k1
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do
      Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
         Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
            Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
               jmed=sggMiEy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEy) then
                       if (i1 < minx(jmed)) minx(jmed)=i1
                       if (j1 < miny(jmed)) miny(jmed)=j1
                       if (k1 < minz(jmed)) minz(jmed)=k1
                       if (i1 > maxx(jmed)) maxx(jmed)=i1
                       if (j1 > maxy(jmed)) maxy(jmed)=j1
                       if (k1 > maxz(jmed)) maxz(jmed)=k1
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
         Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
            Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
               jmed=sggMiEz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEz) then
                       if (i1 < minx(jmed)) minx(jmed)=i1
                       if (j1 < miny(jmed)) miny(jmed)=j1
                       if (k1 < minz(jmed)) minz(jmed)=k1
                       if (i1 > maxx(jmed)) maxx(jmed)=i1
                       if (j1 > maxy(jmed)) maxy(jmed)=j1
                       if (k1 > maxz(jmed)) maxz(jmed)=k1
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do

      !!!!!!!!!!!!!!!!!!!!!!
      ThereArePMLbodies=(conta /=0)
      if (.not.therearePMLbodies) then
         return
      endif
!!!!!
      berpmlE%NumNodes=conta
      allocate (berpmlE%Nodes(1 : berpmlE%NumNodes))

      conta=0
      Do k1=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
         Do j1=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
            Do i1=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
               jmed=sggMiHx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEx) then
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do
      Do k1=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
         Do j1=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
            Do i1=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
               jmed=sggMiHy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEy) then
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
         Do j1=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
            Do i1=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
               jmed=sggMiHz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                   orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                   if (orient/=iEz) then
                       conta=conta+1
                   endif
               endif
            end do
         end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!
      ThereArePMLbodies=ThereArePMLbodies.and.(conta /=0)
      if (.not.therearePMLbodies) then
         WRITE (buff, *)    'Buggy ERROR: In PMLbodies. fields exist withouth Hfields. '
        CALL WarnErrReport (buff,.TRUE.)
      endif
      berpmlH%NumNodes=conta
      allocate (berpmlH%Nodes(1 : berpmlH%NumNodes))
!!!!
      !!!!!!!!
      conta=0
      Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
         Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
            Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
               jmed=sggMiEx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEx) then
                  conta=conta+1
                  PML_ => berpmlE%Nodes(conta)
                  PML_%field    => Ex(i1,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                  case (iEy)
                     PML_%del               =    1.0_RKIND/IDye(   j1     )
                     PML_%transversalDelta  =    1.0_RKIND/IDYh(   j1     )
                     PML_%Plus   => Hz(i1,j1  ,k1)
                     PML_%Minu   => Hz(i1,j1-1,k1)
                     PML_%gx2=G2(jmed)
                     PML_%minTotal=minY(jmed)
                     PML_%maxTotal=maxY(jmed)
                     PML_%posi=j1
                  case (iEz)
                     PML_%del               =   1.0_RKIND/IDze(   j1     )
                     PML_%transversalDelta    = 1.0_RKIND/IDzh(      k1  )
                     PML_%Plus   => Hy(i1,j1  ,k1)
                     PML_%Minu   => Hy(i1,j1  ,k1-1)
                     PML_%gx2=-G2(jmed)
                     PML_%minTotal=minZ(jmed)
                     PML_%maxTotal=maxZ(jmed)
                     PML_%posi=k1
                  case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do
      !!!!!!!!!!
      Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
         Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
            Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
               jmed=sggMiEy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEy) then
                  conta=conta+1
                  PML_ => berpmlE%Nodes(conta)
                  PML_%field    => Ey(i1  ,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                   case (iEz)
                                  PML_%del = 1.0_RKIND/IDze(      k1  )
                     PML_%transversalDelta = 1.0_RKIND/IDzh(      k1  )
                     PML_%Plus   => Hx(i1  ,j1  ,k1  )
                     PML_%Minu   => Hx(i1  ,j1  ,k1-1)
                     PML_%minTotal=minZ(jmed)
                     PML_%maxTotal=maxZ(jmed)
                     PML_%gx2=G2(jmed)
                     PML_%posi=k1
                   case (iEx)
                                  PML_%del = 1.0_RKIND/IDxe(i1        )
                     PML_%transversalDelta = 1.0_RKIND/IDxh(i1        )
                     PML_%Plus   => Hz(i1  ,j1  ,k1  )
                     PML_%Minu   => Hz(i1-1,j1  ,k1  )
                     PML_%gx2=-G2(jmed)
                     PML_%minTotal=minX(jmed)
                     PML_%maxTotal=maxX(jmed)
                     PML_%posi=i1
                   case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
         Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
            Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
               jmed=sggMiEz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEz) then
                  conta=conta+1
                  PML_ => berpmlE%Nodes(conta)
                  PML_%field  => Ez(i1  ,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                   case (iEx)
                                  PML_%del = 1.0_RKIND/IDxe(i1        )
                     PML_%transversalDelta = 1.0_RKIND/IDxh(i1        )
                     PML_%Plus => Hy(i1  ,j1  ,k1  )
                     PML_%Minu => Hy(i1-1,j1  ,k1)
                     PML_%gx2=G2(jmed)
                     PML_%minTotal=minX(jmed)
                     PML_%maxTotal=maxX(jmed)
                     PML_%posi=i1
                   case (iEy)
                                  PML_%del = 1.0_RKIND/IDye(    j1     )
                     PML_%transversalDelta = 1.0_RKIND/IDyh(    j1     )
                     PML_%Plus => Hx(i1  ,j1  ,k1  )
                     PML_%Minu => Hx(i1  ,j1-1,k1  )
                     PML_%gx2=-G2(jmed)
                     PML_%minTotal=minY(jmed)
                     PML_%maxTotal=maxY(jmed)
                     PML_%posi=j1
                   case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do

!magneticos


      !!!!!!!!
      conta=0
      Do k1=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
         Do j1=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
            Do i1=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
               jmed=sggMiHx(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEx) then
                  conta=conta+1
                  PML_ => berpmlH%Nodes(conta)
                  PML_%field    => Hx(i1,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                  case (iEy)
                                   PML_%del = 1.0_RKIND/IDYe(   j1     )
                     PML_%transversalDelta  = 1.0_RKIND/IDYe(   j1     )
                     PML_%Plus   => Ez(i1,j1+1,k1)
                     PML_%Minu   => Ez(i1,j1  ,k1)
                     PML_%gx2=GM2(jmed)
                     PML_%minTotal=minY(jmed)
                     PML_%maxTotal=maxY(jmed)
                     PML_%posi=j1+0.5
                  case (iEz)
                                     PML_%del = 1.0_RKIND/IDze(      k1  )
                     PML_%transversalDelta    = 1.0_RKIND/IDze(      k1  )
                     PML_%Plus   => Ey(i1,j1  ,k1+1)
                     PML_%Minu   => Ey(i1,j1  ,k1  )
                     PML_%gx2=-GM2(jmed)
                     PML_%minTotal=minZ(jmed)
                     PML_%maxTotal=maxZ(jmed)
                     PML_%posi=k1+0.5
                  case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do
      !!!!!!!!!!
      Do k1=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
         Do j1=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
            Do i1=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
               jmed=sggMiHy(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEy) then
                  conta=conta+1
                  PML_ => berpmlH%Nodes(conta)
                  PML_%field    => Hy(i1  ,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                   case (iEz)
                                  PML_%del = 1.0_RKIND/IDze(      k1  )
                     PML_%transversalDelta = 1.0_RKIND/IDze(      k1  )
                     PML_%Plus   => Ex(i1  ,j1  ,k1+1)
                     PML_%Minu   => Ex(i1  ,j1  ,k1  )
                     PML_%gx2=GM2(jmed)
                     PML_%minTotal=minZ(jmed)
                     PML_%maxTotal=maxZ(jmed)
                     PML_%posi=k1+0.5
                   case (iEx)
                                  PML_%del =  1.0_RKIND/IDxe(i1        )
                     PML_%transversalDelta = 1.0_RKIND/IDxe(i1        )
                     PML_%Plus   => Ez(i1+1,j1  ,k1  )
                     PML_%Minu   => Ez(i1  ,j1  ,k1  )
                     PML_%gx2=-GM2(jmed)
                     PML_%minTotal=minX(jmed)
                     PML_%maxTotal=maxX(jmed)
                     PML_%posi=i1+0.5
                   case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do

      Do k1=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
         Do j1=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
            Do i1=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
               jmed=sggMiHz(i1,j1,k1)
               if (SGG%Med(jmed)%Is%PMLbody) then
                orient=abs(SGG%Med(jmed)%PMLbody(1)%orient)
                if (orient/=iEz) then
                  conta=conta+1
                  PML_ => berpmlH%Nodes(conta)
                  PML_%field  => Hz(i1  ,j1  ,k1  )
                  berpmlE%orient=orient
                  select case (orient)
                   case (iEx)
                                  PML_%del = 1.0_RKIND/IDxe(i1        )
                     PML_%transversalDelta = 1.0_RKIND/IDxe(i1        )
                     PML_%Plus => Ey(i1+1  ,j1  ,k1  )
                     PML_%Minu => Ey(i1    ,j1  ,k1)
                     PML_%gx2=GM2(jmed)
                     PML_%minTotal=minX(jmed)
                     PML_%maxTotal=maxX(jmed)
                     PML_%posi=i1+0.5
                   case (iEy)
                                  PML_%del = 1.0_RKIND/IDye(    j1     )
                     PML_%transversalDelta = 1.0_RKIND/IDye(    j1     )
                     PML_%Plus => Ex(i1  ,j1+1  ,k1  )
                     PML_%Minu => Ex(i1  ,j1    ,k1  )
                     PML_%gx2=-GM2(jmed)
                     PML_%minTotal=minY(jmed)
                     PML_%maxTotal=maxY(jmed)
                     PML_%posi=j1+0.5
                   case DEFAULT
                     WRITE (buff, *)    'Buggy ERROR: In PMLbodies. '
                     CALL WarnErrReport (buff,.TRUE.)
                  end select
                  !call calc_pmlbodypar
                endif
               endif  
            end do
         end do
      end do
!!!inicializa las constantes
      call calc_pmlbodypar(sgg,eps00,mu00)
!!!

      !!!!!!!!!resuming
      if (.not.resume) then  
         do conta=1,berpmlE%numnodes
            PML_ => berpmlE%Nodes(conta)
            PML_%Psi=0.0_RKIND
         end do
         do conta=1,berpmlH%numnodes
            PML_ => berpmlH%Nodes(conta)
            PML_%Psi=0.0_RKIND
         end do
      else  
            READ (14) (berpmlE%Nodes(conta)%Psi,conta=1,berpmlE%numnodes)
            READ (14) (berpmlH%Nodes(conta)%Psi,conta=1,berpmlH%numnodes)
      endif
      return

   end subroutine InitPMLbodies



  subroutine calc_pmlbodypar(sgg,eps00,mu00)

        type (SGGFDTDINFO), intent(IN)     ::  sgg  
        REAL (KIND=RKIND)           ::  eps00,mu00
        integer (kind=4)  ::  conta
        type (BerPML__t), pointer :: PML_
        integer (kind=4) :: nn
        real (kind=RKIND) :: sigmamax,sigma
!
        eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
!         
        do conta=1,berpmlE%numnodes
                PML_ => berpmlE%Nodes(conta) 
                sigmamax = -((log(CoeffReflPML)*(PMLorden+1))/ &
                            (2*sqrt(Mu0/eps0)*(PML_%del*(PML_%maxTotal-PML_%minTotal)/2.0_RKIND)))
                if (PML_%posi <= (PML_%maxTotal+PML_%minTotal)/2.0_RKIND) then
                sigma=Sigmamax * (abs(PML_%posi-1.0*PML_%minTotal) /((1.0*PML_%maxTotal-1.0*PML_%minTotal)/2.0_RKIND))**PMLorden 
                else
                sigma=Sigmamax * (abs(1.0*PML_%maxTotal-PML_%posi) /((1.0*PML_%maxTotal-1.0*PML_%minTotal)/2.0_RKIND))**PMLorden
                endif
                PML_%P_be=Exp(-(sigma)*sgg%dt/Eps0)
                PML_%P_ce=  (PML_%P_be-1.0_RKIND)/PML_%transversalDelta
                !!!  write(30,*) jmed,orient,i1,j1,k1,sigma,PML_%P_be,PML_%P_ce
        end do
        do conta=1,berpmlH%numnodes
                PML_ => berpmlH%Nodes(conta) 
                sigmamax = -((log(CoeffReflPML)*(PMLorden+1))/ &
                            (2*sqrt(Mu0/eps0)*(PML_%del*(PML_%maxTotal-PML_%minTotal)/2.0_RKIND)))
                if (PML_%posi <= (PML_%maxTotal+PML_%minTotal)/2.0_RKIND) then
                sigma=Sigmamax * (abs(PML_%posi-1.0*PML_%minTotal) /((1.0*PML_%maxTotal-1.0*PML_%minTotal)/2.0_RKIND))**PMLorden 
                else
                sigma=Sigmamax * (abs(1.0*PML_%maxTotal-PML_%posi) /((1.0*PML_%maxTotal-1.0*PML_%minTotal)/2.0_RKIND))**PMLorden
                endif
                PML_%P_be=Exp(-(sigma)*sgg%dt/Eps0)
                PML_%P_ce=  (PML_%P_be-1.0_RKIND)/PML_%transversalDelta
                !!!  write(30,*) jmed,orient,i1,j1,k1,sigma,PML_%P_be,PML_%P_ce
        end do
        return

     end subroutine calc_pmlbodypar  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the PMLbody: Usual Yee
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine AdvancePMLbodyE
      integer (kind=4)  ::  conta
      type (BerPML__t), pointer :: PML_
      do conta=1,berpmlE%numnodes
         PML_ => berpmlE%Nodes(conta) 
         PML_%Psi=PML_%P_be * PML_%Psi + PML_%P_ce * (PML_%Plus   - PML_%Minu)    
         PML_%field = PML_%field + PML_%gx2 * PML_%Psi 
      end do
      return
   end subroutine AdvancePMLbodyE



   subroutine AdvancePMLbodyH
      integer (kind=4)  ::  conta
      type (BerPML__t), pointer :: PML_
      do conta=1,berpmlH%numnodes
         PML_ => berpmlH%Nodes(conta) 
         PML_%Psi=PML_%P_be * PML_%Psi + PML_%P_ce * (PML_%Plus   - PML_%Minu)    
         PML_%field = PML_%field - PML_%gx2 * PML_%Psi 
      end do
      return
   end subroutine AdvancePMLbodyH




 

   subroutine StorefieldsPMLbodies
         integer (kind=4) :: conta
            write (14,err=634) (berpmlE%Nodes(conta)%Psi,conta=1,berpmlE%numnodes)
            write (14,err=634) (berpmlH%Nodes(conta)%Psi,conta=1,berpmlH%numnodes)

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'PMLBODIES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StorefieldsPMLbodies


   subroutine DestroyPMLbodies(sgg)

      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      integer (kind=4)  ::  i

      !free up memory
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%PMLbody).and.(.not.sgg%Med(i)%Is%PML)) then
             if (associated (sgg%Med(i)%PMLbody)) deallocate (sgg%Med(i)%PMLbody)
         endif
      end do


      if(allocated(berpmlE%nodes)) deallocate (berpmlE%nodes)
      if(allocated(berpmlH%nodes)) deallocate (berpmlH%nodes)

   end subroutine




end module PMLbodies
