
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Creation date Date :  October, 24, 2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CALC_CONSTANTS
   use fdetypes
   use report
    implicit none
    private
    public calc_G1G2Gm1Gm2
    contains
   subroutine calc_G1G2Gm1Gm2(sgg,G1,G2,Gm1,Gm2,eps0,mu0)
        type (SGGFDTDINFO), intent(IN)   ::  sgg
        REAL (KIND=RKIND) , intent(inout) :: Eps0, Mu0
        REAL (KIND=RKIND) , dimension (0 : sgg%NumMedia),intent(inout) ::  g1,g2,gM1,gM2
        integer :: r,i
        REAL (KIND=RKIND) :: Sigmam , Epsilon , Mu , Sigma,width, epr
        character(len=BUFSIZE) :: buff
        
        do r=0,sgg%NumMedia
         Sigmam = sgg%Med(r)%SigmaM
         Epsilon = Eps0*sgg%Med(r)%Epr
         Mu = Mu0*sgg%Med(r)%Mur
         Sigma = sgg%Med(r)%Sigma
         !In case of Multiport set updating Ca, Cbfficients trivially, since the field
         !updating are handled locally in the Multiport module
         if  ((sgg%Med(R)%Is%already_YEEadvanced_byconformal).or.(sgg%Med(R)%Is%split_and_useless)) then
            G1(r)=1.0_RKIND
            G2(r)=0.0_RKIND
            GM1(r)=1.0_RKIND
            GM2(r)=0.0_RKIND
         else
            if ((sgg%Med(R)%Is%multiport).or.(sgg%Med(R)%Is%AnisMultiport)) then
               G1(r)=0.0_RKIND !mull fields on the main procedure (good both for Ian and for me)
               G2(r)=0.0_RKIND
               GM1(r)=0.0_RKIND
               GM2(r)=0.0_RKIND
            elseif ((sgg%Med(R)%Is%pec).or.(r==0)) then
               !Trivially PEC updating Ca, Cbfficients are set to 0.0_RKIND
               G1(r)=0.0_RKIND ;  G2(r)=0.0_RKIND;  GM1(r)=0.0_RKIND;  GM2(r)=0.0_RKIND ;
            elseif (sgg%Med(R)%Is%lumped) then
               !Trivially PEC NOT updating coefficients para el avance en E. They are set to 1.0_RKIND. La rutina propia ya se encargara
               G1(r)=1.0_RKIND ;  G2(r)=0.0_RKIND; 
               GM1(r)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               GM2(r)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               if (gm1(r) < 0.0_RKIND) then !exponential time stepping
                   gm1(r)=exp(- Sigmam * sgg%dt / (Mu ))
                   gm2(r)=(1.0_RKIND-gm1(r))/ Sigmam
               endif
            elseif (sgg%Med(R)%Is%SGBC) then
!!!!! 090519 He quitado todo este calculo que luego hara InitSGBCs para no duplicar codigo propenso a errores. Uso valores absurdos por lo que truene.
!!!!! ojo que los parametros stochastic tambien se obtendran en InitSGBCs, por eso lo he quitado esto de aqui
                      G1(r)=0.0; G2(r)=0.0;
                       !!! G1(r)=1.0e31; G2(r)=-1.0e23 !quitado este sinsentido por si diese lugar a errores de redondeo 170519
!!!!!                       
!!!!!                       !quitado soporte multicapas en SGBC a 260419
!!!!!               !!!!!los intrinsecos suyos promediados si es una multicapa!!. La rutina de avance machaca los campos pero le paso estos coefs a initSGBCs para que los aproveche
!!!!!                       width=0.; sigma=0.; epr=0.;
!!!!!                       do i=1,sgg%Med(r)%multiport(1)%numcapas
!!!!!                          width=width         + sgg%Med(r)%multiport(1)%width(i) 
!!!!!                       end do
!!!!!                       do i=1,sgg%Med(r)%multiport(1)%numcapas
!!!!!                          sigma=sigma + sgg%Med(r)%multiport(1)%sigma(i)*sgg%Med(r)%multiport(1)%width(i) / width
!!!!!                          epr=epr     +  eps0*sgg%Med(r)%multiport(1)%epr(i)  *sgg%Med(r)%multiport(1)%width(i) / width
!!!!!                       end do
!!!!! !!!!!!091215
!!!!!                     G1(r)=(1 -  Sigma * sgg%dt / (2.0_RKIND * Epsilon ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                     G2(r)=sgg%dt /Epsilon                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
!!!!!                     if (g1(r) < 0.0_RKIND) then !exponential time stepping
!!!!!                        g1(r)=exp(- Sigma * sgg%dt / (Epsilon ))
!!!!!                        g2(r)=(1.0_RKIND-g1(r))/ Sigma
!!!!!                     endif
!!!!! hasta aqui lo comentado a 090519. Los calculos de GM1 y GM2 no los hace InitSGBCs y sus incertidumbres son nulas. asi que los hago aqui                       
                          GM1(r)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
                          GM2(r)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
                          if (gm1(r) < 0.0_RKIND) then !exponential time stepping
                             gm1(r)=exp(- Sigmam * sgg%dt / (Mu ))
                             gm2(r)=(1.0_RKIND-gm1(r))/ Sigmam
                          endif
            elseif (sgg%Med(R)%Is%Anisotropic) then
               G1(r)=1.0_RKIND !para que no haga nada en el bucle principal evitando los ifs
               G2(r)=0.0_RKIND
               GM1(r)=1.0_RKIND !para que no haga nada en el bucle principal evitando los ifs
               GM2(r)=0.0_RKIND
            elseif  ((sgg%Med(R)%Is%EDispersive).and.(.not.sgg%Med(R)%Is%MDispersive).and.(.not.sgg%Med(r)%Is%EdispersiveANIS)) then
               !solo cierto para ISOTROPOS
               G1(r)=0.0_RKIND !will be overwritten by own values created by InitEDispersives
               G2(r)=0.0_RKIND !will be overwritten by own values created by InitEDispersives
               GM1(r)=(1.0_RKIND- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) / (1+ SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               GM2(r)=sgg%dt / Mu                   / (1+ SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               if (gm1(r) < 0.0_RKIND) then !exponential time stepping
                  gm1(r)=exp(- Sigmam * sgg%dt / (Mu ))
                  gm2(r)=(1.0_RKIND-gm1(r))/ Sigmam
               endif
            elseif  ((sgg%Med(R)%Is%MDispersive).and.(.not.sgg%Med(R)%Is%EDispersive).and.(.not.sgg%Med(r)%Is%MdispersiveANIS)) then
               !solo cierto para ISOTROPOS
               G1(r)=(1.0_RKIND- Sigma*sgg%dt/(2.0_RKIND * Epsilon)) / (1.0_RKIND + Sigma*sgg%dt/(2.0_RKIND * Epsilon))
               G2(r)=sgg%dt / Epsilon                / (1.0_RKIND + Sigma*sgg%dt/(2.0_RKIND * Epsilon))
               if (g1(r) < 0.0_RKIND) then !exponential time stepping
                  g1(r)=exp(- Sigma * sgg%dt / (Epsilon ))
                  g2(r)=(1.0_RKIND-g1(r))/ Sigma
               endif
               GM1(r)=0.0_RKIND !will be overwritten by own values created by InitMDispersives !when written this routine
               GM2(r)=0.0_RKIND !will be overwritten by own values created by InitMDispersives
            elseif  ((sgg%Med(r)%Is%MdispersiveANIS).OR.(sgg%Med(r)%Is%EdispersiveANIS)) then
               BUFF='ERROR: ANISOTROPIC DISPERSIVE CURRENTLY UNSUPPORTED IN THE ENGINE'
               CALL StopOnError (0,0,buff)  !lo deberia reportar y parar antes SEMBA_FDTD.F90 !quitar algun dia para que no ralentice 170719
            else
               G1(r)=(1 -  Sigma * sgg%dt / (2.0_RKIND * Epsilon ) ) / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
               G2(r)=sgg%dt /Epsilon                        / (1.0_RKIND + Sigma * sgg%dt / (2.0_RKIND * Epsilon ))
               if (g1(r) < 0.0_RKIND) then !exponential time stepping
                  g1(r)=exp(- Sigma * sgg%dt / (Epsilon ))
                  g2(r)=(1.0_RKIND-g1(r))/ Sigma
               endif
               GM1(r)=(1- SigmaM*sgg%dt/(2.0_RKIND *  Mu )) /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               GM2(r)=sgg%dt/ Mu                   /(1.0_RKIND + SigmaM*sgg%dt/(2.0_RKIND *  Mu ))
               if (gm1(r) < 0.0_RKIND) then !exponential time stepping
                  gm1(r)=exp(- Sigmam * sgg%dt / (Mu ))
                  gm2(r)=(1.0_RKIND-gm1(r))/ Sigmam
               endif
            endif !Multiports
         endif
      End do
   end subroutine calc_G1G2Gm1Gm2

end module CALC_CONSTANTS

!!!!!!INNECESARIO!      call InitOtherBorders    ()
!!!!!!DONE       !      call InitCPMLBorders     ()
!!!!!!DONE       !      call InitPMLbodies()
!!!!!!DONE       !      call InitMURBorders      ()
!!!!!!DONE       !!!!            if ((trim(adjustl(wiresflavor))=='holland') .or. (trim(adjustl(wiresflavor))=='transition')) then
!!!!!!DONE       !!!!               call InitWires       ()
!!!!!!UNSUPPORTED       !!!      elseif (trim(adjustl(wiresflavor))=='berenger') then
!!!!!!UNSUPPORTED       !!!         call InitWires_Berenger()
!!!!!!UNSUPPORTED       !!!      elseif((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
!!!!!!UNSUPPORTED       !!!         call InitWires_Slanted()
!!!!!!UNSUPPORTED       !!!      endif
!!!!!!DONE       !               CALL InitLumped()
!!!!!!DONE       !               call InitAnisotropic()
!!!!!!DONE       !               call InitSGBCs()
!!!!!!UNSUPPORTED       !!!      call InitMultiports        ()
!!!!!!UNSUPPORTED       !!!      call InitAnisMultiports  ()
!!!!!!UNSUPPORTED       !        call initialize_memory_FDTD_conf_fields ()
!!!!!!UNSUPPORTED       !        call load_antenna()
!!!!!!UNSUPPORTED       !!!      call InitEDispersives()
!!!!!!UNSUPPORTED       !!!      call InitMDispersives()
!!!!!!DONE              !!!      call InitPlaneWave  
!!!!!!DONE              !!!      call InitNodalSources
!!!!!!DONE EXCEPT DFT & NF2FF   !!!      call InitObservation 

!!call AdvanceAnisotropicE(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
!!call Advance_Ex          (Ex, Hy, Hz, Idyh, Idzh, sggMiEx, b,g1,g2)
!!call Advance_Ey          (Ey, Hz, Hx, Idzh, Idxh, sggMiEy, b,g1,g2)
!!call Advance_Ez          (Ez, Hx, Hy, Idxh, Idyh, sggMiEz, b,g1,g2)
!!! no se ganada nada de tiempo        Call Advance_ExEyEz(Ex,Ey,Ez,Hx,Hy,Hz,Idxh,Idyh,Idzh,sggMiEx,sggMiEy,sggMiEz,b,g1,g2)
!!call conformal_advance_E()
!!call AdvanceWiresEcrank(sgg,n, layoutnumber,wiresflavor)
!!call AdvanceWiresE(sgg,n, layoutnumber,wiresflavor,simu_devia,stochastic)
!!call AdvanceWiresE_Berenger(sgg,n)
!!call AdvanceWiresE_Slanted(sgg,n)
!!call AdvancePMLbodyE
!!call AdvanceelectricCPML          (sgg%NumMedia, b       ,sggMiEx,sggMiEy,sggMiEz,G2,Ex,Ey,Ez,Hx,Hy,Hz)
!!call AdvanceMultiportE(sgg%alloc,Ex, Ey, Ez)
!!call AdvanceAnisMultiportE(sgg%alloc,Ex, Ey, Ez)
!!call AdvanceSGBCE(sgg%dt,SGBCDispersive)
!!call AdvanceLumpedE(sgg)
!!call AdvanceEDispersiveE(sgg)
!!call AdvancePlaneWaveE(sgg,n, b       ,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez)
!!call AdvanceNodalE(sgg,sggMiEx,sggMiEy,sggMiEz,sgg%NumMedia,n, b,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez)
!!
!!call AdvanceAnisotropicH(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
!!call Advance_Hx           (Hx, Ey, Ez, Idye, Idze, sggMiHx, b,gm1,gm2)
!!call Advance_Hy           (Hy, Ez, Ex, Idze, Idxe, sggMiHy, b,gm1,gm2)
!!call Advance_Hz           (Hz, Ex, Ey, Idxe, Idye, sggMiHz, b,gm1,gm2)
!!call FreeSpace_Advance_Hx(Hx, Ey, Ez, Idye, Idze,           b,gm1,gm2)
!!call FreeSpace_Advance_Hy(Hy, Ez, Ex, Idze, Idxe,           b,gm1,gm2)
!!call FreeSpace_Advance_Hz(Hz, Ex, Ey, Idxe, Idye,           b,gm1,gm2)
!!! no se ganada nada de tiempo                 Call Advance_HxHyHz(Hx,Hy,Hz,Ex,Ey,Ez,IdxE,IdyE,IdzE,sggMiHx,sggMiHy,sggMiHz,b,gm1,gm2)
!!call AdvancePMLbodyH
!!call AdvanceMagneticCPML          ( sgg%NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
!!call FreeSpace_AdvanceMagneticCPML( sgg%NumMedia, b,                            gm2, Hx, Hy, Hz, Ex, Ey, Ez)
!!!!INNECESARIO!call MinusCloneMagneticPMC(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
!!!!INNECESARIO!call CloneMagneticPeriodic(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
!!call AdvanceSGBCH
!!call AdvanceMDispersiveH(sgg)
!!call AdvanceMultiportH    (sgg%alloc,Hx,Hy,Hz,Ex,Ey,Ez,Idxe,Idye,Idze,sggMiHx,sggMiHy,sggMiHz,gm2,sgg%nummedia,conformalskin)
!!call AdvanceAnisMultiportH(sgg%alloc,Hx,Hy,Hz,Ex,Ey,Ez,Idxe,Idye,Idze,sggMiHx,sggMiHy,sggMiHz,gm2,sgg%nummedia)
!!call AdvancePlaneWaveH(sgg,n, b        , GM2, Idxe,Idye, Idze, Hx, Hy, Hz)
!!call AdvanceNodalH(sgg,sggMiHx,sggMiHy,sggMiHz,sgg%NumMedia,n, b       ,GM2,Idxe,Idye,Idze,Hx,Hy,Hz)
!!call MinusCloneMagneticPMC(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
!!call CloneMagneticPeriodic(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,layoutnumber,size)
!!call conformal_advance_H()
!!call AdvanceMagneticMUR              (b, sgg,sggMiHx, sggMiHy, sggMiHz, Hx, Hy, Hz,mur_second)
!!call UpdateObservation(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag, n,ini_save, b, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh,wiresflavor,SINPML_FULLSIZE,wirecrank)
!!call FlushObservationFiles(sgg,ini_save,mindum,layoutnumber,size, dxe, dye, dze, dxh, dyh, dzh,b,singlefilewrite,facesNF2FF,.FALSE.) !no se flushean los farfields ahora
!!call flush_and_save_resume(sgg, b, layoutnumber, size, nEntradaroot, nresumeable2, thereare, n,eps0,mu0, everflushed,Ex, Ey, Ez, Hx, Hy, Hz,wiresflavor,simu_devia,stochastic)
!!call PostProcessOnthefly(layoutnumber,size,sgg,nEntradaRoot,at,somethingdone,permitscaling)
!!call createvtkOnTheFly(layoutnumber,size,sgg,vtkindex,somethingdone)
!!call createxdmfOnTheFly(sgg,layoutnumber,size,vtkindex,somethingdone)
!!call unpacksinglefiles(sgg,layoutnumber,size,singlefilewrite,initialtimestep,resume) !dump the remaining to disk


   

