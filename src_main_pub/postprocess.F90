
    
module PostProcessing

   use fdetypes
   use Report
   use Observa

#ifdef CompileWithHDF
   USE HDF5
#endif

   implicit none

   private
   public PostProcess,postprocessonthefly
contains

   !
   !Subrutine to parse the .h5 file in search for outputrequest sections and wrap the ASCII data into .h5
   !
   subroutine PostProcess(layoutnumber,size,sgg,nEntradaRoot,rinstant,somethingdone,niapapostprocess,forceresampled)


      type (SGGFDTDINFO), intent(IN)     ::  sgg
      integer (kind=4), intent(in)  :: layoutnumber,size
      character (len=*), intent(in) :: nEntradaRoot

      type(output_t), pointer, dimension( : )  ::  output
      character (LEN=BUFSIZE) :: cabecera,cabeceraNew,path,path2,path3,path_resampled
      integer (kind=4) :: numComp
      REAL (KIND=RKIND), allocatable, dimension(:,:) :: valores
      complex (kind=CKIND), allocatable, dimension(:,:) :: valoresDF,valoresDF2
      REAL (KIND=RKIND_tiempo), allocatable, dimension(:) :: tiempo,samplingtime
	  integer (kind=4) :: pozi
      logical :: existe,neverprecounted,escribir,escribirBloque,niapapostprocess,forceresampled,somethingdone
      integer (kind=4) :: fqLength,ii,i,i1,j1,field,ns,timesteps,compo,iii,pp,pasadas
      real (kind=RKIND) :: dummy
      real (kind=RKIND_tiempo) :: rinstant
      !
      real (kind=RKIND), allocatable, dimension(:)   :: fqPos,signal
      real (kind=RKIND)   :: fmin,fmax,fstep,value_interp
      complex (kind=CKIND), allocatable, dimension(:) :: fqValues
      character(len=BUFSIZE) :: buff
      character (LEN=BUFSIZE)     ::  dubuf
      !
      character (LEN=BUFSIZE)  ::  whoami,whoamishort
      real (kind=RKIND) :: t_pedido
      integer (kind=4) columna,jjj
      integer :: my_iostat
     !!!!!
      if (niapapostprocess) then  !niapa 200120 ojooo
           print *,'Copiar a mano los .dat en tiempo para que se postrocesen bien...'
           pause
           print *,'Continuing...'
      endif
     !!!!! 
      write(whoamishort,'(i5)') layoutnumber+1
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '

      Output => GetOutput() !get the output private info from observation

      !
      somethingdone=.false.
      open (119,file=trim(adjustl(nEntradaRoot))//'_Outputrequests_'// &
      trim(adjustl(whoamishort))//'.txt')
      cabecera=' '
      do while (trim(adjustl(cabecera)) /= '!END')
         read(119,'(a)') cabecera
      end do
      backspace(119)
      !
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
            field=sgg%Observation(ii)%P(i)%What
            if (field /= nothing) then
            if (.not.sgg%Observation(ii)%Volumic) then

#ifdef CompileWithMPI
               escribirBloque=((field == iBloqueJx).or.(field == iBloqueJy).or. &
               (field == iBloqueMx).or.(field == iBloqueMy)).and. &
               (layoutnumber == output(ii)%item(i)%MPIRoot)
#else
               escribirBloque=.true.
#endif
               escribir= (( ((field /= iBloqueJx).and.(field /= iBloqueJy)  .and. &
                             (field /= iBloqueMx).and.(field /= iBloqueMy)) .or. &
                             escribirBloque) .and. &
                            (SGG%Observation(ii)%FreqDomain .or. SGG%Observation(ii)%Transfer))
               if (escribir) then
                !
                  !
                  if (SGG%Observation(ii)%FreqDomain) then
                     pasadas=1
                        pp=1
                        if (pp==1) path=trim(adjustl(output(ii)%item(i)%path))

                        inquire(file=trim(adjustl(path)),exist=existe)
                        if (.not.existe) then
                           buff='Not processing: Inexistent file '//trim(adjustl(path))
                           call print11(layoutnumber,BUFF)
                        ELSE
                           numComp=output(ii)%item(i)%columnas
                           neverprecounted=.true.
                           if (neverprecounted) then !all should have the same timesteps !only once read
                              neverprecounted=.false.
                              open (output(ii)%item(i)%unit,file=trim(adjustl(path)),form='formatted') !precounting
                              timesteps=0
                              read (output(ii)%item(i)%unit,'(a)') cabecera
                              do
                                 read (output(ii)%item(i)%unit,*,end=178) dummy,dummy
                                 timesteps=timesteps+1
                              end do
178                           continue
                              close (output(ii)%item(i)%unit)
                              timesteps=timesteps-1 !le resto 1 para evitar lo del ultimo time step a tiempos grandes que saca en el .exc !sgg 180516
                              if (allocated(valores)) deallocate (valores,tiempo,signal)
                              allocate (valores(1:timesteps,1:numComp))
                              allocate (tiempo(1:timesteps))
                              allocate (signal(1:timesteps))
                              allocate (samplingtime(1:timesteps))
                           endif

                           open (output(ii)%item(i)%unit,file=trim(adjustl(path)),form='formatted')
                           read (output(ii)%item(i)%unit,'(a)') cabecera
                           do ns=1,timesteps
                              read (output(ii)%item(i)%unit,*) tiempo(ns), (valores(ns,compo),compo=1,numComp)
                           end do
                           close (output(ii)%item(i)%unit)
                           
                           if (niapapostprocess) then         !niapa 200120 ojooo
                             print *,'Correcting in FreqDomain postprocess ',timesteps,trim(adjustl(path))
                             do ns=1,timesteps
                               tiempo(ns)=real(ns*sgg%dt,RKIND_tiempo)
                             end do
                           endif
                           ! 
                           !!buscanoigual2: do ns=3,timesteps
                           !!   if (.not.almostequal(tiempo(ns)-tiempo(ns-1),tiempo(2)-tiempo(1))) then
                           !!      !check that steps are equally spaced
                           !!      call print11(layoutnumber,'WARNING: (DFT) Unequally spaced timesteps in '//trim(adjustl(path)))
                           !!      exit buscanoigual2
                           !!   endif
                           !!end do buscanoigual2
                           
                           
                          
                            if (forceresampled) then !no sacar resampleadas salvo forzando 120123 
                            !!!write resampled data en tiempo re-sampleado con la peticion temporal OLD ruido 200319
                            !  if (output(ii)%Trancos.ne.1) then      !write always a peticion de OLD 200120
                                  path_resampled=trim(adjustl(path(1:index(path,'.dat')-1)))//'_resampled_time.dat'
                                  columna=1
                                  open (output(ii)%item(i)%unit,file=trim(adjustl(path_resampled)),form='formatted')
                                  t_pedido=tiempo(1)                           
      !                            write(output(ii)%item(i)%unit,'(a)') trim(adjustl(cabecera))
                                  write (output(ii)%item(i)%unit,fmt) t_pedido, valores(1,columna) !el primero se escribe sigual
                                  do iii=2,timesteps
                                      do while (t_pedido<=tiempo(iii))
                                          t_pedido=(t_pedido+sgg%Observation(ii)%TimeStep)
                                          buscinterpol: do jjj=iii-1,timesteps-1
                                              if ((t_pedido>=tiempo(jjj).and.t_pedido<tiempo(jjj+1))) then
                                                  value_interp=(valores(jjj+1,columna)-valores(jjj,columna))/(tiempo(jjj+1)-tiempo(jjj)) * (t_pedido-tiempo(jjj)) + valores(jjj,columna)
                                                  write (output(ii)%item(i)%unit,fmt) t_pedido, value_interp
                                                  exit buscinterpol
                                              else
                                                  cycle buscinterpol
                                              endif
                                          enddo buscinterpol
                                      enddo
                                  enddo  
                                  close (output(ii)%item(i)%unit)
                            !  endif 
                            endif
                                                               
!                              !!!fin write resampled data en tiempo sampleado con la peticion temporal 200319

                           fmin=(min(sgg%observation(ii)%FinalFreq,sgg%observation(ii)%InitialFreq))
                           fmax=(max(sgg%observation(ii)%FinalFreq,sgg%observation(ii)%InitialFreq))

                           IF ((sgg%observation(ii)%FreqStep == 0.0_RKIND).or. &
                           (sgg%observation(ii)%FreqStep > fmax-fmin)) then
                              fstep = fmax-fmin
                           else
                              fstep = (sgg%observation(ii)%FreqStep)
                           endif
                           fqLength=int((fmax - fmin)/fstep)+2
                           allocate (fqPos(1:fqLength),fqValues(1:fqLength),valoresDF(1:fqLength,1:numComp))

                           pozi=index(path,'_log_')
                           if (pozi/=0) then
                              fmin=max(1.0_RKIND,log10(fmin))
                              fmax=log10(fmax)
                              fstep=(fmax-fmin)/(fqLength-2.0_RKIND)
                           endif
                           !
                           do i1=1,fqLength
                              fqPos(i1)=fmin+(i1-1.0_RKIND)*fstep
                           end do

                           if (pozi/=0) then
                              do i1=1,fqLength
                                 fqPos(i1)=10.0_RKIND **fqPos(i1)
                              end do
                           endif
                           !
                           ! Computes dtft.
                           signal=0.0_RKIND
                           samplingtime=0.0_RKIND
                           samplingtime(1:timesteps)=(tiempo(1:timesteps))
                           do i1=1,numComp
                              signal(1:timesteps)=(valores(1:timesteps,i1))
                              call dtft(fqValues, fqPos, fqLength, samplingTime, signal, timesteps)
                              valoresDF(1:fqLength,i1)=fqValues(1:fqLength) !overwrite since no longer needed
                           end do
                           path2=trim(adjustl(path(1:index(path,'.dat')-1)))//'_df.dat'
                           open (output(ii)%item(i)%unit,file=trim(adjustl(path2)),form='formatted')
                           !
                           write (119,'(a)') trim(adjustl(path2)) !update the outputlistsfile
                           !
                           call conviertecabecera(cabecera,cabeceraNew,numComp+1,rinstant)
                           write(output(ii)%item(i)%unit,'(a)') trim(adjustl(cabeceraNew))
                           do i1=1,fqLength
                              write(output(ii)%item(i)%unit,fmt) fqPos(i1), &
                              (abs(valoresDF(i1,j1)),atan2(aimag(valoresDF(i1,j1)),real(valoresDF(i1,j1))), j1=1,numComp)
                           enddo
                           close (output(ii)%item(i)%unit)
                           !
                           deallocate(valores,tiempo,signal,samplingtime) !no longer needed
                           !
                           if (sgg%observation(ii)%Transfer) then
                              !
                              path3=trim(adjustl(sgg%observation(ii)%FileNormalize))
                              inquire(file=trim(adjustl(path3)),exist=existe)
                              if (.not.existe) then
                                 buff='Not processing: Inexistent file '//trim(adjustl(path3))
                                 call print11(layoutnumber,BUFF)
                              ELSE
                                 neverprecounted=.true.
                                 if (neverprecounted) then !all should have the same timesteps !only once read
                                    neverprecounted=.false.
                                    open (34,file=trim(adjustl(path3)),form='formatted') !precounting
                                    timesteps=0
                                    do
                                       read (34,*,end=778) dummy,dummy
                                       timesteps=timesteps+1
                                    end do
778                                 continue
                                    close (34)
                                    timesteps=timesteps-1 !le resto 1 para evitar lo del ultimo time step a tiempos grandes que saca en el .exc !sgg 180516
                                    if (allocated(valores)) deallocate (valores,tiempo,signal)
                                    allocate (valores(1:timesteps,1:1))
                                    allocate (tiempo(1:timesteps))
                                    allocate (signal(1:timesteps))
                                    allocate (samplingtime(1:timesteps))
                                 endif

                                 open (34,file=trim(adjustl(sgg%observation(ii)%FileNormalize)),form='formatted')
                                 do ns=1,timesteps
                                    read (34,*) tiempo(ns), valores(ns,1)
                                 end do
                                 close (34)
                                 
                                 if (niapapostprocess) then  !niapa 200120 ojooo
                                    print *,'Correcting in Transfer postprocess ',timesteps,trim(adjustl(path))
                                    do ns=1,timesteps
                                       tiempo(ns)=real(ns*sgg%dt,RKIND_tiempo)
                                    end do
                                 endif
                                 
                          
                                if (forceresampled) then !no sacar resampleadas salvo forzando 120123 
                          !!!write resampled data en tiempo re-sampleado con la peticion temporal OLD ruido 200319
                                ! if (output(ii)%Trancos.ne.1) then  !write always a peticion de OLD 200120
                                    path2=sgg%observation(ii)%FileNormalize
                                    path_resampled=trim(adjustl(path(1:index(path,'.dat')-1)))//'_resampled_time.dat'  !no se, corregido 120920
                                    path2=trim(adjustl(path2))//'_at_'//trim(adjustl(path_resampled)) !no necesariamente .dat sino .exc u otra cosa
                                    columna=1
                                    open (34,file=trim(adjustl(path2)),form='formatted')
                                    t_pedido=tiempo(1)                           
                                    write (34,fmt) t_pedido, valores(1,columna) !el primero se escribe sigual
                                    do iii=2,timesteps
                                        do while (t_pedido<=tiempo(iii))
                                            t_pedido=(t_pedido+sgg%Observation(ii)%TimeStep)
                                            buscinterpol2: do jjj=iii-1,timesteps-1
                                                if ((t_pedido>=tiempo(jjj).and.t_pedido<tiempo(jjj+1))) then
                                                    value_interp=(valores(jjj+1,columna)-valores(jjj,columna))/(tiempo(jjj+1)-tiempo(jjj)) * (t_pedido-tiempo(jjj)) + valores(jjj,columna)
                                                    write (34,fmt) t_pedido, value_interp
                                                    exit buscinterpol2
                                                else
                                                    cycle buscinterpol2
                                                endif
                                            enddo buscinterpol2
                                        enddo
                                    enddo  
                                    close (34)  
                                 ! endif
                                 
                                endif 
                          !!!end write resampled 
                                 
                                 !
!!!                                 buscanoigual3: do ns=3,timesteps
!!!                                    if (.not.almostequal(tiempo(ns)-tiempo(ns-1),tiempo(2)-tiempo(1))) then
!!!                                       !check that steps are equally spaced
!!!                                       call print11(layoutnumber,'WARNING: (DFT) Unequally spaced timesteps in '//trim(adjustl(path3)))
!!!!                                       exit buscanoigual3
!!!                                    endif
!!!                                 end do buscanoigual3
                                 !assume equally spaced!!! !may be different between computed and reference
                                 !the normalization file only has one column
                                 allocate (valoresDF2(1:fqLength,1:numComp))
                                 ! Computes dtft.
                                 signal=0.0_RKIND
                                 samplingtime=0.0_RKIND_tiempo
                                 signal(1:timesteps)=(valores(1:timesteps,1))
                                 samplingtime(1:timesteps)=(tiempo(1:timesteps))
                                 call dtft(fqValues, fqPos, fqLength, samplingTime, signal, timesteps)
                                 !normalization
                                 do iii=1,fqLength
                                    !do j1=1,numComp-1
                                    do j1=1,numComp
                                       valoresDF2(iii,j1)=(valoresDF(iii,j1)/fqValues(iii)) !overwrite since no longer needed
                                    enddo
                                    !do j1=numComp,numComp
                                    !    valoresDF2(iii,j1)=fqValues(iii) !la entrada se deja sin normalizar pero se saca
                                    !enddo
                                 end do
                                 !
                                 path2=trim(adjustl(path(1:index(path,'.dat')-1)))//'_normalization_df.dat'
                                 open (output(ii)%item(i)%unit,file=trim(adjustl(path2)),form='formatted')
                                 write(output(ii)%item(i)%unit,'(a)') 'Freq Mod_Normalization_data Phase_Normalization_data'
                                 !!!!!! write normalization data !bug serio en sondas transfer detectado 20/07/15. Por azar las Z_input se calculaban bien para thin_wires pq. en el voltage se daba el valor de normalizacion (tocate los guevos!!!!)
                                 do i1=1,fqLength
                                    write(output(ii)%item(i)%unit,fmt) fqPos(i1), abs(fqValues(i1)),atan2(aimag(fqValues(i1)),real(fqValues(i1)))
                                 end do
                                 close (output(ii)%item(i)%unit)
                                 !!!!!!end write normalization data                                
                                 
                                 path2=trim(adjustl(path(1:index(path,'.dat')-1)))//'_tr.dat'
                                 open (output(ii)%item(i)%unit,file=trim(adjustl(path2)),form='formatted')
                                 !
                                 write (119,'(a)') trim(adjustl(path2)) !update the outputlistsfile
                                 !
                                 call conviertecabecera(cabecera,cabeceraNew,numComp+1,rinstant)
                                 !
                                 pozi=index(path,'_log_')
                                 if (pozi/=0) cabeceraNew=trim(adjustl(cabeceraNew))//'_(In_dB_(20_log10)_and_radians)'
                                 !
                                 write(output(ii)%item(i)%unit,'(a)') trim(adjustl(cabeceraNew))
                                 do i1=1,fqLength
                                    !ojo no saco la SE sino la funcion de transferencia. Habia antes del 13/07 un signo -1 delante del 20.0_RKIND * log10 de mas abajo. Cambiado!!! por eso salia el coef de reflexion numer mayor que la unidad
                                     do j1=1,numComp
                                        if  (abs(valoresDF2(i1,j1)).lt.1e-30) valoresDF2(i1,j1)=1e-30 !para evitar calcular el log10(0)
                                     end do
                                     write(output(ii)%item(i)%unit,fmt) fqPos(i1), &
                                        (20.0_RKIND * log10(abs(valoresDF2(i1,j1))),&
                                            atan2(aimag(valoresDF2(i1,j1)),real(valoresDF2(i1,j1))), j1=1,numComp)

                                 enddo
                                 close (output(ii)%item(i)%unit)
                                 deallocate(valoresDF2,valores,tiempo,signal,samplingtime)
                              endif !del not existe
                           endif !del transfer
                           !
                           deallocate(fqPos,fqValues,valoresDF)
                        endif !del existe
                     
                     somethingdone=.true.
                  endif !del if fredomain
                  
                  !!!!!!!!!peticion correo 300120 !siempre escribe resammpled
                  !
                  if ((SGG%Observation(ii)%TimeDomain).or.(SGG%Observation(ii)%FreqDomain).or.(SGG%Observation(ii)%Transfer)) then
                        path=trim(adjustl(output(ii)%item(i)%path))
                        inquire(file=trim(adjustl(path)),exist=existe)
                        if (.not.existe) then
                           buff='Not processing: Inexistent file '//trim(adjustl(path))
                           call print11(layoutnumber,BUFF)
                        ELSE
                           numComp=output(ii)%item(i)%columnas
                           neverprecounted=.true.
                           if (neverprecounted) then !all should have the same timesteps !only once read
                              neverprecounted=.false.
                              open (output(ii)%item(i)%unit,file=trim(adjustl(path)),form='formatted') !precounting
                              timesteps=0
                              read (output(ii)%item(i)%unit,'(a)') cabecera
                              do
                                 read (output(ii)%item(i)%unit,*,end=5178) dummy,dummy
                                 timesteps=timesteps+1
                              end do
5178                           continue
                              close (output(ii)%item(i)%unit)
                              timesteps=timesteps-1 !le resto 1 para evitar lo del ultimo time step a tiempos grandes que saca en el .exc !sgg 180516
                              if (allocated(valores)) deallocate (valores,tiempo,signal)
                              allocate (valores(1:timesteps,1:numComp))
                              allocate (tiempo(1:timesteps))
                              allocate (signal(1:timesteps))
                              allocate (samplingtime(1:timesteps))
                           endif

                           open (output(ii)%item(i)%unit,file=trim(adjustl(path)),form='formatted')
                           read (output(ii)%item(i)%unit,'(a)') cabecera
                           do ns=1,timesteps
                              read (output(ii)%item(i)%unit,*) tiempo(ns), (valores(ns,compo),compo=1,numComp)
                           end do
                           close (output(ii)%item(i)%unit)
                           
                           if (niapapostprocess) then         !niapa 200120 ojooo
                             print *,'Correcting in FreqDomain postprocess ',timesteps,trim(adjustl(path))
                             do ns=1,timesteps
                               tiempo(ns)=real(ns*sgg%dt,RKIND_tiempo)
                             end do
                           endif

                          
                          if (forceresampled) then !no sacar resampleadas salvo forzando 120123 
                          !!!write resampled data en tiempo re-sampleado con la peticion temporal OLD ruido 200319
                          !  if (output(ii)%Trancos.ne.1) then      !write always a peticion de OLD 200120
                                path_resampled=trim(adjustl(path(1:index(path,'.dat')-1)))//'_resampled_time.dat'
                                columna=1
                                open (output(ii)%item(i)%unit,file=trim(adjustl(path_resampled)),form='formatted')
                                t_pedido=tiempo(1)                           
    !                            write(output(ii)%item(i)%unit,'(a)') trim(adjustl(cabecera))
                                write (output(ii)%item(i)%unit,fmt) t_pedido, valores(1,columna) !el primero se escribe sigual
                                do iii=2,timesteps
                                    do while (t_pedido<=tiempo(iii))
                                        t_pedido=(t_pedido+sgg%Observation(ii)%TimeStep)
                                        buscinterpoli: do jjj=iii-1,timesteps-1
                                            if ((t_pedido>=tiempo(jjj).and.t_pedido<tiempo(jjj+1))) then
                                                value_interp=(valores(jjj+1,columna)-valores(jjj,columna))/(tiempo(jjj+1)-tiempo(jjj)) * (t_pedido-tiempo(jjj)) + valores(jjj,columna)
                                                write (output(ii)%item(i)%unit,fmt) t_pedido, value_interp
                                                exit buscinterpoli
                                            else
                                                cycle buscinterpoli
                                            endif
                                        enddo buscinterpoli
                                    enddo
                                enddo  
                                close (output(ii)%item(i)%unit)
                            !endif del trancos
                          endif !del resampled      
                          deallocate(valores,tiempo,signal,samplingtime) !no longer needed
                        endif !del existe
                  endif !del if timedomain

                  !!!!!!!!!!!!
               endif !del escribir
            endif !del volumic
            endif !del nothing
         end do
         !
      end do !barrido puntos de observacion

      !no olvidar updatear el outpurequestfile con los nuevos fichero frecuencia 

      write (119,'(a)') '!END '
      close (119)

      return


   end subroutine postprocess

   !
   !
   function almostequal(a,b) result(almost)
      real (kind=RKIND) a,b,ratio
      logical almost
      real (kind=RKIND), parameter :: tolerancia=0.01

      ratio=a/b
      if (abs(ratio)<1.0_RKIND) ratio=1.0_RKIND / ratio
      if ((ratio > 0.0_RKIND).and.(ratio < 1.0_RKIND+tolerancia)) then
         almost=.true.
      else
         almost=.false.
      endif
   end function


   subroutine conviertecabecera(c,cNew,columnas,rinstant)
      integer :: columnas
      character(LEN=BUFSIZE) :: c,cNew,chninstant
      character(LEN=BUFSIZE), dimension(1:columnas) :: c2
      integer :: i,j,k,longi,ii
      real (kind=RKIND_tiempo) rinstant

      write(chninstant,fmt) rinstant


      c2=' '
      longi=len(trim(c))+1

      k=1
      i=1
      buscam: do while(i <= longi)
         if (c(i:i) /= ' ') then
            interno: do j=i,longi
               if (c(j:j) == ' ') then
                  if (k==1) then
                     c2(k)=' f_at_'//trim(adjustl(chninstant))
                  else
                     c2(k)=c(i:j-1)//'_Mod   '//c(i:j-1)//'_Phase   '
                  endif
                  k=k+1
                  if (k > columnas) exit buscam
                  do ii=j,longi
                     i=ii
                     if ((c(ii:ii) /= ' ').or.(ii == longi)) exit interno
                  end do
               endif
            end do interno
         else
            i=i+1
         endif
      end do buscam

      cNew=' '
      do i=1,columnas
         cNew=trim(adjustl(cNew))//'   '//trim(adjustl(c2(i)))
      enddo
      return
   end subroutine



   subroutine postprocessonthefly(layoutnumber,size,sgg,nEntradaRoot,rinstant,somethingdone,niapapostprocess,forceresampled)


      type (SGGFDTDINFO), intent(IN)     ::  sgg
      integer (kind=4), intent(in)  :: layoutnumber,size
      character (len=*), intent(in) :: nEntradaRoot

      type(output_t), pointer, dimension( : )  ::  output
      character (LEN=BUFSIZE) :: path
      logical :: existe,escribir,escribirBloque,somethingdone,niapapostprocess,forceresampled
      integer (kind=4) :: ii,i,field
      character(len=BUFSIZE) :: buff
      real (kind=RKIND_tiempo) :: rinstant

      Output => GetOutput() !get the output private info from observation
      !
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
            field=sgg%Observation(ii)%P(i)%What
#ifdef CompileWithMPI
            escribirBloque=((field == iBloqueJx).or.(field == iBloqueJy).or. &
            (field == iBloqueMx).or.(field == iBloqueMy)).and. &
            (layoutnumber == output(ii)%item(i)%MPIRoot)
#else
            escribirBloque=.true.
#endif
            escribir= (( ((field /= iBloqueJx).and.(field /= iBloqueJy)  .and. &
                          (field /= iBloqueMx).and.(field /= iBloqueMy)) .or. &
                          escribirBloque) .and. &
                         (SGG%Observation(ii)%FreqDomain .or. SGG%Observation(ii)%Transfer))
            if (escribir) then
               if (SGG%Observation(ii)%FreqDomain) then
                  path=trim(adjustl(output(ii)%item(i)%path))
                  inquire(file=trim(adjustl(path)),exist=existe)
                  if (.not.existe) then
                     buff='Not processing: Inexistent file '//trim(adjustl(path))
                     call print11(layoutnumber,BUFF)
                  ELSE
                     close (output(ii)%item(i)%unit)
                  endif !del existe
               endif !del if fredomain
            endif !del escribir
         end do
         !
      end do !barrido puntos de observacion
      !
      call postprocess(layoutnumber,size,sgg,nEntradaRoot,rinstant,somethingdone,niapapostprocess,forceresampled)
      !
      do ii=1,sgg%NumberRequest
         do i=1,sgg%Observation(ii)%nP
            field=sgg%Observation(ii)%P(i)%What
#ifdef CompileWithMPI
            escribirBloque=((field == iBloqueJx).or.(field == iBloqueJy).or. &
            (field == iBloqueMx).or.(field == iBloqueMy)).and. &
            (layoutnumber == output(ii)%item(i)%MPIRoot)
#else
            escribirBloque=.true.
#endif  
            escribir= (( ((field /= iBloqueJx).and.(field /= iBloqueJy)  .and. &
                          (field /= iBloqueMx).and.(field /= iBloqueMy)) .or. &
                          escribirBloque) .and. &
                         (SGG%Observation(ii)%FreqDomain .or. SGG%Observation(ii)%Transfer))
            if (escribir) then
               if (SGG%Observation(ii)%FreqDomain) then
                  path=trim(adjustl(output(ii)%item(i)%path))
                  inquire(file=trim(adjustl(path)),exist=existe)
                  if (.not.existe) then
                     buff='ERROR: Inexistent file. Creating new '//trim(adjustl(path))
                     call print11(layoutnumber,BUFF)
                  ELSE
                     open (output(ii)%item(i)%unit,file=trim(adjustl(path)),form='formatted',position='append')
                  endif !del existe
               endif !del if fredomain
            endif !del escribir
         end do
      end do !barrido puntos de observacion
   end subroutine postprocessonthefly

end module
