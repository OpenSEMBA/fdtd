!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD SOLVER MODULE
!  creation date Date :  June,  6, 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module solver_mod

    use fdetypes
    implicit none


    type, public :: solver_t
        ! real (kind=rkind), pointer, dimension (:,:,:)  ::  Ex,Ey,Ez,Hx,Hy,Hz
        integer(kind=4) :: numberOfSteps
        type (bounds_t)  ::  bounds


    contains
        procedure :: init => solver_init
        procedure :: run => solver_run
        procedure :: end => solver_end
    end type

    interface solver_t
        module procedure solver_ctor 
    end interface

contains

    function solver_ctor(sgg) result(res)
        type(solver_t) :: res
        type(sggfdtdinfo) :: sgg

        call findBounds(sgg, res%bounds)
    end function

    subroutine solver_init(this)
        class(solver_t) :: this
        ! findBounds
        ! readFields
        ! crea_timevector
        ! calc_g1g2gm1gm2
        ! init reporting
        ! init other borders
        ! init PML bodies
        ! init lumped
        ! init wires
        ! init mur abc slanted
        ! init anisotropic
        ! init sgbcs
        ! init multiports
        ! init edispersives
        ! init mdispersives
        ! init planewave
        ! init nodalsources
        ! init hopf
        ! fill mtag
        ! init observation
        ! init timing
    end subroutine

    subroutine solver_run(this)
        class(solver_t) :: this
        integer :: i
        do i = 1, this%numberOfSteps
            ! call step from timeStepping class?
            ! call this%stepper%step()
            ! call this%step()
        end do

    end subroutine

    subroutine solver_end(this)
        class(solver_t) :: this
        ! end confornal simulation
        ! timing
        ! flush observation files
        ! close observation files
        ! postprocess
        ! create vtk
        ! create xdmf
    end subroutine


    subroutine findBounds(sgg)
        type (sggfdtdinfo), intent(in) :: sgg
        type (bounds_t), intent(out) :: b

        b%dxe%XI=sgg%alloc(iHx)%XI
        b%dxe%XE=sgg%alloc(iHx)%XE
        b%dye%YI=sgg%alloc(iHy)%YI
        b%dye%YE=sgg%alloc(iHy)%YE
        b%dze%ZI=sgg%alloc(iHz)%ZI
        b%dze%ZE=sgg%alloc(iHz)%ZE

        b%dxh%XI=sgg%alloc(iEx)%XI
        b%dxh%XE=sgg%alloc(iEx)%XE
        b%dyh%YI=sgg%alloc(iEy)%YI
        b%dyh%YE=sgg%alloc(iEy)%YE
        b%dzh%ZI=sgg%alloc(iEz)%ZI
        b%dzh%ZE=sgg%alloc(iEz)%ZE

        b%Ex%XI=sgg%Alloc(iEx)%XI
        b%Ex%XE=sgg%Alloc(iEx)%XE
        b%Ey%XI=sgg%Alloc(iEy)%XI
        b%Ey%XE=sgg%Alloc(iEy)%XE
        b%Ez%XI=sgg%Alloc(iEz)%XI
        b%Ez%XE=sgg%Alloc(iEz)%XE

        b%Hx%XI=sgg%Alloc(iHx)%XI
        b%Hx%XE=sgg%Alloc(iHx)%XE
        b%Hy%XI=sgg%Alloc(iHy)%XI
        b%Hy%XE=sgg%Alloc(iHy)%XE
        b%Hz%XI=sgg%Alloc(iHz)%XI
        b%Hz%XE=sgg%Alloc(iHz)%XE
        
        b%Ex%YI=sgg%Alloc(iEx)%YI
        b%Ex%YE=sgg%Alloc(iEx)%YE
        b%Ey%YI=sgg%Alloc(iEy)%YI
        b%Ey%YE=sgg%Alloc(iEy)%YE
        b%Ez%YI=sgg%Alloc(iEz)%YI
        b%Ez%YE=sgg%Alloc(iEz)%YE
        
        b%Hx%YI=sgg%Alloc(iHx)%YI
        b%Hx%YE=sgg%Alloc(iHx)%YE
        b%Hy%YI=sgg%Alloc(iHy)%YI
        b%Hy%YE=sgg%Alloc(iHy)%YE
        b%Hz%YI=sgg%Alloc(iHz)%YI
        b%Hz%YE=sgg%Alloc(iHz)%YE
        
        b%Ex%ZI=sgg%Alloc(iEx)%ZI
        b%Ex%ZE=sgg%Alloc(iEx)%ZE
        b%Ey%ZI=sgg%Alloc(iEy)%ZI
        b%Ey%ZE=sgg%Alloc(iEy)%ZE
        b%Ez%ZI=sgg%Alloc(iEz)%ZI
        b%Ez%ZE=sgg%Alloc(iEz)%ZE
        
        b%Hx%ZI=sgg%Alloc(iHx)%ZI
        b%Hx%ZE=sgg%Alloc(iHx)%ZE
        b%Hy%ZI=sgg%Alloc(iHy)%ZI
        b%Hy%ZE=sgg%Alloc(iHy)%ZE
        b%Hz%ZI=sgg%Alloc(iHz)%ZI
        b%Hz%ZE=sgg%Alloc(iHz)%ZE
        
        b%sggMiEx%XI=sgg%Alloc(iEx)%XI
        b%sggMiEx%XE=sgg%Alloc(iEx)%XE
        b%sggMiEy%XI=sgg%Alloc(iEy)%XI
        b%sggMiEy%XE=sgg%Alloc(iEy)%XE
        b%sggMiEz%XI=sgg%Alloc(iEz)%XI
        b%sggMiEz%XE=sgg%Alloc(iEz)%XE
        
        b%sggMiHx%XI=sgg%Alloc(iHx)%XI
        b%sggMiHx%XE=sgg%Alloc(iHx)%XE
        b%sggMiHy%XI=sgg%Alloc(iHy)%XI
        b%sggMiHy%XE=sgg%Alloc(iHy)%XE
        b%sggMiHz%XI=sgg%Alloc(iHz)%XI
        b%sggMiHz%XE=sgg%Alloc(iHz)%XE
        
        b%sggMiEx%YI=sgg%Alloc(iEx)%YI
        b%sggMiEx%YE=sgg%Alloc(iEx)%YE
        b%sggMiEy%YI=sgg%Alloc(iEy)%YI
        b%sggMiEy%YE=sgg%Alloc(iEy)%YE
        b%sggMiEz%YI=sgg%Alloc(iEz)%YI
        b%sggMiEz%YE=sgg%Alloc(iEz)%YE
        
        b%sggMiHx%YI=sgg%Alloc(iHx)%YI
        b%sggMiHx%YE=sgg%Alloc(iHx)%YE
        b%sggMiHy%YI=sgg%Alloc(iHy)%YI
        b%sggMiHy%YE=sgg%Alloc(iHy)%YE
        b%sggMiHz%YI=sgg%Alloc(iHz)%YI
        b%sggMiHz%YE=sgg%Alloc(iHz)%YE
        
        b%sggMiEx%ZI=sgg%Alloc(iEx)%ZI
        b%sggMiEx%ZE=sgg%Alloc(iEx)%ZE
        b%sggMiEy%ZI=sgg%Alloc(iEy)%ZI
        b%sggMiEy%ZE=sgg%Alloc(iEy)%ZE
        b%sggMiEz%ZI=sgg%Alloc(iEz)%ZI
        b%sggMiEz%ZE=sgg%Alloc(iEz)%ZE
        
        b%sggMiHx%ZI=sgg%Alloc(iHx)%ZI
        b%sggMiHx%ZE=sgg%Alloc(iHx)%ZE
        b%sggMiHy%ZI=sgg%Alloc(iHy)%ZI
        b%sggMiHy%ZE=sgg%Alloc(iHy)%ZE
        b%sggMiHz%ZI=sgg%Alloc(iHz)%ZI
        b%sggMiHz%ZE=sgg%Alloc(iHz)%ZE
        
        
        
        b%sweepEx%XI=sgg%Sweep(iEx)%XI
        b%sweepEx%XE=sgg%Sweep(iEx)%XE
        b%sweepEy%XI=sgg%Sweep(iEy)%XI
        b%sweepEy%XE=sgg%Sweep(iEy)%XE
        b%sweepEz%XI=sgg%Sweep(iEz)%XI
        b%sweepEz%XE=sgg%Sweep(iEz)%XE
        
        b%sweepHx%XI=sgg%Sweep(iHx)%XI
        b%sweepHx%XE=sgg%Sweep(iHx)%XE
        b%sweepHy%XI=sgg%Sweep(iHy)%XI
        b%sweepHy%XE=sgg%Sweep(iHy)%XE
        b%sweepHz%XI=sgg%Sweep(iHz)%XI
        b%sweepHz%XE=sgg%Sweep(iHz)%XE
        
        
        b%sweepEx%YI=sgg%Sweep(iEx)%YI
        b%sweepEx%YE=sgg%Sweep(iEx)%YE
        b%sweepEy%YI=sgg%Sweep(iEy)%YI
        b%sweepEy%YE=sgg%Sweep(iEy)%YE
        b%sweepEz%YI=sgg%Sweep(iEz)%YI
        b%sweepEz%YE=sgg%Sweep(iEz)%YE
        
        b%sweepHx%YI=sgg%Sweep(iHx)%YI
        b%sweepHx%YE=sgg%Sweep(iHx)%YE
        b%sweepHy%YI=sgg%Sweep(iHy)%YI
        b%sweepHy%YE=sgg%Sweep(iHy)%YE
        b%sweepHz%YI=sgg%Sweep(iHz)%YI
        b%sweepHz%YE=sgg%Sweep(iHz)%YE
        
        b%sweepEx%ZI=sgg%Sweep(iEx)%ZI
        b%sweepEx%ZE=sgg%Sweep(iEx)%ZE
        b%sweepEy%ZI=sgg%Sweep(iEy)%ZI
        b%sweepEy%ZE=sgg%Sweep(iEy)%ZE
        b%sweepEz%ZI=sgg%Sweep(iEz)%ZI
        b%sweepEz%ZE=sgg%Sweep(iEz)%ZE
        
        b%sweepHx%ZI=sgg%Sweep(iHx)%ZI
        b%sweepHx%ZE=sgg%Sweep(iHx)%ZE
        b%sweepHy%ZI=sgg%Sweep(iHy)%ZI
        b%sweepHy%ZE=sgg%Sweep(iHy)%ZE
        b%sweepHz%ZI=sgg%Sweep(iHz)%ZI
        b%sweepHz%ZE=sgg%Sweep(iHz)%ZE
        
        b%sweepSINPMLEx%XI=sgg%SINPMLSweep(iEx)%XI
        b%sweepSINPMLEy%XI=sgg%SINPMLSweep(iEy)%XI
        b%sweepSINPMLEz%XI=sgg%SINPMLSweep(iEz)%XI
        b%sweepSINPMLHx%XI=sgg%SINPMLSweep(iHx)%XI
        b%sweepSINPMLHy%XI=sgg%SINPMLSweep(iHy)%XI
        b%sweepSINPMLHz%XI=sgg%SINPMLSweep(iHz)%XI
        
        b%sweepSINPMLEx%XE=sgg%SINPMLSweep(iEx)%XE
        b%sweepSINPMLEy%XE=sgg%SINPMLSweep(iEy)%XE
        b%sweepSINPMLEz%XE=sgg%SINPMLSweep(iEz)%XE
        b%sweepSINPMLHx%XE=sgg%SINPMLSweep(iHx)%XE
        b%sweepSINPMLHy%XE=sgg%SINPMLSweep(iHy)%XE
        b%sweepSINPMLHz%XE=sgg%SINPMLSweep(iHz)%XE
        
        b%sweepSINPMLEx%YI=sgg%SINPMLSweep(iEx)%YI
        b%sweepSINPMLEy%YI=sgg%SINPMLSweep(iEy)%YI
        b%sweepSINPMLEz%YI=sgg%SINPMLSweep(iEz)%YI
        b%sweepSINPMLHx%YI=sgg%SINPMLSweep(iHx)%YI
        b%sweepSINPMLHy%YI=sgg%SINPMLSweep(iHy)%YI
        b%sweepSINPMLHz%YI=sgg%SINPMLSweep(iHz)%YI
        
        b%sweepSINPMLEx%YE=sgg%SINPMLSweep(iEx)%YE
        b%sweepSINPMLEy%YE=sgg%SINPMLSweep(iEy)%YE
        b%sweepSINPMLEz%YE=sgg%SINPMLSweep(iEz)%YE
        b%sweepSINPMLHx%YE=sgg%SINPMLSweep(iHx)%YE
        b%sweepSINPMLHy%YE=sgg%SINPMLSweep(iHy)%YE
        b%sweepSINPMLHz%YE=sgg%SINPMLSweep(iHz)%YE
        
        b%sweepSINPMLEx%ZI=sgg%SINPMLSweep(iEx)%ZI
        b%sweepSINPMLEy%ZI=sgg%SINPMLSweep(iEy)%ZI
        b%sweepSINPMLEz%ZI=sgg%SINPMLSweep(iEz)%ZI
        b%sweepSINPMLHx%ZI=sgg%SINPMLSweep(iHx)%ZI
        b%sweepSINPMLHy%ZI=sgg%SINPMLSweep(iHy)%ZI
        b%sweepSINPMLHz%ZI=sgg%SINPMLSweep(iHz)%ZI
        
        b%sweepSINPMLEx%ZE=sgg%SINPMLSweep(iEx)%ZE
        b%sweepSINPMLEy%ZE=sgg%SINPMLSweep(iEy)%ZE
        b%sweepSINPMLEz%ZE=sgg%SINPMLSweep(iEz)%ZE
        b%sweepSINPMLHx%ZE=sgg%SINPMLSweep(iHx)%ZE
        b%sweepSINPMLHy%ZE=sgg%SINPMLSweep(iHy)%ZE
        b%sweepSINPMLHz%ZE=sgg%SINPMLSweep(iHz)%ZE

        b%Ex%NX=b%Ex%XE-b%Ex%XI+1
        b%Ex%NY=b%Ex%YE-b%Ex%YI+1
        b%Ex%NZ=b%Ex%ZE-b%Ex%ZI+1

        b%Ey%NX=b%Ey%XE-b%Ey%XI+1
        b%Ey%NY=b%Ey%YE-b%Ey%YI+1
        b%Ey%NZ=b%Ey%ZE-b%Ey%ZI+1

        b%Ez%NX=b%Ez%XE-b%Ez%XI+1
        b%Ez%NY=b%Ez%YE-b%Ez%YI+1
        b%Ez%NZ=b%Ez%ZE-b%Ez%ZI+1

        b%Hx%NX=b%Hx%XE-b%Hx%XI+1
        b%Hx%NY=b%Hx%YE-b%Hx%YI+1
        b%Hx%NZ=b%Hx%ZE-b%Hx%ZI+1

        b%Hy%NX=b%Hy%XE-b%Hy%XI+1
        b%Hy%NY=b%Hy%YE-b%Hy%YI+1
        b%Hy%NZ=b%Hy%ZE-b%Hy%ZI+1

        b%Hz%NX=b%Hz%XE-b%Hz%XI+1
        b%Hz%NY=b%Hz%YE-b%Hz%YI+1
        b%Hz%NZ=b%Hz%ZE-b%Hz%ZI+1

        b%sweepEx%NX=b%sweepEx%XE-b%sweepEx%XI+1
        b%sweepEx%NY=b%sweepEx%YE-b%sweepEx%YI+1
        b%sweepEx%NZ=b%sweepEx%ZE-b%sweepEx%ZI+1

        b%sweepEy%NX=b%sweepEy%XE-b%sweepEy%XI+1
        b%sweepEy%NY=b%sweepEy%YE-b%sweepEy%YI+1
        b%sweepEy%NZ=b%sweepEy%ZE-b%sweepEy%ZI+1

        b%sweepEz%NX=b%sweepEz%XE-b%sweepEz%XI+1
        b%sweepEz%NY=b%sweepEz%YE-b%sweepEz%YI+1
        b%sweepEz%NZ=b%sweepEz%ZE-b%sweepEz%ZI+1

        b%sweepHx%NX=b%sweepHx%XE-b%sweepHx%XI+1
        b%sweepHx%NY=b%sweepHx%YE-b%sweepHx%YI+1
        b%sweepHx%NZ=b%sweepHx%ZE-b%sweepHx%ZI+1

        b%sweepHy%NX=b%sweepHy%XE-b%sweepHy%XI+1
        b%sweepHy%NY=b%sweepHy%YE-b%sweepHy%YI+1
        b%sweepHy%NZ=b%sweepHy%ZE-b%sweepHy%ZI+1

        b%sweepHz%NX=b%sweepHz%XE-b%sweepHz%XI+1
        b%sweepHz%NY=b%sweepHz%YE-b%sweepHz%YI+1
        b%sweepHz%NZ=b%sweepHz%ZE-b%sweepHz%ZI+1

        b%sggMiEx%NX=b%sggMiEx%XE-b%sggMiEx%XI+1
        b%sggMiEx%NY=b%sggMiEx%YE-b%sggMiEx%YI+1
        b%sggMiEx%NZ=b%sggMiEx%ZE-b%sggMiEx%ZI+1
        b%sggMiEy%NX=b%sggMiEy%XE-b%sggMiEy%XI+1
        b%sggMiEy%NY=b%sggMiEy%YE-b%sggMiEy%YI+1
        b%sggMiEy%NZ=b%sggMiEy%ZE-b%sggMiEy%ZI+1
        b%sggMiEz%NX=b%sggMiEz%XE-b%sggMiEz%XI+1
        b%sggMiEz%NY=b%sggMiEz%YE-b%sggMiEz%YI+1
        b%sggMiEz%NZ=b%sggMiEz%ZE-b%sggMiEz%ZI+1

        b%sggMiHx%NX=b%sggMiHx%XE-b%sggMiHx%XI+1
        b%sggMiHx%NY=b%sggMiHx%YE-b%sggMiHx%YI+1
        b%sggMiHx%NZ=b%sggMiHx%ZE-b%sggMiHx%ZI+1
        b%sggMiHy%NX=b%sggMiHy%XE-b%sggMiHy%XI+1
        b%sggMiHy%NY=b%sggMiHy%YE-b%sggMiHy%YI+1
        b%sggMiHy%NZ=b%sggMiHy%ZE-b%sggMiHy%ZI+1
        b%sggMiHz%NX=b%sggMiHz%XE-b%sggMiHz%XI+1
        b%sggMiHz%NY=b%sggMiHz%YE-b%sggMiHz%YI+1
        b%sggMiHz%NZ=b%sggMiHz%ZE-b%sggMiHz%ZI+1

        b%dxe%NX=b%dxe%XE-b%dxe%XI+1
        b%dye%NY=b%dye%YE-b%dye%YI+1
        b%dze%NZ=b%dze%ZE-b%dze%ZI+1

        b%dxh%NX=b%dxh%XE-b%dxh%XI+1
        b%dyh%NY=b%dyh%YE-b%dyh%YI+1
        b%dzh%NZ=b%dzh%ZE-b%dzh%ZI+1

    end subroutine

end module solver_mod