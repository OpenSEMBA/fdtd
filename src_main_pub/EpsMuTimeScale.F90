
 
!//===========================================================================//
!// HYSTORY & VERSION:                                                        //
!//    DATE beginning: 2018-12-18 17:17:09                                   //
!//===========================================================================//
module EpsMuTimeScale_m

USE FDETYPES
private :: new_input_, checkError_

type EpsMuTimeScale_input_parameters_t
    real (kind=RKind) :: tini, tend, alpha_max
    logical :: electric, magnetic
    logical :: are_there
contains
    procedure, pass, public  :: get_slope  => get_slope_
    procedure, pass, public  :: init0      => new_input_
    procedure, pass, public  :: checkError => checkError_
end type EpsMuTimeScale_input_parameters_t

contains

function get_slope_ (this) result (slope)
    class (EpsMuTimeScale_input_parameters_t) :: this
    real (kind=RKind) :: slope
    slope = (this%alpha_max-1.0_Rkind)/(this%tend-this%tini)
end function get_slope_

subroutine new_input_ (this)
    class (EpsMuTimeScale_input_parameters_t) :: this
    this%alpha_max = 1.0_Rkind
    this%tini        = 1e20_Rkind
    this%tend      = 1e20_Rkind
    this%electric = .false.
    this%magnetic  = .false.
    this%are_there = .false.
end subroutine new_input_

function checkError_ (this) result (res)
    class (EpsMuTimeScale_input_parameters_t) :: this
    integer :: res
    res = 0
    if (this%alpha_max<=0.0 .or. this%tini<0.0) then
        res = -1
    end if

end function checkError_

end module EpsMuTimeScale_m