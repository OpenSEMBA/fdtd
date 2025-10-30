MODULE vtk
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: vtkWriter

    !---------------------------------
    ! Define interface for updater
    !---------------------------------
    ABSTRACT INTERFACE
        SUBROUTINE updater_interface(this)
            IMPORT :: vtkWriter
            CLASS(vtkWriter), INTENT(INOUT) :: this
        END SUBROUTINE updater_interface
    END INTERFACE

    !---------------------------------
    ! VTK Writer type
    !---------------------------------
    TYPE, PUBLIC :: vtkWriter
        CHARACTER(LEN=100) :: fileName = ''
        LOGICAL :: isOpen = .FALSE.
        INTEGER :: timesWritten = 0
        PROCEDURE(updater_interface), POINTER :: updater => NULL()
    CONTAINS
        PROCEDURE :: create
        PROCEDURE :: append => appendData
    END TYPE vtkWriter

CONTAINS

    !---------------------------------
    ! Create and initialize VTK file
    !---------------------------------
    SUBROUTINE create(this, fname, dataFlag)
        CLASS(vtkWriter), INTENT(INOUT) :: this
        CHARACTER(LEN=*), INTENT(IN) :: fname
        INTEGER, INTENT(IN) :: dataFlag
        INTEGER :: unitNum, ios

        ! Store file name
        this%fileName = fname

        ! Assign updater procedure based on flag
        SELECT CASE(dataFlag)
        CASE(1)
            this%updater => appendCurrent
        CASE(2)
            this%updater => appendEField
        CASE(3)
            this%updater => appendHField
        CASE DEFAULT
            PRINT *, 'Unknown dataFlag. No updater assigned.'
            this%updater => NULL()
        END SELECT

        ! Initialize VTK file
        OPEN(NEWUNIT=unitNum, FILE=TRIM(this%fileName), STATUS='REPLACE', FORM='FORMATTED', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, 'Error opening file: ', TRIM(this%fileName)
            this%isOpen = .FALSE.
            RETURN
        END IF

        WRITE(unitNum,'(A)') '# vtk DataFile Version 3.0'
        if (what==mapvtk) THEN
         write(unitNum,'(a)') 'PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, CONF=5/6, OTHER=-1 (ADD +0.5 for borders)'
      else                                      
         if (TimeDomain) then
            write(unitNum,'(a,e21.12e3)') 'Time= ',time
         else
            write(unitNum,'(a,e21.12e3)') 'Frequency= ',time   
         endif
      endif
        WRITE(unitNum,'(A)') 'VTK file initialized by vtkWriter'
        WRITE(unitNum,'(A)') 'ASCII'
        WRITE(unitNum,'(A)') ' '
        WRITE(unitNum,'(A)') 'DATASET UNSTRUCTURED_GRID'

        CLOSE(unitNum)
        this%isOpen = .TRUE.
        this%timesWritten = 0
    END SUBROUTINE create

    !---------------------------------
    ! Append data (calls assigned updater)
    !---------------------------------
    SUBROUTINE appendData(this)
        CLASS(vtkWriter), INTENT(INOUT) :: this

        IF (.NOT. this%isOpen) THEN
            PRINT *, 'VTK file not initialized.'
            RETURN
        END IF

        IF (ASSOCIATED(this%updater)) THEN
            CALL this%updater()
            this%timesWritten = this%timesWritten + 1
        ELSE
            PRINT *, 'No updater assigned.'
        END IF
    END SUBROUTINE appendData

    !---------------------------------
    ! Example updater: Current
    !---------------------------------
    SUBROUTINE appendCurrent(this)
        CLASS(vtkWriter), INTENT(INOUT) :: this
        ! Open file for appending and write your current data
        ! Placeholder
        PRINT *, 'Appending CURRENT data to: ', TRIM(this%fileName)
    END SUBROUTINE appendCurrent

    !---------------------------------
    ! Example updater: EField
    !---------------------------------
    SUBROUTINE appendEField(this)
        CLASS(vtkWriter), INTENT(INOUT) :: this
        ! Placeholder for E-field data writing
        PRINT *, 'Appending E-FIELD data to: ', TRIM(this%fileName)
    END SUBROUTINE appendEField

    !---------------------------------
    ! Example updater: HField
    !---------------------------------
    SUBROUTINE appendHField(this)
        CLASS(vtkWriter), INTENT(INOUT) :: this
        ! Placeholder for H-field data writing
        PRINT *, 'Appending H-FIELD data to: ', TRIM(this%fileName)
    END SUBROUTINE appendHField

END MODULE vtk
