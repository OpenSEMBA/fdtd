module test_nfde_rotate_m
    use nfde_rotate_m
    use NFDETypes
    use test_utils_m
    implicit none
    private
    public :: test_nfde_rotate

contains

    subroutine test_nfde_rotate()
        call test_rotate_generateSpaceSteps()
        call test_rotate_generateCurrent_Field_Sources()
        call test_rotate_generatePlaneWaves()
        call test_rotate_generateBoxSources()
        call test_rotate_generateFronteras()
        call test_rotate_generatePECs()
        call test_rotate_generatePMCs()
        call test_rotate_generateNONMetals()
        call test_rotate_generateANISOTROPICs()
        call test_rotate_generateThinWires()
        call test_rotate_generateSlantedWires()
        call test_rotate_generateThinSlots()
        call test_rotate_generateLossyThinSurface()
        call test_rotate_generateFDMs()
        call test_rotate_generateSONDAs()
        call test_rotate_generateMasSondas()
        call test_rotate_generateBloqueProbes()
        call test_rotate_generateVolumicProbes()
    end subroutine test_nfde_rotate


    

    

    

    
    

    

    

    subroutine test_rotate_generateThinWires()
        type(Parseador) :: this
        integer(kind=4) :: mpidir, i
        
        ! Test case 1: Rotate thin wires in mpidir=2 direction
        mpidir = 2
        allocate(this%THINWIRES)
        
        ! Initialize test data with some thin wires
        this%THINWIRES%nwires = 2
        
        ! Allocate and initialize wire data
        allocate(this%THINWIRES%x1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%y1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%z1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%x2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%y2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%z2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%dirx(this%THINWIRES%nwires))
        allocate(this%THINWIRES%diry(this%THINWIRES%nwires))
        allocate(this%THINWIRES%dirz(this%THINWIRES%nwires))
        
        ! Wire 1: Horizontal wire along x-axis
        this%THINWIRES%x1(1) = 1.0d0
        this%THINWIRES%y1(1) = 2.0d0
        this%THINWIRES%z1(1) = 3.0d0
        this%THINWIRES%x2(1) = 4.0d0
        this%THINWIRES%y2(1) = 2.0d0
        this%THINWIRES%z2(1) = 3.0d0
        this%THINWIRES%dirx(1) = 1.0d0
        this%THINWIRES%diry(1) = 0.0d0
        this%THINWIRES%dirz(1) = 0.0d0
        
        ! Wire 2: Vertical wire along z-axis
        this%THINWIRES%x1(2) = 5.0d0
        this%THINWIRES%y1(2) = 6.0d0
        this%THINWIRES%z1(2) = 7.0d0
        this%THINWIRES%x2(2) = 5.0d0
        this%THINWIRES%y2(2) = 6.0d0
        this%THINWIRES%z2(2) = 10.0d0
        this%THINWIRES%dirx(2) = 0.0d0
        this%THINWIRES%diry(2) = 0.0d0
        this%THINWIRES%dirz(2) = 1.0d0
        
        ! Call the rotation routine
        call rotate_generateThinWires(this, mpidir)
        
        ! Verify the results for wire 1 (should rotate around y-axis)
        call assert_equal(this%THINWIRES%x1(1), -3.0d0, "rotate_generateThinWires: x1(1) should be rotated")
        call assert_equal(this%THINWIRES%y1(1), 2.0d0, "rotate_generateThinWires: y1(1) should remain unchanged")
        call assert_equal(this%THINWIRES%z1(1), -1.0d0, "rotate_generateThinWires: z1(1) should be rotated")
        call assert_equal(this%THINWIRES%x2(1), -3.0d0, "rotate_generateThinWires: x2(1) should be rotated")
        call assert_equal(this%THINWIRES%y2(1), 2.0d0, "rotate_generateThinWires: y2(1) should remain unchanged")
        call assert_equal(this%THINWIRES%z2(1), -4.0d0, "rotate_generateThinWires: z2(1) should be rotated")
        call assert_equal(this%THINWIRES%dirx(1), -1.0d0, "rotate_generateThinWires: dirx(1) should be rotated")
        call assert_equal(this%THINWIRES%diry(1), 0.0d0, "rotate_generateThinWires: diry(1) should remain unchanged")
        call assert_equal(this%THINWIRES%dirz(1), 0.0d0, "rotate_generateThinWires: dirz(1) should be rotated")
        
        ! Verify the results for wire 2 (should rotate around y-axis)
        call assert_equal(this%THINWIRES%x1(2), -7.0d0, "rotate_generateThinWires: x1(2) should be rotated")
        call assert_equal(this%THINWIRES%y1(2), 6.0d0, "rotate_generateThinWires: y1(2) should remain unchanged")
        call assert_equal(this%THINWIRES%z1(2), -5.0d0, "rotate_generateThinWires: z1(2) should be rotated")
        call assert_equal(this%THINWIRES%x2(2), -10.0d0, "rotate_generateThinWires: x2(2) should be rotated")
        call assert_equal(this%THINWIRES%y2(2), 6.0d0, "rotate_generateThinWires: y2(2) should remain unchanged")
        call assert_equal(this%THINWIRES%z2(2), -5.0d0, "rotate_generateThinWires: z2(2) should be rotated")
        call assert_equal(this%THINWIRES%dirx(2), 0.0d0, "rotate_generateThinWires: dirx(2) should be rotated")
        call assert_equal(this%THINWIRES%diry(2), 0.0d0, "rotate_generateThinWires: diry(2) should remain unchanged")
        call assert_equal(this%THINWIRES%dirz(2), -1.0d0, "rotate_generateThinWires: dirz(2) should be rotated")
        
        deallocate(this%THINWIRES%x1, this%THINWIRES%y1, this%THINWIRES%z1)
        deallocate(this%THINWIRES%x2, this%THINWIRES%y2, this%THINWIRES%z2)
        deallocate(this%THINWIRES%dirx, this%THINWIRES%diry, this%THINWIRES%dirz)
        deallocate(this%THINWIRES)
        
        ! Test case 2: Rotate thin wires in mpidir=1 direction
        mpidir = 1
        allocate(this%THINWIRES)
        
        ! Initialize test data with some thin wires
        this%THINWIRES%nwires = 2
        
        ! Allocate and initialize wire data
        allocate(this%THINWIRES%x1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%y1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%z1(this%THINWIRES%nwires))
        allocate(this%THINWIRES%x2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%y2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%z2(this%THINWIRES%nwires))
        allocate(this%THINWIRES%dirx(this%THINWIRES%nwires))
        allocate(this%THINWIRES%diry(this%THINWIRES%nwires))
        allocate(this%THINWIRES%dirz(this%THINWIRES%nwires))
        
        ! Wire 1: Horizontal wire along x-axis
        this%THINWIRES%x1(1) = 1.0d0
        this%THINWIRES%y1(1) = 2.0d0
        this%THINWIRES%z1(1) = 3.0d0
        this%THINWIRES%x2(1) = 4.0d0
        this%THINWIRES%y2(1) = 2.0d0
        this%THINWIRES%z2(1) = 3.0d0
        this%THINWIRES%dirx(1) = 1.0d0
        this%THINWIRES%diry(1) = 0.0d0
        this%THINWIRES%dirz(1) = 0.0d0
        
        ! Wire 2: Vertical wire along z-axis
        this%THINWIRES%x1(2) = 5.0d0
        this%THINWIRES%y1(2) = 6.0d0
        this%THINWIRES%z1(2) = 7.0d0
        this%THINWIRES%x2(2) = 5.0d0
        this%THINWIRES%y2(2) = 6.0d0
        this%THINWIRES%z2(2) = 10.0d0
        this%THINWIRES%dirx(2) = 0.0d0
        this%THINWIRES%diry(2) = 0.0d0
        this%THINWIRES%dirz(2) = 1.0d0
        
        ! Call the rotation routine
        call rotate_generateThinWires(this, mpidir)
        
        ! Verify the results for wire 1 (should rotate around x-axis)
        call assert_equal(this%THINWIRES%x1(1), 1.0d0, "rotate_generateThinWires: x1(1) should remain unchanged")
        call assert_equal(this%THINWIRES%y1(1), -3.0d0, "rotate_generateThinWires: y1(1) should be rotated")
        call assert_equal(this%THINWIRES%z1(1), 2.0d0, "rotate_generateThinWires: z1(1) should be rotated")
        call assert_equal(this%THINWIRES%x2(1), 4.0d0, "rotate_generateThinWires: x2(1) should remain unchanged")
        call assert_equal(this%THINWIRES%y2(1), -3.0d0, "rotate_generateThinWires: y2(1) should be rotated")
        call assert_equal(this%THINWIRES%z2(1), 2.0d0, "rotate_generateThinWires: z2(1) should be rotated")
        call assert_equal(this%THINWIRES%dirx(1), 1.0d0, "rotate_generateThinWires: dirx(1) should remain unchanged")
        call assert_equal(this%THINWIRES%diry(1), 0.0d0, "rotate_generateThinWires: diry(1) should be rotated")
        call assert_equal(this%THINWIRES%dirz(1), 0.0d0, "rotate_generateThinWires: dirz(1) should be rotated")
        
        ! Verify the results for wire 2 (should rotate around x-axis)
        call assert_equal(this%THINWIRES%x1(2), 5.0d0, "rotate_generateThinWires: x1(2) should remain unchanged")
        call assert_equal(this%THINWIRES%y1(2), -7.0d0, "rotate_generateThinWires: y1(2) should be rotated")
        call assert_equal(this%THINWIRES%z1(2), 6.0d0, "rotate_generateThinWires: z1(2) should be rotated")
        call assert_equal(this%THINWIRES%x2(2), 5.0d0, "rotate_generateThinWires: x2(2) should remain unchanged")
        call assert_equal(this%THINWIRES%y2(2), -10.0d0, "rotate_generateThinWires: y2(2) should be rotated")
        call assert_equal(this%THINWIRES%z2(2), 6.0d0, "rotate_generateThinWires: z2(2) should be rotated")
        call assert_equal(this%THINWIRES%dirx(2), 0.0d0, "rotate_generateThinWires: dirx(2) should remain unchanged")
        call assert_equal(this%THINWIRES%diry(2), 0.0d0, "rotate_generateThinWires: diry(2) should be rotated")
        call assert_equal(this%THINWIRES%dirz(2), -1.0d0, "rotate_generateThinWires: dirz(2) should be rotated")
        
        deallocate(this%THINWIRES%x1, this%THINWIRES%y1, this%THINWIRES%z1)
        deallocate(this%THINWIRES%x2, this%THINWIRES%y2, this%THINWIRES%z2)
        deallocate(this%THINWIRES%dirx, this%THINWIRES%diry, this%THINWIRES%dirz)
        deallocate(this%THINWIRES)
        
        ! Test case 3: Verify behavior with no wires
        mpidir = 2
        allocate(this%THINWIRES)
        this%THINWIRES%nwires = 0
        
        ! Call the rotation routine - should handle empty case gracefully
        call rotate_generateThinWires(this, mpidir)
        
        ! Verify that nwires remains 0
        call assert_equal(this%THINWIRES%nwires, 0, "rotate_generateThinWires: nwires should remain 0")
        
        deallocate(this%THINWIRES)
    end subroutine test_rotate_generateThinWires

    subroutine test_rotate_generateSlantedWires()
        type(Parseador) :: this
        integer(kind=4) :: mpidir, i
        
        ! Test case 1: Rotate slanted wires in mpidir=2 direction
        mpidir = 2
        allocate(this%SLANTEDWIRES)
        
        ! Initialize test data with some slanted wires
        this%SLANTEDWIRES%nwires = 2
        
        ! Allocate and initialize wire data
        allocate(this%SLANTEDWIRES%x1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%y1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%z1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%x2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%y2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%z2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%dirx(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%diry(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%dirz(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%radio(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angx(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angy(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angz(this%SLANTEDWIRES%nwires))
        
        ! Wire 1: Slanted wire with 45-degree angle in x-z plane
        this%SLANTEDWIRES%x1(1) = 1.0d0
        this%SLANTEDWIRES%y1(1) = 2.0d0
        this%SLANTEDWIRES%z1(1) = 3.0d0
        this%SLANTEDWIRES%x2(1) = 4.0d0
        this%SLANTEDWIRES%y2(1) = 2.0d0
        this%SLANTEDWIRES%z2(1) = 6.0d0
        this%SLANTEDWIRES%dirx(1) = 0.7071067811865475d0  ! cos(45°)
        this%SLANTEDWIRES%diry(1) = 0.0d0
        this%SLANTEDWIRES%dirz(1) = 0.7071067811865475d0  ! sin(45°)
        this%SLANTEDWIRES%radio(1) = 0.1d0
        this%SLANTEDWIRES%angx(1) = 45.0d0
        this%SLANTEDWIRES%angy(1) = 0.0d0
        this%SLANTEDWIRES%angz(1) = 45.0d0
        
        ! Wire 2: Slanted wire with 30-degree angle in y-z plane
        this%SLANTEDWIRES%x1(2) = 5.0d0
        this%SLANTEDWIRES%y1(2) = 6.0d0
        this%SLANTEDWIRES%z1(2) = 7.0d0
        this%SLANTEDWIRES%x2(2) = 5.0d0
        this%SLANTEDWIRES%y2(2) = 8.0d0
        this%SLANTEDWIRES%z2(2) = 10.0d0
        this%SLANTEDWIRES%dirx(2) = 0.0d0
        this%SLANTEDWIRES%diry(2) = 0.5d0  ! cos(60°)
        this%SLANTEDWIRES%dirz(2) = 0.8660254037844386d0  ! sin(60°)
        this%SLANTEDWIRES%radio(2) = 0.2d0
        this%SLANTEDWIRES%angx(2) = 0.0d0
        this%SLANTEDWIRES%angy(2) = 30.0d0
        this%SLANTEDWIRES%angz(2) = 60.0d0
        
        ! Call the rotation routine
        call rotate_generateSlantedWires(this, mpidir)
        
        ! Verify the results for wire 1 (should rotate around y-axis)
        call assert_equal(this%SLANTEDWIRES%x1(1), -3.0d0, "rotate_generateSlantedWires: x1(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%y1(1), 2.0d0, "rotate_generateSlantedWires: y1(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%z1(1), -1.0d0, "rotate_generateSlantedWires: z1(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%x2(1), -6.0d0, "rotate_generateSlantedWires: x2(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%y2(1), 2.0d0, "rotate_generateSlantedWires: y2(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%z2(1), -4.0d0, "rotate_generateSlantedWires: z2(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirx(1), -0.7071067811865475d0, "rotate_generateSlantedWires: dirx(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%diry(1), 0.0d0, "rotate_generateSlantedWires: diry(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%dirz(1), -0.7071067811865475d0, "rotate_generateSlantedWires: dirz(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%radio(1), 0.1d0, "rotate_generateSlantedWires: radio(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angx(1), -45.0d0, "rotate_generateSlantedWires: angx(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%angy(1), 0.0d0, "rotate_generateSlantedWires: angy(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angz(1), -45.0d0, "rotate_generateSlantedWires: angz(1) should be rotated")
        
        ! Verify the results for wire 2 (should rotate around y-axis)
        call assert_equal(this%SLANTEDWIRES%x1(2), -7.0d0, "rotate_generateSlantedWires: x1(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%y1(2), 6.0d0, "rotate_generateSlantedWires: y1(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%z1(2), -5.0d0, "rotate_generateSlantedWires: z1(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%x2(2), -10.0d0, "rotate_generateSlantedWires: x2(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%y2(2), 8.0d0, "rotate_generateSlantedWires: y2(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%z2(2), -5.0d0, "rotate_generateSlantedWires: z2(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirx(2), -0.8660254037844386d0, "rotate_generateSlantedWires: dirx(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%diry(2), 0.5d0, "rotate_generateSlantedWires: diry(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%dirz(2), -0.0d0, "rotate_generateSlantedWires: dirz(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%radio(2), 0.2d0, "rotate_generateSlantedWires: radio(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angx(2), -60.0d0, "rotate_generateSlantedWires: angx(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%angy(2), 30.0d0, "rotate_generateSlantedWires: angy(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angz(2), -90.0d0, "rotate_generateSlantedWires: angz(2) should be rotated")
        
        deallocate(this%SLANTEDWIRES%x1, this%SLANTEDWIRES%y1, this%SLANTEDWIRES%z1)
        deallocate(this%SLANTEDWIRES%x2, this%SLANTEDWIRES%y2, this%SLANTEDWIRES%z2)
        deallocate(this%SLANTEDWIRES%dirx, this%SLANTEDWIRES%diry, this%SLANTEDWIRES%dirz)
        deallocate(this%SLANTEDWIRES%radio)
        deallocate(this%SLANTEDWIRES%angx, this%SLANTEDWIRES%angy, this%SLANTEDWIRES%angz)
        deallocate(this%SLANTEDWIRES)
        
        ! Test case 2: Rotate slanted wires in mpidir=1 direction
        mpidir = 1
        allocate(this%SLANTEDWIRES)
        
        ! Initialize test data with some slanted wires
        this%SLANTEDWIRES%nwires = 2
        
        ! Allocate and initialize wire data (same as test case 1)
        allocate(this%SLANTEDWIRES%x1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%y1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%z1(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%x2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%y2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%z2(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%dirx(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%diry(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%dirz(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%radio(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angx(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angy(this%SLANTEDWIRES%nwires))
        allocate(this%SLANTEDWIRES%angz(this%SLANTEDWIRES%nwires))
        
        ! Initialize with same values as test case 1
        ! Wire 1: Slanted wire with 45-degree angle in x-z plane
        this%SLANTEDWIRES%x1(1) = 1.0d0
        this%SLANTEDWIRES%y1(1) = 2.0d0
        this%SLANTEDWIRES%z1(1) = 3.0d0
        this%SLANTEDWIRES%x2(1) = 4.0d0
        this%SLANTEDWIRES%y2(1) = 2.0d0
        this%SLANTEDWIRES%z2(1) = 6.0d0
        this%SLANTEDWIRES%dirx(1) = 0.7071067811865475d0
        this%SLANTEDWIRES%diry(1) = 0.0d0
        this%SLANTEDWIRES%dirz(1) = 0.7071067811865475d0
        this%SLANTEDWIRES%radio(1) = 0.1d0
        this%SLANTEDWIRES%angx(1) = 45.0d0
        this%SLANTEDWIRES%angy(1) = 0.0d0
        this%SLANTEDWIRES%angz(1) = 45.0d0
        
        ! Wire 2: Slanted wire with 30-degree angle in y-z plane
        this%SLANTEDWIRES%x1(2) = 5.0d0
        this%SLANTEDWIRES%y1(2) = 6.0d0
        this%SLANTEDWIRES%z1(2) = 7.0d0
        this%SLANTEDWIRES%x2(2) = 5.0d0
        this%SLANTEDWIRES%y2(2) = 8.0d0
        this%SLANTEDWIRES%z2(2) = 10.0d0
        this%SLANTEDWIRES%dirx(2) = 0.0d0
        this%SLANTEDWIRES%diry(2) = 0.5d0
        this%SLANTEDWIRES%dirz(2) = 0.8660254037844386d0
        this%SLANTEDWIRES%radio(2) = 0.2d0
        this%SLANTEDWIRES%angx(2) = 0.0d0
        this%SLANTEDWIRES%angy(2) = 30.0d0
        this%SLANTEDWIRES%angz(2) = 60.0d0
        
        ! Call the rotation routine
        call rotate_generateSlantedWires(this, mpidir)
        
        ! Verify the results for wire 1 (should rotate around x-axis)
        call assert_equal(this%SLANTEDWIRES%x1(1), 1.0d0, "rotate_generateSlantedWires: x1(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%y1(1), -3.0d0, "rotate_generateSlantedWires: y1(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%z1(1), 2.0d0, "rotate_generateSlantedWires: z1(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%x2(1), 4.0d0, "rotate_generateSlantedWires: x2(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%y2(1), -6.0d0, "rotate_generateSlantedWires: y2(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%z2(1), 2.0d0, "rotate_generateSlantedWires: z2(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirx(1), 0.7071067811865475d0, "rotate_generateSlantedWires: dirx(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%diry(1), -0.7071067811865475d0, "rotate_generateSlantedWires: diry(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirz(1), 0.0d0, "rotate_generateSlantedWires: dirz(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%radio(1), 0.1d0, "rotate_generateSlantedWires: radio(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angx(1), 45.0d0, "rotate_generateSlantedWires: angx(1) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angy(1), -45.0d0, "rotate_generateSlantedWires: angy(1) should be rotated")
        call assert_equal(this%SLANTEDWIRES%angz(1), 0.0d0, "rotate_generateSlantedWires: angz(1) should be rotated")
        
        ! Verify the results for wire 2 (should rotate around x-axis)
        call assert_equal(this%SLANTEDWIRES%x1(2), 5.0d0, "rotate_generateSlantedWires: x1(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%y1(2), -7.0d0, "rotate_generateSlantedWires: y1(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%z1(2), 6.0d0, "rotate_generateSlantedWires: z1(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%x2(2), 5.0d0, "rotate_generateSlantedWires: x2(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%y2(2), -10.0d0, "rotate_generateSlantedWires: y2(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%z2(2), 8.0d0, "rotate_generateSlantedWires: z2(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirx(2), 0.0d0, "rotate_generateSlantedWires: dirx(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%diry(2), -0.8660254037844386d0, "rotate_generateSlantedWires: diry(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%dirz(2), 0.5d0, "rotate_generateSlantedWires: dirz(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%radio(2), 0.2d0, "rotate_generateSlantedWires: radio(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angx(2), 0.0d0, "rotate_generateSlantedWires: angx(2) should remain unchanged")
        call assert_equal(this%SLANTEDWIRES%angy(2), -60.0d0, "rotate_generateSlantedWires: angy(2) should be rotated")
        call assert_equal(this%SLANTEDWIRES%angz(2), 30.0d0, "rotate_generateSlantedWires: angz(2) should be rotated")
        
        deallocate(this%SLANTEDWIRES%x1, this%SLANTEDWIRES%y1, this%SLANTEDWIRES%z1)
        deallocate(this%SLANTEDWIRES%x2, this%SLANTEDWIRES%y2, this%SLANTEDWIRES%z2)
        deallocate(this%SLANTEDWIRES%dirx, this%SLANTEDWIRES%diry, this%SLANTEDWIRES%dirz)
        deallocate(this%SLANTEDWIRES%radio)
        deallocate(this%SLANTEDWIRES%angx, this%SLANTEDWIRES%angy, this%SLANTEDWIRES%angz)
        deallocate(this%SLANTEDWIRES)
        
        ! Test case 3: Verify behavior with no wires
        mpidir = 2
        allocate(this%SLANTEDWIRES)
        this%SLANTEDWIRES%nwires = 0
        
        ! Call the rotation routine - should handle empty case gracefully
        call rotate_generateSlantedWires(this, mpidir)
        
        ! Verify that nwires remains 0
        call assert_equal(this%SLANTEDWIRES%nwires, 0, "rotate_generateSlantedWires: nwires should remain 0")
        
        deallocate(this%SLANTEDWIRES)
    end subroutine test_rotate_generateSlantedWires

    subroutine test_rotate_generateThinSlots()
        type(Parseador) :: this
        integer(kind=4) :: mpidir, i
        
        ! Test case 1: Rotate thin slots in mpidir=2 direction
        mpidir = 2
        allocate(this%THINSLOTS)
        
        ! Initialize test data with some thin slots
        this%THINSLOTS%nslots = 2
        
        ! Allocate and initialize slot data
        allocate(this%THINSLOTS%x1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%y1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%z1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%x2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%y2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%z2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%dirx(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%diry(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%dirz(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%ancho(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angx(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angy(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angz(this%THINSLOTS%nslots))
        
        ! Slot 1: Horizontal slot along x-axis
        this%THINSLOTS%x1(1) = 1.0d0
        this%THINSLOTS%y1(1) = 2.0d0
        this%THINSLOTS%z1(1) = 3.0d0
        this%THINSLOTS%x2(1) = 4.0d0
        this%THINSLOTS%y2(1) = 2.0d0
        this%THINSLOTS%z2(1) = 3.0d0
        this%THINSLOTS%dirx(1) = 1.0d0
        this%THINSLOTS%diry(1) = 0.0d0
        this%THINSLOTS%dirz(1) = 0.0d0
        this%THINSLOTS%ancho(1) = 0.1d0
        this%THINSLOTS%angx(1) = 0.0d0
        this%THINSLOTS%angy(1) = 0.0d0
        this%THINSLOTS%angz(1) = 0.0d0
        
        ! Slot 2: Vertical slot along z-axis
        this%THINSLOTS%x1(2) = 5.0d0
        this%THINSLOTS%y1(2) = 6.0d0
        this%THINSLOTS%z1(2) = 7.0d0
        this%THINSLOTS%x2(2) = 5.0d0
        this%THINSLOTS%y2(2) = 6.0d0
        this%THINSLOTS%z2(2) = 10.0d0
        this%THINSLOTS%dirx(2) = 0.0d0
        this%THINSLOTS%diry(2) = 0.0d0
        this%THINSLOTS%dirz(2) = 1.0d0
        this%THINSLOTS%ancho(2) = 0.2d0
        this%THINSLOTS%angx(2) = 90.0d0
        this%THINSLOTS%angy(2) = 0.0d0
        this%THINSLOTS%angz(2) = 0.0d0
        
        ! Call the rotation routine
        call rotate_generateThinSlots(this, mpidir)
        
        ! Verify the results for slot 1 (should rotate around y-axis)
        call assert_equal(this%THINSLOTS%x1(1), -3.0d0, "rotate_generateThinSlots: x1(1) should be rotated")
        call assert_equal(this%THINSLOTS%y1(1), 2.0d0, "rotate_generateThinSlots: y1(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%z1(1), -1.0d0, "rotate_generateThinSlots: z1(1) should be rotated")
        call assert_equal(this%THINSLOTS%x2(1), -3.0d0, "rotate_generateThinSlots: x2(1) should be rotated")
        call assert_equal(this%THINSLOTS%y2(1), 2.0d0, "rotate_generateThinSlots: y2(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%z2(1), -4.0d0, "rotate_generateThinSlots: z2(1) should be rotated")
        call assert_equal(this%THINSLOTS%dirx(1), -1.0d0, "rotate_generateThinSlots: dirx(1) should be rotated")
        call assert_equal(this%THINSLOTS%diry(1), 0.0d0, "rotate_generateThinSlots: diry(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%dirz(1), 0.0d0, "rotate_generateThinSlots: dirz(1) should be rotated")
        call assert_equal(this%THINSLOTS%ancho(1), 0.1d0, "rotate_generateThinSlots: ancho(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%angx(1), 180.0d0, "rotate_generateThinSlots: angx(1) should be rotated")
        call assert_equal(this%THINSLOTS%angy(1), 0.0d0, "rotate_generateThinSlots: angy(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%angz(1), 0.0d0, "rotate_generateThinSlots: angz(1) should be rotated")
        
        ! Verify the results for slot 2 (should rotate around y-axis)
        call assert_equal(this%THINSLOTS%x1(2), -7.0d0, "rotate_generateThinSlots: x1(2) should be rotated")
        call assert_equal(this%THINSLOTS%y1(2), 6.0d0, "rotate_generateThinSlots: y1(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%z1(2), -5.0d0, "rotate_generateThinSlots: z1(2) should be rotated")
        call assert_equal(this%THINSLOTS%x2(2), -10.0d0, "rotate_generateThinSlots: x2(2) should be rotated")
        call assert_equal(this%THINSLOTS%y2(2), 6.0d0, "rotate_generateThinSlots: y2(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%z2(2), -5.0d0, "rotate_generateThinSlots: z2(2) should be rotated")
        call assert_equal(this%THINSLOTS%dirx(2), 0.0d0, "rotate_generateThinSlots: dirx(2) should be rotated")
        call assert_equal(this%THINSLOTS%diry(2), 0.0d0, "rotate_generateThinSlots: diry(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%dirz(2), -1.0d0, "rotate_generateThinSlots: dirz(2) should be rotated")
        call assert_equal(this%THINSLOTS%ancho(2), 0.2d0, "rotate_generateThinSlots: ancho(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%angx(2), 90.0d0, "rotate_generateThinSlots: angx(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%angy(2), 0.0d0, "rotate_generateThinSlots: angy(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%angz(2), 180.0d0, "rotate_generateThinSlots: angz(2) should be rotated")
        
        deallocate(this%THINSLOTS%x1, this%THINSLOTS%y1, this%THINSLOTS%z1)
        deallocate(this%THINSLOTS%x2, this%THINSLOTS%y2, this%THINSLOTS%z2)
        deallocate(this%THINSLOTS%dirx, this%THINSLOTS%diry, this%THINSLOTS%dirz)
        deallocate(this%THINSLOTS%ancho)
        deallocate(this%THINSLOTS%angx, this%THINSLOTS%angy, this%THINSLOTS%angz)
        deallocate(this%THINSLOTS)
        
        ! Test case 2: Rotate thin slots in mpidir=1 direction
        mpidir = 1
        allocate(this%THINSLOTS)
        
        ! Initialize test data with some thin slots
        this%THINSLOTS%nslots = 2
        
        ! Allocate and initialize slot data (same as test case 1)
        allocate(this%THINSLOTS%x1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%y1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%z1(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%x2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%y2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%z2(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%dirx(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%diry(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%dirz(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%ancho(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angx(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angy(this%THINSLOTS%nslots))
        allocate(this%THINSLOTS%angz(this%THINSLOTS%nslots))
        
        ! Initialize with same values as test case 1
        ! Slot 1: Horizontal slot along x-axis
        this%THINSLOTS%x1(1) = 1.0d0
        this%THINSLOTS%y1(1) = 2.0d0
        this%THINSLOTS%z1(1) = 3.0d0
        this%THINSLOTS%x2(1) = 4.0d0
        this%THINSLOTS%y2(1) = 2.0d0
        this%THINSLOTS%z2(1) = 3.0d0
        this%THINSLOTS%dirx(1) = 1.0d0
        this%THINSLOTS%diry(1) = 0.0d0
        this%THINSLOTS%dirz(1) = 0.0d0
        this%THINSLOTS%ancho(1) = 0.1d0
        this%THINSLOTS%angx(1) = 0.0d0
        this%THINSLOTS%angy(1) = 0.0d0
        this%THINSLOTS%angz(1) = 0.0d0
        
        ! Slot 2: Vertical slot along z-axis
        this%THINSLOTS%x1(2) = 5.0d0
        this%THINSLOTS%y1(2) = 6.0d0
        this%THINSLOTS%z1(2) = 7.0d0
        this%THINSLOTS%x2(2) = 5.0d0
        this%THINSLOTS%y2(2) = 6.0d0
        this%THINSLOTS%z2(2) = 10.0d0
        this%THINSLOTS%dirx(2) = 0.0d0
        this%THINSLOTS%diry(2) = 0.0d0
        this%THINSLOTS%dirz(2) = 1.0d0
        this%THINSLOTS%ancho(2) = 0.2d0
        this%THINSLOTS%angx(2) = 90.0d0
        this%THINSLOTS%angy(2) = 0.0d0
        this%THINSLOTS%angz(2) = 0.0d0
        
        ! Call the rotation routine
        call rotate_generateThinSlots(this, mpidir)
        
        ! Verify the results for slot 1 (should rotate around x-axis)
        call assert_equal(this%THINSLOTS%x1(1), 1.0d0, "rotate_generateThinSlots: x1(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%y1(1), -3.0d0, "rotate_generateThinSlots: y1(1) should be rotated")
        call assert_equal(this%THINSLOTS%z1(1), 2.0d0, "rotate_generateThinSlots: z1(1) should be rotated")
        call assert_equal(this%THINSLOTS%x2(1), 4.0d0, "rotate_generateThinSlots: x2(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%y2(1), -3.0d0, "rotate_generateThinSlots: y2(1) should be rotated")
        call assert_equal(this%THINSLOTS%z2(1), 2.0d0, "rotate_generateThinSlots: z2(1) should be rotated")
        call assert_equal(this%THINSLOTS%dirx(1), 1.0d0, "rotate_generateThinSlots: dirx(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%diry(1), 0.0d0, "rotate_generateThinSlots: diry(1) should be rotated")
        call assert_equal(this%THINSLOTS%dirz(1), 0.0d0, "rotate_generateThinSlots: dirz(1) should be rotated")
        call assert_equal(this%THINSLOTS%ancho(1), 0.1d0, "rotate_generateThinSlots: ancho(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%angx(1), 0.0d0, "rotate_generateThinSlots: angx(1) should remain unchanged")
        call assert_equal(this%THINSLOTS%angy(1), 180.0d0, "rotate_generateThinSlots: angy(1) should be rotated")
        call assert_equal(this%THINSLOTS%angz(1), 0.0d0, "rotate_generateThinSlots: angz(1) should be rotated")
        
        ! Verify the results for slot 2 (should rotate around x-axis)
        call assert_equal(this%THINSLOTS%x1(2), 5.0d0, "rotate_generateThinSlots: x1(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%y1(2), -7.0d0, "rotate_generateThinSlots: y1(2) should be rotated")
        call assert_equal(this%THINSLOTS%z1(2), 6.0d0, "rotate_generateThinSlots: z1(2) should be rotated")
        call assert_equal(this%THINSLOTS%x2(2), 5.0d0, "rotate_generateThinSlots: x2(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%y2(2), -10.0d0, "rotate_generateThinSlots: y2(2) should be rotated")
        call assert_equal(this%THINSLOTS%z2(2), 6.0d0, "rotate_generateThinSlots: z2(2) should be rotated")
        call assert_equal(this%THINSLOTS%dirx(2), 0.0d0, "rotate_generateThinSlots: dirx(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%diry(2), 0.0d0, "rotate_generateThinSlots: diry(2) should be rotated")
        call assert_equal(this%THINSLOTS%dirz(2), -1.0d0, "rotate_generateThinSlots: dirz(2) should be rotated")
        call assert_equal(this%THINSLOTS%ancho(2), 0.2d0, "rotate_generateThinSlots: ancho(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%angx(2), 90.0d0, "rotate_generateThinSlots: angx(2) should remain unchanged")
        call assert_equal(this%THINSLOTS%angy(2), 180.0d0, "rotate_generateThinSlots: angy(2) should be rotated")
        call assert_equal(this%THINSLOTS%angz(2), 90.0d0, "rotate_generateThinSlots: angz(2) should be rotated")
        
        deallocate(this%THINSLOTS%x1, this%THINSLOTS%y1, this%THINSLOTS%z1)
        deallocate(this%THINSLOTS%x2, this%THINSLOTS%y2, this%THINSLOTS%z2)
        deallocate(this%THINSLOTS%dirx, this%THINSLOTS%diry, this%THINSLOTS%dirz)
        deallocate(this%THINSLOTS%ancho)
        deallocate(this%THINSLOTS%angx, this%THINSLOTS%angy, this%THINSLOTS%angz)
        deallocate(this%THINSLOTS)
        
        ! Test case 3: Verify behavior with no slots
        mpidir = 2
        allocate(this%THINSLOTS)
        this%THINSLOTS%nslots = 0
        
        ! Call the rotation routine - should handle empty case gracefully
        call rotate_generateThinSlots(this, mpidir)
        
        ! Verify that nslots remains 0
        call assert_equal(this%THINSLOTS%nslots, 0, "rotate_generateThinSlots: nslots should remain 0")
        
        deallocate(this%THINSLOTS)
    end subroutine test_rotate_generateThinSlots

end module test_nfde_rotate_m 