module mtln_types_mod
   use fdetypes, ONLY: RKIND   !sggmtln
   implicit none

   integer, parameter :: TERMINATION_UNDEFINED  = -1
   integer, parameter :: TERMINATION_SHORT      =  1
   integer, parameter :: TERMINATION_OPEN       =  2
   integer, parameter :: TERMINATION_SERIES     =  3
   integer, parameter :: TERMINATION_LCpRs      =  4
   integer, parameter :: TERMINATION_RLsCp      =  5
   integer, parameter :: TERMINATION_CIRCUIT      =  6

   integer, parameter :: TERMINAL_NODE_SIDE_UNDEFINED = -1
   integer, parameter :: TERMINAL_NODE_SIDE_INI       =  1
   integer, parameter :: TERMINAL_NODE_SIDE_END       =  2

   integer, parameter :: TRANSFER_IMPEDANCE_DIRECTION_INWARDS   =  1
   integer, parameter :: TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS  =  2
   integer, parameter :: TRANSFER_IMPEDANCE_DIRECTION_BOTH      =  3

   integer, parameter :: PROBE_TYPE_UNDEFINED = -1
   integer, parameter :: PROBE_TYPE_VOLTAGE   =  1
   integer, parameter :: PROBE_TYPE_CURRENT   =  2

   integer, parameter :: SOURCE_TYPE_UNDEFINED = -1
   integer, parameter :: SOURCE_TYPE_VOLTAGE   =  1
   integer, parameter :: SOURCE_TYPE_CURRENT   =  2

   integer, parameter :: DIRECTION_X_POS   =  1
   integer, parameter :: DIRECTION_X_NEG   =  -1
   integer, parameter :: DIRECTION_Y_POS   =  2
   integer, parameter :: DIRECTION_Y_NEG   =  -2
   integer, parameter :: DIRECTION_Z_POS   =  3
   integer, parameter :: DIRECTION_Z_NEG   =  -3

   type :: external_field_segment_t
      integer, dimension(3) ::position
      integer :: direction
      real :: radius, dielectricRadius, dielectricRelativePermittivity
      real :: effectiveRelativePermittivity
      real (kind=rkind) , pointer  ::  field => null()
   contains
      private
      procedure :: external_field_segments_eq 
      generic, public :: operator(==) => external_field_segments_eq
   end type

   type node_source_t
      character(len=256) :: path_to_excitation = ""
      integer :: source_type = SOURCE_TYPE_UNDEFINED
   end type

   type terminal_circuit_t
      character(len=256) :: file = ""
      character(len=256) :: model_name = ""
   end type

   type, public :: termination_t
      integer :: termination_type = TERMINATION_UNDEFINED
      real :: resistance = 0.0
      real :: inductance = 0.0
      real :: capacitance = 1e22
      type(node_source_t) :: source
      type(terminal_circuit_t) :: model
      integer :: subcircuitPort = -1
   contains
      private
      procedure :: termination_eq
      generic, public :: operator(==) => termination_eq
   end type


   type :: terminal_node_t
      type(cable_t), pointer :: belongs_to_cable => null()
      integer :: conductor_in_cable
      integer :: side = TERMINAL_NODE_SIDE_UNDEFINED
      type(termination_t) :: termination
   contains
      private
      procedure :: terminal_node_eq
      generic, public :: operator(==) => terminal_node_eq
   end type


   type :: subcircuit_t
      character(len=256) :: model_file = ""
      character(len=256) :: model_name = ""
      character(len=256) :: subcircuit_name = ""
      integer :: numberOfPorts
      integer :: nodeId
   end type
      
   type :: terminal_connection_t
      type(terminal_node_t), dimension(:), allocatable :: nodes
      type(subcircuit_t) :: subcircuit
      logical :: has_subcircuit = .false.
   contains
      private
      procedure :: terminal_connection_eq
      generic, public :: operator(==) => terminal_connection_eq
      procedure, public :: add_node => terminal_connection_add_node
   end type

   type :: terminal_network_t
      type(terminal_connection_t), dimension(:), allocatable :: connections
   contains
      private
      procedure :: terminal_network_eq
      generic, public :: operator(==) => terminal_network_eq
      procedure, public :: add_connection => terminal_network_add_connection
   end type

   type, public :: transfer_impedance_per_meter_t
      real :: inductive_term = 0.0
      real :: resistive_term = 0.0
      complex, dimension(:), allocatable :: poles, residues ! poles and residues
      integer :: direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS
   contains
      private
      procedure :: transfer_impedance_per_meter_eq
      generic, public :: operator(==) => transfer_impedance_per_meter_eq

   end type

   type :: connector_t
      integer :: id
      real, dimension(:), allocatable :: resistances
      type(transfer_impedance_per_meter_t) :: transfer_impedance_per_meter
   contains
      private
      procedure :: connector_eq
      generic, public :: operator(==) => connector_eq
   end type

   type, public :: cable_t
      character (len=:), allocatable :: name
      real, allocatable, dimension(:,:) :: resistance_per_meter
      real, allocatable, dimension(:,:) :: capacitance_per_meter
      real, allocatable, dimension(:,:) :: inductance_per_meter
      real, allocatable, dimension(:,:) :: conductance_per_meter
      real, allocatable, dimension(:) :: step_size
      type(transfer_impedance_per_meter_t) :: transfer_impedance
      type(cable_t), pointer :: parent_cable => null()
      integer :: conductor_in_parent = -1
      type(connector_t), pointer :: initial_connector => null()
      type(connector_t), pointer :: end_connector => null()
      type(external_field_segment_t), allocatable, dimension(:) :: external_field_segments
      logical :: isPassthrough = .false.

   contains
      private
      procedure :: cable_eq
      generic, public :: operator(==) => cable_eq
   end type

   type :: probe_t
      type(cable_t), pointer :: attached_to_cable => null()
      integer :: index
      integer :: probe_type = PROBE_TYPE_UNDEFINED
      character (len=:), allocatable :: probe_name
      real, dimension(3) :: probe_position
   contains
      private
      procedure :: probe_eq
      generic, public :: operator(==) => probe_eq
   end type

   type, public :: mtln_t
      type(cable_t), dimension(:), pointer :: cables
      type(terminal_network_t), dimension(:), allocatable :: networks
      type(probe_t), dimension(:), allocatable :: probes
      type(connector_t), dimension(:), pointer :: connectors
      real :: time_step
      integer :: number_of_steps
      logical :: has_multiwires
   contains
      private
      procedure :: mtln_eq
      generic, public :: operator(==) => mtln_eq
   end type


contains

   logical function mtln_eq(a,b)
      class(mtln_t), intent(in) :: a,b
      integer :: i
         
      if (size(a%cables) /= size(b%cables)) then
         mtln_eq = .false.
         return
      end if
      do i = 1, size(a%cables)
         if (.not. a%cables(i) == b%cables(i)) then
            mtln_eq = .false.
            return
         end if
      end do

      if (size(a%probes) /= size(b%probes)) then
         mtln_eq = .false.
         return
      end if
      do i = 1, size(a%probes)
         if (.not. a%probes(i) == b%probes(i)) then
            mtln_eq = .false.
            return
         end if
      end do

      if (size(a%networks) /= size(b%networks)) then
         mtln_eq = .false.
         return
      end if
      do i = 1, size(a%networks)
         if (.not. a%networks(i) == b%networks(i)) then
            mtln_eq = .false.
            return
         end if
      end do

      mtln_eq = .true.
   end function

   elemental logical function transfer_impedance_per_meter_eq(a,b)
      class(transfer_impedance_per_meter_t), intent(in) :: a, b
      transfer_impedance_per_meter_eq = &
         (a%inductive_term == b%inductive_term) .and. &
         (a%resistive_term == b%resistive_term) .and. &
         all(a%poles == b%poles) .and. &
         all(a%residues == b%residues) .and. &
         (a%direction == b%direction)
   end function

   recursive logical function cable_eq(a,b)
      class(cable_t), intent(in) :: a, b
      cable_eq = .true.
      cable_eq = cable_eq .and.  (a%name == b%name) 
      cable_eq = cable_eq .and.  all(a%inductance_per_meter == b%inductance_per_meter)
      cable_eq = cable_eq .and.  all(a%capacitance_per_meter == b%capacitance_per_meter)
      cable_eq = cable_eq .and.  all(a%resistance_per_meter == b%resistance_per_meter)
      cable_eq = cable_eq .and.  all(a%conductance_per_meter == b%conductance_per_meter)
      cable_eq = cable_eq .and.  all(a%step_size == b%step_size)
      cable_eq = cable_eq .and.  (a%transfer_impedance == b%transfer_impedance)
      cable_eq = cable_eq .and.  (a%conductor_in_parent == b%conductor_in_parent)
      cable_eq = cable_eq .and.  all(a%external_field_segments == b%external_field_segments)


      if (.not. cable_eq) then
         cable_eq = .false.
      end if

      if (.not. associated(a%parent_cable) .and. .not. associated(b%parent_cable)) then
         cable_eq = cable_eq .and. .true.
      else if ((associated(a%parent_cable) .and. .not. associated(b%parent_cable)) .or. &
         (.not. associated(a%parent_cable) .and. associated(b%parent_cable))) then
         cable_eq = cable_eq .and. .false.
      else
         cable_eq = cable_eq .and. (a%parent_cable == b%parent_cable)
      end if

      if (.not. cable_eq) then
         cable_eq = .false.
      end if

      if (.not. associated(a%initial_connector) .and. .not. associated(b%initial_connector)) then
         cable_eq = cable_eq .and. .true.
      else if ((associated(a%initial_connector) .and. .not. associated(b%initial_connector)) .or. &
         (.not. associated(a%initial_connector) .and. associated(b%initial_connector))) then
         cable_eq = cable_eq .and. .false.
      else
         cable_eq = cable_eq .and. (a%initial_connector == b%initial_connector)
      end if
      if (.not. cable_eq) then
         cable_eq = .false.
      end if

      if (.not. associated(a%end_connector) .and. .not. associated(b%end_connector)) then
         cable_eq = cable_eq .and. .true.
      else if ((associated(a%end_connector) .and. .not. associated(b%end_connector)) .or. &
         (.not. associated(a%end_connector) .and. associated(b%end_connector))) then
         cable_eq = cable_eq .and. .false.
      else
         cable_eq = cable_eq .and. (a%end_connector == b%end_connector)
      end if
      if (.not. cable_eq) then
         cable_eq = .false.
      end if

   end function

   elemental logical function connector_eq(a,b)
      class(connector_t), intent(in) :: a, b
      logical :: l
      connector_eq = &
         (a%id == b%id) .and. &
         (all(a%resistances == b%resistances)) .and. &
         (a%transfer_impedance_per_meter == b%transfer_impedance_per_meter)
   end function

   elemental logical function termination_eq(a, b)
      class(termination_t), intent(in) :: a
      type(termination_t), intent(in) :: b
      termination_eq = &
         (a%termination_type == b%termination_type) .and. &
         (a%resistance == b%resistance) .and. &
         (a%inductance == b%inductance) .and. &
         (a%capacitance == b%capacitance) .and. &
         a%source%path_to_excitation == b%source%path_to_excitation .and. &
         a%source%source_type == b%source%source_type
   end function

   logical function probe_eq(a,b)
      class(probe_t), intent(in) :: a,b
      probe_eq = &
         (a%index == b%index) .and. &
         (a%probe_type == b%probe_type) .and. &
         (a%probe_name == b%probe_name) .and. &
         all(a%probe_position == b%probe_position)


      if (.not. associated(a%attached_to_cable) .and. .not. associated(b%attached_to_cable)) then
         probe_eq = probe_eq .and. .true.
      else if ((associated(a%attached_to_cable) .and. .not. associated(b%attached_to_cable)) .or. &
         (.not. associated(a%attached_to_cable) .and. associated(b%attached_to_cable))) then
         probe_eq = probe_eq .and. .false.
      else
         probe_eq = probe_eq .and. (a%attached_to_cable == b%attached_to_cable)
      end if
      if (probe_eq .eqv. .false.) then
         probe_eq = .false.
      end if
   end function

   logical function terminal_node_eq(a, b)
      class(terminal_node_t), intent(in) :: a, b

      terminal_node_eq = &
         (a%conductor_in_cable == b%conductor_in_cable) .and. &
         (a%side == b%side) .and. &
         (a%termination == b%termination)

      if (.not. associated(a%belongs_to_cable) .and. .not. associated(b%belongs_to_cable)) then
         terminal_node_eq = terminal_node_eq .and. .true.
      else if ((associated(a%belongs_to_cable) .and. .not. associated(b%belongs_to_cable)) .or. &
         (.not. associated(a%belongs_to_cable) .and. associated(b%belongs_to_cable))) then
         terminal_node_eq = terminal_node_eq .and. .false.
      else
         terminal_node_eq = terminal_node_eq .and. (a%belongs_to_cable == b%belongs_to_cable)
      end if

   end function

   logical function terminal_connection_eq(a,b)
      class(terminal_connection_t), intent(in) :: a,b
      integer :: i
      if (size(a%nodes) /= size(b%nodes)) then
         terminal_connection_eq = .false.
         return
      end if
      do i = 1, size(a%nodes)
         if (.not. (a%nodes(i) == b%nodes(i))) then
            terminal_connection_eq = .false.
            return
         end if
      end do
      terminal_connection_eq = .true.
   end function

   logical function terminal_network_eq(a,b)
      class(terminal_network_t), intent(in) :: a,b
      integer :: i
      if (size(a%connections) /= size(b%connections)) then
         terminal_network_eq = .false.
         return
      end if
      do i = 1, size(a%connections)
         if (.not. (a%connections(i) == b%connections(i))) then
            terminal_network_eq = .false.
            return
         end if
      end do
      terminal_network_eq = .true.
   end function

   elemental logical function external_field_segments_eq(a,b)
      class(external_field_segment_t), intent(in) :: a,b
      external_field_segments_eq = &
         all(a%position == b%position) .and. &
         a%direction == b%direction

      if (.not. associated(a%field) .and. .not. associated(b%field)) then
         external_field_segments_eq = external_field_segments_eq .and. .true.
      else if ((associated(a%field) .and. .not. associated(b%field)) .or. &
         (.not. associated(a%field) .and. associated(b%field))) then
            external_field_segments_eq = external_field_segments_eq .and. .false.
      else
         external_field_segments_eq = external_field_segments_eq .and. (a%field == b%field)
      end if


   end function

   subroutine terminal_connection_add_node(this, node)
      class(terminal_connection_t) :: this
      type(terminal_node_t) :: node
      type(terminal_node_t), dimension(:), allocatable :: newNodes
      integer :: newNodesSize

      if (.not. allocated(this%nodes))  allocate(this%nodes(0))

      allocate(newNodes( size(this%nodes) + 1 ) )
      newNodesSize = size(newNodes)
      newNodes(1:newNodesSize-1) = this%nodes
      newNodes(newNodesSize) = node
      call MOVE_ALLOC(from=newNodes, to=this%nodes)

   end subroutine

   subroutine terminal_network_add_connection(this, connection)
      class(terminal_network_t) :: this
      type(terminal_connection_t) :: connection
      type(terminal_connection_t), dimension(:), allocatable :: newConnections
      integer :: newConnectionsSize
      
      if (.not. allocated(this%connections))  allocate(this%connections(0))
      
      allocate(newConnections( size(this%connections) + 1 ) )
      newConnectionsSize = size(newConnections)
      newConnections(1:newConnectionsSize-1) = this%connections
      newConnections(newConnectionsSize) = connection
      call MOVE_ALLOC(from=newConnections, to=this%connections)
   end subroutine

end module
