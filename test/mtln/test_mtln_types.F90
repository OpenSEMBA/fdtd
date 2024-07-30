integer function test_derived_type_submodule() bind(C) result(err)

   use mtln_solver_mod
   
   implicit none
   
   type(terminal_node_t) :: node
   type(terminal_connection_t) :: connection
   type(terminal_network_t) :: network
   type(parsed_mtln_t) :: parsed 
   character(len=256), parameter :: square_excitation = trim('coaxial_line_paul_8_6_0.25_square.exc')
   type(termination_t) :: t
   type(node_source_t) :: source
   err = 0 

   node%conductor_in_cable = 1
   node%side = TERMINAL_NODE_SIDE_INI

   source%path_to_excitation = trim(square_excitation)
   source%source_type = SOURCE_TYPE_VOLTAGE
   node%termination = termination_t(source = source, &
                                    termination_type = TERMINATION_SERIES, &
                                    resistance = 150, &
                                    inductance = 0.0, &
                                    capacitance = 1e22)

   allocate(connection%nodes(1))                                                   
   connection%nodes(1) = node
   allocate(network%connections(1))
   network%connections(1) = connection
   allocate(parsed%networks(1))
   parsed%networks(1) = network

   if (parsed%networks(1)%connections(1)%nodes(1)%termination%source%path_to_excitation &
      /= trim(square_excitation)) then
      err = err + 1
   end if

end function

integer function test_mtln_types() bind(C) result(err)

   use mtln_types_mod
   implicit none

   
   type(terminal_node_t) :: t
   type(node_source_t) :: source

   err = 0

   source%path_to_excitation = "path"
   source%source_type = SOURCE_TYPE_VOLTAGE

   t%termination = termination_t(source = source, &
                                 termination_type = TERMINATION_SERIES, &
                                 resistance = 150, &
                                 inductance = 0.0, &
                                 capacitance = 1e22)


   if (t%termination%source%path_to_excitation /= "path") then 
         err = err + 1
   end if

end function

