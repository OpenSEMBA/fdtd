module smbjson_labels_mod

#ifdef CompileWithSMBJSON    
   ! LABELS
   ! -- common labels
   character (len=*), parameter :: J_NAME = "name"
   character (len=*), parameter :: J_ID = "id"
   character (len=*), parameter :: J_TYPE = "type"
   character (len=*), parameter :: J_ELEMENTIDS = "elementIds"

   character (len=*), parameter :: J_DIR_X = "x"
   character (len=*), parameter :: J_DIR_Y = "y"
   character (len=*), parameter :: J_DIR_Z = "z"
   character (len=*), parameter :: J_DIR_M = 'magnitude'

   character (len=*), parameter :: J_FIELD = "field"
   character (len=*), parameter :: J_FIELD_ELECTRIC = "electric"
   character (len=*), parameter :: J_FIELD_MAGNETIC = "magnetic"
   character (len=*), parameter :: J_FIELD_VOLTAGE = "voltage"
   character (len=*), parameter :: J_FIELD_CURRENT = "current"
   character (len=*), parameter :: J_FIELD_CURRENT_DENSITY = "currentDensity"
   character (len=*), parameter :: J_FIELD_CHARGE = "charge"

   ! -- materials
   character (len=*), parameter :: J_MATERIALS = "materials"
   character (len=*), parameter :: J_MAT_REL_PERMITTIVITY = "relativePermittivity"
   character (len=*), parameter :: J_MAT_REL_PERMEABILITY = "relativePermeability"
   character (len=*), parameter :: J_MAT_ELECTRIC_CONDUCTIVITY = "electricConductivity"
   character (len=*), parameter :: J_MAT_MAGNETIC_CONDUCTIVITY = "magneticConductivity"
   
   character (len=*), parameter :: J_MAT_TYPE_PEC = "pec"
   character (len=*), parameter :: J_MAT_TYPE_PMC = "pmc"
   character (len=*), parameter :: J_MAT_TYPE_ISOTROPIC = "isotropic"
   character (len=*), parameter :: J_MAT_TYPE_LUMPED = "lumped"
   character (len=*), parameter :: J_MAT_TYPE_MULTILAYERED_SURFACE = "multilayeredSurface"
   character (len=*), parameter :: J_MAT_TYPE_SLOT = "thinSlot"
   character (len=*), parameter :: J_MAT_TYPE_WIRE = "wire"
   character (len=*), parameter :: J_MAT_TYPE_MULTIWIRE = "multiwire"
   character (len=*), parameter :: J_MAT_TYPE_TERMINAL = "terminal"
   character (len=*), parameter :: J_MAT_TYPE_CONNECTOR = "connector"
   
   character (len=*), parameter :: J_MAT_WIRE_RADIUS = "radius"
   character (len=*), parameter :: J_MAT_WIRE_RESISTANCE = "resistancePerMeter"
   character (len=*), parameter :: J_MAT_WIRE_INDUCTANCE = "inductancePermeter"
   character (len=*), parameter :: J_MAT_WIRE_REF_CAPACITANCE = "__referenceCapacitancePerMeter"
   character (len=*), parameter :: J_MAT_WIRE_REF_INDUCTANCE = "__referenceInductancePerMeter"
   character (len=*), parameter :: J_MAT_WIRE_DIELECTRIC = "dielectric"
   character (len=*), parameter :: J_MAT_WIRE_DIELECTRIC_RADIUS = "radius"
   character (len=*), parameter :: J_MAT_WIRE_DIELECTRIC_PERMITTIVITY = "relativePermittivity"
   character (len=*), parameter :: J_MAT_WIRE_PASS = "isPassthrough"

   character (len=*), parameter :: J_MAT_LUMPED_MODEL = "model"
   character (len=*), parameter :: J_MAT_LUMPED_MODEL_RESISTOR = "resistor"
   character (len=*), parameter :: J_MAT_LUMPED_MODEL_INDUCTOR = "inductor"
   character (len=*), parameter :: J_MAT_LUMPED_MODEL_CAPACITOR = "capacitor"
   character (len=*), parameter :: J_MAT_LUMPED_RESISTANCE = "resistance"
   character (len=*), parameter :: J_MAT_LUMPED_STARTING_TIME = "startingTime"
   character (len=*), parameter :: J_MAT_LUMPED_END_TIME = "endTime"
   character (len=*), parameter :: J_MAT_LUMPED_INDUCTANCE = "inductance"
   character (len=*), parameter :: J_MAT_LUMPED_CAPACITANCE = "capacitance" 
      
   character (len=*), parameter :: J_MAT_TERM_TERMINATIONS = "terminations"
   character (len=*), parameter :: J_MAT_TERM_TYPE_OPEN = "open"
   character (len=*), parameter :: J_MAT_TERM_TYPE_SHORT = "short"
   character (len=*), parameter :: J_MAT_TERM_TYPE_SERIES = "series"
   character (len=*), parameter :: J_MAT_TERM_TYPE_PARALLEL = "parallel"
   character (len=*), parameter :: J_MAT_TERM_TYPE_LsRCp = "LsRCp"
   character (len=*), parameter :: J_MAT_TERM_TYPE_CsLRp = "CsLRp"
   character (len=*), parameter :: J_MAT_TERM_TYPE_RCsLp = "RCsLp"
   character (len=*), parameter :: J_MAT_TERM_TYPE_LCsRp = "LCsRp"

   character (len=*), parameter :: J_MAT_TERM_TYPE_RsLCp = "RsLCp"
   character (len=*), parameter :: J_MAT_TERM_TYPE_RLsCp = "RLsCp"
   character (len=*), parameter :: J_MAT_TERM_TYPE_CIRCUIT = "circuit"

   character (len=*), parameter :: J_MAT_TERM_RESISTANCE = "resistance"
   character (len=*), parameter :: J_MAT_TERM_INDUCTANCE = "inductance"
   character (len=*), parameter :: J_MAT_TERM_CAPACITANCE = "capacitance"
   character (len=*), parameter :: J_MAT_TERM_EXCITATION = "path_to_excitation"
   character (len=*), parameter :: J_MAT_TERM_MODEL_FILE = "file"
   character (len=*), parameter :: J_MAT_TERM_MODEL_NAME = "name"
   character (len=*), parameter :: J_MAT_TERM_MODEL_PORT = "subcircuitPort"

   character (len=*), parameter :: J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE = "transferImpedancePerMeter"
   character (len=*), parameter :: J_MAT_MULTIWIRE_CAPACITANCE = "capacitancePerMeter"
   character (len=*), parameter :: J_MAT_MULTIWIRE_INDUCTANCE = "inductancePerMeter"
   character (len=*), parameter :: J_MAT_MULTIWIRE_RESISTANCE = "resistancePerMeter"
   character (len=*), parameter :: J_MAT_MULTIWIRE_CONDUCTANCE = "conductancePerMeter"

   character (len=*), parameter :: J_MAT_MULTILAYERED_SURF_LAYERS = "layers"
   character (len=*), parameter :: J_MAT_MULTILAYERED_SURF_THICKNESS = "thickness"
   
   character (len=*), parameter :: J_MAT_THINSLOT_WIDTH = "width"

   ! -- materialAssociations
   character (len=*), parameter :: J_MATERIAL_ASSOCIATIONS = "materialAssociations"
   character (len=*), parameter :: J_MATERIAL_ID = "materialId"
   character (len=*), parameter :: J_MATERIAL_PASS = "materialId"

   character (len=*), parameter :: J_MAT_ASS_CAB_INI_TERM_ID = "initialTerminalId"
   character (len=*), parameter :: J_MAT_ASS_CAB_END_TERM_ID = "endTerminalId"
   character (len=*), parameter :: J_MAT_ASS_CAB_INI_CONN_ID = "initialConnectorId"
   character (len=*), parameter :: J_MAT_ASS_CAB_END_CONN_ID = "endConnectorId"
   character (len=*), parameter :: J_MAT_ASS_CAB_CONTAINED_WITHIN_ID = "containedWithinElementId"
   
   ! -- connector
   character (len=*), parameter :: J_MAT_CONN_RESISTANCES = "resistances"
   character (len=*), parameter :: J_MAT_CONN_TRANSFER_IMPEDANCE = "transferImpedancePerMeter"

   ! -- transferImpedancePerMeter
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_RESISTANCE = "resistiveTerm"
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE = "inductiveTerm"
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_DIRECTION = "direction"
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_POLES = "poles"
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_RESIDUES = "residues"
   character (len=*), parameter :: J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES = "numberOfPoles"

   ! --  SPICE subcircuits
   character (len=*), parameter :: J_SUBCIRCUITS  = "subcircuits"
   character (len=*), parameter :: J_SUBCKT_NAME  = "name"
   character (len=*), parameter :: J_SUBCKT_PORTS = "numberOfPorts"
   character (len=*), parameter :: J_SUBCKT_FILE  = "file"
   
   ! -- Mesh and geometry.
   character (len=*), parameter :: J_MESH = "mesh"
   
   character (len=*), parameter :: J_GRID = "grid"
   character (len=*), parameter :: J_COORDINATES = "coordinates"
   character (len=*), parameter :: J_ELEMENTS = "elements"
   
   character (len=*), parameter :: J_GRID_NUMBER_OF_CELLS = "numberOfCells"
   character (len=*), parameter :: J_GRID_STEPS = "steps"
   
   character (len=*), parameter :: J_COORDINATE_IDS = "coordinateIds"
   character (len=*), parameter :: J_COORDINATE_POS = "relativePosition"
   
   character (len=*), parameter :: J_ELEM_TYPE_NODE = "node"
   character (len=*), parameter :: J_ELEM_TYPE_POLYLINE = "polyline"
   character (len=*), parameter :: J_ELEM_TYPE_CELL = "cell"
   character (len=*), parameter :: J_CELL_INTERVALS = "intervals"

   ! type(NFDEGeneral)
   character (len=*), parameter :: J_GENERAL = "general"
   character (len=*), parameter :: J_GEN_TIME_STEP = "timeStep"
   character (len=*), parameter :: J_GEN_NUMBER_OF_STEPS = "numberOfSteps"
   character (len=*), parameter :: J_GEN_MTLN_PROBLEM = "mtlnProblem"
   character (len=*), parameter :: J_GEN_ADDITIONAL_ARGUMENTS = "additionalArguments"


   ! type(Frontera)
   character (len=*), parameter :: J_BOUNDARY = "boundary"
   character (len=*), parameter :: J_BND_ALL = "all"
   character (len=*), parameter :: J_BND_XL = "xLower"
   character (len=*), parameter :: J_BND_XU = "xUpper"
   character (len=*), parameter :: J_BND_YL = "yLower"
   character (len=*), parameter :: J_BND_YU = "yUpper"
   character (len=*), parameter :: J_BND_ZL = "zLower"
   character (len=*), parameter :: J_BND_ZU = "zUpper"


   character (len=*), parameter :: J_BND_TYPE_PEC = "pec"
   character (len=*), parameter :: J_BND_TYPE_PMC = "pmc"
   character (len=*), parameter :: J_BND_TYPE_PERIODIC = "periodic"
   character (len=*), parameter :: J_BND_TYPE_MUR = "mur"
   character (len=*), parameter :: J_BND_TYPE_PML = "pml"
   character (len=*), parameter :: J_BND_PML_LAYERS = "layers"
   character (len=*), parameter :: J_BND_PML_ORDER = "order"
   character (len=*), parameter :: J_BND_PML_REFLECTION = "reflection"

   ! -- source types
   character (len=*), parameter :: J_SOURCES = "sources"
   character (len=*), parameter :: J_SRC_MAGNITUDE_FILE = "magnitudeFile"
   
   character (len=*), parameter :: J_SRC_TYPE_PW = "planewave"
   character (len=*), parameter :: J_SRC_TYPE_NS = "nodalSource"
   character (len=*), parameter :: J_SRC_TYPE_GEN = "generator"
   character (len=*), parameter :: J_SRC_ATTACHED_ID = "attachedToLineId"

   ! type(Planewave)
   character (len=*), parameter :: J_SRC_PW_DIRECTION = "direction"
   character (len=*), parameter :: J_SRC_PW_POLARIZATION = "polarization"
   character (len=*), parameter :: J_SRC_PW_THETA = "theta"
   character (len=*), parameter :: J_SRC_PW_PHI = "phi"

   ! type(NodalSource)
   character (len=*), parameter :: J_SRC_NS_HARDNESS = "hardness"
   character (len=*), parameter :: J_SRC_NS_HARDNESS_SOFT = "soft"
   character (len=*), parameter :: J_SRC_NS_HARDNESS_HARD = "hard"

   ! --- probe types
   character (len=*), parameter :: J_PROBES = "probes"
   
   character (len=*), parameter :: J_PR_TYPE_POINT = "point"
   character (len=*), parameter :: J_PR_TYPE_WIRE = "wire"
   character (len=*), parameter :: J_PR_TYPE_BULK_CURRENT = "bulkCurrent"
   character (len=*), parameter :: J_PR_TYPE_FARFIELD = "farField"
   character (len=*), parameter :: J_PR_TYPE_MOVIE = "movie"
   character (len=*), parameter :: J_PR_TYPE_LINE = "line"
   
   character (len=*), parameter :: J_PR_POINT_DIRECTIONS = "directions"

   character (len=*), parameter :: J_PR_MOVIE_COMPONENT = "component"

   character (len=*), parameter :: J_PR_FAR_FIELD_THETA = "theta"
   character (len=*), parameter :: J_PR_FAR_FIELD_PHI = "phi"
   character (len=*), parameter :: J_PR_FAR_FIELD_DIR_INITIAL = "initial"
   character (len=*), parameter :: J_PR_FAR_FIELD_DIR_FINAL = "final"
   character (len=*), parameter :: J_PR_FAR_FIELD_DIR_STEP = "step"

   ! domain stuff
   character (len=*), parameter :: J_PR_DOMAIN = "domain"
   character (len=*), parameter :: J_PR_DOMAIN_TYPE = "type"
   
   character (len=*), parameter :: J_PR_DOMAIN_TYPE_TIME = "time"
   character (len=*), parameter :: J_PR_DOMAIN_TYPE_FREQ = "frequency"
   character (len=*), parameter :: J_PR_DOMAIN_TYPE_TIMEFREQ = "timeFrequency"

   character (len=*), parameter :: J_PR_DOMAIN_MAGNITUDE_FILE = "magnitudeFile"

   character (len=*), parameter :: J_PR_DOMAIN_TIME_START = "initialTime"
   character (len=*), parameter :: J_PR_DOMAIN_TIME_STOP   = "finalTime"
   character (len=*), parameter :: J_PR_DOMAIN_TIME_STEP  = "samplingPeriod"
   
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_START = "initialFrequency"
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_STOP   = "finalFrequency"
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_NUMBER  = "numberOfFrequencies"
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_SPACING  = "frequencySpacing"
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_SPACING_LINEAR  = "linear"
   character (len=*), parameter :: J_PR_DOMAIN_FREQ_SPACING_LOGARITHMIC  = "logarithmic"
#endif
end module
