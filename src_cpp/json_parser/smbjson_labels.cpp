// smbjson_labels_m.cpp
// This file contains the constants defined in the Fortran module smbjson_labels_m.
// Since the module only contains parameters (constants), it is best represented
// as a namespace with constexpr or const std::string members in C++.

#include <string>

namespace smbjson_labels_m {

#ifdef CompileWithSMBJSON

    // LABELS
    // -- common labels
    constexpr const char* J_NAME = "name";
    constexpr const char* J_ID = "id";
    constexpr const char* J_TYPE = "type";
    constexpr const char* J_SUBTYPE = "subtype";
    constexpr const char* J_ELEMENTIDS = "elementIds";

    constexpr const char* J_DIR = "direction";
    constexpr const char* J_DIR_X = "x";
    constexpr const char* J_DIR_Y = "y";
    constexpr const char* J_DIR_Z = "z";
    constexpr const char* J_DIR_M = "magnitude";

    constexpr const char* J_FIELD = "field";
    constexpr const char* J_FIELD_ELECTRIC = "electric";
    constexpr const char* J_FIELD_MAGNETIC = "magnetic";
    constexpr const char* J_FIELD_VOLTAGE = "voltage";
    constexpr const char* J_FIELD_CURRENT = "current";
    constexpr const char* J_FIELD_CURRENT_DENSITY = "currentDensity";
    constexpr const char* J_FIELD_CHARGE = "charge";

    // -- materials
    constexpr const char* J_MATERIALS = "materials";
    constexpr const char* J_MAT_ABS_PERMITTIVITY = "absolutePermittivity";
    constexpr const char* J_MAT_ABS_PERMEABILITY = "absolutePermeability";
    constexpr const char* J_MAT_REL_PERMITTIVITY = "relativePermittivity";
    constexpr const char* J_MAT_REL_PERMEABILITY = "relativePermeability";
    constexpr const char* J_MAT_ELECTRIC_CONDUCTIVITY = "electricConductivity";
    constexpr const char* J_MAT_MAGNETIC_CONDUCTIVITY = "magneticConductivity";
    
    constexpr const char* J_MAT_TYPE_PEC = "pec";
    constexpr const char* J_MAT_TYPE_PMC = "pmc";
    constexpr const char* J_MAT_TYPE_ISOTROPIC = "isotropic";
    constexpr const char* J_MAT_TYPE_LUMPED = "lumped";
    constexpr const char* J_MAT_TYPE_MULTILAYERED_SURFACE = "multilayeredSurface";
    constexpr const char* J_MAT_TYPE_SLOT = "thinSlot";
    constexpr const char* J_MAT_TYPE_WIRE = "wire";
    constexpr const char* J_MAT_TYPE_SHIELDED_MULTIWIRE = "shieldedMultiwire";
    constexpr const char* J_MAT_TYPE_UNSHIELDED_MULTIWIRE = "unshieldedMultiwire";
    constexpr const char* J_MAT_TYPE_TERMINAL = "terminal";
    constexpr const char* J_MAT_TYPE_CONNECTOR = "connector";
    
    constexpr const char* J_MAT_WIRE_RADIUS = "radius";
    constexpr const char* J_MAT_WIRE_RESISTANCE = "resistancePerMeter";
    constexpr const char* J_MAT_WIRE_INDUCTANCE = "inductancePermeter";

    constexpr const char* J_MAT_LUMPED_MODEL = "model";
    constexpr const char* J_MAT_LUMPED_MODEL_RESISTOR = "resistor";
    constexpr const char* J_MAT_LUMPED_MODEL_INDUCTOR = "inductor";
    constexpr const char* J_MAT_LUMPED_MODEL_CAPACITOR = "capacitor";
    constexpr const char* J_MAT_LUMPED_RESISTANCE = "resistance";
    constexpr const char* J_MAT_LUMPED_STARTING_TIME = "startingTime";
    constexpr const char* J_MAT_LUMPED_END_TIME = "endTime";
    constexpr const char* J_MAT_LUMPED_INDUCTANCE = "inductance";
    constexpr const char* J_MAT_LUMPED_CAPACITANCE = "capacitance"; 
      
    constexpr const char* J_MAT_TERM_TERMINATIONS = "terminations";
    constexpr const char* J_MAT_TERM_TYPE_OPEN = "open";
    constexpr const char* J_MAT_TERM_TYPE_SHORT = "short";
    constexpr const char* J_MAT_TERM_TYPE_SERIES = "series";
    constexpr const char* J_MAT_TERM_TYPE_PARALLEL = "parallel";
    constexpr const char* J_MAT_TERM_TYPE_LsRCp = "LsRCp";
    constexpr const char* J_MAT_TERM_TYPE_CsLRp = "CsLRp";
    constexpr const char* J_MAT_TERM_TYPE_RCsLp = "RCsLp";
    constexpr const char* J_MAT_TERM_TYPE_LCsRp = "LCsRp";
    constexpr const char* J_MAT_TERM_TYPE_RsLCp = "RsLCp";
    constexpr const char* J_MAT_TERM_TYPE_RLsCp = "RLsCp";
    constexpr const char* J_MAT_TERM_TYPE_CIRCUIT = "circuit";
    constexpr const char* J_MAT_TERM_TYPE_NETWORK = "network";

    constexpr const char* J_MAT_TERM_RESISTANCE = "resistance";
    constexpr const char* J_MAT_TERM_INDUCTANCE = "inductance";
    constexpr const char* J_MAT_TERM_CAPACITANCE = "capacitance";
    constexpr const char* J_MAT_TERM_EXCITATION = "path_to_excitation";
    constexpr const char* J_MAT_TERM_MODEL_FILE = "file";
    constexpr const char* J_MAT_TERM_MODEL_NAME = "name";
    constexpr const char* J_MAT_TERM_MODEL_NODE = "node";

    constexpr const char* J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE = "transferImpedancePerMeter";
    constexpr const char* J_MAT_MULTIWIRE_CAPACITANCE = "capacitancePerMeter";
    constexpr const char* J_MAT_MULTIWIRE_INDUCTANCE = "inductancePerMeter";
    constexpr const char* J_MAT_MULTIWIRE_RESISTANCE = "resistancePerMeter";
    constexpr const char* J_MAT_MULTIWIRE_CONDUCTANCE = "conductancePerMeter";
    
    constexpr const char* J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION = "multipolarExpansion";
    // ME = Multipolar Expansion
    constexpr const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX = "innerRegionBox";
    constexpr const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MAX = "max";
    constexpr const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MIN = "min";
    constexpr const char* J_MAT_MULTIWIRE_ME_ELECTRIC = "electric";
    constexpr const char* J_MAT_MULTIWIRE_ME_MAGNETIC = "magnetic";
    // MEFR = Multipolar Expansion Field Reconstruction
    constexpr const char* J_MAT_MULTIWIRE_MEFR_AB = "ab";
    constexpr const char* J_MAT_MULTIWIRE_MEFR_CONDUCTOR_POTENTIALS = "conductorPotentials";
    constexpr const char* J_MAT_MULTIWIRE_MEFR_EXPANSION_CENTER = "expansionCenter";
    constexpr const char* J_MAT_MULTIWIRE_MEFR_INNER_REGION_AVERAGE_POTENTIAL = "innerRegionAveragePotential";


    constexpr const char* J_MAT_MULTILAYERED_SURF_LAYERS = "layers";
    constexpr const char* J_MAT_MULTILAYERED_SURF_THICKNESS = "thickness";
    
    constexpr const char* J_MAT_THINSLOT_WIDTH = "width";

    // -- materialAssociations
    constexpr const char* J_MATERIAL_ASSOCIATIONS = "materialAssociations";
    constexpr const char* J_MATERIAL_ID = "materialId";
    constexpr const char* J_MATERIAL_PASS = "materialId";

    constexpr const char* J_MAT_ASS_CAB_INI_TERM_ID = "initialTerminalId";
    constexpr const char* J_MAT_ASS_CAB_END_TERM_ID = "endTerminalId";
    constexpr const char* J_MAT_ASS_CAB_INI_CONN_ID = "initialConnectorId";
    constexpr const char* J_MAT_ASS_CAB_END_CONN_ID = "endConnectorId";
    constexpr const char* J_MAT_ASS_CAB_CONTAINED_WITHIN_ID = "containedWithinElementId";
    constexpr const char* J_MAT_ASS_TOTAL_RESISTANCE = "totalResistance";
    
    // -- connector
    constexpr const char* J_MAT_CONN_RESISTANCES = "resistances";
    constexpr const char* J_MAT_CONN_TRANSFER_IMPEDANCES = "transferImpedancesPerMeter";

    // -- transferImpedancePerMeter
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_RESISTANCE = "resistiveTerm";
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE = "inductiveTerm";
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_DIRECTION = "direction";
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_POLES = "poles";
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_RESIDUES = "residues";
    constexpr const char* J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES = "numberOfPoles";
    
    // -- Mesh and geometry.
    constexpr const char* J_MESH = "mesh";
    
    constexpr const char* J_GRID = "grid";
    constexpr const char* J_COORDINATES = "coordinates";
    constexpr const char* J_ELEMENTS = "elements";
    
    constexpr const char* J_GRID_NUMBER_OF_CELLS = "numberOfCells";
    constexpr const char* J_GRID_STEPS = "steps";
    constexpr const char* J_GRID_ORIGIN = "origin";
    
    constexpr const char* J_COORDINATE_IDS = "coordinateIds";
    constexpr const char* J_COORDINATE_POS = "relativePosition";
    
    constexpr const char* J_ELEM_TYPE_NODE = "node";
    constexpr const char* J_ELEM_TYPE_POLYLINE = "polyline";
    constexpr const char* J_ELEM_TYPE_CELL = "cell";
    constexpr const char* J_CELL_INTERVALS = "intervals";
    constexpr const char* J_CONF_VOLUME_TRIANGLES = "triangles";
    constexpr const char* J_CONF_VOLUME_INTERVALS = "intervals";
    
    constexpr const char* J_CONF_SUBTYPE_VOLUME  = "volume";
    constexpr const char* J_CONF_SUBTYPE_SURFACE = "surface";

    // type(NFDEGeneral_t)
    constexpr const char* J_GENERAL = "general";
    constexpr const char* J_GEN_TIME_STEP = "timeStep";
    constexpr const char* J_GEN_NUMBER_OF_STEPS = "numberOfSteps";
    constexpr const char* J_GEN_MTLN_PROBLEM = "mtlnProblem";
    constexpr const char* J_GEN_ADDITIONAL_ARGUMENTS = "additionalArguments";

    // background
    constexpr const char* J_BACKGROUND = "background";
    constexpr const char* J_BKG_ABS_PERMITTIVITY = "absolutePermittivity";
    constexpr const char* J_BKG_ABS_PERMEABILITY = "absolutePermeability";

    // type(Frontera_t)
    constexpr const char* J_BOUNDARY = "boundary";
    constexpr const char* J_BND_ALL = "all";
    constexpr const char* J_BND_XL = "xLower";
    constexpr const char* J_BND_XU = "xUpper";
    constexpr const char* J_BND_YL = "yLower";
    constexpr const char* J_BND_YU = "yUpper";
    constexpr const char* J_BND_ZL = "zLower";
    constexpr const char* J_BND_ZU = "zUpper";


    constexpr const char* J_BND_TYPE_PEC = "pec";
    constexpr const char* J_BND_TYPE_PMC = "pmc";
    constexpr const char* J_BND_TYPE_PERIODIC = "periodic";
    constexpr const char* J_BND_TYPE_MUR = "mur";
    constexpr const char* J_BND_TYPE_PML = "pml";
    constexpr const char* J_BND_PML_LAYERS = "layers";
    constexpr const char* J_BND_PML_ORDER = "order";
    constexpr const char* J_BND_PML_REFLECTION = "reflection";

    // -- source types
    constexpr const char* J_SOURCES = "sources";
    constexpr const char* J_SRC_MAGNITUDE_FILE = "magnitudeFile";
    
    constexpr const char* J_SRC_TYPE_PW = "planewave";
    constexpr const char* J_SRC_TYPE_NS = "nodalSource";
    constexpr const char* J_SRC_TYPE_GEN = "generator";
    constexpr const char* J_SRC_ATTACHED_ID = "attachedToLineId";
    constexpr const char* J_SRC_RESISTANCE_GEN = "resistance";

    // type(PlaneWave_t)
    constexpr const char* J_SRC_PW_DIRECTION = "direction";
    constexpr const char* J_SRC_PW_POLARIZATION = "polarization";
    constexpr const char* J_SRC_PW_THETA = "theta";
    constexpr const char* J_SRC_PW_PHI = "phi";

    // type(NodalSource)
    constexpr const char* J_SRC_NS_HARDNESS = "hardness";
    constexpr const char* J_SRC_NS_HARDNESS_SOFT = "soft";
    constexpr const char* J_SRC_NS_HARDNESS_HARD = "hard";

    // --- probe types
    constexpr const char* J_PROBES = "probes";
    
    constexpr const char* J_PR_TYPE_POINT = "point";
    constexpr const char* J_PR_TYPE_WIRE = "wire";
    constexpr const char* J_PR_TYPE_BULK_CURRENT = "bulkCurrent";
    constexpr const char* J_PR_TYPE_FARFIELD = "farField";
    constexpr const char* J_PR_TYPE_MOVIE = "movie";
    constexpr const char* J_PR_TYPE_LINE = "line";
    
    constexpr const char* J_PR_POINT_DIRECTIONS = "directions";

    constexpr const char* J_PR_MOVIE_COMPONENT = "component";

    constexpr const char* J_PR_FAR_FIELD_THETA = "theta";
    constexpr const char* J_PR_FAR_FIELD_PHI = "phi";
    constexpr const char* J_PR_FAR_FIELD_DIR_INITIAL = "initial";
    constexpr const char* J_PR_FAR_FIELD_DIR_FINAL = "final";
    constexpr const char* J_PR_FAR_FIELD_DIR_STEP = "step";

    // domain stuff
    constexpr const char* J_PR_DOMAIN = "domain";
    constexpr const char* J_PR_DOMAIN_TYPE = "type";
    
    constexpr const char* J_PR_DOMAIN_TYPE_TIME = "time";
    constexpr const char* J_PR_DOMAIN_TYPE_FREQ = "frequency";
    constexpr const char* J_PR_DOMAIN_TYPE_TIMEFREQ = "timeFrequency";

    constexpr const char* J_PR_DOMAIN_MAGNITUDE_FILE = "magnitudeFile";

    constexpr const char* J_PR_DOMAIN_TIME_START = "initialTime";
    constexpr const char* J_PR_DOMAIN_TIME_STOP   = "finalTime";
    constexpr const char* J_PR_DOMAIN_TIME_STEP  = "samplingPeriod";
    
    constexpr const char* J_PR_DOMAIN_FREQ_START = "initialFrequency";
    constexpr const char* J_PR_DOMAIN_FREQ_STOP   = "finalFrequency";
    constexpr const char* J_PR_DOMAIN_FREQ_NUMBER  = "numberOfFrequencies";
    constexpr const char* J_PR_DOMAIN_FREQ_SPACING  = "frequencySpacing";
    constexpr const char* J_PR_DOMAIN_FREQ_SPACING_LINEAR  = "linear";
    constexpr const char* J_PR_DOMAIN_FREQ_SPACING_LOGARITHMIC  = "logarithmic";

#endif

} // namespace smbjson_labels_m