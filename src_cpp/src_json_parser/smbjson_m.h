#ifndef SMBJSON_M_H
#define SMBJSON_M_H

#include "../src_main_pub/nfde_types.h"
#include "smbjson_labels_m.h"
#include "mesh_m.h"
#include "parser_tools_m.h"
#include "idchildtable_m.h"
#include "Report_m.h"
#include "json_module.h"
#include "json_kinds.h"
#include "conformal_types.h"

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <numeric>

#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

namespace jmod = json_module;
namespace jlbl = smbjson_labels_m;
namespace NFDE = NFDETypes_m;
namespace Mesh = mesh_m;
namespace Cell = cells_m;
namespace IdTable = idchildtable_m;
namespace Pt = parser_tools_m;
namespace Report = Report_m;
namespace Conf = conformal_types_m;

namespace smbjson {

class parser_t {
private:
    std::string filename;
    jmod::json_file* jsonfile = nullptr;
    jmod::json_core* core = nullptr;
    jmod::json_value* root = nullptr;
    Mesh::mesh_t mesh;
    IdTable::IdChildTable_t matTable, elementTable;

    struct thinwiretermination_t {
        int terminationType = 0;
        double r = 0.0, l = 0.0, c = 0.0;
    };
    struct generator_description_t {
        std::string srctype, srcfile;
        double multiplier = 1.0;
    };
    struct materialAssociation_t {
        std::string name;
        int materialId = 0;
        std::vector<int> elementIds;
        std::string matAssType;
        int initialTerminalId = -1;
        int endTerminalId = -1;
        int initialConnectorId = -1;
        int endConnectorId = -1;
        int containedWithinElementId = -1;
        std::vector<double> totalResistance;
        bool hasTotalResistance = false;
    };
    struct domain_t {
        double tstart = 0.0, tstop = 0.0, tstep = 0.0;
        double fstart = 0.0, fstop = 0.0, fstep = 0.0;
        std::string filename;
        int type1 = NFDE::NP_T1_PLAIN;
        int type2 = NFDE::NP_T2_TIME;
        bool isLogarithmicFrequencySpacing = false;
    };

public:
    bool isInitialized = false;
    parser_t(const std::string& filename);
    NFDE::Parseador_t readProblemDescription();
    Mesh::mesh_t readMesh();

private:
    std::string readAdditionalArguments();
    NFDE::NFDEGeneral_t readGeneral();
    NFDE::MatrizMedios_t readMediaMatrix();
    NFDE::Desplazamiento_t readGrid();
    NFDE::Frontera_t readBoundary();
    void readBackgroundMaterial(NFDE::Materials_t& mats);
    NFDE::PECRegions_t readPECRegions();
    NFDE::PECRegions_t readPMCRegions();
    NFDE::PECRegions_t buildPECPMCRegions(const std::string& matType);
    NFDE::ConformalPECRegions_t readConformalRegions();
    NFDE::DielectricRegions_t readDielectricRegions();
    NFDE::LossyThinSurfaces_t readLossyThinSurfaces();
    NFDE::PlaneWaves_t readPlanewaves();
    NFDE::NodSource_t readNodalSources();
    NFDE::Sondas_t readProbes();
    NFDE::MasSondas_t readMoreProbes();
    NFDE::BloqueProbes_t readBlockProbes();
    NFDE::VolProbes_t readVolumicProbes();
    void readThinWires(NFDE::ThinWires_t& res, NFDE::MasSondas_t& sonda);
    NFDE::ThinSlots_t readThinSlots();
    bool getLogicalAt(const jmod::json_value* val, const std::string& key, bool default_val, bool* foundOut = nullptr);
    int getIntAt(const jmod::json_value* val, const std::string& key, int default_val, bool* foundOut = nullptr);
    std::vector<int> getIntsAt(const jmod::json_value* val, const std::string& key, bool* foundOut = nullptr);
    double getRealAt(const jmod::json_value* val, const std::string& key, double default_val, bool* foundOut = nullptr);
    std::vector<double> getRealsAt(const jmod::json_value* val, const std::string& key, bool* foundOut = nullptr);
    std::vector<std::vector<double>> getMatrixAt(const jmod::json_value* val, const std::string& key, bool* foundOut = nullptr);
    std::string getStrAt(const jmod::json_value* val, const std::string& key, const std::string& default_val = "", bool* foundOut = nullptr);
    bool existsAt(const jmod::json_value* val, const std::string& key);
    int dimensionAt(const jmod::json_value* val, const std::string& key);
    domain_t getDomain(const jmod::json_value* place, const std::string& path);
    std::vector<materialAssociation_t> getMaterialAssociations(
        const std::vector<std::string>& materialTypes,
        const std::vector<std::string>& elementLabels);
    materialAssociation_t parseMaterialAssociation(const jmod::json_value* matAss);
    void matAssToCoords(const materialAssociation_t& mA, std::vector<NFDE::coords_t>& res, int cellType);
    std::string buildTagName(int matId, int elementId);
    std::vector<jmod::json_value*> jsonValueFilterByKeyValue(
        const jmod::json_value* place, const std::string& key, const std::string& value);
    std::vector<jmod::json_value*> jsonValueFilterByKeyValues(
        const jmod::json_value* place, const std::string& key, const std::vector<std::string>& values);
    void addCoordinates(Mesh::mesh_t& mesh);
    void addElements(Mesh::mesh_t& mesh);
    std::vector<Cell::cell_interval_t> readCellIntervals(const jmod::json_value* place, const std::string& path);
    std::vector<Conf::triangle_t> readTriangles(const jmod::json_value* place, const std::string& path);
    std::vector<Cell::cell_interval_t> getSingleVolumeInElementsIds(const jmod::json_value* pw);
    NFDE::FronteraPML_t readPMLProperties(const std::string& p);
    int labelToBoundaryPlace(const std::string& str);
    int labelToBoundaryType(const std::string& str);
    void fillDielectricsOfCellType(std::vector<NFDE::Dielectric_t>& res, int cellType);
    NFDE::Dielectric_t readDielectric(const materialAssociation_t& mA, int cellType);
    NFDE::Dielectric_t readLumped(const materialAssociation_t& mA, int cellType);
    bool containsCellRegionsWithType(const materialAssociation_t& mA, int cellType);
    NFDE::LossyThinSurface_t readLossyThinSurface(const materialAssociation_t& mA);
    NFDE::LossyThinSurfaces_t emptyLossyThinSurfaces();
    NFDE::PlaneWave_t readPlanewave(const jmod::json_value* pw);
    NFDE::Curr_Field_Src_t readField(const jmod::json_value* jns);
    NFDE::abstractSonda_t readFarFieldProbe(const jmod::json_value* p);
    void readDirection(const jmod::json_value* p, const std::string& label, double& initial, double& final, double& step);
    bool isMoreProbe(const jmod::json_value* p);
    bool isLineProbe(const jmod::json_value* p);
    bool isPointProbe(const jmod::json_value* p);
    NFDE::MasSonda_t readLineProbe(const jmod::json_value* p);
    NFDE::MasSonda_t readPointProbe(const jmod::json_value* p);
    std::vector<char> buildDirLabels(const jmod::json_value* dirLabelsPtr);
    void setDomain(NFDE::MasSonda_t& res, const domain_t& domain);
    int strToFieldType(const std::string& fieldLabel, char dirLabel);
    NFDE::BloqueProbe_t readBlockProbe(const jmod::json_value* bp);
    void setDomainBlock(NFDE::BloqueProbe_t& res, const domain_t& domain);
    NFDE::VolProbe_t readVolProbe(const jmod::json_value* p);
    int buildVolProbeType(const std::string& fieldType, const std::string& component);
    void setDomainVol(NFDE::VolProbe_t& res, const domain_t& domain);
    NFDE::VolProbes_t buildNoVolProbes();
    NFDE::ThinSlot_t readThinSlot(const materialAssociation_t& mA);
    void coordsToThinSlotComp(const std::vector<NFDE::coords_t>& cs, std::vector<NFDE::ThinSlotComp_t>& tc);
    NFDE::ThinSlotComp_t buildBaseThinSlotComponent(const NFDE::coords_t& cs);
    NFDE::ThinWire_t readThinWire(const materialAssociation_t& cable);
    int getOrAssignNodeIndex(int coordId);
    std::vector<generator_description_t> readGeneratorOnThinWire(
        const std::vector<Cell::linel_t>& linels, const std::vector<int>& plineElemIds);
    int orientFieldFromGenerator(const std::vector<Cell::linel_t>& linels, int position);
    int findSourcePositionInLinels(const std::vector<int>& srcElemIds, const std::vector<Cell::linel_t>& linels);
    thinwiretermination_t readThinWireTermination(const jmod::json_value* terminal);
    int strToTerminationType(const std::string& label);
    bool isThinWire(const materialAssociation_t& mA);
    NFDE::MasSonda_t readWireProbe(const jmod::json_value* p);
    void setDomainOfWireProbe(NFDE::MasSonda_t& res, const domain_t& domain);
    int getSegmentNdWhichMatchesCoord(int coordId, const Mesh::coordinate_t& probe_coord);
    void appendLogSufix(std::string& fn);
    int getNPDomainType(const std::string& typeLabel, bool hasTransferFunction);
    void showLabelNotFoundError(const std::string& label);
    bool isMaterialIdOfType(int matId, const std::string& matType);
    bool isAssociatedWithMaterial(const jmod::json_value* mAPtr, const std::string& materialType);
    bool isAssociatedWithElementLabel(const jmod::json_value* mAPtr, const std::vector<std::string>& elementLabels);
    std::string adaptName(const std::string& str);
    void checkIsValidName(const std::string& str);
};

} // namespace smbjson

#endif // SMBJSON_M_H
