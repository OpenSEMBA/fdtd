#ifndef SMBJSON_M_H
#define SMBJSON_M_H

#include "../main/nfde_types.h"
#include "smbjson_labels_m.h"
#include "mesh_m.h"
#include "parser_tools_m.h"
#include "id_map_m.h"
#include "Report_m.h"
#include <nlohmann/json.hpp>
#include "conformal_types.h"
#ifdef CompileWithMTLN
#include "mtln_types.h"
#endif

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <numeric>

#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

namespace jlbl = smbjson_labels_m;
namespace NFDE = NFDETypes_m;
namespace Mesh = mesh_m;
namespace Cell = cells_m;
namespace Pt = parser_tools_m;
namespace Report = Report_m;
namespace Conf = conformal_types_m;

namespace smbjson {

class parser_t {
public:
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

private:
    std::string filename;
    nlohmann::json rootJson;
    Mesh::mesh_t mesh;
    id_map_m::id_map_t matTable, elementTable;

    struct thinwiretermination_t {
        int terminationType = 0;
        double r = 0.0, l = 0.0, c = 0.0;
    };
    struct generator_description_t {
        std::string srctype, srcfile;
        double multiplier = 1.0;
    };
    struct domain_t {
        double tstart = 0.0, tstop = 0.0, tstep = 0.0;
        double fstart = 0.0, fstop = 0.0, fstep = 0.0;
        std::string filename;
        bool hasTransferFunction = false;
        int type1 = NFDE::NP_T1_PLAIN;
        int type2 = NFDE::NP_T2_TIME;
        bool isLogarithmicFrequencySpacing = false;
    };

    // Transient state used while parsing thin wires
    int wireNGlobal_ = 0;
    int wireNNodes_ = 0;
    std::vector<int> wireNodeCoordIds_;
    std::vector<int> wireNodeNodeIdx_;
    std::vector<materialAssociation_t> wireMAs_;
    NFDE::ThinWires_t* wireRes_ = nullptr;

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
    bool getLogicalAt(const nlohmann::json* val, const std::string& key, bool default_val, bool* foundOut = nullptr);
    int getIntAt(const nlohmann::json* val, const std::string& key, int default_val, bool* foundOut = nullptr);
    std::vector<int> getIntsAt(const nlohmann::json* val, const std::string& key, bool* foundOut = nullptr);
    double getRealAt(const nlohmann::json* val, const std::string& key, double default_val, bool* foundOut = nullptr);
    std::vector<double> getRealsAt(const nlohmann::json* val, const std::string& key, bool* foundOut = nullptr);
    std::vector<std::vector<double>> getMatrixAt(const nlohmann::json* val, const std::string& key, bool* foundOut = nullptr);
    std::string getStrAt(const nlohmann::json* val, const std::string& key, const std::string& default_val = "", bool* foundOut = nullptr);
    domain_t getDomain(const nlohmann::json* place, const std::string& path);
    std::vector<materialAssociation_t> getMaterialAssociations(
        const std::vector<std::string>& materialTypes,
        const std::vector<std::string>& elementLabels);
    materialAssociation_t parseMaterialAssociation(const nlohmann::json* matAss);
    void matAssToCoords(const materialAssociation_t& mA, std::vector<NFDE::coords_t>& res, int cellType);
    std::string buildTagName(int matId, int elementId);
    std::vector<const nlohmann::json*> jsonValueFilterByKeyValue(
        const nlohmann::json* place, const std::string& key, const std::string& value);
    std::vector<const nlohmann::json*> jsonValueFilterByKeyValues(
        const nlohmann::json* place, const std::string& key, const std::vector<std::string>& values);
    void addCoordinates(Mesh::mesh_t& mesh);
    void addElements(Mesh::mesh_t& mesh);
    std::vector<Cell::cell_interval_t> readCellIntervals(const nlohmann::json* place, const std::string& path);
    std::vector<Conf::triangle_t> readTriangles(const nlohmann::json* place, const std::string& path);
    std::vector<Cell::cell_interval_t> getSingleVolumeInElementsIds(const nlohmann::json* pw);
    NFDE::FronteraPML_t readPMLProperties(const std::string& p);
    int labelToBoundaryPlace(const std::string& str);
    int labelToBoundaryType(const std::string& str);
    void fillDielectricsOfCellType(std::vector<NFDE::Dielectric_t>& res, int cellType);
    NFDE::Dielectric_t readDielectric(const materialAssociation_t& mA, int cellType);
    NFDE::Dielectric_t readLumped(const materialAssociation_t& mA, int cellType);
    bool containsCellRegionsWithType(const materialAssociation_t& mA, int cellType);
    NFDE::LossyThinSurface_t readLossyThinSurface(const materialAssociation_t& mA);
    NFDE::LossyThinSurfaces_t emptyLossyThinSurfaces();
    NFDE::PlaneWave_t readPlanewave(const nlohmann::json* pw);
    NFDE::Curr_Field_Src_t readField(const nlohmann::json* jns);
    NFDE::abstractSonda_t readFarFieldProbe(const nlohmann::json* p);
    void readDirection(const nlohmann::json* p, const std::string& label, double& initial, double& final, double& step);
    bool isMoreProbe(const nlohmann::json* p);
    bool isLineProbe(const nlohmann::json* p);
    bool isPointProbe(const nlohmann::json* p);
    NFDE::MasSonda_t readLineProbe(const nlohmann::json* p);
    NFDE::MasSonda_t readPointProbe(const nlohmann::json* p);
    std::vector<char> buildDirLabels(const nlohmann::json* dirLabelsPtr);
    void setDomain(NFDE::MasSonda_t& res, const domain_t& domain);
    int strToFieldType(const std::string& fieldLabel, char dirLabel);
    NFDE::BloqueProbe_t readBlockProbe(const nlohmann::json* bp);
    void setDomainBlock(NFDE::BloqueProbe_t& res, const domain_t& domain);
    NFDE::VolProbe_t readVolProbe(const nlohmann::json* p);
    int buildVolProbeType(const std::string& fieldType, const std::string& component);
    void setDomainVol(NFDE::VolProbe_t& res, const domain_t& domain);
    NFDE::VolProbes_t buildNoVolProbes();
    NFDE::ThinSlot_t readThinSlot(const materialAssociation_t& mA);
    void coordsToThinSlotComp(const std::vector<NFDE::coords_t>& cs, std::vector<NFDE::ThinSlotComp_t>& tc);
    NFDE::ThinSlotComp_t buildBaseThinSlotComponent(const NFDE::coords_t& cs);
    NFDE::ThinWire_t readThinWire(const materialAssociation_t& cable);
    int getOrAssignNodeIndex(int coordId);
    std::vector<generator_description_t> readGeneratorOnThinWire(
        const std::vector<Mesh::linel_t>& linels, const std::vector<int>& plineElemIds);
    int orientFieldFromGenerator(const std::vector<Mesh::linel_t>& linels, int position);
    int findSourcePositionInLinels(const std::vector<int>& srcElemIds, const std::vector<Mesh::linel_t>& linels);
    thinwiretermination_t readThinWireTermination(const nlohmann::json* terminal);
    int strToTerminationType(const std::string& label);
    bool isThinWire(const materialAssociation_t& mA);
    NFDE::MasSonda_t readWireProbe(const nlohmann::json* p);
    void setDomainOfWireProbe(NFDE::MasSonda_t& res, const domain_t& domain);
    int getSegmentNdWhichMatchesCoord(int coordId, const Mesh::coordinate_t& probe_coord);
    void appendLogSufix(std::string& fn);
    int getNPDomainType(const std::string& typeLabel, bool hasTransferFunction);
    void showLabelNotFoundError(const std::string& label);
    bool isMaterialIdOfType(int matId, const std::string& matType);
    bool isAssociatedWithMaterial(const nlohmann::json* mAPtr, const std::string& materialType);
    bool isAssociatedWithElementLabel(const nlohmann::json* mAPtr, const std::vector<std::string>& elementLabels);
    std::string adaptName(const std::string& str);
    void checkIsValidName(const std::string& str);

#ifdef CompileWithMTLN
    class MtlnReader;
    mtln_types_m::mtln_t readMTLN();
#endif
};

} // namespace smbjson

#endif // SMBJSON_M_H
