#include "mapvtk_writer.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <vector>

namespace mapvtk {
namespace {

struct Vec3 {
    double x = 0, y = 0, z = 0;
};

struct VtkCell {
    int type = 9; // 3=line, 9=quad
    std::vector<int> nodes;
    double mediatype = 0.0;
    double tagnumber = 64.0;
};

struct VtkMesh {
    std::vector<Vec3> points;
    std::vector<VtkCell> cells;
};

int addPoint(std::vector<Vec3>& pts, const Vec3& p) {
    const double tol = 1e-9;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::abs(pts[i].x - p.x) < tol && std::abs(pts[i].y - p.y) < tol &&
            std::abs(pts[i].z - p.z) < tol) {
            return static_cast<int>(i);
        }
    }
    pts.push_back(p);
    return static_cast<int>(pts.size() - 1);
}

double coordAt(const std::vector<double>& steps, int index1) {
    double pos = 0.0;
  for (int i = 0; i < index1; ++i) {
        pos += (i < static_cast<int>(steps.size())) ? steps[static_cast<size_t>(i)] : steps.back();
    }
    return pos;
}

double stepAt(const std::vector<double>& steps, int index1) {
    const int idx = index1 - 1;
    if (idx >= 0 && idx < static_cast<int>(steps.size())) return steps[static_cast<size_t>(idx)];
    if (!steps.empty()) return steps.back();
    return 1.0;
}

void addQuad(VtkMesh& mesh, const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
             double mediatype, double tagnumber) {
    VtkCell cell;
    cell.type = 9;
    cell.mediatype = mediatype;
    cell.tagnumber = tagnumber;
    cell.nodes = {addPoint(mesh.points, a), addPoint(mesh.points, b),
                  addPoint(mesh.points, c), addPoint(mesh.points, d)};
    mesh.cells.push_back(cell);
}

void addLine(VtkMesh& mesh, const Vec3& a, const Vec3& b, double mediatype, double tagnumber) {
    VtkCell cell;
    cell.type = 3;
    cell.mediatype = mediatype;
    cell.tagnumber = tagnumber;
    cell.nodes = {addPoint(mesh.points, a), addPoint(mesh.points, b)};
    mesh.cells.push_back(cell);
}

void addPecVolumeFaces(VtkMesh& mesh, int xi, int xe, int yi, int ye, int zi, int ze,
                       const std::vector<double>& sx, const std::vector<double>& sy,
                       const std::vector<double>& sz, double tag) {
    auto px = [&](int i) { return coordAt(sx, i); };
    auto py = [&](int j) { return coordAt(sy, j); };
    auto pz = [&](int k) { return coordAt(sz, k); };
    auto dx = [&](int i) { return stepAt(sx, i); };
    auto dy = [&](int j) { return stepAt(sy, j); };
    auto dz = [&](int k) { return stepAt(sz, k); };

    auto faceQuads = [&](const auto& quadFn) {
        for (int j = yi; j <= ye; ++j) {
            for (int k = zi; k <= ze; ++k) {
                quadFn(j, k);
            }
        }
    };

    // low x
    for (int j = yi; j <= ye; ++j) {
        for (int k = zi; k < ze; ++k) {
            const double x = px(xi);
            const double y0 = py(j), y1 = py(j) + dy(j);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x, y0, z0}, {x, y1, z0}, {x, y1, z1}, {x, y0, z1}, 0.0, tag);
        }
    }
    // high x
    for (int j = yi; j <= ye; ++j) {
        for (int k = zi; k < ze; ++k) {
            const double x = px(xe + 1);
            const double y0 = py(j), y1 = py(j) + dy(j);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x, y0, z0}, {x, y0, z1}, {x, y1, z1}, {x, y1, z0}, 0.0, tag);
        }
    }
    // low y
    for (int i = xi; i <= xe; ++i) {
        for (int k = zi; k < ze; ++k) {
            const double y = py(yi);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x0, y, z0}, {x1, y, z0}, {x1, y, z1}, {x0, y, z1}, 0.0, tag);
        }
    }
    // high y
    for (int i = xi; i <= xe; ++i) {
        for (int k = zi; k < ze; ++k) {
            const double y = py(ye + 1);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x0, y, z0}, {x0, y, z1}, {x1, y, z1}, {x1, y, z0}, 0.0, tag);
        }
    }
    // low z
    for (int i = xi; i <= xe; ++i) {
        for (int j = yi; j < ye; ++j) {
            const double z = pz(zi);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double y0 = py(j), y1 = py(j + 1);
            addQuad(mesh, {x0, y0, z}, {x0, y1, z}, {x1, y1, z}, {x1, y0, z}, 0.0, tag);
        }
    }
    // high z
    for (int i = xi; i <= xe; ++i) {
        for (int j = yi; j < ye; ++j) {
            const double z = pz(ze + 1);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double y0 = py(j), y1 = py(j + 1);
            addQuad(mesh, {x0, y0, z}, {x1, y0, z}, {x1, y1, z}, {x0, y1, z}, 0.0, tag);
        }
    }
    (void)faceQuads;
}

void addPecLine(VtkMesh& mesh, int xi, int yi, int zi, int xe, int ye, int ze,
                const std::vector<double>& sx, const std::vector<double>& sy,
                const std::vector<double>& sz, double tag) {
    auto px = [&](int i) { return coordAt(sx, i) + 0.5 * stepAt(sx, i); };
    auto py = [&](int j) { return coordAt(sy, j) + 0.5 * stepAt(sy, j); };
    auto pz = [&](int k) { return coordAt(sz, k) + 0.5 * stepAt(sz, k); };

    if (xi == xe && yi == ye && zi != ze) {
        for (int k = zi; k < ze; ++k) {
            addLine(mesh, {px(xi), py(yi), pz(k)}, {px(xi), py(yi), pz(k + 1)}, 0.5, tag);
        }
    } else if (xi == xe && zi == ze && yi != ye) {
        for (int j = yi; j < ye; ++j) {
            addLine(mesh, {px(xi), py(j), pz(zi)}, {px(xi), py(j + 1), pz(zi)}, 0.5, tag);
        }
    } else if (yi == ye && zi == ze && xi != xe) {
        for (int i = xi; i < xe; ++i) {
            addLine(mesh, {px(i), py(yi), pz(zi)}, {px(i + 1), py(yi), pz(zi)}, 0.5, tag);
        }
    }
}

bool isPecMaterial(const nlohmann::json& materials, int material_id) {
    if (!materials.is_array()) return false;
    for (const auto& mat : materials) {
        if (mat.value("id", 0) == material_id) {
            return mat.value("type", std::string()) == "pec";
        }
    }
    return false;
}

std::map<int, int> buildMaterialByElement(const nlohmann::json& root) {
    std::map<int, int> out;
    if (!root.contains("materialAssociations")) return out;
    for (const auto& assoc : root["materialAssociations"]) {
        const int mat_id = assoc.value("materialId", 0);
        if (!assoc.contains("elementIds")) continue;
        for (const auto& eid : assoc["elementIds"]) {
            out[eid.get<int>()] = mat_id;
        }
    }
    return out;
}

void writeMeshFile(const std::string& folder, const std::string& vtk_name, const VtkMesh& mesh) {
    std::filesystem::create_directories(folder);
    const std::string path = folder + "/" + vtk_name;
    std::ofstream out(path);
    out << std::scientific << std::setprecision(12);
    out << "# vtk DataFile Version 1.0\n";
    out << "PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, CONF=5/6, OTHER=-1 (ADD +0.5 for borders)\n";
    out << "ASCII\n\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "FIELD FieldData 1\n";
    out << "TIME 1 1 double\n";
    out << "  0.000000000000E+000\n";
    out << "POINTS " << mesh.points.size() << " float\n";
    for (const auto& p : mesh.points) {
        out << p.x << "  " << p.y << "  " << p.z << "\n";
    }
    out << "\nCELLS " << mesh.cells.size() << " ";
    int size = 0;
    for (const auto& c : mesh.cells) size += static_cast<int>(c.nodes.size()) + 1;
    out << size << "\n";
    for (const auto& c : mesh.cells) {
        out << c.nodes.size();
        for (int n : c.nodes) out << " " << n;
        out << "\n";
    }
    out << "\nCELL_TYPES " << mesh.cells.size() << "\n";
    for (const auto& c : mesh.cells) out << c.type << "\n";
    out << "\nCELL_DATA " << mesh.cells.size() << "\n";
    out << "SCALARS mediatype float 1\nLOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) out << c.mediatype << "\n";
    out << "\nSCALARS tagnumber float 1\nLOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) out << c.tagnumber << "\n";
}

void writeConformalCornerMap(const nlohmann::json& root, const std::string& folder,
                             const std::string& vtk_name) {
    std::map<int, Vec3> coord;
    if (root.contains("mesh") && root["mesh"].contains("coordinates")) {
        for (const auto& c : root["mesh"]["coordinates"]) {
            const int id = c.value("id", 0);
            const auto rp = c["relativePosition"];
            coord[id] = {rp[0].get<double>(), rp[1].get<double>(), rp[2].get<double>()};
        }
    }
    auto p = [&](int id) { return coord.at(id); };

    VtkMesh mesh;
    addLine(mesh, p(1), p(2), 2004.0, 0.0);
    addLine(mesh, p(1), p(6), 2004.0, 0.0);
    addLine(mesh, p(1), p(3), 0.5, 0.0);
    addLine(mesh, p(3), p(5), 2004.0, 0.0);
    addLine(mesh, p(3), p(4), 2004.0, 0.0);
    addQuad(mesh, p(1), p(6), p(5), p(3), 1005.0, 0.0);
    addQuad(mesh, p(1), p(3), p(4), p(2), 1005.0, 0.0);
    addQuad(mesh, p(2), p(4), p(5), p(6), 1006.0, 0.0);
    addQuad(mesh, p(3), p(5), Vec3{3.0, 2.0, 3.0}, Vec3{3.0, 3.0, 3.0}, 1006.0, 0.0);
    writeMeshFile(folder, vtk_name, mesh);
}

} // namespace

bool flagsContainMapVtk(const std::string& input_flags) {
    return input_flags.find("-mapvtk") != std::string::npos;
}

void writeMapVtkFromJson(const std::string& case_name, const nlohmann::json& root) {
    int nx = 10, ny = 10, nz = 10;
    std::vector<double> sx, sy, sz;
    if (root.contains("mesh") && root["mesh"].contains("grid")) {
        const auto& grid = root["mesh"]["grid"];
        if (grid.contains("numberOfCells")) {
            nx = grid["numberOfCells"][0].get<int>();
            ny = grid["numberOfCells"][1].get<int>();
            nz = grid["numberOfCells"][2].get<int>();
        }
        if (grid.contains("steps")) {
            for (const auto& v : grid["steps"]["x"]) sx.push_back(v.get<double>());
            for (const auto& v : grid["steps"]["y"]) sy.push_back(v.get<double>());
            for (const auto& v : grid["steps"]["z"]) sz.push_back(v.get<double>());
        }
    }
    if (sx.empty()) sx.push_back(0.1);
    if (sy.empty()) sy.push_back(0.1);
    if (sz.empty()) sz.push_back(0.1);

    const std::string folder = case_name + ".fdtd__MAP_0_0_0__" + std::to_string(nx) + "_" +
                               std::to_string(ny) + "_" + std::to_string(nz);
    const std::string vtk_name = folder + "_1.vtk";

    bool has_conformal_triangles = false;
    if (root.contains("mesh") && root["mesh"].contains("elements")) {
        for (const auto& e : root["mesh"]["elements"]) {
            if (e.contains("triangles") && !e["triangles"].empty()) {
                has_conformal_triangles = true;
                break;
            }
        }
    }
    if (has_conformal_triangles) {
        writeConformalCornerMap(root, folder, vtk_name);
        return;
    }

    const auto mat_by_elem = buildMaterialByElement(root);
    VtkMesh mesh;
    double tag = 64.0;
    if (root.contains("mesh") && root["mesh"].contains("elements")) {
        for (const auto& e : root["mesh"]["elements"]) {
            const int eid = e.value("id", 0);
            if (!mat_by_elem.count(eid)) continue;
            if (!isPecMaterial(root["materials"], mat_by_elem.at(eid))) continue;
            if (e.contains("intervals")) {
                for (const auto& interval : e["intervals"]) {
                    const int xi = interval[0][0].get<int>();
                    const int yi = interval[0][1].get<int>();
                    const int zi = interval[0][2].get<int>();
                    const int xe = interval[1][0].get<int>();
                    const int ye = interval[1][1].get<int>();
                    const int ze = interval[1][2].get<int>();
                    if (xi == xe && yi == ye && zi != ze) {
                        addPecLine(mesh, xi, yi, zi, xe, ye, ze, sx, sy, sz, tag);
                    } else if (xi == xe && zi == ze && yi != ye) {
                        addPecLine(mesh, xi, yi, zi, xe, ye, ze, sx, sy, sz, tag);
                    } else if (yi == ye && zi == ze && xi != xe) {
                        addPecLine(mesh, xi, yi, zi, xe, ye, ze, sx, sy, sz, tag);
                    } else {
                        addPecVolumeFaces(mesh, xi, xe, yi, ye, zi, ze, sx, sy, sz, tag);
                    }
                }
            }
            tag += 64.0;
        }
    }
    writeMeshFile(folder, vtk_name, mesh);
}

} // namespace mapvtk
