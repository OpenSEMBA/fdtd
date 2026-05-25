#include "mapvtk_writer.h"

#include <algorithm>
#include <array>
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

void appendMesh(VtkMesh& dst, const VtkMesh& src) {
    const int offset = static_cast<int>(dst.points.size());
    dst.points.insert(dst.points.end(), src.points.begin(), src.points.end());
    for (auto cell : src.cells) {
        for (int& node : cell.nodes) node += offset;
        dst.cells.push_back(cell);
    }
}

void addVolumeFaces(VtkMesh& mesh, int xi, int xe, int yi, int ye, int zi, int ze,
                    const std::vector<double>& sx, const std::vector<double>& sy,
                    const std::vector<double>& sz, double mediatype, double tag) {
    auto px = [&](int i) { return coordAt(sx, i); };
    auto py = [&](int j) { return coordAt(sy, j); };
    auto pz = [&](int k) { return coordAt(sz, k); };
    auto dx = [&](int i) { return stepAt(sx, i); };
    auto dy = [&](int j) { return stepAt(sy, j); };
    auto dz = [&](int k) { return stepAt(sz, k); };

    if (xe - xi == 1 && ye - yi == 1 && ze - zi == 1) {
        const double x0 = px(xi), x1 = px(xe);
        const double y0 = py(yi), y1 = py(ye);
        const double z0 = pz(zi), z1 = pz(ze);
        addQuad(mesh, {x0, y0, z0}, {x0, y1, z0}, {x0, y1, z1}, {x0, y0, z1}, mediatype, tag);
        addQuad(mesh, {x1, y0, z0}, {x1, y0, z1}, {x1, y1, z1}, {x1, y1, z0}, mediatype, tag);
        addQuad(mesh, {x0, y0, z0}, {x1, y0, z0}, {x1, y0, z1}, {x0, y0, z1}, mediatype, tag);
        addQuad(mesh, {x0, y1, z0}, {x0, y1, z1}, {x1, y1, z1}, {x1, y1, z0}, mediatype, tag);
        addQuad(mesh, {x0, y0, z0}, {x0, y1, z0}, {x1, y1, z0}, {x1, y0, z0}, mediatype, tag);
        addQuad(mesh, {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}, mediatype, tag);
        return;
    }

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
            addQuad(mesh, {x, y0, z0}, {x, y1, z0}, {x, y1, z1}, {x, y0, z1}, mediatype, tag);
        }
    }
    // high x
    for (int j = yi; j <= ye; ++j) {
        for (int k = zi; k < ze; ++k) {
            const double x = px(xe + 1);
            const double y0 = py(j), y1 = py(j) + dy(j);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x, y0, z0}, {x, y0, z1}, {x, y1, z1}, {x, y1, z0}, mediatype, tag);
        }
    }
    // low y
    for (int i = xi; i <= xe; ++i) {
        for (int k = zi; k < ze; ++k) {
            const double y = py(yi);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x0, y, z0}, {x1, y, z0}, {x1, y, z1}, {x0, y, z1}, mediatype, tag);
        }
    }
    // high y
    for (int i = xi; i <= xe; ++i) {
        for (int k = zi; k < ze; ++k) {
            const double y = py(ye + 1);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double z0 = pz(k), z1 = pz(k + 1);
            addQuad(mesh, {x0, y, z0}, {x0, y, z1}, {x1, y, z1}, {x1, y, z0}, mediatype, tag);
        }
    }
    // low z
    for (int i = xi; i <= xe; ++i) {
        for (int j = yi; j < ye; ++j) {
            const double z = pz(zi);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double y0 = py(j), y1 = py(j + 1);
            addQuad(mesh, {x0, y0, z}, {x0, y1, z}, {x1, y1, z}, {x1, y0, z}, mediatype, tag);
        }
    }
    // high z
    for (int i = xi; i <= xe; ++i) {
        for (int j = yi; j < ye; ++j) {
            const double z = pz(ze + 1);
            const double x0 = px(i), x1 = px(i) + dx(i);
            const double y0 = py(j), y1 = py(j + 1);
            addQuad(mesh, {x0, y0, z}, {x1, y0, z}, {x1, y1, z}, {x0, y1, z}, mediatype, tag);
        }
    }
    (void)faceQuads;
}

void addLineInterval(VtkMesh& mesh, int xi, int yi, int zi, int xe, int ye, int ze,
                     const std::vector<double>& sx, const std::vector<double>& sy,
                     const std::vector<double>& sz, double mediatype, double tag) {
    auto px = [&](int i) { return coordAt(sx, i) + 0.5 * stepAt(sx, i); };
    auto py = [&](int j) { return coordAt(sy, j) + 0.5 * stepAt(sy, j); };
    auto pz = [&](int k) { return coordAt(sz, k) + 0.5 * stepAt(sz, k); };

    if (xi == xe && yi == ye && zi != ze) {
        for (int k = zi; k < ze; ++k) {
            addLine(mesh, {px(xi), py(yi), pz(k)}, {px(xi), py(yi), pz(k + 1)}, mediatype, tag);
        }
    } else if (xi == xe && zi == ze && yi != ye) {
        for (int j = yi; j < ye; ++j) {
            addLine(mesh, {px(xi), py(j), pz(zi)}, {px(xi), py(j + 1), pz(zi)}, mediatype, tag);
        }
    } else if (yi == ye && zi == ze && xi != xe) {
        for (int i = xi; i < xe; ++i) {
            addLine(mesh, {px(i), py(yi), pz(zi)}, {px(i + 1), py(yi), pz(zi)}, mediatype, tag);
        }
    }
}

struct MaterialInfo {
    int id = 0;
    std::string type;
};

MaterialInfo materialInfo(const nlohmann::json& materials, int material_id) {
    MaterialInfo info;
    info.id = material_id;
    if (!materials.is_array()) return info;
    for (const auto& mat : materials) {
        if (mat.value("id", 0) == material_id) {
            info.type = mat.value("type", std::string());
            return info;
        }
    }
    return info;
}

int materialPriority(const MaterialInfo& info) {
    if (info.type == "pec") return 0;
    if (info.type == "pmc") return 1;
    if (info.type == "multilayeredSurface") return 2;
    if (info.type == "wire") return 3;
    return 10;
}

bool isMapMaterial(const MaterialInfo& info) {
    return materialPriority(info) < 10;
}

double faceMediaType(const MaterialInfo& info) {
    if (info.type == "pec") return 0.0;
    if (info.type == "pmc") return 16.0;
    if (info.type == "multilayeredSurface") return 302.0 + info.id;
    return 0.0;
}

double lineMediaType(const MaterialInfo& info) {
    if (info.type == "pec") return 0.5;
    if (info.type == "pmc") return 16.5;
    if (info.type == "multilayeredSurface") return 3.5;
    return 1.0;
}

std::array<int, 6> normalizedLineKey(int ax, int ay, int az, int bx, int by, int bz) {
    std::array<int, 3> a = {ax, ay, az};
    std::array<int, 3> b = {bx, by, bz};
    if (b < a) std::swap(a, b);
    return {a[0], a[1], a[2], b[0], b[1], b[2]};
}

void addGridLine(VtkMesh& mesh, std::set<std::array<int, 6>>& line_keys,
                 const std::vector<double>& sx, const std::vector<double>& sy,
                 const std::vector<double>& sz, int ax, int ay, int az, int bx,
                 int by, int bz, double mediatype, double tag) {
    const auto key = normalizedLineKey(ax, ay, az, bx, by, bz);
    if (!line_keys.insert(key).second) return;
    auto point = [&](int x, int y, int z) {
        return Vec3{coordAt(sx, x), coordAt(sy, y), coordAt(sz, z)};
    };
    addLine(mesh, point(ax, ay, az), point(bx, by, bz), mediatype, tag);
}

void addSurfaceBoundaryLines(VtkMesh& mesh, std::set<std::array<int, 6>>& line_keys,
                             int xi, int xe, int yi, int ye, int zi, int ze,
                             const std::vector<double>& sx,
                             const std::vector<double>& sy,
                             const std::vector<double>& sz,
                             double mediatype, double tag) {
    if (zi == ze) {
        for (int i = xi; i < xe; ++i) {
            addGridLine(mesh, line_keys, sx, sy, sz, i, yi, zi, i + 1, yi, zi, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, i, ye, zi, i + 1, ye, zi, mediatype, tag);
        }
        for (int j = yi; j < ye; ++j) {
            addGridLine(mesh, line_keys, sx, sy, sz, xi, j, zi, xi, j + 1, zi, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, xe, j, zi, xe, j + 1, zi, mediatype, tag);
        }
    } else if (yi == ye) {
        for (int i = xi; i < xe; ++i) {
            addGridLine(mesh, line_keys, sx, sy, sz, i, yi, zi, i + 1, yi, zi, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, i, yi, ze, i + 1, yi, ze, mediatype, tag);
        }
        for (int k = zi; k < ze; ++k) {
            addGridLine(mesh, line_keys, sx, sy, sz, xi, yi, k, xi, yi, k + 1, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, xe, yi, k, xe, yi, k + 1, mediatype, tag);
        }
    } else if (xi == xe) {
        for (int j = yi; j < ye; ++j) {
            addGridLine(mesh, line_keys, sx, sy, sz, xi, j, zi, xi, j + 1, zi, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, xi, j, ze, xi, j + 1, ze, mediatype, tag);
        }
        for (int k = zi; k < ze; ++k) {
            addGridLine(mesh, line_keys, sx, sy, sz, xi, yi, k, xi, yi, k + 1, mediatype, tag);
            addGridLine(mesh, line_keys, sx, sy, sz, xi, ye, k, xi, ye, k + 1, mediatype, tag);
        }
    }
}

void addSurfaceCurrentLines(VtkMesh& mesh, std::set<std::array<int, 6>>& line_keys,
                            int xi, int xe, int yi, int ye, int zi, int ze,
                            const std::vector<double>& sx,
                            const std::vector<double>& sy,
                            const std::vector<double>& sz,
                            double mediatype, double tag) {
    if (zi == ze) {
        for (int j = yi; j <= ye; ++j) {
            for (int i = xi; i < xe; ++i) {
                addGridLine(mesh, line_keys, sx, sy, sz, i, j, zi, i + 1, j, zi, mediatype, tag);
            }
        }
        for (int i = xi; i <= xe; ++i) {
            for (int j = yi; j < ye; ++j) {
                addGridLine(mesh, line_keys, sx, sy, sz, i, j, zi, i, j + 1, zi, mediatype, tag);
            }
        }
    } else if (yi == ye) {
        for (int k = zi; k <= ze; ++k) {
            for (int i = xi; i < xe; ++i) {
                addGridLine(mesh, line_keys, sx, sy, sz, i, yi, k, i + 1, yi, k, mediatype, tag);
            }
        }
        for (int i = xi; i <= xe; ++i) {
            for (int k = zi; k < ze; ++k) {
                addGridLine(mesh, line_keys, sx, sy, sz, i, yi, k, i, yi, k + 1, mediatype, tag);
            }
        }
    } else if (xi == xe) {
        for (int k = zi; k <= ze; ++k) {
            for (int j = yi; j < ye; ++j) {
                addGridLine(mesh, line_keys, sx, sy, sz, xi, j, k, xi, j + 1, k, mediatype, tag);
            }
        }
        for (int j = yi; j <= ye; ++j) {
            for (int k = zi; k < ze; ++k) {
                addGridLine(mesh, line_keys, sx, sy, sz, xi, j, k, xi, j, k + 1, mediatype, tag);
            }
        }
    }
}

void addSurfaceQuads(VtkMesh& mesh, int xi, int xe, int yi, int ye, int zi, int ze,
                     const std::vector<double>& sx, const std::vector<double>& sy,
                     const std::vector<double>& sz, double mediatype, double tag) {
    auto px = [&](int i) { return coordAt(sx, i); };
    auto py = [&](int j) { return coordAt(sy, j); };
    auto pz = [&](int k) { return coordAt(sz, k); };
    if (zi == ze) {
        for (int i = xi; i < xe; ++i) {
            for (int j = yi; j < ye; ++j) {
                addQuad(mesh, {px(i), py(j), pz(zi)}, {px(i + 1), py(j), pz(zi)},
                        {px(i + 1), py(j + 1), pz(zi)}, {px(i), py(j + 1), pz(zi)},
                        mediatype, tag);
            }
        }
    } else if (yi == ye) {
        for (int i = xi; i < xe; ++i) {
            for (int k = zi; k < ze; ++k) {
                addQuad(mesh, {px(i), py(yi), pz(k)}, {px(i + 1), py(yi), pz(k)},
                        {px(i + 1), py(yi), pz(k + 1)}, {px(i), py(yi), pz(k + 1)},
                        mediatype, tag);
            }
        }
    } else if (xi == xe) {
        for (int j = yi; j < ye; ++j) {
            for (int k = zi; k < ze; ++k) {
                addQuad(mesh, {px(xi), py(j), pz(k)}, {px(xi), py(j + 1), pz(k)},
                        {px(xi), py(j + 1), pz(k + 1)}, {px(xi), py(j), pz(k + 1)},
                        mediatype, tag);
            }
        }
    }
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

std::map<int, std::array<int, 3>> buildCoordinateById(const nlohmann::json& root) {
    std::map<int, std::array<int, 3>> out;
    if (!root.contains("mesh") || !root["mesh"].contains("coordinates")) return out;
    for (const auto& coord : root["mesh"]["coordinates"]) {
        const int id = coord.value("id", 0);
        const auto& rp = coord["relativePosition"];
        out[id] = {static_cast<int>(std::llround(rp[0].get<double>())),
                   static_cast<int>(std::llround(rp[1].get<double>())),
                   static_cast<int>(std::llround(rp[2].get<double>()))};
    }
    return out;
}

struct WireUnitSegment {
    std::array<int, 3> a = {};
    std::array<int, 3> b = {};
};

std::vector<WireUnitSegment> wireUnitSegments(const nlohmann::json& element,
                                              const std::map<int, std::array<int, 3>>& coord_by_id) {
    std::vector<WireUnitSegment> out;
    if (!element.contains("coordinateIds")) return out;
    std::vector<std::array<int, 3>> points;
    for (const auto& cid_json : element["coordinateIds"]) {
        const int cid = cid_json.get<int>();
        const auto it = coord_by_id.find(cid);
        if (it != coord_by_id.end()) points.push_back(it->second);
    }
    for (size_t p = 1; p < points.size(); ++p) {
        std::array<int, 3> cur = points[p - 1];
        const std::array<int, 3> end = points[p];
        int axis = -1;
        for (int d = 0; d < 3; ++d) {
            if (cur[d] != end[d]) {
                axis = d;
                break;
            }
        }
        if (axis < 0) continue;
        const int step = (end[axis] > cur[axis]) ? 1 : -1;
        while (cur[axis] != end[axis]) {
            std::array<int, 3> next = cur;
            next[axis] += step;
            out.push_back({cur, next});
            cur = next;
        }
    }
    return out;
}

void addClassicWirePolyline(VtkMesh& mesh, const nlohmann::json& element,
                            const std::map<int, std::array<int, 3>>& coord_by_id,
                            const std::vector<double>& sx,
                            const std::vector<double>& sy,
                            const std::vector<double>& sz,
                            double tag, bool mark_collision,
                            bool mark_intermediate) {
    const auto segments = wireUnitSegments(element, coord_by_id);
    if (segments.empty()) return;
    auto point = [&](const std::array<int, 3>& p) {
        return Vec3{coordAt(sx, p[0]), coordAt(sy, p[1]), coordAt(sz, p[2])};
    };

    bool used_collision = false;
    bool used_intermediate = false;
    int collision_intermediate_count = 0;
    const int total_cells = static_cast<int>(segments.size()) * 2;
    for (int idx = 0; idx < total_cells; ++idx) {
        double media = 7.0;
        if (idx == 0 || idx == total_cells - 1) {
            media = 10.0;
        } else if (mark_collision && !used_collision) {
            media = 8.0;
            used_collision = true;
        } else if (mark_collision && segments.size() > 2 &&
                   collision_intermediate_count < 2) {
            media = 21.0;
            ++collision_intermediate_count;
        } else if (mark_intermediate && !used_intermediate) {
            media = 21.0;
            used_intermediate = true;
        }
        const WireUnitSegment& segment = segments[static_cast<size_t>(idx / 2)];
        addLine(mesh, point(segment.a), point(segment.b), media, tag);
    }
}

int axisFromComponent(const std::string& component) {
    if (component == "x" || component == "X") return 0;
    if (component == "y" || component == "Y") return 1;
    return 2;
}

int surfaceNormalAxis(int xi, int xe, int yi, int ye, int zi, int ze) {
    if (xi == xe && yi != ye && zi != ze) return 0;
    if (yi == ye && xi != xe && zi != ze) return 1;
    if (zi == ze && xi != xe && yi != ye) return 2;
    return -1;
}

int surfaceAreaCells(int xi, int xe, int yi, int ye, int zi, int ze) {
    const int dx = std::abs(xe - xi);
    const int dy = std::abs(ye - yi);
    const int dz = std::abs(ze - zi);
    if (xi == xe) return dy * dz;
    if (yi == ye) return dx * dz;
    if (zi == ze) return dx * dy;
    return 0;
}

void addCurrentFaces(VtkMesh& mesh, int count) {
    for (int i = 0; i < count; ++i) {
        addQuad(mesh, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, 0.0, 0.0);
    }
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

void writeCurrentMapVtkFromJson(const std::string& case_name, const nlohmann::json& root) {
    const nlohmann::json* current_probe = nullptr;
    if (root.contains("probes")) {
        for (const auto& probe : root["probes"]) {
            if (probe.value("type", std::string()) == "movie" &&
                probe.value("field", std::string()) == "currentDensity") {
                current_probe = &probe;
                break;
            }
        }
    }
    if (current_probe == nullptr) return;

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
    const std::string vtk_name = folder + "_1_current.vtk";

    const auto mat_by_elem = buildMaterialByElement(root);
    const auto coord_by_id = buildCoordinateById(root);
    std::vector<const nlohmann::json*> ordered_elements;
    bool has_pec_line = false;
    if (root.contains("mesh") && root["mesh"].contains("elements") &&
        root.contains("materials")) {
        for (const auto& e : root["mesh"]["elements"]) {
            const int eid = e.value("id", 0);
            if (!mat_by_elem.count(eid)) continue;
            const MaterialInfo info = materialInfo(root["materials"], mat_by_elem.at(eid));
            if (!isMapMaterial(info)) continue;
            ordered_elements.push_back(&e);
            if (info.type == "pec" && e.contains("intervals")) {
                for (const auto& interval : e["intervals"]) {
                    const int x0 = interval[0][0].get<int>();
                    const int y0 = interval[0][1].get<int>();
                    const int z0 = interval[0][2].get<int>();
                    const int x1 = interval[1][0].get<int>();
                    const int y1 = interval[1][1].get<int>();
                    const int z1 = interval[1][2].get<int>();
                    const int num_same = static_cast<int>(x0 == x1) +
                        static_cast<int>(y0 == y1) + static_cast<int>(z0 == z1);
                    if (num_same == 2) has_pec_line = true;
                }
            }
        }
        std::stable_sort(ordered_elements.begin(), ordered_elements.end(),
                         [&](const nlohmann::json* lhs, const nlohmann::json* rhs) {
                             const int lmat = mat_by_elem.at(lhs->value("id", 0));
                             const int rmat = mat_by_elem.at(rhs->value("id", 0));
                             const MaterialInfo li = materialInfo(root["materials"], lmat);
                             const MaterialInfo ri = materialInfo(root["materials"], rmat);
                             return materialPriority(li) < materialPriority(ri);
                         });
    }

    int face_count = 0;
    if (current_probe->contains("elementIds") && root.contains("mesh") &&
        root["mesh"].contains("elements")) {
        std::set<int> probe_elements;
        for (const auto& eid : (*current_probe)["elementIds"]) {
            probe_elements.insert(eid.get<int>());
        }
        for (const auto& e : root["mesh"]["elements"]) {
            if (!probe_elements.count(e.value("id", 0)) || !e.contains("intervals")) {
                continue;
            }
            for (const auto& interval : e["intervals"]) {
                const int x0 = interval[0][0].get<int>();
                const int y0 = interval[0][1].get<int>();
                const int z0 = interval[0][2].get<int>();
                const int x1 = interval[1][0].get<int>();
                const int y1 = interval[1][1].get<int>();
                const int z1 = interval[1][2].get<int>();
                face_count += std::abs(x1 - x0) * std::abs(y1 - y0) * std::abs(z1 - z0);
            }
        }
    }

    const int component_axis = axisFromComponent(current_probe->value("component", std::string("z")));
    for (const auto* eptr : ordered_elements) {
        const auto& e = *eptr;
        if (!e.contains("intervals")) continue;
        for (const auto& interval : e["intervals"]) {
            const int x0 = interval[0][0].get<int>();
            const int y0 = interval[0][1].get<int>();
            const int z0 = interval[0][2].get<int>();
            const int x1 = interval[1][0].get<int>();
            const int y1 = interval[1][1].get<int>();
            const int z1 = interval[1][2].get<int>();
            const int xi = std::min(x0, x1);
            const int yi = std::min(y0, y1);
            const int zi = std::min(z0, z1);
            const int xe = std::max(x0, x1);
            const int ye = std::max(y0, y1);
            const int ze = std::max(z0, z1);
            if (surfaceNormalAxis(xi, xe, yi, ye, zi, ze) == component_axis) {
                face_count -= surfaceAreaCells(xi, xe, yi, ye, zi, ze);
            }
        }
    }
    if (face_count < 0) face_count = 0;

    VtkMesh mesh;
    VtkMesh lineMesh;
    std::set<std::array<int, 6>> currentLineKeys;
    addCurrentFaces(mesh, face_count);

    double tag = 64.0;
    bool wire_collision_marker_available = has_pec_line;
    bool wire_intermediate_marker_available = true;
    for (const auto* eptr : ordered_elements) {
        const auto& e = *eptr;
        const int eid = e.value("id", 0);
        const MaterialInfo info = materialInfo(root["materials"], mat_by_elem.at(eid));
        if (info.type == "wire" && e.value("type", std::string()) == "polyline") {
            const auto segments = wireUnitSegments(e, coord_by_id);
            const bool mark_collision = wire_collision_marker_available && !segments.empty();
            if (mark_collision) wire_collision_marker_available = false;
            const bool mark_intermediate = !mark_collision &&
                wire_intermediate_marker_available && segments.size() > 2;
            if (mark_intermediate) wire_intermediate_marker_available = false;
            addClassicWirePolyline(lineMesh, e, coord_by_id, sx, sy, sz, tag,
                                   mark_collision, mark_intermediate);
        } else if (e.contains("intervals")) {
            for (const auto& interval : e["intervals"]) {
                const int x0 = interval[0][0].get<int>();
                const int y0 = interval[0][1].get<int>();
                const int z0 = interval[0][2].get<int>();
                const int x1 = interval[1][0].get<int>();
                const int y1 = interval[1][1].get<int>();
                const int z1 = interval[1][2].get<int>();
                const int xi = std::min(x0, x1);
                const int yi = std::min(y0, y1);
                const int zi = std::min(z0, z1);
                const int xe = std::max(x0, x1);
                const int ye = std::max(y0, y1);
                const int ze = std::max(z0, z1);
                const int num_same = static_cast<int>(xi == xe) +
                    static_cast<int>(yi == ye) + static_cast<int>(zi == ze);
                if (num_same == 2) {
                    addLineInterval(lineMesh, xi, yi, zi, xe, ye, ze, sx, sy, sz,
                                    lineMediaType(info), tag);
                } else if (num_same == 1) {
                    addSurfaceCurrentLines(lineMesh, currentLineKeys, xi, xe, yi, ye, zi, ze,
                                           sx, sy, sz, lineMediaType(info), tag);
                }
            }
        }
        tag += 64.0;
    }

    appendMesh(mesh, lineMesh);
    writeMeshFile(folder, vtk_name, mesh);
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
    const auto coord_by_id = buildCoordinateById(root);
    VtkMesh mesh;
    VtkMesh lineMesh;
    std::set<std::array<int, 6>> surfaceLineKeys;
    std::vector<const nlohmann::json*> ordered_elements;
    bool has_pec_line = false;
    if (root.contains("mesh") && root["mesh"].contains("elements") &&
        root.contains("materials")) {
        for (const auto& e : root["mesh"]["elements"]) {
            const int eid = e.value("id", 0);
            if (!mat_by_elem.count(eid)) continue;
            const MaterialInfo info = materialInfo(root["materials"], mat_by_elem.at(eid));
            if (!isMapMaterial(info)) continue;
            ordered_elements.push_back(&e);
            if (info.type == "pec" && e.contains("intervals")) {
                for (const auto& interval : e["intervals"]) {
                    const int x0 = interval[0][0].get<int>();
                    const int y0 = interval[0][1].get<int>();
                    const int z0 = interval[0][2].get<int>();
                    const int x1 = interval[1][0].get<int>();
                    const int y1 = interval[1][1].get<int>();
                    const int z1 = interval[1][2].get<int>();
                    const int num_same = static_cast<int>(x0 == x1) +
                        static_cast<int>(y0 == y1) + static_cast<int>(z0 == z1);
                    if (num_same == 2) has_pec_line = true;
                }
            }
        }
        std::stable_sort(ordered_elements.begin(), ordered_elements.end(),
                         [&](const nlohmann::json* lhs, const nlohmann::json* rhs) {
                             const int lmat = mat_by_elem.at(lhs->value("id", 0));
                             const int rmat = mat_by_elem.at(rhs->value("id", 0));
                             const MaterialInfo li = materialInfo(root["materials"], lmat);
                             const MaterialInfo ri = materialInfo(root["materials"], rmat);
                             return materialPriority(li) < materialPriority(ri);
                         });
    }
    double tag = 64.0;
    bool wire_collision_marker_available = has_pec_line;
    bool wire_intermediate_marker_available = true;
    if (!ordered_elements.empty()) {
        for (const auto* eptr : ordered_elements) {
            const auto& e = *eptr;
            const int eid = e.value("id", 0);
            const MaterialInfo info = materialInfo(root["materials"], mat_by_elem.at(eid));
            if (info.type == "wire" && e.value("type", std::string()) == "polyline") {
                const auto segments = wireUnitSegments(e, coord_by_id);
                const bool mark_collision = wire_collision_marker_available && !segments.empty();
                if (mark_collision) wire_collision_marker_available = false;
                const bool mark_intermediate = !mark_collision &&
                    wire_intermediate_marker_available && segments.size() > 2;
                if (mark_intermediate) wire_intermediate_marker_available = false;
                addClassicWirePolyline(lineMesh, e, coord_by_id, sx, sy, sz, tag,
                                       mark_collision, mark_intermediate);
            } else if (e.contains("intervals")) {
                for (const auto& interval : e["intervals"]) {
                    const int x0 = interval[0][0].get<int>();
                    const int y0 = interval[0][1].get<int>();
                    const int z0 = interval[0][2].get<int>();
                    const int x1 = interval[1][0].get<int>();
                    const int y1 = interval[1][1].get<int>();
                    const int z1 = interval[1][2].get<int>();
                    const int xi = std::min(x0, x1);
                    const int yi = std::min(y0, y1);
                    const int zi = std::min(z0, z1);
                    const int xe = std::max(x0, x1);
                    const int ye = std::max(y0, y1);
                    const int ze = std::max(z0, z1);
                    const int num_same = static_cast<int>(xi == xe) +
                        static_cast<int>(yi == ye) + static_cast<int>(zi == ze);
                    if (num_same == 2) {
                        addLineInterval(lineMesh, xi, yi, zi, xe, ye, ze, sx, sy, sz,
                                        lineMediaType(info), tag);
                    } else if (num_same == 1) {
                        addSurfaceQuads(mesh, xi, xe, yi, ye, zi, ze, sx, sy, sz,
                                        faceMediaType(info), tag);
                        addSurfaceBoundaryLines(lineMesh, surfaceLineKeys, xi, xe, yi, ye, zi, ze,
                                                sx, sy, sz, lineMediaType(info), tag);
                    } else {
                        addVolumeFaces(mesh, xi, xe, yi, ye, zi, ze, sx, sy, sz,
                                       faceMediaType(info), tag);
                        if (info.type == "pec" && ordered_elements.size() > 1 &&
                            xe - xi == 1 && ye - yi == 1 && ze - zi == 1) {
                            addGridLine(lineMesh, surfaceLineKeys, sx, sy, sz,
                                        xi, ye, zi, xe, ye, zi, lineMediaType(info), tag);
                        }
                    }
                }
            }
            tag += 64.0;
        }
    }
    appendMesh(mesh, lineMesh);
    writeMeshFile(folder, vtk_name, mesh);
}

} // namespace mapvtk
