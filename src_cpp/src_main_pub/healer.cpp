```cpp
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>

// Forward declarations and includes for external types used in the module
// These would typically be defined in Report_m and FDETYPES_m
// struct Shared_t;
// struct MediaData_t;
// struct XYZlimit_t;
// struct taglist_t;
// struct limit_t;
// struct SGGFDTDINFO_t;
// struct Border_t;
// struct MedioExtra_t;
// enum field_index { iEx, iEy, iEz, iHx, iHy, iHz };
// enum coord_index { icoord, jcoord, kcoord };
// enum side_index { fine, comi };
// constexpr int BUFSIZE = 256; // Example value, adjust as needed
// constexpr int IKINDMTAG = 4;
// constexpr int INTEGERSIZEOFMEDIAMATRICES = 4;
// using RKIND = double;
// constexpr RKIND prior_BV = -1.0;

namespace CreateMatrices_m {

    struct crosscheck_t {
        int32_t actual;
        int32_t NewActual;
        int32_t NewActual2;
        int32_t tent[4];
    };

    // Constants in
    const int32_t in[6][3][2] = {
        {
            {0, 1}, {1, 1}, {1, 0}
        },
        {
            {0, 0}, {1, 0}, {1, 0}
        },
        {
            {1, 1}, {0, 0}, {0, 0}
        },
        {
            {1, 1}, {0, 0}, {0, 0}
        },
        {
            {-1, -1}, {-1, -1}, {-1, -1}
        },
        {
            {-1, -1}, {-1, -1}, {-1, -1}
        }
    };

    // Constants on
    const int32_t on[6][3][2] = {
        {
            {0, 0}, {0, 0}, {0, 0}
        },
        {
            {0, 0}, {0, 0}, {0, 0}
        },
        {
            {0, 0}, {0, 0}, {0, 0}
        },
        {
            {0, 0}, {0, 0}, {0, 0}
        },
        {
            {-1, 0}, {0, 0}, {-1, -1}
        },
        {
            {0, -1}, {0, 0}, {0, -1}
        }
    };

    void SortInitEndWithIncreasingOrder(XYZlimit_t& p) {
        int32_t aux;
        if (p.XI > p.XE) {
            aux = p.XI;
            p.XI = p.XE;
            p.XE = aux;
        }
        if (p.YI > p.YE) {
            aux = p.YI;
            p.YI = p.YE;
            p.YE = aux;
        }
        if (p.ZI > p.ZE) {
            aux = p.ZI;
            p.ZI = p.ZE;
            p.ZE = aux;
        }
    }

    // Helper functions for CreateConformalPECVolume
    namespace internal {
        
        // These functions need access to MMi arrays, med, tags, Mtag, numertag, indicemedio
        // We will pass them as references or encapsulate in a context struct if needed.
        // For simplicity, we'll assume they are passed or accessible via a context.
        // However, to keep the signature close to the original, we will define them as static helpers
        // that take the necessary arrays.

        void fillBoundaryFaceIfAllEdgesPEC(int32_t i, int32_t j, int32_t k, int32_t face,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiEx,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiEy,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiEz,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiHx,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiHy,
                                           const std::vector<std::vector<std::vector<int32_t>>>& MMiHz,
                                           const std::vector<MediaData_t>& med,
                                           int32_t indicemedio,
                                           taglist_t& tags,
                                           std::vector<std::vector<std::vector<int32_t>>>& Mtag,
                                           int32_t numertag) {
            int32_t m1, m2, m3, m4;
            bool on_boundary = false;
            switch (face) {
                case 1: // FACE_X
                    m1 = MMiEy[i][j][k];
                    m2 = MMiEz[i][j][k];
                    m3 = MMiEy[i][j][k + 1];
                    m4 = MMiEz[i][j + 1][k];
                    break;
                case 2: // FACE_Y
                    m1 = MMiEx[i][j][k];
                    m2 = MMiEz[i][j][k];
                    m3 = MMiEx[i][j][k + 1];
                    m4 = MMiEz[i + 1][j][k];
                    break;
                case 3: // FACE_Z
                    m1 = MMiEy[i][j][k];
                    m2 = MMiEx[i][j][k];
                    m3 = MMiEy[i + 1][j][k];
                    m4 = MMiEx[i][j + 1][k];
                    break;
                default:
                    return;
            }
            on_boundary = (med[m1].Is.PEC) && (med[m2].Is.PEC) && (med[m3].Is.PEC) && (med[m4].Is.PEC);
            if (on_boundary) {
                Mtag[i][j][k] = 64 * numertag;
                switch (face) {
                    case 1: // FACE_X
                        MMiHx[i][j][k] = indicemedio;
                        tags.face.x[i][j][k] = 64 * numertag;
                        break;
                    case 2: // FACE_Y
                        MMiHy[i][j][k] = indicemedio;
                        tags.face.y[i][j][k] = 64 * numertag;
                        break;
                    case 3: // FACE_Z
                        MMiHz[i][j][k] = indicemedio;
                        tags.face.z[i][j][k] = 64 * numertag;
                        break;
                }
            }
        }

        bool hasCrossedPEC(int32_t m1, int32_t m2, int32_t m3, int32_t m4, const std::vector<MediaData_t>& med) {
            return (med[m1].Is.ConformalPEC || med[m1].Is.PEC) &&
                   (med[m2].Is.ConformalPEC || med[m2].Is.PEC) &&
                   (med[m3].Is.ConformalPEC || med[m3].Is.PEC) &&
                   (med[m4].Is.ConformalPEC || med[m4].Is.PEC);
        }

        bool hasCrossedPECOrConformalPEC(int32_t m1, int32_t m2, int32_t m3, int32_t m4, const std::vector<MediaData_t>& med) {
            return (med[m1].Is.ConformalPEC || med[m1].Is.PEC) ||
                   (med[m2].Is.ConformalPEC || med[m2].Is.PEC) ||
                   (med[m3].Is.ConformalPEC || med[m3].Is.PEC) ||
                   (med[m4].Is.ConformalPEC || med[m4].Is.PEC);
        }

        void fillPECFaceInsideVolume(int32_t i, int32_t j, int32_t k, int32_t face,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiEx,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiEy,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiEz,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiHx,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiHy,
                                     const std::vector<std::vector<std::vector<int32_t>>>& MMiHz,
                                     const std::vector<MediaData_t>& med,
                                     int32_t indicemedio,
                                     taglist_t& tags,
                                     std::vector<std::vector<std::vector<int32_t>>>& Mtag,
                                     int32_t numertag) {
            int32_t m1, m2, m3, m4, m;
            bool on_boundary;
            switch (face) {
                case 1: // FACE_X
                    m1 = MMiEy[i][j][k];
                    m2 = MMiEz[i][j][k];
                    m3 = MMiEy[i][j][k + 1];
                    m4 = MMiEz[i][j + 1][k];
                    m = MMiHx[i][j][k];
                    break;
                case 2: // FACE_Y
                    m1 = MMiEx[i][j][k];
                    m2 = MMiEz[i][j][k];
                    m3 = MMiEx[i][j][k + 1];
                    m4 = MMiEz[i + 1][j][k];
                    m = MMiHy[i][j][k];
                    break;
                case 3: // FACE_Z
                    m1 = MMiEy[i][j][k];
                    m2 = MMiEx[i][j][k];
                    m3 = MMiEy[i + 1][j][k];
                    m4 = MMiEx[i][j + 1][k];
                    m = MMiHz[i][j][k];
                    break;
                default:
                    return;
            }

            on_boundary = (med[m1].Is.PEC || med[m1].Is.conformalPEC) &&
                          (med[m2].Is.PEC || med[m2].Is.conformalPEC) &&
                          (med[m3].Is.PEC || med[m3].Is.conformalPEC) &&
                          (med[m4].Is.PEC || med[m4].Is.conformalPEC);

            if (on_boundary && !(med[m].Is.PEC || med[m].Is.ConformalPEC)) {
                Mtag[i][j][k] = 64 * numertag;
                switch (face) {
                    case 1: // FACE_X
                        MMiHx[i][j][k] = indicemedio;
                        tags.face.x[i][j][k] = 64 * numertag;
                        break;
                    case 2: // FACE_Y
                        MMiHy[i][j][k] = indicemedio;
                        tags.face.y[i][j][k] = 64 * numertag;
                        break;
                    case 3: // FACE_Z
                        MMiHz[i][j][k] = indicemedio;
                        tags.face.z[i][j][k] = 64 * numertag;
                        break;
                }
            } else if (on_boundary && (med[m].Is.PEC || med[m].Is.ConformalPEC)) {
                Mtag[i][j][k] = 64 * numertag;
                switch (face) {
                    case 1: // FACE_X
                        tags.face.x[i][j][k] = 64 * numertag;
                        break;
                    case 2: // FACE_Y
                        tags.face.y[i][j][k] = 64 * numertag;
                        break;
                    case 3: // FACE_Z
                        tags.face.z[i][j][k] = 64 * numertag;
                        break;
                }
            }
        }

        int32_t countCrossesX(int32_t j, int32_t k, int32_t xi, int32_t xe,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiEx,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiHx,
                              const std::vector<MediaData_t>& med) {
            int32_t res = 0;
            bool crossed = false;
            for (int32_t i = xi; i <= xe + 1; ++i) {
                int32_t mE = MMiEx[i][j][k];
                int32_t mEPrev = MMiEx[i - 1][j][k];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHx[i - 1][j][k], MMiHx[i - 1][j - 1][k], MMiHx[i - 1][j][k - 1], MMiHx[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHx[i + 1][j][k], MMiHx[i + 1][j - 1][k], MMiHx[i + 1][j][k - 1], MMiHx[i + 1][j - 1][k - 1], med);
                    }
                }
                if (crossed) res++;
            }
            if (res != 0) {
                if (res % 2 != 0) {
                    throw std::runtime_error("uneven number of crosses");
                }
            }
            return res;
        }

        void fillEdgesInsideVolumeX(int32_t j, int32_t k, int32_t xi, int32_t xe,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiEx,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiHx,
                                    const std::vector<MediaData_t>& med,
                                    int32_t indicemedio,
                                    taglist_t& tags,
                                    std::vector<std::vector<std::vector<int32_t>>>& Mtag,
                                    int32_t numertag) {
            bool inside_volume = false;
            bool crossed = false;
            int32_t n_crosses = countCrossesX(j, k, xi, xe, MMiEx, MMiHx, med);
            std::vector<int32_t> idx_in(n_crosses / 2, 0);
            std::vector<int32_t> idx_out(n_crosses / 2, 0);

            for (int32_t i = xi; i <= xe + 1; ++i) {
                int32_t mE = MMiEx[i][j][k];
                int32_t mEPrev = MMiEx[i - 1][j][k];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHx[i - 1][j][k], MMiHx[i - 1][j - 1][k], MMiHx[i - 1][j][k - 1], MMiHx[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHx[i][j][k], MMiHx[i][j - 1][k], MMiHx[i][j][k - 1], MMiHx[i][j - 1][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHx[i + 1][j][k], MMiHx[i + 1][j - 1][k], MMiHx[i + 1][j][k - 1], MMiHx[i + 1][j - 1][k - 1], med);
                    }
                }
                if (crossed) inside_volume = !inside_volume;
                if (crossed && inside_volume) {
                    // Find first zero in idx_in and assign i
                    for (auto& val : idx_in) {
                        if (val == 0) {
                            val = i;
                            break;
                        }
                    }
                }
                if (crossed && !inside_volume) {
                    // Find first zero in idx_out and assign i-1
                    for (auto& val : idx_out) {
                        if (val == 0) {
                            val = i - 1;
                            break;
                        }
                    }
                }
            }

            for (size_t ii = 0; ii < idx_in.size(); ++ii) {
                if (idx_in[ii] != 0 && idx_out[ii] != 0) {
                    for (int32_t i = idx_in[ii]; i < idx_out[ii]; ++i) {
                        if (MMiEx[i][j][k] == 1) {
                            MMiEx[i][j][k] = indicemedio;
                            Mtag[i][j][k] = 64 * numertag;
                            tags.edge.x[i][j][k] = 64 * numertag;
                        }
                    }
                }
            }
        }

        int32_t countCrossesY(int32_t i, int32_t k, int32_t yi, int32_t ye,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiEy,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiHy,
                              const std::vector<MediaData_t>& med) {
            int32_t res = 0;
            bool crossed = false;
            for (int32_t j = yi; j <= ye + 1; ++j) {
                int32_t mE = MMiEy[i][j][k];
                int32_t mEPrev = MMiEy[i][j - 1][k];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHy[i][j][k], MMiHy[i - 1][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k - 1], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHy[i][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k], MMiHy[i - 1][j][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHy[i][j - 1][k], MMiHy[i][j - 1][k - 1], MMiHy[i - 1][j - 1][k], MMiHy[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHy[i][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k], MMiHy[i - 1][j][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHy[i][j + 1][k], MMiHy[i][j + 1][k - 1], MMiHy[i - 1][j + 1][k], MMiHy[i - 1][j + 1][k - 1], med);
                    }
                }
                if (crossed) res++;
            }
            if (res != 0) {
                if (res % 2 != 0) {
                    throw std::runtime_error("uneven number of crosses");
                }
            }
            return res;
        }

        void fillEdgesInsideVolumeY(int32_t i, int32_t k, int32_t yi, int32_t ye,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiEy,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiHy,
                                    const std::vector<MediaData_t>& med,
                                    int32_t indicemedio,
                                    taglist_t& tags,
                                    std::vector<std::vector<std::vector<int32_t>>>& Mtag,
                                    int32_t numertag) {
            bool inside_volume = false;
            bool crossed = false;
            int32_t n_crosses = countCrossesY(i, k, yi, ye, MMiEy, MMiHy, med);
            std::vector<int32_t> idx_in(n_crosses / 2, 0);
            std::vector<int32_t> idx_out(n_crosses / 2, 0);

            for (int32_t j = yi; j <= ye + 1; ++j) {
                int32_t mE = MMiEy[i][j][k];
                int32_t mEPrev = MMiEy[i][j - 1][k];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHy[i][j][k], MMiHy[i - 1][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k - 1], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHy[i][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k], MMiHy[i - 1][j][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHy[i][j - 1][k], MMiHy[i][j - 1][k - 1], MMiHy[i - 1][j - 1][k], MMiHy[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHy[i][j][k], MMiHy[i][j][k - 1], MMiHy[i - 1][j][k], MMiHy[i - 1][j][k - 1], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHy[i][j + 1][k], MMiHy[i][j + 1][k - 1], MMiHy[i - 1][j + 1][k], MMiHy[i - 1][j + 1][k - 1], med);
                    }
                }
                if (crossed) inside_volume = !inside_volume;
                if (crossed && inside_volume) {
                    for (auto& val : idx_in) {
                        if (val == 0) {
                            val = j;
                            break;
                        }
                    }
                }
                if (crossed && !inside_volume) {
                    for (auto& val : idx_out) {
                        if (val == 0) {
                            val = j - 1;
                            break;
                        }
                    }
                }
            }

            for (size_t jj = 0; jj < idx_in.size(); ++jj) {
                if (idx_in[jj] != 0 && idx_out[jj] != 0) {
                    for (int32_t j = idx_in[jj]; j < idx_out[jj]; ++j) {
                        if (MMiEy[i][j][k] == 1) {
                            MMiEy[i][j][k] = indicemedio;
                            Mtag[i][j][k] = 64 * numertag;
                            tags.edge.y[i][j][k] = 64 * numertag;
                        }
                    }
                }
            }
        }

        int32_t countCrossesZ(int32_t i, int32_t j, int32_t zi, int32_t ze,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiEz,
                              const std::vector<std::vector<std::vector<int32_t>>>& MMiHz,
                              const std::vector<MediaData_t>& med) {
            int32_t res = 0;
            bool crossed = false;
            for (int32_t k = zi; k <= ze + 1; ++k) {
                int32_t mE = MMiEz[i][j][k];
                int32_t mEPrev = MMiEz[i][j][k - 1];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHz[i][j][k - 1], MMiHz[i - 1][j][k - 1], MMiHz[i][j - 1][k - 1], MMiHz[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHz[i][j][k + 1], MMiHz[i - 1][j][k + 1], MMiHz[i][j - 1][k + 1], MMiHz[i - 1][j - 1][k + 1], med);
                    }
                }
                if (crossed) res++;
            }
            if (res != 0) {
                if (res % 2 != 0) {
                    throw std::runtime_error("uneven number of crosses");
                }
            }
            return res;
        }

        void fillEdgesInsideVolumeZ(int32_t i, int32_t j, int32_t zi, int32_t ze,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiEz,
                                    const std::vector<std::vector<std::vector<int32_t>>>& MMiHz,
                                    const std::vector<MediaData_t>& med,
                                    int32_t indicemedio,
                                    taglist_t& tags,
                                    std::vector<std::vector<std::vector<int32_t>>>& Mtag,
                                    int32_t numertag) {
            bool inside_volume = false;
            bool crossed = false;
            int32_t n_crosses = countCrossesZ(i, j, zi, ze, MMiEz, MMiHz, med);
            std::vector<int32_t> idx_in(n_crosses / 2, 0);
            std::vector<int32_t> idx_out(n_crosses / 2, 0);

            for (int32_t k = zi; k <= ze + 1; ++k) {
                int32_t mE = MMiEz[i][j][k];
                int32_t mEPrev = MMiEz[i][j][k - 1];
                crossed = false;
                if (!(med[mE].Is.ConformalPEC || med[mE].Is.PEC)) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                    } else if (med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHz[i][j][k - 1], MMiHz[i - 1][j][k - 1], MMiHz[i][j - 1][k - 1], MMiHz[i - 1][j - 1][k - 1], med);
                    }
                } else if (med[mE].Is.ConformalPEC || med[mE].Is.PEC) {
                    if (!(med[mEPrev].Is.ConformalPEC || med[mEPrev].Is.PEC)) {
                        crossed = hasCrossedPECOrConformalPEC(MMiHz[i][j][k], MMiHz[i - 1][j][k], MMiHz[i][j - 1][k], MMiHz[i - 1][j - 1][k], med);
                        crossed = crossed || hasCrossedPECOrConformalPEC(MMiHz[i][j][k + 1], MMiHz[i - 1][j][k + 1], MMiHz[i][j - 1][k + 1], MMiHz[i - 1][j - 1][k + 1], med);
                    }
                }
                if (crossed) inside_volume = !inside_volume;
                if (crossed && inside_volume) {
                    for (auto& val : idx_in) {
                        if (val == 0) {
                            val = k;
                            break;
                        }
                    }
                }
                if (crossed && !inside_volume) {
                    for (auto& val : idx_out) {
                        if (val == 0) {
                            val = k - 1;
                            break;
                        }
                    }
                }
            }

            for (size_t kk = 0; kk < idx_in.size(); ++kk) {
                if (idx_in[kk] != 0 && idx_out[kk] != 0) {
                    for (int32