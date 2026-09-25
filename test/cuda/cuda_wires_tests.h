#ifndef CUDA_WIRES_TESTS_H
#define CUDA_WIRES_TESTS_H

#ifndef CompileWithCUDA
#error "test/cuda is only built when SEMBA_FDTD_ENABLE_CUDA=ON (CompileWithCUDA)"
#endif

#include "cuda_test_utils.h"

using namespace cuda_test;

/*
 * Host AdvanceWiresE (holland, thickness 1), in fdtd_real:
 *   Q = CteProp*Qpast - CtePlain*(Iplus-Iminus)
 *   terminal: the missing side mirrors with a minus (periodic uses a plus)
 *   E -= cte5*I   (skipped when coupled=0)
 *   I = cte1*I - cte3*(fPlus*Qplus - fMinus*Qminus) + cte2*E
 *   IsPMC sets I=0 after the E injection. Soft V adds vscale*evolucion(t).
 */

static fdtd_wire_node make_node()
{
    fdtd_wire_node n{};
    n.exists = 1;
    n.node_inside = -1;
    for (int k = 0; k < FDTD_WIRE_NEIGH; ++k) {
        n.plus_seg[k] = -1;
        n.minus_seg[k] = -1;
    }
    n.deltaevol = (fdtd_real)1;
    return n;
}

static fdtd_wire_seg make_seg()
{
    fdtd_wire_seg s{};
    s.coupled = 1;
    s.field_comp = 2;
    s.charge_plus = -1;
    s.charge_minus = -1;
    s.deltaevol = (fdtd_real)1;
    s.fraction_plus = (fdtd_real)1;
    s.fraction_minus = (fdtd_real)1;
    return s;
}

static int upload_wire_grid(fdtd_cuda_ctx *ctx, fdtd_dims3 d)
{
    if (fdtd_cuda_alloc_fields(ctx, d, d, d, d, d, d) != 1) return 0;
    SixFields f = make_six_zero(d);
    return upload_six(ctx, f);
}

TEST(cuda, wire_straight_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    const int i = 2, j = 3, k = 4;
    SixFields f = make_six_zero(d);
    f.Ez[idx3(i, j, k, d.nx, d.ny)] = (fdtd_real)5;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_wire_node nodes[2] = {make_node(), make_node()};
    nodes[0].n_minus = 1;
    nodes[0].minus_seg[0] = 0;
    nodes[0].cte_prop = (fdtd_real)0.8;
    nodes[0].cte_plain = (fdtd_real)0.4;
    nodes[1].n_plus = 1;
    nodes[1].plus_seg[0] = 0;
    nodes[1].cte_prop = (fdtd_real)0.8;
    nodes[1].cte_plain = (fdtd_real)0.4;

    fdtd_wire_seg seg = make_seg();
    seg.i = i; seg.j = j; seg.k = k;
    seg.charge_minus = 0;
    seg.charge_plus = 1;
    seg.cte1 = (fdtd_real)0.5;
    seg.cte2 = (fdtd_real)0.25;
    seg.cte3 = (fdtd_real)0.1;
    seg.cte5 = (fdtd_real)0.2;

    fdtd_real current[1] = {(fdtd_real)3};
    fdtd_real charge[2] = {(fdtd_real)1, (fdtd_real)2};
    fdtd_real qpast[2] = {(fdtd_real)0, (fdtd_real)0};
    fdtd_real samples[1] = {(fdtd_real)0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, &seg, 1, nodes, 2, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);

    const fdtd_real I = current[0];
    const fdtd_real q0 = nodes[0].cte_prop * charge[0] - nodes[0].cte_plain * ((-I) - I);
    const fdtd_real q1 = nodes[1].cte_prop * charge[1] - nodes[1].cte_plain * (I - (-I));
    const fdtd_real Ez = (fdtd_real)5 - seg.cte5 * I;
    const fdtd_real expect_I = seg.cte1 * I - seg.cte3 * (q1 - q0) + seg.cte2 * Ez;

    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, (fdtd_real)1, (fdtd_real)0.5), 1);
    fdtd_real Ic[1], Ip[1], Qc[2], Qp[2];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    EXPECT_NEAR(Qp[0], charge[0], 1e-5);
    EXPECT_NEAR(Qp[1], charge[1], 1e-5);
    EXPECT_NEAR(Qc[0], q0, 1e-5);
    EXPECT_NEAR(Qc[1], q1, 1e-5);
    EXPECT_NEAR(Ip[0], I, 1e-5);
    EXPECT_NEAR(Ic[0], expect_I, 1e-4);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(i, j, k, d.nx, d.ny)], Ez, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_terminal_mirror_and_junction)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);

    fdtd_wire_node nodes[2] = {make_node(), make_node()};
    nodes[0].n_minus = 1;
    nodes[0].minus_seg[0] = 0;
    nodes[0].is_periodic = 1;
    nodes[0].cte_prop = (fdtd_real)1;
    nodes[0].cte_plain = (fdtd_real)1;
    nodes[1].n_plus = 2;
    nodes[1].plus_seg[0] = 0;
    nodes[1].plus_seg[1] = 1;
    nodes[1].cte_prop = (fdtd_real)1;
    nodes[1].cte_plain = (fdtd_real)1;

    fdtd_wire_seg segs[2] = {make_seg(), make_seg()};
    segs[0].coupled = 0;
    segs[0].charge_minus = 0;
    segs[0].charge_plus = 1;
    segs[1].coupled = 0;
    segs[1].charge_minus = 1;
    segs[1].charge_plus = 1;
    segs[0].cte1 = segs[1].cte1 = (fdtd_real)0;
    segs[0].cte3 = segs[1].cte3 = (fdtd_real)0;

    fdtd_real current[2] = {(fdtd_real)4, (fdtd_real)1};
    fdtd_real charge[2] = {(fdtd_real)0, (fdtd_real)0};
    fdtd_real qpast[2] = {0, 0};
    fdtd_real samples[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, segs, 2, nodes, 2, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, 0, 0), 1);
    fdtd_real Ic[2], Ip[2], Qc[2], Qp[2];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    /* Periodic terminal: Iplus = +Iminus = 4, Q = -(4-4) = 0. */
    EXPECT_NEAR(Qc[0], (fdtd_real)0, 1e-5);
    /* Junction sums both plus currents and does not mirror (n_plus != 1): Q = -(5-0). */
    EXPECT_NEAR(Qc[1], (fdtd_real)-5, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_pmc_zeros_current_after_inject)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    SixFields f = make_six_zero(d);
    f.Ez[idx3(2, 2, 2, d.nx, d.ny)] = (fdtd_real)1;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_wire_node node = make_node();
    node.n_plus = 1;
    node.plus_seg[0] = 0;
    node.cte_prop = (fdtd_real)1;
    fdtd_wire_seg seg = make_seg();
    seg.i = seg.j = seg.k = 2;
    seg.is_pmc = 1;
    seg.charge_plus = 0;
    seg.charge_minus = 0;
    seg.cte5 = (fdtd_real)0.5;
    seg.cte1 = (fdtd_real)9;
    fdtd_real current[1] = {(fdtd_real)4};
    fdtd_real charge[1] = {0};
    fdtd_real qpast[1] = {0};
    fdtd_real samples[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, &seg, 1, &node, 1, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, 0, 0), 1);
    fdtd_real Ic[1], Ip[1], Qc[1], Qp[1];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    EXPECT_NEAR(Ic[0], (fdtd_real)0, 1e-5);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(2, 2, 2, d.nx, d.ny)], (fdtd_real)1 - (fdtd_real)0.5 * 4, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_shielded_skips_field)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    SixFields f = make_six_zero(d);
    f.Ez[idx3(2, 2, 2, d.nx, d.ny)] = (fdtd_real)7;
    ASSERT_EQ(upload_six(ctx, f), 1);

    fdtd_wire_node node = make_node();
    node.cte_prop = (fdtd_real)1;
    node.n_plus = 1;
    node.plus_seg[0] = 0;
    fdtd_wire_seg seg = make_seg();
    seg.coupled = 0;
    seg.i = seg.j = seg.k = 2;
    seg.charge_plus = 0;
    seg.charge_minus = 0;
    seg.cte1 = (fdtd_real)0.5;
    seg.cte2 = (fdtd_real)3;
    seg.cte5 = (fdtd_real)4;
    fdtd_real current[1] = {(fdtd_real)2};
    fdtd_real charge[1] = {0}, qpast[1] = {0}, samples[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, &seg, 1, &node, 1, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, 0, 0), 1);
    fdtd_real Ic[1], Ip[1], Qc[1], Qp[1];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    EXPECT_NEAR(Ic[0], (fdtd_real)1, 1e-5);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(2, 2, 2, d.nx, d.ny)], (fdtd_real)7, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_mur_node_matches_host)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    fdtd_wire_node nodes[2] = {make_node(), make_node()};
    nodes[0].cte_prop = (fdtd_real)1;
    nodes[0].cte_plain = (fdtd_real)0;
    nodes[1].is_mur = 1;
    nodes[1].node_inside = 0;
    nodes[1].cte_mur = (fdtd_real)0.5;
    nodes[1].cte_prop = (fdtd_real)9;
    fdtd_wire_seg seg = make_seg();
    seg.coupled = 0;
    seg.charge_plus = 0;
    seg.charge_minus = 1;
    seg.cte1 = (fdtd_real)0;
    seg.cte3 = (fdtd_real)0;
    fdtd_real current[1] = {0};
    fdtd_real charge[2] = {(fdtd_real)4, (fdtd_real)1};
    fdtd_real qpast[2] = {(fdtd_real)2, (fdtd_real)8};
    fdtd_real samples[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, &seg, 1, nodes, 2, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, 0, 0), 1);
    fdtd_real Ic[1], Ip[1], Qc[2], Qp[2];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    /* Inside saves past=4 and stays at 4 (cte_plain 0). Mur: past_in + cte*(present_in - mur_past). */
    const fdtd_real expect = (fdtd_real)4 + (fdtd_real)0.5 * ((fdtd_real)4 - (fdtd_real)1);
    EXPECT_NEAR(Qp[1], (fdtd_real)1, 1e-5);
    EXPECT_NEAR(Qc[1], expect, 1e-5);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_two_segments_share_cell)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    SixFields f = make_six_zero(d);
    f.Ez[idx3(3, 3, 3, d.nx, d.ny)] = (fdtd_real)10;
    ASSERT_EQ(upload_six(ctx, f), 1);
    fdtd_wire_node node = make_node();
    node.cte_prop = (fdtd_real)1;
    fdtd_wire_seg segs[2] = {make_seg(), make_seg()};
    segs[0].i = segs[0].j = segs[0].k = 3;
    segs[1].i = segs[1].j = segs[1].k = 3;
    segs[0].cte5 = (fdtd_real)1;
    segs[1].cte5 = (fdtd_real)2;
    segs[0].charge_plus = segs[0].charge_minus = 0;
    segs[1].charge_plus = segs[1].charge_minus = 0;
    segs[0].cte1 = segs[1].cte1 = (fdtd_real)0;
    segs[0].cte2 = segs[1].cte2 = (fdtd_real)0;
    segs[0].cte3 = segs[1].cte3 = (fdtd_real)0;
    fdtd_real current[2] = {(fdtd_real)1, (fdtd_real)3};
    fdtd_real charge[1] = {0}, qpast[1] = {0}, samples[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, segs, 2, &node, 1, samples, 1, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, 0, 0), 1);
    SixFields out = make_six_zero(d);
    ASSERT_EQ(download_six(ctx, out), 1);
    EXPECT_NEAR(out.Ez[idx3(3, 3, 3, d.nx, d.ny)], (fdtd_real)10 - (fdtd_real)1 * 1 - (fdtd_real)2 * 3, 1e-4);
    fdtd_cuda_destroy(ctx);
}

TEST(cuda, wire_voltage_out_of_range_is_zero)
{
    fdtd_cuda_ctx *ctx = nullptr;
    CUDA_REQUIRE_CTX(ctx);
    fdtd_dims3 d = kTiny;
    ASSERT_EQ(upload_wire_grid(ctx, d), 1);
    fdtd_wire_node node = make_node();
    node.cte_prop = (fdtd_real)1;
    fdtd_wire_seg seg = make_seg();
    seg.coupled = 0;
    seg.charge_plus = 0;
    seg.charge_minus = 0;
    seg.has_v = 1;
    seg.numus = 1;
    seg.vscale = (fdtd_real)5;
    seg.cte1 = (fdtd_real)1;
    fdtd_real samples[2] = {(fdtd_real)1, (fdtd_real)2};
    fdtd_real current[1] = {(fdtd_real)3};
    fdtd_real charge[1] = {0}, qpast[1] = {0};
    ASSERT_EQ(fdtd_cuda_upload_wires(ctx, &seg, 1, &node, 1, samples, 2, current, charge, qpast,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0), 1);
    ASSERT_EQ(fdtd_cuda_advance_wires_e(ctx, (fdtd_real)20, (fdtd_real)0), 1);
    fdtd_real Ic[1], Ip[1], Qc[1], Qp[1];
    ASSERT_EQ(fdtd_cuda_download_wires(ctx, Ic, Ip, Qc, Qp), 1);
    EXPECT_NEAR(Ic[0], (fdtd_real)3, 1e-5);
    fdtd_cuda_destroy(ctx);
}

#endif
