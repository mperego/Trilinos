// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file  test_static_rank_views.cpp
    \brief PROTOTYPE.  Checks the static-rank Kokkos::View overloads against a hand-computed reference, and
           checks that the Kokkos::DynRankView overloads (which now convert their arguments and forward to
           the static-rank implementation) give the same answer, for double and Sacado Fad scalars,
           LayoutLeft/LayoutRight/LayoutStride containers, and subviews.
*/
#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Sacado.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include <iostream>
#include <cmath>

using Dev = Kokkos::DefaultExecutionSpace::device_type;
using at  = Intrepid2::ArrayTools<Dev>;
using fst = Intrepid2::FunctionSpaceTools<Dev>;
using ct  = Intrepid2::CellTools<Dev>;

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; std::cout << "  FAILED: " #cond " (line " << __LINE__ << ")\n"; } } while (0)

// ---- compile-time checks: the static overloads only accept the documented ranks, and never a DynRankView
template<class O, class L, class R> concept CanFieldProduct = requires(O o, L l, R r) { at::matvecProductDataField(o, l, r, 'N'); };
template<class O, class J, class I> concept CanTransformGrad = requires(O o, J j, I i) { fst::HGRADtransformGRAD(o, j, i); };
using V2 = Kokkos::View<double**, Dev>;   using V3 = Kokkos::View<double***, Dev>;
using V4 = Kokkos::View<double****, Dev>; using V5 = Kokkos::View<double*****, Dev>;
using DRV = Kokkos::DynRankView<double, Dev>;
static_assert( CanFieldProduct<V4, V4, V3>);   // (C,F,P,D) = (C,P,D,D) x (F,P,D)
static_assert( CanFieldProduct<V4, V2, V4>);   // (C,F,P,D) = (C,P)     x (C,F,P,D)
static_assert(!CanFieldProduct<V3, V4, V3>);   // output must have rank 4
static_assert(!CanFieldProduct<V4, V5, V3>);   // data rank 5 not allowed
static_assert( CanFieldProduct<DRV, DRV, DRV>);// DynRankView overload still exists
static_assert(!CanFieldProduct<V4, DRV, V3>);  // mixing static and dynamic rank arguments is not offered by this prototype

template<class Scalar> double val(const Scalar& s) { return Sacado::scalarValue(s); }
template<class Scalar> double dx(const Scalar& s, int k) { if constexpr (std::is_arithmetic_v<Scalar>) return 0.0; else return s.dx(k); }
template<class Scalar> Scalar mkScalar(double v, int k, int nderiv) { if constexpr (std::is_arithmetic_v<Scalar>) return v; else return Scalar(nderiv, k % nderiv, v); }
// allocates a view of the given logical dimensions, adding the hidden derivative dimension for Fad types
template<class V, class... Dims> V alloc(const std::string& label, int nderiv, Dims... dims) {
  if constexpr (std::is_arithmetic_v<typename V::value_type>) return V(label, dims...); else return V(label, dims..., nderiv + 1);
}

template<class A, class B>
double maxDiff(const A& a, const B& b, int C, int F, int P, int D, int nderiv) {
  double e = 0;
  for (int c = 0; c < C; ++c) for (int f = 0; f < F; ++f) for (int p = 0; p < P; ++p) for (int d = 0; d < D; ++d) {
    const typename A::value_type x = a(c, f, p, d); const typename B::value_type y = b(c, f, p, d);
    e = std::max(e, std::abs(val(x) - val(y)));
    for (int k = 0; k < nderiv; ++k) e = std::max(e, std::abs(dx(x, k) - dx(y, k)));
  }
  return e;
}

template<class Scalar, class Layout>
void runCase(const char* name, int nderiv) {
  const int C = 5, F = 4, P = 3, D = 3;
  using Stat4 = Kokkos::View<Scalar****, Layout, Dev>; using Stat3 = Kokkos::View<Scalar***, Layout, Dev>;
  using Dyn   = Kokkos::DynRankView<Scalar, Layout, Dev>;
  auto make4 = [&](const char* l, int a, int b, int c, int d) { return alloc<Stat4>(l, nderiv, a, b, c, d); };
  auto make3 = [&](const char* l, int a, int b, int c)        { return alloc<Stat3>(l, nderiv, a, b, c); };
  auto makeD = [&](const char* l, int a, int b, int c, int d) { return d ? alloc<Dyn>(l, nderiv, a, b, c, d) : alloc<Dyn>(l, nderiv, a, b, c); };

  Stat4 sJacInv = make4("sJacInv", C, P, D, D), sOut = make4("sOut", C, F, P, D), sRef = make4("sRef", C, F, P, D);
  Stat3 sGrads  = make3("sGrads", F, P, D);
  Dyn   dJacInv = makeD("dJacInv", C, P, D, D), dOut = makeD("dOut", C, F, P, D), dGrads = makeD("dGrads", F, P, D, 0);

  auto mk = [&](double v, int k) { return mkScalar<Scalar>(v, k, nderiv); };
  for (int c = 0; c < C; ++c) for (int p = 0; p < P; ++p) for (int i = 0; i < D; ++i) for (int j = 0; j < D; ++j) {
    Scalar v = mk(std::sin(1.0 + c + 2 * p + 3 * i + 5 * j), c + p + i + j); sJacInv(c, p, i, j) = v; dJacInv(c, p, i, j) = v; }
  for (int f = 0; f < F; ++f) for (int p = 0; p < P; ++p) for (int j = 0; j < D; ++j) {
    Scalar v = mk(std::cos(0.5 + f + 3 * p + 7 * j), f + p + j); sGrads(f, p, j) = v; dGrads(f, p, j) = v; }
  // reference: out(c,f,p,i) = sum_j JacInv(c,p,j,i) * grads(f,p,j)   (transpose)
  for (int c = 0; c < C; ++c) for (int f = 0; f < F; ++f) for (int p = 0; p < P; ++p) for (int i = 0; i < D; ++i) {
    Scalar t = 0.0; for (int j = 0; j < D; ++j) t += sJacInv(c, p, j, i) * sGrads(f, p, j); sRef(c, f, p, i) = t; }

  fst::HGRADtransformGRAD(sOut, sJacInv, sGrads);     // static-rank Views
  fst::HGRADtransformGRAD(dOut, dJacInv, dGrads);     // DynRankViews: shim -> static implementation
  const double eS = maxDiff(sOut, sRef, C, F, P, D, nderiv), eD = maxDiff(dOut, sRef, C, F, P, D, nderiv);
  CHECK(eS < 1e-13); CHECK(eD < 1e-13);

  // LayoutStride arguments (what subviews produce, and what Basis::OutputViewType uses today)
  Kokkos::DynRankView<Scalar, Kokkos::LayoutStride, Dev> dOutS = dOut, dJacInvS = dJacInv, dGradsS = dGrads;
  Kokkos::deep_copy(dOut, Scalar(0.0));
  fst::HGRADtransformGRAD(dOutS, dJacInvS, dGradsS);
  const double eStr = maxDiff(dOut, sRef, C, F, P, D, nderiv);
  CHECK(eStr < 1e-13);

  // a strided subview in the cell dimension (cells 1..3) of the DynRankViews
  auto dOutSub = Kokkos::subview(dOut, std::make_pair(1, 4), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  auto dJacSub = Kokkos::subview(dJacInv, std::make_pair(1, 4), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  Kokkos::deep_copy(dOut, Scalar(0.0));
  fst::HGRADtransformGRAD(dOutSub, dJacSub, dGrads);
  double eSub = 0;
  for (int c = 1; c < 4; ++c) for (int f = 0; f < F; ++f) for (int p = 0; p < P; ++p) for (int i = 0; i < D; ++i) {
    const Scalar a = dOut(c, f, p, i), b = sRef(c, f, p, i);
    eSub = std::max(eSub, std::abs(val(a) - val(b)));
    for (int k = 0; k < nderiv; ++k) eSub = std::max(eSub, std::abs(dx(a, k) - dx(b, k)));
  }
  CHECK(eSub < 1e-13); CHECK(val(Scalar(dOut(0, 0, 0, 0))) == 0.0); CHECK(val(Scalar(dOut(4, 0, 0, 0))) == 0.0);

  // static subviews work too
  auto sOutSub = Kokkos::subview(sOut, std::make_pair(1, 4), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  auto sJacSub = Kokkos::subview(sJacInv, std::make_pair(1, 4), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
  Kokkos::deep_copy(sOut, Scalar(0.0));
  fst::HGRADtransformGRAD(sOutSub, sJacSub, sGrads);
  double eSSub = 0;
  for (int c = 1; c < 4; ++c) for (int f = 0; f < F; ++f) for (int p = 0; p < P; ++p) for (int i = 0; i < D; ++i)
    eSSub = std::max(eSSub, std::abs(val(Scalar(sOut(c, f, p, i))) - val(Scalar(sRef(c, f, p, i)))));
  CHECK(eSSub < 1e-13);

  std::cout << name << ": static " << eS << ", dyn shim " << eD << ", dyn LayoutStride " << eStr << ", dyn subview " << eSub << ", static subview " << eSSub << "\n";
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    runCase<double, Kokkos::LayoutRight>("double  LayoutRight", 0);
    runCase<double, Kokkos::LayoutLeft >("double  LayoutLeft ", 0);
    runCase<Sacado::Fad::DFad<double>, Kokkos::LayoutRight>("DFad    LayoutRight", 4);
    runCase<Sacado::Fad::DFad<double>, Kokkos::LayoutLeft >("DFad    LayoutLeft ", 4);
    runCase<Sacado::Fad::SFad<double, 4>, Kokkos::LayoutRight>("SFad<4> LayoutRight", 4);

    // unsupported run-time rank combination through the DynRankView shim -> exception, as before
    {
      Kokkos::DynRankView<double, Dev> out("out", 2, 3, 4), jac("jac", 2, 4, 3, 3), in("in", 3, 4, 3); // output rank 3 is invalid here
      bool threw = false;
      try { fst::HGRADtransformGRAD(out, jac, in); } catch (const std::invalid_argument&) { threw = true; }
      CHECK(threw);
    }

    // setJacobian with static-rank arguments, both (F,P,D) and (C,F,P,D) gradients, against the DynRankView result
    {
      const int C = 6, N = 8, P = 4, D = 3;
      Kokkos::View<double***, Dev> sNodes("n", C, N, D), sGrads("g", N, P, D); Kokkos::View<double****, Dev> sGradsC("gc", C, N, P, D), sJ("J", C, P, D, D), sJc("Jc", C, P, D, D);
      Kokkos::DynRankView<double, Dev> dNodes("n", C, N, D), dGrads("g", N, P, D), dGradsC("gc", C, N, P, D), dJ("J", C, P, D, D), dJc("Jc", C, P, D, D);
      for (int c = 0; c < C; ++c) for (int n = 0; n < N; ++n) for (int d = 0; d < D; ++d) { double v = std::sin(1. + c + 2 * n + 3 * d); sNodes(c, n, d) = v; dNodes(c, n, d) = v; }
      for (int n = 0; n < N; ++n) for (int p = 0; p < P; ++p) for (int d = 0; d < D; ++d) { double v = std::cos(2. + n + 3 * p + 5 * d); sGrads(n, p, d) = v; dGrads(n, p, d) = v;
        for (int c = 0; c < C; ++c) { sGradsC(c, n, p, d) = v * (1 + c); dGradsC(c, n, p, d) = v * (1 + c); } }
      ct::setJacobian(sJ, sNodes, sGrads);   ct::setJacobian(dJ, dNodes, dGrads);
      ct::setJacobian(sJc, sNodes, sGradsC); ct::setJacobian(dJc, dNodes, dGradsC);
      double e1 = 0, e2 = 0;
      for (int c = 0; c < C; ++c) for (int p = 0; p < P; ++p) for (int i = 0; i < D; ++i) for (int j = 0; j < D; ++j) {
        e1 = std::max(e1, std::abs(sJ(c, p, i, j) - dJ(c, p, i, j))); e2 = std::max(e2, std::abs(sJc(c, p, i, j) - dJc(c, p, i, j))); }
      CHECK(e1 == 0.0); CHECK(e2 == 0.0);
      std::cout << "setJacobian static vs dynamic: (F,P,D) grads " << e1 << ", (C,F,P,D) grads " << e2 << "\n";
    }
  }
  Kokkos::finalize();
  std::cout << (failures == 0 ? "End Result: TEST PASSED" : "End Result: TEST FAILED") << "\n";
  return failures;
}
