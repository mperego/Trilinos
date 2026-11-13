// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file  bench_static_rank_views.cpp
    \brief PROTOTYPE.  Times a few Intrepid2 kernels called with Kokkos::DynRankView arguments and with
           static-rank Kokkos::View arguments, for double and for Sacado Fad scalars.

    To obtain the "before" numbers, build this same file against an unpatched Intrepid2 with
    -DINTREPID2_BENCH_UNPATCHED: it then skips the calls that only exist (or only compile) with the patch.
*/

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_Timer.hpp"
#include "Sacado.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

  using Dev = Kokkos::DefaultExecutionSpace::device_type;

  template<class F>
  double bestOf(const int reps, F &&f) {
    double best = std::numeric_limits<double>::max();
    for (int r=0;r<reps;++r) {
      Kokkos::Timer timer;
      f();
      Kokkos::fence();
      best = std::min(best, timer.seconds());
    }
    return best;
  }

  // allocates a view with the given logical dimensions, adding the hidden derivative dimension for Fad scalars
  template<class ViewType, class... Dims>
  ViewType alloc(const std::string &label, const int numDerivs, Dims... dims) {
    if constexpr (std::is_arithmetic_v<typename ViewType::value_type>)
      return ViewType(label, dims...);
    else
      return ViewType(label, dims..., numDerivs + 1);
  }

  void report(const std::string &scalar, const std::string &kernel, const double tDyn, const double tStatic) {
    std::cout << std::left << std::setw(9) << scalar << std::setw(30) << kernel << std::right << std::fixed << std::setprecision(3)
              << std::setw(10) << tDyn*1e3;
    if (tStatic >= 0) std::cout << std::setw(10) << tStatic*1e3 << std::setw(8) << std::setprecision(2) << tDyn/tStatic << "x";
    else              std::cout << std::setw(10) << "n/a";
    std::cout << "\n";
  }

  template<class Scalar>
  void run(const std::string &scalar, const int C, const int reps, const int numDerivs) {
    using rst = Intrepid2::RealSpaceTools<Dev>;
    using ct  = Intrepid2::CellTools<Dev>;
    using fst = Intrepid2::FunctionSpaceTools<Dev>;

    const int N = 8, F = 8, P = 8, D = 3; // hexahedra, 2x2x2 points
    using Dyn = Kokkos::DynRankView<Scalar,Dev>;
    using S2  = Kokkos::View<Scalar**,  Dev>;
    using S3  = Kokkos::View<Scalar***, Dev>;
    using S4  = Kokkos::View<Scalar****,Dev>;

    Dyn dNodes = alloc<Dyn>("nodes", numDerivs, C, N, D),     dRefGrads = alloc<Dyn>("refGrads", numDerivs, F, P, D),
        dJac   = alloc<Dyn>("jac",   numDerivs, C, P, D, D),  dJacInv   = alloc<Dyn>("jacInv",   numDerivs, C, P, D, D),
        dDet   = alloc<Dyn>("det",   numDerivs, C, P),        dGrads    = alloc<Dyn>("grads",    numDerivs, C, F, P, D),
        dPts   = alloc<Dyn>("pts",   numDerivs, P*C, D),      dVals     = alloc<Dyn>("vals",     numDerivs, F, P*C),
        dBasisGrads = alloc<Dyn>("basisGrads", numDerivs, F, P*C, D);
    S3  sNodes = alloc<S3>("nodes",  numDerivs, C, N, D),     sRefGrads = alloc<S3>("refGrads",  numDerivs, F, P, D);
    S4  sJac   = alloc<S4>("jac",    numDerivs, C, P, D, D),  sJacInv   = alloc<S4>("jacInv",    numDerivs, C, P, D, D),
        sGrads = alloc<S4>("grads",  numDerivs, C, F, P, D);
    S2  sDet   = alloc<S2>("det",    numDerivs, C, P),        sPts      = alloc<S2>("pts",       numDerivs, P*C, D),
        sVals  = alloc<S2>("vals",   numDerivs, F, P*C);
    S3  sBasisGrads = alloc<S3>("basisGrads", numDerivs, F, P*C, D);

    // perturbed reference hexahedra; trilinear reference gradients at pseudo-random points
    const double ref[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    for (int c=0;c<C;++c) for (int n=0;n<N;++n) for (int d=0;d<D;++d) {
      const double v = ref[n][d] + 0.05*std::sin(1.0 + c + 3*n + 7*d);
      dNodes(c,n,d) = v; sNodes(c,n,d) = v;
    }
    for (int p=0;p<P;++p) {
      const double x[3] = {0.6*std::sin(1.0 + p), 0.6*std::cos(2.0 + p), 0.6*std::sin(3.0 + 2*p)};
      for (int n=0;n<N;++n) for (int d=0;d<D;++d) {
        double g = 0.125*ref[n][d];
        for (int e=0;e<D;++e) if (e != d) g *= (1.0 + ref[n][e]*x[e]);
        dRefGrads(n,p,d) = g; sRefGrads(n,p,d) = g;
      }
    }
    for (int p=0;p<P*C;++p) for (int d=0;d<D;++d) {
      const double v = 0.9*std::sin(1.0 + p + 3*d);
      dPts(p,d) = v; sPts(p,d) = v;
    }

    // ---- DynRankView arguments
    const double tJacDyn = bestOf(reps, [&] { ct::setJacobian(dJac, dNodes, dRefGrads); });
    const double tDetDyn = bestOf(reps, [&] { rst::det(dDet, dJac); });
    const double tInvDyn = bestOf(reps, [&] { rst::inverse(dJacInv, dJac); });
    const double tTrfDyn = bestOf(reps, [&] { fst::HGRADtransformGRAD(dGrads, dJacInv, dRefGrads); });
    Intrepid2::Basis_HGRAD_HEX_C1_FEM<Dev,Scalar,Scalar> hexBasis;
    const Intrepid2::Basis<Dev,Scalar,Scalar> &basis = hexBasis; // call through the base class, as applications do
    const double tValDyn = bestOf(reps, [&] { basis.getValues(dVals,       dPts, Intrepid2::OPERATOR_VALUE); });
    const double tGrdDyn = bestOf(reps, [&] { basis.getValues(dBasisGrads, dPts, Intrepid2::OPERATOR_GRAD); });

    // ---- static-rank View arguments
    double tJacSta = -1, tTrfSta = -1, tValSta = -1, tGrdSta = -1;
#ifdef INTREPID2_BENCH_UNPATCHED
    Kokkos::deep_copy(sJac, Kokkos::Impl::as_view_of_rank_n<4>(dJac));
#else
    tJacSta = bestOf(reps, [&] { ct::setJacobian(sJac, sNodes, sRefGrads); });
#endif
    const double tDetSta = bestOf(reps, [&] { rst::det(sDet, sJac); });
    const double tInvSta = bestOf(reps, [&] { rst::inverse(sJacInv, sJac); });
#ifndef INTREPID2_BENCH_UNPATCHED
    tTrfSta = bestOf(reps, [&] { fst::HGRADtransformGRAD(sGrads, sJacInv, sRefGrads); });
    tValSta = bestOf(reps, [&] { basis.getValues(sVals,       sPts, Intrepid2::OPERATOR_VALUE); });
    tGrdSta = bestOf(reps, [&] { basis.getValues(sBasisGrads, sPts, Intrepid2::OPERATOR_GRAD); });
#endif

    report(scalar, "CellTools::setJacobian",        tJacDyn, tJacSta);
    report(scalar, "RealSpaceTools::det",           tDetDyn, tDetSta);
    report(scalar, "RealSpaceTools::inverse",       tInvDyn, tInvSta);
    report(scalar, "FST::HGRADtransformGRAD",       tTrfDyn, tTrfSta);
    report(scalar, "HGRAD_HEX_C1 getValues VALUE",  tValDyn, tValSta);
    report(scalar, "HGRAD_HEX_C1 getValues GRAD",   tGrdDyn, tGrdSta);

    // the two containers must hold the same numbers
    double maxDiff = 0;
    for (int c=0;c<C;++c) for (int p=0;p<P;++p) {
      maxDiff = std::max(maxDiff, std::abs(Sacado::scalarValue(Scalar(dDet(c,p))) - Sacado::scalarValue(Scalar(sDet(c,p)))));
      for (int i=0;i<D;++i) for (int j=0;j<D;++j)
        maxDiff = std::max(maxDiff, std::abs(Sacado::scalarValue(Scalar(dJacInv(c,p,i,j))) - Sacado::scalarValue(Scalar(sJacInv(c,p,i,j)))));
#ifndef INTREPID2_BENCH_UNPATCHED
      for (int f=0;f<F;++f) for (int i=0;i<D;++i)
        maxDiff = std::max(maxDiff, std::abs(Sacado::scalarValue(Scalar(dGrads(c,f,p,i))) - Sacado::scalarValue(Scalar(sGrads(c,f,p,i)))));
#endif
    }
    std::cout << "         max |DynRankView result - View result| = " << std::scientific << std::setprecision(1) << maxDiff << "\n";
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const int C    = argc > 1 ? std::atoi(argv[1]) : 8000;
    const int reps = argc > 2 ? std::atoi(argv[2]) : 21;
#ifdef INTREPID2_BENCH_UNPATCHED
    std::cout << "Intrepid2 WITHOUT the static-rank prototype";
#else
    std::cout << "Intrepid2 WITH the static-rank prototype";
#endif
    std::cout << "; " << C << " hexahedra, 8 points per cell; best of " << reps << " runs; times in ms\n";
    std::cout << std::left << std::setw(9) << "scalar" << std::setw(30) << "kernel" << std::right
              << std::setw(10) << "DynRank" << std::setw(10) << "View" << std::setw(9) << "ratio" << "\n";
    run<double>("double", C, reps, 0);
    run<Sacado::Fad::SFad<double,8>>("SFad<8>", C/4, reps, 8);
    run<Sacado::Fad::DFad<double>>("DFad(8)", C/4, reps, 8);
  }
  Kokkos::finalize();
  return 0;
}
