// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ViewConcepts.hpp
    \brief  PROTOTYPE.  Concepts that distinguish static-rank Kokkos::View arguments from
            Kokkos::DynRankView arguments, and the machinery used by the (deprecated) DynRankView
            overloads to forward to the static-rank implementations:

              - toStaticView<N>(dynRankView): aliasing conversion to a Kokkos::View of rank N;
              - dispatchStaticRank(name, f, views...): converts every DynRankView in views... to
                the static-rank View matching its run-time rank and calls f with the result.  Only
                the rank combinations f accepts (i.e. that satisfy the constraints of the static
                overload f forwards to) are instantiated; any other combination found at run time
                throws std::invalid_argument.
*/

#ifndef __INTREPID2_VIEW_CONCEPTS_HPP__
#define __INTREPID2_VIEW_CONCEPTS_HPP__

#include <concepts>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Utils.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

/** Attribute placed on the DynRankView overloads that have a static-rank replacement.  Warnings are
    opt-in for now, so that the deprecation can be switched on package by package once the callers
    inside Trilinos have been moved over. */
#ifdef INTREPID2_SHOW_DYNRANKVIEW_DEPRECATION_WARNINGS
#define INTREPID2_DEPRECATED_DYNRANKVIEW \
  [[deprecated("Kokkos::DynRankView arguments are deprecated here; pass static-rank Kokkos::View arguments instead.")]]
#else
#define INTREPID2_DEPRECATED_DYNRANKVIEW
#endif

namespace Intrepid2 {

  //! A Kokkos::View: its rank is known at compile time.
  template<class V>
  concept StaticRankView = Kokkos::is_view_v<std::remove_cvref_t<V>>;

  //! A Kokkos::DynRankView: its rank is known only at run time.
  template<class V>
  concept DynamicRankView = Kokkos::is_dyn_rank_view<std::remove_cvref_t<V>>::value;

  //! A Kokkos::View whose rank is one of \c Ranks.
  template<class V, std::size_t... Ranks>
  concept ViewOfRank = StaticRankView<V> && ((std::remove_cvref_t<V>::rank() == Ranks) || ...);

  namespace Impl {

    template<class V>
    struct StaticRank : std::integral_constant<int, -1> {};

    template<StaticRankView V>
    struct StaticRank<V> : std::integral_constant<int, static_cast<int>(V::rank())> {};

    //! The rank of \c V when it is known at compile time (Kokkos::View), -1 otherwise (DynRankView, functors, Data, ...).
    template<class V>
    inline constexpr int static_rank_v = StaticRank<std::remove_cvref_t<V>>::value;

    template<class View7, std::size_t... Kept, std::size_t... Dropped>
    KOKKOS_INLINE_FUNCTION
    auto sliceLeadingDims(const View7 &view, std::index_sequence<Kept...>, std::index_sequence<Dropped...>) {
      return Kokkos::subview(view, ((void)Kept, Kokkos::ALL())..., ((void)Dropped, 0)...);
    }

    /** \brief Returns a Kokkos::View of rank \c N that aliases the data of \c view.  Aborts if \c view.rank() != N.

        For LayoutLeft/LayoutRight this defers to Kokkos, and the result keeps the (contiguous) layout.
        For LayoutStride we slice the rank-7 View that underlies the DynRankView instead: in the Kokkos
        version this prototype was written against (5.2.2), Kokkos::Impl::as_view_of_rank_n places the
        stride of the hidden derivative dimension of Sacado Fad types incorrectly for LayoutStride.
    */
    template<std::size_t N, class T, class... P>
    KOKKOS_INLINE_FUNCTION
    auto toStaticView(const Kokkos::DynRankView<T, P...> &view) {
      using dyn_view_type = Kokkos::DynRankView<T, P...>;
      if constexpr (std::is_same_v<typename dyn_view_type::array_layout, Kokkos::LayoutStride>) {
        if (view.rank() != N) Kokkos::abort("Intrepid2::Impl::toStaticView: rank of DynRankView does not match the requested static rank");
        return sliceLeadingDims(view.DownCast(), std::make_index_sequence<N>(), std::make_index_sequence<7 - N>());
      }
      else {
        return Kokkos::Impl::as_view_of_rank_n<N>(view);
      }
    }

    /** \brief Returns a DynRankView of type \c DynView that aliases the data of the static-rank \c view.

        Goes through a DynRankView with the layout of \c view first: in the Kokkos version this prototype
        was written against (5.2.2), constructing a LayoutStride DynRankView directly from a
        LayoutLeft/LayoutRight View does not compile.
    */
    template<DynamicRankView DynView, StaticRankView StaticView>
    DynView toDynRankView(const StaticView &view) {
      using same_layout_dyn_view = Kokkos::DynRankView<typename StaticView::value_type,
                                                       typename StaticView::array_layout,
                                                       typename StaticView::device_type,
                                                       typename StaticView::memory_traits>;
      return DynView(same_layout_dyn_view(view));
    }

    //! Restricts the ranks dispatchStaticRank() tries for a DynRankView argument to [Lo,Hi]; see ranks().
    template<class DynView, std::size_t Lo, std::size_t Hi>
    struct RankBounded {
      const DynView &view;
    };

    /** \brief Declares that the DynRankView \c view may only have a rank in [Lo,Hi].  Use this for every
               DynRankView handed to dispatchStaticRank(): it bounds the number of rank combinations the
               compiler has to consider. */
    template<std::size_t Lo, std::size_t Hi, DynamicRankView DynView>
    RankBounded<DynView, Lo, Hi> ranks(const DynView &view) { return {view}; }

    template<class T>
    struct IsRankBounded : std::false_type {};
    template<class DynView, std::size_t Lo, std::size_t Hi>
    struct IsRankBounded<RankBounded<DynView, Lo, Hi>> : std::true_type {};

    // all arguments converted: call f if it accepts this combination of static ranks
    template<class F, class... Converted>
    bool dispatchStaticRankImpl(F &f, std::tuple<Converted...> &&converted) {
      if constexpr (std::invocable<F &, const Converted &...>) {
        std::apply(f, std::as_const(converted));
        return true;
      }
      else {
        return false;
      }
    }

    template<class F, class... Converted, class Head, class... Tail>
    bool dispatchStaticRankImpl(F &f, std::tuple<Converted...> &&converted, const Head &head, const Tail &...tail);

    // tries rank N, N+1, ..., Hi for the DynRankView head
    template<std::size_t N, std::size_t Hi, class F, class... Converted, class DynView, class... Tail>
    bool dispatchRankRange(F &f, std::tuple<Converted...> &&converted, const DynView &head, const Tail &...tail) {
      if (head.rank() == N)
        return dispatchStaticRankImpl(f, std::tuple_cat(std::move(converted), std::make_tuple(toStaticView<N>(head))), tail...);
      if constexpr (N < Hi)
        return dispatchRankRange<N + 1, Hi>(f, std::move(converted), head, tail...);
      else
        return false;
    }

    template<class F, class... Converted, class Head, class... Tail>
    bool dispatchStaticRankImpl(F &f, std::tuple<Converted...> &&converted, const Head &head, const Tail &...tail) {
      if constexpr (IsRankBounded<Head>::value) {
        return [&]<class DynView, std::size_t Lo, std::size_t Hi>(const RankBounded<DynView, Lo, Hi> &bounded) {
          return dispatchRankRange<Lo, Hi>(f, std::move(converted), bounded.view, tail...);
        }(head);
      }
      else if constexpr (DynamicRankView<Head>) {
        return dispatchRankRange<0, 7>(f, std::move(converted), head, tail...);
      }
      else { // already a static-rank View, or not a view at all: pass through unchanged
        return dispatchStaticRankImpl(f, std::tuple_cat(std::move(converted), std::tuple<Head>(head)), tail...);
      }
    }

    /** \brief Converts each DynRankView in \c views to the Kokkos::View of matching rank and calls \c f.

        \param functionName [in] - used in the exception message
        \param f            [in] - callable; must be constrained (e.g. a generic lambda with a requires-clause),
                                   so that it is only invocable for the supported rank combinations
        \param views        [in] - DynRankViews (optionally wrapped in ranks<Lo,Hi>()), static Views, or other arguments
    */
    template<class F, class... Views>
    void dispatchStaticRank(const char *functionName, F &&f, const Views &...views) {
      const bool called = dispatchStaticRankImpl(f, std::tuple<>(), views...);
      if (!called) {
        const std::string message = std::string(">>> ERROR (") + functionName + "): unsupported combination of DynRankView ranks.";
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument, message.c_str() );
      }
    }

  } // namespace Impl
} // namespace Intrepid2

#endif
