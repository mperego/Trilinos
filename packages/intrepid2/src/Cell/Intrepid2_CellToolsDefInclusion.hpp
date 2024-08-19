// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefInclusion.hpp
    \brief  Definition file for point inclusion functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__
#define __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                                        Inclusion tests                                     //
  //                                                                                            //
  //============================================================================================//

  
  template<typename DeviceType>
  template<typename pointValueType, class ...pointProperties>
  bool 
  CellTools<DeviceType>::
  checkPointInclusion( const Kokkos::DynRankView<pointValueType,pointProperties...> point,
                       const shards::CellTopology cellTopo,
                       const double               threshold) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( point.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point must have rank 1. ");
    INTREPID2_TEST_FOR_EXCEPTION( point.extent(0) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point and cell dimensions do not match. ");
#endif
    bool testResult = true;
  
    const double t = threshold; //(threshold < 0 ? threshold() : threshold); 

    // Using these values in the tests effectievly inflates the reference element to a larger one
    const double minus_one =  -1.0 - t;
    const double plus_one  =   1.0 + t;
    const double minus_zero = -t;

    // A cell with extended topology has the same reference cell as a cell with base topology. 
    // => testing for inclusion in a reference Triangle<> and a reference Triangle<6> relies on 
    // on the same set of inequalities. To eliminate unnecessary cases we switch on the base topology
    const auto key = cellTopo.getBaseKey();
    switch (key) {
    
    case shards::Line<>::key :
      if( !(minus_one <= point(0) && point(0) <= plus_one))  testResult = false;
      break;
      
    case shards::Triangle<>::key : {
      const auto distance = Util<pointValueType>::max( std::max( -point(0), -point(1) ), point(0) + point(1) - 1.0 );
      if( distance > threshold ) testResult = false;
      break;
    }
      
    case shards::Quadrilateral<>::key :
      if(!( (minus_one <= point(0) && point(0) <= plus_one) &&          
            (minus_one <= point(1) && point(1) <= plus_one) ) ) testResult = false;   
      break;
      
    case shards::Tetrahedron<>::key : {
      const auto distance = Util<pointValueType>::max(  Util<pointValueType>::max(-point(0),-point(1)), 
                                        Util<pointValueType>::max(-point(2), point(0) + point(1) + point(2) - 1)  );
      if( distance > threshold ) testResult = false;
      break;
    }
      
    case shards::Hexahedron<>::key :
      if(!((minus_one <= point(0) && point(0) <= plus_one) && 
           (minus_one <= point(1) && point(1) <= plus_one) && 
           (minus_one <= point(2) && point(2) <= plus_one)))  
        testResult = false;
      break;
      
      // The base of the reference prism is the same as the reference triangle => apply triangle test
      // to X and Y coordinates and test whether Z is in [-1,1]
    case shards::Wedge<>::key : {
      const auto distance = Util<pointValueType>::max( Util<pointValueType>::max( -point(0), -point(1) ), point(0) + point(1) - 1 );
      if( distance > threshold ||                     
          point(2) < minus_one || point(2) > plus_one) 
        testResult = false;
      break;
    }

      // The base of the reference pyramid is the same as the reference quad cell => a horizontal plane
      // through a point P(x,y,z) intersects the pyramid at a quadrilateral that equals the base quad 
      // scaled by (1-z) => P(x,y,z) is inside the pyramid <=> (x,y) is in [-1+z,1-z]^2 && 0 <= Z <= 1 
    case shards::Pyramid<>::key : {
      const auto left  = minus_one + point(2);
      const auto right = plus_one  - point(2);
      if(!( (left       <= point(0) && point(0) <= right) && \
            (left       <= point(1) && point(1) <= right) && 
            (minus_zero <= point(2) && point(2) <= plus_one) ) )  \
        testResult = false;  
      break;
    }
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( !( (key == shards::Line<>::key ) ||
                                       (key == shards::Triangle<>::key)  ||
                                       (key == shards::Quadrilateral<>::key) ||
                                       (key == shards::Tetrahedron<>::key)  ||
                                       (key == shards::Hexahedron<>::key)  ||
                                       (key == shards::Wedge<>::key)  ||
                                       (key == shards::Pyramid<>::key) ),
                                    std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Invalid cell topology. ");
    }
    return testResult;
  }



  template<typename cellTopologyTagType,
             typename OutputViewType,
             typename inputViewType>
    struct checkPointInclusionFunctor {
      OutputViewType output_;
      inputViewType input_;
      double threshold_;

      KOKKOS_INLINE_FUNCTION
      checkPointInclusionFunctor(OutputViewType output,
                            inputViewType input,
                            double threshold)
        : output_(output), 
          input_(input),
          threshold_(threshold) {}

      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type i) const {
        const auto in = Kokkos::subview(input_,i,Kokkos::ALL());
        const auto check = cellTopologyTagType::checkPointInclusion(in, threshold_);
        output_(i) = check;        
      }
      
      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type i, const ordinal_type j) const {
        const auto in = Kokkos::subview(input_,i,j,Kokkos::ALL());
        const auto check = cellTopologyTagType::checkPointInclusion(in, threshold_);
        output_(i,j) = check;        
      }
    };

  template<typename DeviceType>
  template<typename cellTopologyTagType,
           typename OutputViewType,
           typename inputViewType>
  void CellTools<DeviceType>::
  checkPointwiseInclusion(OutputViewType inCell, 
                                      const inputViewType points,
                                      const double threshold) {
      
    using FunctorType = checkPointInclusionFunctor<cellTopologyTagType,decltype(inCell),decltype(points)>; 
    if (points.rank() == 2) {
      Kokkos::RangePolicy<typename DeviceType::execution_space> policy(0, points.extent(0));
      Kokkos::parallel_for(policy, FunctorType(inCell, points, threshold));
    } else {  //points.rank() == 3
      Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<2>> policy({0,0},{points.extent(0),points.extent(1)});
      Kokkos::parallel_for(policy, FunctorType(inCell, points, threshold));
    }
  }


  template<typename DeviceType>
  template<typename inCellValueType, class ...inCellProperties,
           typename pointValueType, class ...pointProperties>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,
                           const Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const shards::CellTopology cellTopo,
                           const double threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inCell.rank() != (points.rank()-1), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank difference between inCell and points is 1.");  
      const ordinal_type iend = inCell.rank();
      for (ordinal_type i=0;i<iend;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inCell.extent(i) != points.extent(i), std::invalid_argument, 
                                      ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dimension mismatch between inCell and points.");  
      }
    }
#endif

   const auto key = cellTopo.getBaseKey();
   switch (key) {
    
    case shards::Line<>::key :
      checkPointwiseInclusion<Impl::Line<2>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Triangle<>::key :
      checkPointwiseInclusion<Impl::Triangle<3>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Quadrilateral<>::key :
      checkPointwiseInclusion<Impl::Quadrilateral<4>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Hexahedron<>::key :
      checkPointwiseInclusion<Impl::Hexahedron<8>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Wedge<>::key :
      checkPointwiseInclusion<Impl::Wedge<6>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Pyramid<>::key :
      checkPointwiseInclusion<Impl::Pyramid<5>,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( !( (key == shards::Line<>::key ) ||
                                       (key == shards::Triangle<>::key)  ||
                                       (key == shards::Quadrilateral<>::key) ||
                                       (key == shards::Tetrahedron<>::key)  ||
                                       (key == shards::Hexahedron<>::key)  ||
                                       (key == shards::Wedge<>::key)  ||
                                       (key == shards::Pyramid<>::key) ),
                                    std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Invalid cell topology. ");
    }
  }

  template<typename DeviceType>
  template<typename inCellValueType, class ...inCellProperties,
           typename pointValueType, class ...pointProperties,
           typename cellWorksetValueType, class ...cellWorksetProperties>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,
                           const Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const Kokkos::DynRankView<cellWorksetValueType,cellWorksetProperties...> cellWorkset,
                           const shards::CellTopology cellTopo,
                           const double threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      const auto key = cellTopo.getBaseKey();
      INTREPID2_TEST_FOR_EXCEPTION( key != shards::Line<>::key &&
                                    key != shards::Triangle<>::key &&
                                    key != shards::Quadrilateral<>::key &&
                                    key != shards::Tetrahedron<>::key &&                                                                                                                                               
                                    key != shards::Hexahedron<>::key &&                                                                                                                                                
                                    key != shards::Wedge<>::key &&                                                                                                                                                     
                                    key != shards::Pyramid<>::key, 
                                    std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cell topology not supported");

      INTREPID2_TEST_FOR_EXCEPTION( points.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): Points must have rank 3. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cellWorkset must have rank 3. ");
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Points and cell dimensions do not match. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and cell dimensions do not match. ");
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(0) != cellWorkset.extent(0) , std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and points dimension(0) does not match. ");
    }
#endif    
    const ordinal_type 
      numCells = cellWorkset.extent(0),
      numPoints = points.extent(1), 
      spaceDim = cellTopo.getDimension();

    using result_layout = typename DeduceLayout< decltype(points) >::result_layout;
    auto vcprop = Kokkos::common_view_alloc_prop(points);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, result_layout, DeviceType > refPoints ( Kokkos::view_alloc("CellTools::checkPointwiseInclusion::refPoints", vcprop), numCells, numPoints, spaceDim);
    
    // expect refPoints(CPD), points(CPD), cellWorkset(CND) 
    mapToReferenceFrame(refPoints, points, cellWorkset, cellTopo);
    checkPointwiseInclusion(inCell, refPoints, cellTopo, threshold);
  }

}

#endif
