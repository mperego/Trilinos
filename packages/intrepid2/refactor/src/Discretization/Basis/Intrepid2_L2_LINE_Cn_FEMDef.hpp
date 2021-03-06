// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_L2_LINE_Cn_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on LINE.
    \author Created by M. Perego, based on the Intrepid2::HGRAD_LINE_Cn_FEM class
*/

#ifndef __INTREPID2_L2_LINE_CN_FEM_DEF_HPP__
#define __INTREPID2_L2_LINE_CN_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {
    
    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType,
             typename workViewType,
             typename vinvViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_L2_LINE_Cn_FEM::Serial<opType>::
    getValues(       outputViewType output,
               const inputViewType  input,
                     workViewType   work,
               const vinvViewType   vinv,
               const ordinal_type   operatorDn ) {    
      ordinal_type opDn = operatorDn;

      const ordinal_type card = vinv.dimension(0);
      const ordinal_type npts = input.dimension(0);

      const ordinal_type order = card - 1;
      const double alpha = 0.0, beta = 0.0;

      switch (opType) {
      case OPERATOR_VALUE: {
        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> phis(work.data(), card, npts);

        Impl::Basis_HGRAD_LINE_Cn_FEM_JACOBI::
          Serial<opType>::getValues(phis, input, order, alpha, beta);

        for (ordinal_type i=0;i<card;++i) 
          for (ordinal_type j=0;j<npts;++j) {
            output(i,j) = 0.0;
            for (ordinal_type k=0;k<card;++k)
              output(i,j) += vinv(k,i)*phis(k,j);
          }
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1:
      case OPERATOR_D2:
      case OPERATOR_D3:
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10: 
        opDn = getOperatorOrder(opType);
      case OPERATOR_Dn: {
        // dkcard is always 1 for 1D element
        const ordinal_type dkcard = 1;
        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> phis(work.data(), card, npts, dkcard);
        
        Impl::Basis_HGRAD_LINE_Cn_FEM_JACOBI::
          Serial<opType>::getValues(phis, input, order, alpha, beta, opDn);

        for (ordinal_type i=0;i<card;++i) 
          for (ordinal_type j=0;j<npts;++j) 
            for (ordinal_type k=0;k<dkcard;++k) {
              output(i,j,k) = 0.0;
              for (ordinal_type l=0;l<card;++l)
                output(i,j,k) += vinv(l,i)*phis(l,j,k);
            }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_L2_LINE_Cn_FEM::Serial::getValues) operator is not supported." );
      }
      }
    }
    

    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_L2_LINE_Cn_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef          Kokkos::DynRankView<vinvValueType,       vinvProperties...>                vinvViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // loopSize corresponds to cardinality
      const auto loopSizeTmp1 = (inputPoints.dimension(0)/numPtsPerEval);
      const auto loopSizeTmp2 = (inputPoints.dimension(0)%numPtsPerEval != 0);
      const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
      
      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv) );
        break;
      }
      case OPERATOR_GRAD:
      case OPERATOR_D1:
      case OPERATOR_D2:
      case OPERATOR_D3:
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,
            OPERATOR_Dn,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv,
                                                  getOperatorOrder(operatorType)) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_L2_LINE_Cn_FEM): Operator type not implemented" );
        //break; commented out because this always throws
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_L2_LINE_Cn_FEM<SpT,OT,PT>::
  Basis_L2_LINE_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType ) {
    this->basisCardinality_  = order+1;
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this->basisType_         = BASIS_FEM_FIAT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    const ordinal_type card = this->basisCardinality_;
    
    // points are computed in the host and will be copied 
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("L2::Line::Cn::dofCoords", card, 1);


    switch (pointType) {
    case POINTTYPE_EQUISPACED:
    case POINTTYPE_WARPBLEND: {
      // lattice ordering 
      {
        const ordinal_type offset = 1;
        PointTools::getLattice( dofCoords,
                                this->basisCellTopology_, 
                                order+1+offset, offset,
                                pointType );
        
      }
      // topological order
      // { 

      //   // internal points
      //   typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      //   auto pts = Kokkos::subview(dofCoords, range_type(2, card), Kokkos::ALL());
        
      //   const auto offset = 1;
      //   PointTools::getLattice( pts,
      //                           this->basisCellTopology_, 
      //                           order, offset, 
      //                           pointType );
      // }
      break;
    }
    case POINTTYPE_GAUSS: {
      // internal points only
      PointTools::getGaussPoints( dofCoords, 
                                  order );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( !isValidPointType(pointType),
                                    std::invalid_argument , 
                                    ">>> ERROR: (Intrepid2::Basis_L2_LINE_Cn_FEM) invalid pointType." );
    }
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
    
    // form Vandermonde matrix; actually, this is the transpose of the VDM,
    // this matrix is used in LAPACK so it should be column major and left layout
    const ordinal_type lwork = card*card;
    Kokkos::DynRankView<typename scalarViewType::value_type,Kokkos::LayoutLeft,Kokkos::HostSpace>
      vmat("L2::Line::Cn::vmat", card, card),
      work("L2::Line::Cn::work", lwork),
      ipiv("L2::Line::Cn::ipiv", card);

    const double alpha = 0.0, beta = 0.0;
    Impl::Basis_HGRAD_LINE_Cn_FEM_JACOBI::
      getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>
      (vmat, dofCoords, order, alpha, beta, OPERATOR_VALUE);

    ordinal_type info = 0;
    Teuchos::LAPACK<ordinal_type,typename scalarViewType::value_type> lapack;

    lapack.GETRF(card, card, 
                 vmat.data(), vmat.stride_1(),
                 (ordinal_type*)ipiv.data(),
                 &info);

    INTREPID2_TEST_FOR_EXCEPTION( info != 0,
                                  std::runtime_error , 
                                  ">>> ERROR: (Intrepid2::Basis_L2_LINE_Cn_FEM) lapack.GETRF returns nonzero info." );

    lapack.GETRI(card, 
                 vmat.data(), vmat.stride_1(),
                 (ordinal_type*)ipiv.data(),
                 work.data(), lwork,
                 &info);

    INTREPID2_TEST_FOR_EXCEPTION( info != 0,
                                  std::runtime_error , 
                                  ">>> ERROR: (Intrepid2::Basis_L2_LINE_Cn_FEM) lapack.GETRI returns nonzero info." );
    
    // create host mirror 
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      vinv("L2::Line::Cn::vinv", card, card);

    for (ordinal_type i=0;i<card;++i) 
      for (ordinal_type j=0;j<card;++j) 
        vinv(i,j) = vmat(j,i);

    this->vinv_ = Kokkos::create_mirror_view(typename SpT::memory_space(), vinv);
    Kokkos::deep_copy(this->vinv_ , vinv);

    // initialize tags
    {
      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
      

      ordinal_type tags[Parameters::MaxOrder+1][4];

      for (ordinal_type i=0;i<card;++i) {
        tags[i][0] = 1;    // edge dof
        tags[i][1] = 0;    // edge id
        tags[i][2] = i;    // local dof id
        tags[i][3] = card; // total number of dofs in this edge
      }

      ordinal_type_array_1d_host tagView(&tags[0][0], card*4);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      // tags are constructed on host
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }  
  }
  
}// namespace Intrepid2

#endif















