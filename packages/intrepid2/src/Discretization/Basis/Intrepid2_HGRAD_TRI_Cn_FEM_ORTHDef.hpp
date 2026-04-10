// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TRI_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of arbitrary degree
    for H(grad) functions on TRI.
    \author Created by R. Kirby
    Kokkorized by Kyungjoo Kim and Mauro Perego
 */

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_ORTH_DEF_HPP__

namespace Intrepid2 {
// -------------------------------------------------------------------------------------

namespace Impl {

template<typename OutputViewType,
typename InputViewType,
typename WorkViewType,
bool hasDeriv>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTri<OutputViewType,InputViewType,WorkViewType,hasDeriv,0>::generate(
    OutputViewType output,
    const InputViewType input,
    WorkViewType  /*work*/,
    const ordinal_type order ) {

  constexpr ordinal_type spaceDim = 2;

  typedef typename InputViewType::value_type value_type;

  auto output0 = (hasDeriv) ? Kokkos::subview(output, Kokkos::ALL(),0) : Kokkos::subview(output, Kokkos::ALL());

  // each point needs to be transformed from Pavel's element
  // z(0) --> (2.0 * input(0) - 1.0)
  // z(1) --> (2.0 * input(1) - 1.0)


  // set D^{0,0} = 1.0
  {
    const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(0,0);
    output0(loc) = 1.0;
    if(hasDeriv) {
      output.access(loc,1) = 0;
      output.access(loc,2) = 0;
    }
  }

  if (order > 0) {
    const value_type z0 = 2.0*input(0)-1.0;
    const value_type z1 = 2.0*input(1)-1.0;
    const value_type f1 = z0 + input(1);   // z_0 + (1 + z_1)/2
    const value_type f2 = pow(input(1)-1,2);  //( (1 - z_1)/2 )^2
    
    const value_type df1_0(2.0), df1_1(1.0), df2_1(z1);    

    // set D^{1,0} = f1
    {
      const ordinal_type loc = Intrepid2::getPnEnumeration<spaceDim>(1,0);
      output0(loc) = f1;
      if(hasDeriv) {
        output(loc,1) = df1_0;
        output(loc,2) = df1_1;
      }
    }

    // recurrence in p
    for (ordinal_type p=1;p<order;p++) {
      const ordinal_type
      loc = Intrepid2::getPnEnumeration<spaceDim>(p,0),
      loc_p1 = Intrepid2::getPnEnumeration<spaceDim>(p+1,0),
      loc_m1 = Intrepid2::getPnEnumeration<spaceDim>(p-1,0);

      const value_type a = (2.0*p+1.0)/(1.0+p);
      const value_type b = p / (p+1.0);

      output0(loc_p1) = a * f1 * output0(loc) - b * f2 * output0(loc_m1);
      if(hasDeriv) {
        output(loc_p1,1) =  a * (f1 * output(loc,1) + df1_0 * output0(loc))  -
            b * f2 * output(loc_m1,1) ;
        output(loc_p1,2) =  a * (f1 * output(loc,2) + df1_1 * output0(loc))  -
            b * (df2_1 * output0(loc_m1) + f2 * output(loc_m1,2)) ;
      }
    }

    // D^{p,1}
    for (ordinal_type p=0;p<order;++p) {
      const ordinal_type
      loc_p_0 = Intrepid2::getPnEnumeration<spaceDim>(p,0),
      loc_p_1 = Intrepid2::getPnEnumeration<spaceDim>(p,1);

      output0(loc_p_1) = output0(loc_p_0)*0.5*(1.0+2.0*p+(3.0+2.0*p)*z1);
      if(hasDeriv) {
        output(loc_p_1,1) = output(loc_p_0,1)*0.5*(1.0+2.0*p+(3.0+2.0*p)*z1);
        output(loc_p_1,2) = output(loc_p_0,2)*0.5*(1.0+2.0*p+(3.0+2.0*p)*z1) + output0(loc_p_0)*(3.0+2.0*p);
      }
    }


    // recurrence in q
    for (ordinal_type p=0;p<order-1;++p)
      for (ordinal_type q=1;q<order-p;++q) {
        const ordinal_type
        loc_p_qp1 = Intrepid2::getPnEnumeration<spaceDim>(p,q+1),
        loc_p_q = Intrepid2::getPnEnumeration<spaceDim>(p,q),
        loc_p_qm1 = Intrepid2::getPnEnumeration<spaceDim>(p,q-1);

        value_type a,b,c;
        Intrepid2::getJacobyRecurrenceCoeffs(a,b,c, 2*p+1,0,q);
        output0(loc_p_qp1) =  (a*z1+b)*output0(loc_p_q)
                  - c*output0(loc_p_qm1) ;
        if(hasDeriv) {
          output(loc_p_qp1,1) =  (a*z1+b)*output(loc_p_q,1) - c*output(loc_p_qm1,1) ;
          output(loc_p_qp1,2) =  (a*z1+b)*output(loc_p_q,2) +2*a*output0(loc_p_q) - c*output(loc_p_qm1 ,2) ;
        }
      }
  }

  // orthogonalize
  for (ordinal_type p=0;p<=order;++p)
    for (ordinal_type q=0;q<=order-p;++q){
      output0(Intrepid2::getPnEnumeration<spaceDim>(p,q)) *= std::sqrt( (p+0.5)*(p+q+1.0));
      if(hasDeriv) {
        output.access(Intrepid2::getPnEnumeration<spaceDim>(p,q),1) *= std::sqrt( (p+0.5)*(p+q+1.0));
        output.access(Intrepid2::getPnEnumeration<spaceDim>(p,q),2) *= std::sqrt( (p+0.5)*(p+q+1.0));
      }
    }
}

template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTri<OutputViewType,inputViewType,workViewType,hasDeriv,1>::generate(
    OutputViewType output,
    const inputViewType  input,
    workViewType   work,
    const ordinal_type   order ) {
  constexpr ordinal_type spaceDim = 2;
  const ordinal_type  card = output.extent(0);

  workViewType dummyView;
  OrthPolynomialTri<workViewType,inputViewType,workViewType,hasDeriv,0>::generate(work, input, dummyView, order);
  for (ordinal_type i=0;i<card;++i)
    for (ordinal_type k=0;k<spaceDim;++k)
      output.access(i,k) = work(i,k+1);
}


// when n >= 2, use recursion
template<typename OutputViewType,
typename inputViewType,
typename workViewType,
bool hasDeriv,
ordinal_type n>
KOKKOS_INLINE_FUNCTION
void OrthPolynomialTri<OutputViewType,inputViewType,workViewType,hasDeriv,n>::generate(
    OutputViewType /* output */,
    const inputViewType  /* input */,
    workViewType   /* work */,
    const ordinal_type   /* order */ ) {
INTREPID2_TEST_FOR_ABORT( true,
    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::OrthPolynomialTri) Computing of second and higher-order derivatives is not currently supported");
}



template<EOperator opType>
template<typename OutputViewType,
typename inputViewType,
typename workViewType>
KOKKOS_INLINE_FUNCTION
void
Basis_HGRAD_TRI_Cn_FEM_ORTH::Serial<opType>::
getValues( OutputViewType output,
    const inputViewType  input,
    workViewType   work,
    const ordinal_type   order) {
  switch (opType) {
  case OPERATOR_VALUE: {
    OrthPolynomialTri<OutputViewType,inputViewType,workViewType,false,0>::generate( output, input, work, order );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    OrthPolynomialTri<OutputViewType,inputViewType,workViewType,true,1>::generate( output, input, work, order );
    break;
  }
  case OPERATOR_D2: {
    OrthPolynomialTri<OutputViewType,inputViewType,workViewType,true,2>::generate( output, input, work, order );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::Serial::getValues) operator is not supported");
  }
  }
}

template<typename DT,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties>
void
Basis_HGRAD_TRI_Cn_FEM_ORTH::
getValues(
    const typename DT::execution_space& space,
          Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
    const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
    const ordinal_type order,
    const EOperator operatorType ) {
  typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
  typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
  typedef typename DT::execution_space ExecSpaceType;

  // loopSize corresponds to the # of points
  const auto loopSize = inputPoints.extent(0);
  Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

  const ordinal_type cardinality = outputValues.extent(0);
  const ordinal_type spaceDim = 2;

  typedef typename DeduceDynRankView<inputPointViewType>::type workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType  dummyWorkView;
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_VALUE> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, dummyWorkView, order) );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    workViewType work = createMatchingView<workViewType>(inputPoints, "Basis_HGRAD_TRI_In_FEM_ORTH::getValues::work", cardinality, inputPoints.extent(0), spaceDim+1);
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_D1> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, work, order) );
    break;
  }
  case OPERATOR_D2:{
    workViewType  dummyWorkView;
    typedef Functor<outputValueViewType,inputPointViewType,workViewType,OPERATOR_D2> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints ,dummyWorkView, order) );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( !Intrepid2::isValidOperator(operatorType), std::invalid_argument,
        ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM_ORTH): invalid operator type");
  }
  }
}
}


// -------------------------------------------------------------------------------------

template<typename DT, typename OT, typename PT>
Basis_HGRAD_TRI_Cn_FEM_ORTH<DT,OT,PT>::
Basis_HGRAD_TRI_Cn_FEM_ORTH( const ordinal_type order ) {

  constexpr ordinal_type spaceDim = 2;
  this->basisCardinality_     = Intrepid2::getPnCardinality<spaceDim>(order);
  this->basisDegree_          = order;
  this->basisCellTopologyKey_ = shards::Triangle<3>::key;
  this->basisType_            = BASIS_FEM_HIERARCHICAL;
  this->basisCoordinates_     = COORDINATES_CARTESIAN;
  this->functionSpace_        = FUNCTION_SPACE_HGRAD;

  // initialize tags
  {
    // Basis-dependent initializations
    constexpr ordinal_type tagSize  = 4;    // size of DoF tag, i.e., number of fields in the tag
    const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
    const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
    ordinal_type tags[maxCard][tagSize];
    const ordinal_type card = this->basisCardinality_;
    for (ordinal_type i=0;i<card;++i) {
      tags[i][0] = 2;     // these are all "internal" i.e. "volume" DoFs
      tags[i][1] = 0;     // there is only one line
      tags[i][2] = i;     // local DoF id
      tags[i][3] = card;  // total number of DoFs
    }

    OrdinalTypeArray1DHost tagView(&tags[0][0], card*tagSize);

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

  // dof coords is not applicable to hierarchical functions
}

}
#endif
