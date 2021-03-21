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


/** \file   Intrepid2_CellToolsDefParametrization.hpp
    \brief  Definition file for the parameterization functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_PARAMETRIZATION_HPP__
#define __INTREPID2_CELLTOOLS_DEF_PARAMETRIZATION_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //          Parametrization and hgrad objects caching                                         //
  //                                                                                            //
  //============================================================================================//


  template<typename DeviceType>
  void
  Impl::CellTools<DeviceType>::
  setSubcellParametrization() {
    if(isSubcellParametrizationSet_)
      return;

    ordinal_type subcellDim;
    { 
      const auto tet = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >());

      subcellDim = 2;
      subcellParamData_.tetFaces = subcellParamViewType("CellTools::SubcellParametrization::tetFaces", tet.getSubcellCount(subcellDim), tet.getDimension(), subcellDim+1);
      auto subcell2dParamHost = Kokkos::create_mirror_view(subcellParamData_.tetFaces);
      setSubcellParametrization( subcell2dParamHost, subcellDim, tet );
      deep_copy(subcellParamData_.tetFaces,subcell2dParamHost);

      subcellDim = 1;
      subcellParamData_.tetEdges = subcellParamViewType("CellTools::SubcellParametrization::tetEdges", tet.getSubcellCount(subcellDim), tet.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.tetEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, tet );
      deep_copy(subcellParamData_.tetEdges,subcellParamHost);
    }
    {
      const auto hex = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());

      subcellDim = 2;
      subcellParamData_.hexFaces = subcellParamViewType("CellTools::SubcellParametrization::hexFaces", hex.getSubcellCount(subcellDim), hex.getDimension(), subcellDim+1);
      auto subcell2dParamHost = Kokkos::create_mirror_view(subcellParamData_.hexFaces);
      setSubcellParametrization( subcell2dParamHost, subcellDim, hex );
      deep_copy(subcellParamData_.hexFaces,subcell2dParamHost);

      subcellDim = 1;
      subcellParamData_.hexEdges = subcellParamViewType("CellTools::SubcellParametrization::hexEdges", hex.getSubcellCount(subcellDim), hex.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.hexEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, hex );
      deep_copy(subcellParamData_.hexEdges,subcellParamHost);
    }
    {
      const auto pyr = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >());

      subcellDim = 2;
      subcellParamData_.pyrFaces = subcellParamViewType("CellTools::SubcellParametrization::pyrFaces", pyr.getSubcellCount(subcellDim), pyr.getDimension(), subcellDim+1);
      auto subcell2dParamHost = Kokkos::create_mirror_view(subcellParamData_.pyrFaces);
      setSubcellParametrization( subcell2dParamHost, subcellDim, pyr );
      deep_copy(subcellParamData_.pyrFaces,subcell2dParamHost);

      subcellDim = 1;
      subcellParamData_.pyrEdges = subcellParamViewType("CellTools::SubcellParametrization::pyrEdges", pyr.getSubcellCount(subcellDim), pyr.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.pyrEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, pyr );
      deep_copy(subcellParamData_.pyrEdges,subcellParamHost);
    }
    {
      const auto wedge = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >());

      subcellDim = 2;
      subcellParamData_.wedgeFaces = subcellParamViewType("CellTools::SubcellParametrization::wedgeFaces", wedge.getSubcellCount(subcellDim), wedge.getDimension(), subcellDim+1);
      auto subcell2dParamHost = Kokkos::create_mirror_view(subcellParamData_.wedgeFaces);
      setSubcellParametrization( subcell2dParamHost, subcellDim, wedge );
      deep_copy(subcellParamData_.wedgeFaces,subcell2dParamHost);

      subcellDim = 1;
      subcellParamData_.wedgeEdges = subcellParamViewType("CellTools::SubcellParametrization::wedgeEdges", wedge.getSubcellCount(subcellDim), wedge.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.wedgeEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, wedge );
      deep_copy(subcellParamData_.wedgeEdges,subcellParamHost);
    }
    {
      const auto tri = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >());

      subcellDim = 1;
      subcellParamData_.triEdges = subcellParamViewType("CellTools::SubcellParametrization::triEdges", tri.getSubcellCount(subcellDim), tri.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.triEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, tri );
      deep_copy(subcellParamData_.triEdges,subcellParamHost);
    }
    {
      const auto quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());

      subcellDim = 1;
      subcellParamData_.quadEdges = subcellParamViewType("CellTools::SubcellParametrization::quadEdges", quad.getSubcellCount(subcellDim), quad.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.quadEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, quad );
      deep_copy(subcellParamData_.quadEdges,subcellParamHost);

    }
    {
      const auto line = shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<2> >());

      subcellDim = 1;
      subcellParamData_.lineEdges = subcellParamViewType("CellTools::SubcellParametrization::lineEdges", line.getSubcellCount(subcellDim), line.getDimension(), subcellDim+1);
      auto subcellParamHost = Kokkos::create_mirror_view(subcellParamData_.lineEdges);
      setSubcellParametrization( subcellParamHost, subcellDim, line );
      deep_copy(subcellParamData_.lineEdges,subcellParamHost);

    }

    Kokkos::push_finalize_hook( [=] {
      subcellParamData_.dummy = subcellParamViewType();
      subcellParamData_.lineEdges = subcellParamViewType();
      subcellParamData_.triEdges = subcellParamViewType();
      subcellParamData_.quadEdges = subcellParamViewType();
      subcellParamData_.shellTriEdges = subcellParamViewType();
      subcellParamData_.shellQuadEdges = subcellParamViewType();
      subcellParamData_.tetEdges = subcellParamViewType();
      subcellParamData_.hexEdges = subcellParamViewType();
      subcellParamData_.pyrEdges = subcellParamViewType();
      subcellParamData_.wedgeEdges = subcellParamViewType();
      subcellParamData_.shellTriFaces = subcellParamViewType();
      subcellParamData_.shellQuadFaces = subcellParamViewType();
      subcellParamData_.tetFaces = subcellParamViewType();
      subcellParamData_.hexFaces = subcellParamViewType();
      subcellParamData_.pyrFaces = subcellParamViewType();
      subcellParamData_.wedgeFaces = subcellParamViewType();
    });

    isSubcellParametrizationSet_= true;
  }
  

  template<typename DeviceType>
  void
  CellTools<DeviceType>::
  getSubcellParametrization(       subcellParamViewType &subcellParam,
                             const ordinal_type         subcellDim,
                             const shards::CellTopology parentCell ) {
    Impl::CellTools<DeviceType>::setSubcellParametrization();
    subcellParam = Impl::CellTools<DeviceType>::getSubcellParametrization(subcellDim, parentCell.getKey());
  }

  template<typename DeviceType>
  inline
  typename Impl::CellTools<DeviceType>::subcellParamViewConstType
  Impl::CellTools<DeviceType>::
  getSubcellParametrization( const ordinal_type          subcellDim,
                             const unsigned              parentCellKey ) {

    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( !isSubcellParametrizationSet_, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::getSubcellParametrization): subcell parametrization is not set.");

    subcellParamViewType subcellParam;

    switch (parentCellKey ) {
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<8>::key:
    case shards::Tetrahedron<10>::key:
    case shards::Tetrahedron<11>::key:       subcellParam = ( subcellDim == 2 ? subcellParamData_.tetFaces : subcellParamData_.tetEdges ); break;

    case shards::Hexahedron<8>::key:
    case shards::Hexahedron<20>::key:
    case shards::Hexahedron<27>::key:        subcellParam = ( subcellDim == 2 ? subcellParamData_.hexFaces : subcellParamData_.hexEdges ); break;

    case shards::Pyramid<5>::key:
    case shards::Pyramid<13>::key:
    case shards::Pyramid<14>::key:           subcellParam = ( subcellDim == 2 ? subcellParamData_.pyrFaces : subcellParamData_.pyrEdges ); break;

    case shards::Wedge<6>::key:
    case shards::Wedge<15>::key:
    case shards::Wedge<18>::key:             subcellParam = ( subcellDim == 2 ? subcellParamData_.wedgeFaces : subcellParamData_.wedgeEdges ); break;

    case shards::Triangle<3>::key:
    case shards::Triangle<4>::key:
    case shards::Triangle<6>::key:           subcellParam = subcellParamData_.triEdges; break;

    case shards::Quadrilateral<4>::key:
    case shards::Quadrilateral<8>::key:
    case shards::Quadrilateral<9>::key:      subcellParam = subcellParamData_.quadEdges; break;

    // case shards::ShellTriangle<3>::key:
    // case shards::ShellTriangle<6>::key:      subcellParam = ( subcellDim == 2 ? subcellParamData_.shellTriFaces : subcellParamData_.shellTriEdges ); break;

    // case shards::ShellQuadrilateral<4>::key:
    // case shards::ShellQuadrilateral<8>::key:
    // case shards::ShellQuadrilateral<9>::key: subcellParam = ( subcellDim == 2 ? subcellParamData_.shellQuadFaces : subcellParamData_.shellQuadEdges ); break;

    case shards::ShellLine<2>::key:
    case shards::ShellLine<3>::key:
    case shards::Beam<2>::key:
    case shards::Beam<3>::key:               subcellParam = subcellParamData_.lineEdges; break;
    default: {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( true, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::getSubcellParametrization): invalid cell topology.");
    }
    }
    return subcellParam;
  }
  
  template<typename DeviceType>
  template <typename HostViewType>
  void
  Impl::CellTools<DeviceType>::
  setSubcellParametrization( HostViewType               subcellParam,
                             const ordinal_type         subcellDim,
                             const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell.getKey()), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setSubcellParametrization): the specified cell topology does not have a reference cell.");
#endif
    // subcellParametrization is rank-3 FieldContainer with dimensions (SC, PCD, COEF) where:
    //  - SC    is the subcell count of subcells with the specified dimension in the parent cell
    //  - PCD   is Parent Cell Dimension, which gives the number of coordinate functions in the map
    //          PCD = 2 for standard 2D cells and non-standard 2D cells: shell line and beam
    //          PCD = 3 for standard 3D cells and non-standard 3D cells: shell Tri and Quad
    //  - COEF  is number of coefficients needed to specify a coordinate function:
    //          COEFF = 2 for edge parametrizations
    //          COEFF = 3 for both Quad and Tri face parametrizations. Because all Quad reference faces
    //          are affine, the coefficient of the bilinear term u*v is zero and is not stored, i.e.,
    //          3 coefficients are sufficient to store Quad face parameterization maps.
    //  
    // Edge parametrization maps [-1,1] to edge defined by (v0, v1)
    // Face parametrization maps [-1,1]^2 to quadrilateral face (v0, v1, v2, v3), or
    // standard 2-simplex  {(0,0),(1,0),(0,1)} to traingle face (v0, v1, v2).
    // This defines orientation-preserving parametrizations with respect to reference edge and
    // face orientations induced by their vertex order. 

    // get subcellParametrization dimensions: (sc, pcd, coeff)
    const auto sc    = parentCell.getSubcellCount(subcellDim);
    const auto pcd   = parentCell.getDimension();

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim < 1 || subcellDim > static_cast<ordinal_type>(pcd-1), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::setSubcellParametrization): Parametrizations defined in a range between 1 and (dim-1)");

    Impl::CellTools<HostSpaceType>::setReferenceNodeData();
    const auto refNodes = Impl::CellTools<HostSpaceType>::getReferenceNodes(parentCell.getKey());

    if (subcellDim == 1) {
      // Edge parametrizations of 2D and 3D cells (shell lines and beams are 2D cells with edges)
      for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {
        // vertexK[0] = x_k; vertexK[1] = y_k; vertexK[2] = z_k; z_k = 0 for 2D cells
        // Note that ShellLine and Beam are 2D cells!
        const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
        const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);

        const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
        const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());
        
        // x(t) = (x0 + x1)/2 + t*(x1 - x0)/2 
        subcellParam(subcellOrd, 0, 0) = (v0(0) + v1(0))/2.0;
        subcellParam(subcellOrd, 0, 1) = (v1(0) - v0(0))/2.0;
        
        // y(t) = (y0 + y1)/2 + t*(y1 - y0)/2 
        subcellParam(subcellOrd, 1, 0) = (v0(1) + v1(1))/2.0;
        subcellParam(subcellOrd, 1, 1) = (v1(1) - v0(1))/2.0;
        
        if( pcd == 3 ) {
          // z(t) = (z0 + z1)/2 + t*(z1 - z0)/2 
          subcellParam(subcellOrd, 2, 0) = (v0(2) + v1(2))/2.0;
          subcellParam(subcellOrd, 2, 1) = (v1(2) - v0(2))/2.0;
        }
      }
    }
    else if (subcellDim == 2) {
      // Face parametrizations of 3D cells: (shell Tri and Quad are 3D cells with faces)
      // A 3D cell can have both Tri and Quad faces, but because they are affine images of the
      // parametrization domain, 3 coefficients are enough to store them in both cases.
      for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {
        
        switch (parentCell.getKey(subcellDim,subcellOrd)) {
          
        case shards::Triangle<3>::key:
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key: {
          const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
          const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
          const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);

          //getReferenceVertex(v0, parentCell, v0ord);
          //getReferenceVertex(v1, parentCell, v1ord);
          //getReferenceVertex(v2, parentCell, v2ord);
          const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
          const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());
          const auto v2 = Kokkos::subview(refNodes, v2ord, Kokkos::ALL());
          
          // x(u,v) = x0 + (x1 - x0)*u + (x2 - x0)*v
          subcellParam(subcellOrd, 0, 0) = v0(0);
          subcellParam(subcellOrd, 0, 1) = v1(0) - v0(0);
          subcellParam(subcellOrd, 0, 2) = v2(0) - v0(0);
          
          // y(u,v) = y0 + (y1 - y0)*u + (y2 - y0)*v
          subcellParam(subcellOrd, 1, 0) = v0(1);
          subcellParam(subcellOrd, 1, 1) = v1(1) - v0(1);
          subcellParam(subcellOrd, 1, 2) = v2(1) - v0(1);
          
          // z(u,v) = z0 + (z1 - z0)*u + (z2 - z0)*v
          subcellParam(subcellOrd, 2, 0) = v0(2);
          subcellParam(subcellOrd, 2, 1) = v1(2) - v0(2);
          subcellParam(subcellOrd, 2, 2) = v2(2) - v0(2);
          break;
        }
        case shards::Quadrilateral<4>::key:
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key: {
          const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
          const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
          const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
          const auto v3ord = parentCell.getNodeMap(subcellDim, subcellOrd, 3);

          const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
          const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());
          const auto v2 = Kokkos::subview(refNodes, v2ord, Kokkos::ALL());
          const auto v3 = Kokkos::subview(refNodes, v3ord, Kokkos::ALL());
                
          // x(u,v) = (x0+x1+x2+x3)/4+u*(-x0+x1+x2-x3)/4+v*(-x0-x1+x2+x3)/4+uv*(0=x0-x1+x2-x3)/4 
          subcellParam(subcellOrd, 0, 0) = ( v0(0) + v1(0) + v2(0) + v3(0))/4.0;
          subcellParam(subcellOrd, 0, 1) = (-v0(0) + v1(0) + v2(0) - v3(0))/4.0;
          subcellParam(subcellOrd, 0, 2) = (-v0(0) - v1(0) + v2(0) + v3(0))/4.0;
          
          // y(u,v) = (y0+y1+y2+y3)/4+u*(-y0+y1+y2-y3)/4+v*(-y0-y1+y2+y3)/4+uv*(0=y0-y1+y2-y3)/4 
          subcellParam(subcellOrd, 1, 0) = ( v0(1) + v1(1) + v2(1) + v3(1))/4.0;
          subcellParam(subcellOrd, 1, 1) = (-v0(1) + v1(1) + v2(1) - v3(1))/4.0;
          subcellParam(subcellOrd, 1, 2) = (-v0(1) - v1(1) + v2(1) + v3(1))/4.0;
          
          // z(u,v) = (z0+z1+z2+z3)/4+u*(-z0+z1+z2-z3)/4+v*(-z0-z1+z2+z3)/4+uv*(0=z0-z1+z2-z3)/4 
          subcellParam(subcellOrd, 2, 0) = ( v0(2) + v1(2) + v2(2) + v3(2))/4.0;
          subcellParam(subcellOrd, 2, 1) = (-v0(2) + v1(2) + v2(2) - v3(2))/4.0;
          subcellParam(subcellOrd, 2, 2) = (-v0(2) - v1(2) + v2(2) + v3(2))/4.0;
          break;
        }
        default: {
          INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                        ">>> ERROR (Intrepid2::CellTools::setSubcellParametrization): parametrization not defined for the specified face topology.");
        }
        }
      }
    }
  }

  template<typename DeviceType>
  bool
  Impl::CellTools<DeviceType>::
  isSubcellParametrizationSet_ = false;

  template<typename DeviceType>
  typename Impl::CellTools<DeviceType>::SubcellParamData
  Impl::CellTools<DeviceType>::
  subcellParamData_ = typename Impl::CellTools<DeviceType>::SubcellParamData();
 

}

#endif
