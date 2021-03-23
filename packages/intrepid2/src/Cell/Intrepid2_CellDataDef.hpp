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

/** \file   Intrepid2_CellData.hpp
    \brief  Definition file for the Intrepid2::Impl::CellData class.
    \author Kyungjoo Kim
    \author Mauro Perego
*/

#ifndef __INTREPID2_CELLDATA_DEF_HPP__
#define __INTREPID2_CELLDATA_DEF_HPP__

namespace Intrepid2 {

namespace Impl {

template<typename DeviceType>
inline
typename CellData<DeviceType>::subcellParamViewConstType
CellData<DeviceType>::
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
CellData<DeviceType>::
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

  CellData<HostSpaceType>::setReferenceNodeData();
  const auto refNodes = CellData<HostSpaceType>::getReferenceNodes(parentCell.getKey());

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
inline
typename CellData<DeviceType>::referenceNodeDataViewConstType
CellData<DeviceType>::
getReferenceNodes(const unsigned      cellTopoKey){
#ifdef HAVE_INTREPID2_DEBUG
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( !hasReferenceCell(cellTopoKey), std::invalid_argument,
      ">>> ERROR (Intrepid2::CellTools::getReferenceNodes): the specified cell topology does not have a reference cell." );
#endif

  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( !isReferenceNodeDataSet_, std::invalid_argument,
                                            ">>> ERROR (Intrepid2::CellTools::getReferenceNodes): referenceNodeData not set. \nCall setReferenceNodeData() first.");

  referenceNodeDataViewType refNodes;

  switch (cellTopoKey ) {
  case shards::Line<2>::key:
  case shards::ShellLine<2>::key:
  case shards::Beam<2>::key:               refNodes = refNodeData_.line; break;
  case shards::Line<3>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<3>::key:               refNodes = refNodeData_.line_3; break;

  case shards::Triangle<3>::key:
  case shards::ShellTriangle<3>::key:      refNodes = refNodeData_.triangle; break;
  case shards::Triangle<4>::key:           refNodes = refNodeData_.triangle_4; break;
  case shards::Triangle<6>::key:
  case shards::ShellTriangle<6>::key:      refNodes = refNodeData_.triangle_6; break;

  case shards::Quadrilateral<4>::key:
  case shards::ShellQuadrilateral<4>::key: refNodes = refNodeData_.quadrilateral; break;
  case shards::Quadrilateral<8>::key:
  case shards::ShellQuadrilateral<8>::key: refNodes = refNodeData_.quadrilateral_8; break;
  case shards::Quadrilateral<9>::key:
  case shards::ShellQuadrilateral<9>::key: refNodes = refNodeData_.quadrilateral_9; break;

  case shards::Tetrahedron<4>::key:        refNodes = refNodeData_.tetrahedron; break;
  case shards::Tetrahedron<8>::key:        refNodes = refNodeData_.tetrahedron_8; break;
  case shards::Tetrahedron<10>::key:       refNodes = refNodeData_.tetrahedron_10; break;
  case shards::Tetrahedron<11>::key:       refNodes = refNodeData_.tetrahedron_11; break;

  case shards::Hexahedron<8>::key:         refNodes = refNodeData_.hexahedron; break;
  case shards::Hexahedron<20>::key:        refNodes = refNodeData_.hexahedron_20; break;
  case shards::Hexahedron<27>::key:        refNodes = refNodeData_.hexahedron_27; break;

  case shards::Pyramid<5>::key:            refNodes = refNodeData_.pyramid; break;
  case shards::Pyramid<13>::key:           refNodes = refNodeData_.pyramid_13; break;
  case shards::Pyramid<14>::key:           refNodes = refNodeData_.pyramid_14; break;

  case shards::Wedge<6>::key:              refNodes = refNodeData_.wedge; break;
  case shards::Wedge<15>::key:             refNodes = refNodeData_.wedge_15; break;
  case shards::Wedge<18>::key:             refNodes = refNodeData_.wedge_18; break;

  default: {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( true, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid cell topology.");
  }
  }
  return refNodes;
}

template<typename DeviceType>
inline
typename CellData<DeviceType>::referenceNodeDataViewConstType
CellData<DeviceType>::
getReferenceCellCenter(const unsigned      cellTopoKey){
#ifdef HAVE_INTREPID2_DEBUG
  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( !hasReferenceCell(cellTopoKey), std::invalid_argument,
      ">>> ERROR (Intrepid2::CellTools::getReferenceNodes): the specified cell topology does not have a reference cell." );
#endif

  INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( !isReferenceCellCenterDataSet_, std::invalid_argument,
                                            ">>> ERROR (Intrepid2::CellTools::getReferenceNodes): referenceNodeData not set. \nCall setReferenceNodeData() first.");

  referenceNodeDataViewType cellCenter;

  switch (cellTopoKey ) {
  case shards::Line<2>::key:
  case shards::ShellLine<2>::key:
  case shards::Beam<2>::key:
  case shards::Line<3>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<3>::key:               cellCenter = refCenterData_.line; break;

  case shards::Triangle<3>::key:
  case shards::ShellTriangle<3>::key:
  case shards::Triangle<4>::key:
  case shards::Triangle<6>::key:
  case shards::ShellTriangle<6>::key:      cellCenter = refCenterData_.triangle; break;

  case shards::Quadrilateral<4>::key:
  case shards::ShellQuadrilateral<4>::key:
  case shards::Quadrilateral<8>::key:
  case shards::ShellQuadrilateral<8>::key:
  case shards::Quadrilateral<9>::key:
  case shards::ShellQuadrilateral<9>::key: cellCenter = refCenterData_.quadrilateral; break;

  case shards::Tetrahedron<4>::key:
  case shards::Tetrahedron<8>::key:
  case shards::Tetrahedron<10>::key:
  case shards::Tetrahedron<11>::key:       cellCenter = refCenterData_.tetrahedron; break;

  case shards::Hexahedron<8>::key:
  case shards::Hexahedron<20>::key:
  case shards::Hexahedron<27>::key:        cellCenter = refCenterData_.hexahedron; break;

  case shards::Pyramid<5>::key:
  case shards::Pyramid<13>::key:
  case shards::Pyramid<14>::key:           cellCenter = refCenterData_.pyramid; break;

  case shards::Wedge<6>::key:
  case shards::Wedge<15>::key:
  case shards::Wedge<18>::key:             cellCenter = refCenterData_.wedge; break;

  default: {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( true, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): invalid cell topology.");
  }
  }
  return cellCenter;
}


template<typename DeviceType>
bool
CellData<DeviceType>::
isSubcellParametrizationSet_ = false;

template<typename DeviceType>
typename CellData<DeviceType>::SubcellParamData
CellData<DeviceType>::
subcellParamData_ = typename CellData<DeviceType>::SubcellParamData();


template<typename DeviceType>
bool
CellData<DeviceType>::
isReferenceNodeDataSet_ = false;

template<typename DeviceType>
typename CellData<DeviceType>::ReferenceNodeData
CellData<DeviceType>::
refNodeData_ = typename CellData<DeviceType>::ReferenceNodeData();

template<typename DeviceType>
const typename CellData<DeviceType>::ReferenceNodeDataStatic
CellData<DeviceType>::
refNodeDataStatic_ = {
  // line
  { // 2
    {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}
  },
  { // 3
    {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
  },
  // triangle
  { // 3
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}
  },
  { // 4
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 1.0/3.0, 1.0/3.0, 0.0}
  },
  { // 6
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
    { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
  },
  // quad
  { // 4
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}
  },
  { // 8
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
    { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}
  },
  { // 9
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
    { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
  },
  // tet
  { // 4
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
  },
  { // 8
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
    { 1.0/3.0, 0.0, 1.0/3.0}, { 1.0/3.0, 1.0/3.0, 1.0/3.0}, { 1.0/3.0, 1.0/3.0, 0.0}, { 0.0, 1.0/3.0, 1.0/3.0}
  },
  { // 10
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
    { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
  },
  { // 11
    { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
    { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
  },
  // hex
  { // 8
    {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
    {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
  },
  { // 20
    {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
    {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
    { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
    { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0}
  },
  { // 27
    {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
    {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
    { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
    { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
    { 0.0, 0.0, 0.0},
    { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0}
  },
  // pyramid
  { // 5
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
  },
  { // 13
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
    { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
    {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}
  },
  { // 14
    {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
    { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
    {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0}
  },
  // wedge
  { // 6
    { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}
  },
  { // 15
    { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
    { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
    { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0}
  },
  { // 18
    { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
    { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
    { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
    { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
  }
};

template<typename DeviceType>
bool
CellData<DeviceType>::
isReferenceCellCenterDataSet_ = false;

template<typename DeviceType>
typename CellData<DeviceType>::ReferenceCellCenterData
CellData<DeviceType>::
refCenterData_ = typename CellData<DeviceType>::ReferenceCellCenterData();

template<typename DeviceType>
const typename CellData<DeviceType>::ReferenceCenterDataStatic
CellData<DeviceType>::
refCenterDataStatic_ = {
  // line
  {0.0, 0.0, 0.0},
  // triangle
  { 1.0/3.0, 1.0/3.0, 0.0},
  // quad
  {0.0, 0.0, 0.0},
  // tet
  { 0.25, 0.25, 0.25},
  // hex
  { 0.0, 0.0, 0.0},
  // pyramid
  { 0.0, 0.0, 0.25},
  // wedge
  { 1.0/3.0, 1.0/3.0, 0.0},
};

}
}

#endif

