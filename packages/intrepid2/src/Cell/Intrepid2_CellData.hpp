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

#ifndef __INTREPID2_CELLDATA_HPP__
#define __INTREPID2_CELLDATA_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Kernels.hpp"

namespace Intrepid2 {

namespace Impl {



//============================================================================================//
//                                                                                            //
//          Parametrization coefficients of edges and faces of reference cells                //
//                                                                                            //
//============================================================================================//

template<typename DeviceType>
class RefCellParametrization {
  using ExecSpaceType = typename DeviceType::execution_space;
  using HostSpaceType = typename Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space;
public:
  typedef Kokkos::DynRankView<double,DeviceType> subcellParamViewType;
  typedef Kokkos::DynRankView<const double,DeviceType> subcellParamViewConstType;


  KOKKOS_INLINE_FUNCTION
  static bool
  hasReferenceParametrization( const unsigned cellTopoKey );

  /** \brief  Default constructor.
   */
  RefCellParametrization() = default;

  /** \brief  Destructor
   */
  ~RefCellParametrization() = default;

  /** \brief  Defines orientation-preserving parametrizations of reference edges and faces of cell
      topologies with reference cells.

      Given an edge {V0, V1} of some reference cell, its parametrization is a mapping from
      [-1,1] onto the edge. Parametrization of a triangular face {V0,V1,V2} is mapping from
      the standard 2-simplex {(0,0,0), (1,0,0), (0,1,0)}, embedded in 3D onto that face.
      Parametrization of a quadrilateral face {V0,V1,V2,V3} is mapping from the standard
      2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)}, embedded in 3D, onto that face.

      This method computes the coefficients of edge and face parametrization maps.
      All mappings are affine and orientation-preserving, i.e., they preserve the tangent
      and normal directions implied by the vertex order of the edge or the face relative to
      the reference cell:

      \li     the tangent on [-1,1] from -1 in the direction of 1 is mapped to a tangent on edge {V0,V1}
      from V0 in the direction of V1  (the forward direction of the edge determined by its
      start and end vertices)

      \li     the normal in the direction of (0,0,1) to the standard 2-simplex {(0,0,0),(1,0,0),(0,1,0)}
      and the standard 2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)} is mapped to a normal
      on {V0,V1,V2} and {V0,V1,V2,V3}, determined according to the right-hand rule
      (see http://mathworld.wolfram.com/Right-HandRule.html for definition of right-hand rule
      and Section \ref Section sec_cell_topology_subcell_map for further details).

      Because faces of all reference cells supported in Intrepid are affine images of either
      the standard 2-simplex or the standard 2-cube, the coordinate functions of the respective
      parmetrization maps are linear polynomials in the parameter variables (u,v), i.e., they
      are of the form \c F_i(u,v)=C_0(i)+C_1(i)u+C_2(i)v;  \c 0<=i<3 (face parametrizations
      are supported only for 3D cells, thus parametrization maps have 3 coordinate functions).
      As a result, application of these maps is independent of the face type which is convenient
      for cells such as Wedge or Pyramid that have both types of faces. Also, coefficients of
      coordinate functions for all faces can be stored together in the same array.

   */
  static void setSubcellParametrization();

  static inline
  subcellParamViewConstType
  getSubcellParametrization( const ordinal_type          subcellDim,
      const unsigned              parentCellKey );

private:
  /** \brief  Sets orientation-preserving parametrizations of reference edges and faces of cell
      topologies with reference cells. Used to populate Intrepid2::Impl::CellData::SubcellParamData.

      See Intrepid2::Impl::CellData::setSubcellParametrization and Section \ref sec_cell_topology_subcell_map
      more information about parametrization maps.

      \param  subcellParam           [out]  - array with the coefficients of the parametrization map
      \param  subcellDim             [in]   - dimension of the subcells being parametrized (1 or 2)
      \param  parentCell             [in]   - topology of the parent cell owning the subcells.
   */
  template <typename HostViewType>
  static void
  setSubcellParametrization(       HostViewType          subcellParam,
      const ordinal_type          subcellDim,
      const shards::CellTopology  parentCell );

  // static variables
  /** \struct Intrepid2::Impl::CellData::SubcellParamData
      \brief Parametrization coefficients of edges and faces of reference cells
   */

  static subcellParamViewType dummyParam;
  static subcellParamViewType lineEdgesParam;  // edge maps for 2d non-standard cells; shell line and beam
  static subcellParamViewType triEdgesParam, quadEdgesParam; // edge maps for 2d standard cells
  static subcellParamViewType shellTriEdgesParam, shellQuadEdgesParam; // edge maps for 3d non-standard cells; shell tri and quad
  static subcellParamViewType tetEdgesParam, hexEdgesParam, pyrEdgesParam, wedgeEdgesParam; // edge maps for 3d standard cells
  static subcellParamViewType shellTriFacesParam, shellQuadFacesParam; // face maps for 3d non-standard cells
  static subcellParamViewType tetFacesParam, hexFacesParam, pyrFacesParam, wedgeFacesParam; // face maps for 3d standard cells

  static bool isSubcellParametrizationSet_;


};

template<typename DeviceType>
class RefCellNodes {
  using ExecSpaceType = typename DeviceType::execution_space;
  using HostSpaceType = typename Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space;
public:
  typedef Kokkos::DynRankView<double,DeviceType> referenceNodeDataViewType;
  typedef Kokkos::DynRankView<const double,DeviceType> referenceNodeDataViewConstType;
  typedef Kokkos::DynRankView<const double,HostSpaceType> referenceNodeDataConstViewHostType;


  /** \brief  Default constructor.
   */
  RefCellNodes() = default;

  /** \brief  Destructor
   */
  ~RefCellNodes() = default;

  /** \brief Set reference node coordinates for supported topologies.
   */
  static void setReferenceNodeData();

  /** \brief  Retrieves the Cartesian coordinates of a reference cell node.

    Returns Cartesian coordinates of a reference cell node. Requires cell topology
    with a reference cell. Node coordinates are always returned as an (x,y,z)-triple
    regardlesss of the actual topological cell dimension. The unused coordinates are
    set to zero, e.g., node 0 of Line<2> is returned as {-1,0,0}.

    \remark
    Because the nodes of a cell with a base topology coincide with its vertices, for cells
    with base topology this method is equivalent to Intrepid2::Impl::CellData::getReferenceNodes.

    \param  cell              [in]  - key of the cell topology
   */
  static inline
  referenceNodeDataViewConstType
  getReferenceNodes(const unsigned      cellTopoKey);

private:

  /** \struct Intrepid2::Impl::CellData::ReferenceNodeData
    \brief Reference node data for each supported topology
   */
    static referenceNodeDataViewType lineNodes, line3Nodes;
    static referenceNodeDataViewType triangleNodes, triangle4Nodes, triangle6Nodes;
    static referenceNodeDataViewType quadrilateralNodes, quadrilateral8Nodes, quadrilateral9Nodes;
    static referenceNodeDataViewType tetrahedronNodes, tetrahedron8Nodes, tetrahedron10Nodes, tetrahedron11Nodes;
    static referenceNodeDataViewType hexahedronNodes, hexahedron20Nodes, hexahedron27Nodes;
    static referenceNodeDataViewType pyramidNodes, pyramid13Nodes, pyramid14Nodes;
    static referenceNodeDataViewType wedgeNodes, wedge15Nodes, wedge18Nodes;

  // reference nodes initialized
  /** \struct Intrepid2::CellTools::ReferenceNodeDataStatic
    \brief Reference node containers for each supported topology
   */
  struct ReferenceNodeDataStatic {
    double line[2][3], line_3[3][3];
    double triangle[3][3], triangle_4[4][3], triangle_6[6][3];
    double quadrilateral[4][3], quadrilateral_8[8][3], quadrilateral_9[9][3];
    double tetrahedron[4][3], tetrahedron_8[8][3], tetrahedron_10[10][3], tetrahedron_11[10][3];
    double hexahedron[8][3], hexahedron_20[20][3], hexahedron_27[27][3];
    double pyramid[5][3], pyramid_13[13][3], pyramid_14[14][3];
    double wedge[6][3], wedge_15[15][3], wedge_18[18][3];
  };

  static const ReferenceNodeDataStatic refNodeDataStatic_;
  static bool isReferenceNodeDataSet_;

};


template<typename DeviceType>
class RefCellCenter {
  using ExecSpaceType = typename DeviceType::execution_space;
  using HostSpaceType = typename Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space;
public:
  typedef Kokkos::DynRankView<double,DeviceType> referenceNodeDataViewType;
  typedef Kokkos::DynRankView<const double,DeviceType> referenceNodeDataViewConstType;
  typedef Kokkos::DynRankView<const double,HostSpaceType> referenceNodeDataConstViewHostType;


  /** \brief  Default constructor.
   */
  RefCellCenter() = default;

  /** \brief  Destructor
   */
  ~RefCellCenter() = default;


  /** \brief Set center coordinates of reference cell for supported topologies.
   */
  static void setReferenceCellCenterData();

  /** \brief  Retrieves the Cartesian coordinates of a reference cell node.

    Returns Cartesian coordinates of a reference cell node. Requires cell topology
    with a reference cell. Node coordinates are always returned as an (x,y,z)-triple
    regardlesss of the actual topological cell dimension. The unused coordinates are
    set to zero, e.g., node 0 of Line<2> is returned as {-1,0,0}.

    \remark
    Because the nodes of a cell with a base topology coincide with its vertices, for cells
    with base topology this method is equivalent to Intrepid2::Impl::CellData::getReferenceCellCenter.

    \param  cell              [in]  - key of the cell topology
   */
  static inline
  referenceNodeDataViewConstType
  getReferenceCellCenter(const unsigned      cellTopoKey);

private:

  static referenceNodeDataViewType lineCenter;
  static referenceNodeDataViewType triangleCenter;
  static referenceNodeDataViewType quadrilateralCenter;
  static referenceNodeDataViewType tetrahedronCenter;
  static referenceNodeDataViewType hexahedronCenter;
  static referenceNodeDataViewType pyramidCenter;
  static referenceNodeDataViewType wedgeCenter;

  struct ReferenceCenterDataStatic {
    double line[3];
    double triangle[3];
    double quadrilateral[3];
    double tetrahedron[3];
    double hexahedron[3];
    double pyramid[3];
    double wedge[3];
  };

  static  const ReferenceCenterDataStatic refCenterDataStatic_;
  static bool isReferenceCellCenterDataSet_;
};
}
}

#include "Intrepid2_CellDataDef.hpp"

#endif

