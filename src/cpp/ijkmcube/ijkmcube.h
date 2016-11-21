/// \file ijkmcube.h
/// Generate isosurface in arbitrary dimensions

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006,2007,2009 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*!
  \mainpage IJKMCUBE: IJK MARCHING CUBES

  IJKMCUBE is a program for generating isosurfaces and interval volumes
  using the Marching Cubes algorithm or one of its variants.  It computes
  isosurfaces from two, three or four dimensional volumetric grid data.
  Options include using regular subdivisions or an octree
  for faster isosurface extraction, using negative-equals-positive 
  isosurface lookup tables, using multi-resolution meshes, snapping
  isosurface vertices to grid vertices to eliminate small angled
  isosurface triangles and using the asymptotic decider to determine
  topology of isosurface patches within cubes.
*/

#ifndef _IJKMCUBE_
#define _IJKMCUBE_

#include <string>

#include "ijk.txx"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"
#include "ijkmcube_sub.h"

/// ijkmcube classes and routines.
namespace IJKMCUBE {

// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

  /// Marching Cubes Algorithm.
  void marching_cubes
    (const MC_DATA & mc_data, const SCALAR_TYPE isovalue, 
     MC_ISOSURFACE & mc_isosurface, MCUBE_INFO & mcube_info, const MC_DATA & mc_data_1,
     MC_ISOSURFACE & mc_isosurface_1, MCUBE_INFO & mcube_info_1, float* &f_color,  float* numv);

  /// Interval Volume Algorithm.
  void MCVol
    (const MC_DATA & mc_data, 
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     MC_ISOSURFACE & mc_ivol, MCUBE_INFO & mcube_info);

  /// Multiresolution isosurface extraction using cubes, pyramids and tetrahedra.
  void multires_cubes
    (const MC_DATA & mc_data, const SCALAR_TYPE isovale,
     MC_ISOSURFACE & mc_isosurface, MCUBE_INFO & mcube_info);
     
// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

  /// Marching Cubes Algorithm.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord);

  /// Marching Cubes Algorithm.  Faster, but memory intensive version.
  /// Represents each grid edge by a single integer.
  /// @param merge_data = Data structure for merging edges.  
  /// Requires memory of size(MERGE_INDEX) for each grid edge.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm.  Faster, but memory intensive version.
  /// Represents each grid edge by a single integer.
  /// @param cube_list[i] = Primary vertex of cube containing i'th simplex.
  /// @param merge_data = Data structure for merging edges.  
  /// Requires memory of size(MERGE_INDEX) for each grid edge.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    std::vector<VERTEX_INDEX> & cube_list,
    MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm.  Slower, but less memory intensive version.
  /// Represents edges containing isosurface vertices as pairs of grid vertices.
  /// @param merge_edges_parameter = Parameters for merging edges.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    const MERGE_EDGES_PARAMETERS & merge_edges_parameter, 
    MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm with specified isosurface topology.
  /// @param isosurface_topology = Isosurface topology.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const POLY_ISOTABLE & poly_isotable,
    const SCALAR_TYPE isovalue,
    const ISOSURFACE_TOPOLOGY isosurface_topology,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);

  /// Marching Cubes Algorithm with specified isosurface topology.
  /// @param isosurface_topology = Isosurface topology.
  /// @param interpolation_type = Interpolation type.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const POLY_ISOTABLE & poly_isotable,
    const SCALAR_TYPE isovalue,
    const ISOSURFACE_TOPOLOGY isosurface_topology,
    const INTERPOLATION_TYPE interpolation_type,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);


  // newly added
  /// Marching Cubes Algorithm with specified isosurface topology.
  /// @param isosurface_topology = Isosurface topology.
  /// @param interpolation_type = Interpolation type.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const POLY_ISOTABLE & poly_isotable,
    const SCALAR_TYPE isovalue,
    const ISOSURFACE_TOPOLOGY isosurface_topology,
    const INTERPOLATION_TYPE interpolation_type,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &,const MC_SCALAR_GRID_BASE & scalar_grid_1, 
    const ISOSURFACE_TOPOLOGY isosurface_topology_1,
    std::vector<VERTEX_INDEX> & simplex_vert_1, 
    std::vector<COORD_TYPE> & vertex_coord_1,
    const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &, float* &f_color, float* numvert);	






  /// Marching Cubes Algorithm with specified isosurface topology.
  /// @param isosurface_topology = Isosurface topology.
  /// @pre Topologies cannot require a pyramid isotable (i.e., LINEAR_TOPOLOGY for dimension > 2 is not a valid.)
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & cube_isotable,
    const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info,
    const SCALAR_TYPE isovalue,
    const ISOSURFACE_TOPOLOGY isosurface_topology,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA &, MCUBE_INFO &);

  /// Marching Cubes Algorithm with specified isosurface topology.
  /// @param isosurface_topology = Isosurface topology.
  /// @pre Topologies cannot require a pyramid isotable (i.e., LINEAR_TOPOLOGY for dimension > 2 is not a valid.)
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & cube_isotable,
    const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info,
    const SCALAR_TYPE isovalue,
    const ISOSURFACE_TOPOLOGY isosurface_topology,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);

  /// Marching Cubes Algorithm using partition into uniform regions.
  /// Minimum and maximum scalar value are precomputed for each region.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
    MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm using octree.
  /// Minimum and maximum scalar value are precomputed for each octree node.
  void marching_cubes
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue, 
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord,
    MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree, 
    MCUBE_INFO & mcube_info);

// **************************************************
// NEP MARCHING CUBES (NEP: NEGATIVE-EQUALS-POSITIVE)
// **************************************************

  /// Marching Cubes Algorithm using NEP (negative, equals, positive) isosurface lookup table.
  void marching_cubes_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm using NEP (negative, equals, positive) isosurface lookup table and partition into uniform regions.
  void marching_cubes_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
     MCUBE_INFO & mcube_info);

  /// Marching Cubes Algorithm using NEP (negative, equals, positive) isosurface lookup table and octree.
  void marching_cubes_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree,
     MCUBE_INFO & mcube_info);

// **************************************************
// MCVOL: MARCHING CUBES INTERVAL VOLUME
// **************************************************

  /// Interval Volume Algorithm.
  void MCVol
   (const MC_SCALAR_GRID_BASE & scalar_grid, 
    const ISOSURFACE_TABLE & isotable,
    const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
    std::vector<VERTEX_INDEX> & simplex_vert, 
    std::vector<COORD_TYPE> & vertex_coord, MERGE_DATA & merge_data, 
    MCUBE_INFO & mcube_info);

// **************************************************
// MULTIRES CUBES
// **************************************************

  /// Multiresolution isosurface extraction using cubes, pyramids and tetrahedra.
  void multires_cubes
    (const MC_SCALAR_GRID_BASE & scalar_grid,  
     const MULTIRES_GRID & multires_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue, 
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);

// **************************************************
// CHECK FUNCTIONS
// **************************************************

  /// Check nep (negative/equals/positive) data structures.
  bool check_nep(const ISOSURFACE_TABLE & isotable, 
		 const MERGE_DATA & merge_data,
		 const int num_dup, IJK::ERROR & error);
};

#endif
