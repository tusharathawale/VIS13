/// \file ijkmcube_extract.h
/// Subroutines for extracting isosurface mesh.

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

#ifndef _IJKMCUBE_EXTRACT_
#define _IJKMCUBE_EXTRACT_

#include <string>

#include "ijk.txx"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"


namespace IJKMCUBE {

// **************************************************
// EXTRACT ISOSURFACE MESH
// **************************************************

  /// Extract isosurface simplices from grid cubes.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices from grid cubes.
  void extract_iso_simplices_and_cube_info
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & cube_list,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices from grid cubes.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);

  // Newly added!	
  /// Extract isosurface simplices from grid cubes.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info, const MC_SCALAR_GRID_BASE & scalar_grid_1, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices_1,
     std::vector<VERTEX_INDEX> & endpoint_1,
     MCUBE_INFO & mcube_info_1);


  /// Extract isosurface simplices with specified isosurface_topology.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     const ISOSURFACE_TOPOLOGY isosurface_topology,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

    // newly added
    /// Extract isosurface simplices with specified isosurface_topology.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     const ISOSURFACE_TOPOLOGY isosurface_topology,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info, const MC_SCALAR_GRID_BASE & scalar_grid_1, 
     const ISOSURFACE_TOPOLOGY isosurface_topology_1,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices_1,
     std::vector<VERTEX_INDEX> & endpoint_1, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices_1,
     MCUBE_INFO & mcube_info_1);



  /// Extract isosurface simplices using the asymptotic decider.
  void extract_iso_simplices_adecider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices.
  /// Break cubes at saddle points to approximate linear topology.
  void extract_iso_simplices_saddle
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices resolving some cube ambiguities.
  /// In ambiguous cubes with no ambiguous facets, resolve ambiguity
  ///   based on average scalar values of cube vertices.
  void extract_iso_simplices_cube_decider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, 
     const SCALAR_TYPE isovalue, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices resolving some cube ambiguities.
  /// In ambiguous cubes with no ambiguous facets, resolve ambiguity
  ///   based on average scalar values of cube vertices.
  void extract_iso_simplices_cube_decider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, 
     const SCALAR_TYPE isovalue, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices from grid cubes in list.
  void extract_iso_simplices_from_list
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, 
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices using octree.
  void extract_iso_simplices_from_octree
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info,  const IJKOCTREE::OCTREE & octree);

  /// Extract isosurface simplices using partition into uniform regions.
  void extract_iso_simplices_from_minmax
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info, const MINMAX_REGIONS & minmax);

// **************************************************
// EXTRACT ISOSURFACE MESH USING NEP ISOTABLE
// **************************************************

  /// Extract isosurface simplices from grid cubes.
  /// Use negative, equals, positive (nep) isosurface lookup table.
  void extract_iso_simplices_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info);

  /// Extract nep isosurface vertices lying on the grid boundary.
  void extract_iso_simplices_nep_boundary
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int num_cube_facet_vertices,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_non_empty_cubes);

  /// Extract isosurface simplices using octree.
  /// Use negative, equals, positive (nep) isosurface lookup table.
  void extract_iso_simplices_from_octree_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
     const IJKOCTREE::OCTREE & octree);

  /// Extract isosurface simplices using partition into uniform regions.
  /// Use negative, equals, positive (nep) isosurface lookup table.
  void extract_iso_simplices_from_minmax_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
     const MINMAX_REGIONS & minmax);

  /// Extract isosurface simplices from list of cubes.
  /// Use negative, equals, positive (nep) isosurface lookup table.
  void extract_iso_simplices_from_list_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     const bool extract_from_boundary,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

// **************************************************
// EXTRACT INTERVAL VOLUME MESH
// **************************************************

  /// Extract interval volume simplices.
  void extract_ivol_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_non_empty_cubes);

  /// Extract interval volume simplices from list of cubes.
  void extract_ivol_simplices_from_list
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_non_empty_cubes);

// **************************************************
// EXTRACT ISOSURFACE MESH FROM MULTIRESOLUTION GRID
// **************************************************

  /// Extract isosurface simplices from multiresolution cubes, pyramids and tetrahedra.
  void multires_extract
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const MULTIRES_GRID & multires_grid,
     const POLY_ISOTABLE & poly_isotable, const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);
};

#endif
