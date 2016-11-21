/// \file ijkmcube_sub.h
/// Subroutines for ijkmcube.

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

#ifndef _IJKMCUBE_SUB_
#define _IJKMCUBE_SUB_

#include <string>

#include "ijk.txx"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"

namespace IJKMCUBE {

// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

  /// Merge identical isosurface vertices.
  void merge_identical_vertices
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const std::vector<ISO_VERTEX_INDEX> & vlist0,
     std::vector<ISO_VERTEX_INDEX> & vlist1, 
     std::vector<ISO_VERTEX_INDEX> & vlist0_map, MERGE_DATA & merge_data);

  /// Merge identical edges in list.
  void merge_identical_edges
    (const std::vector<VERTEX_INDEX> & elist0,
     std::vector<VERTEX_INDEX> & elist1, 
     std::vector<ISO_VERTEX_INDEX> & eindex,
     const MERGE_EDGES_PARAMETERS &);

  /// Compute position of isosurface vertices using linear interpolation.
  void position_iso_vertices_linear
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

  /// Compute position of isosurface vertices using linear interpolation.
  void position_iso_vertices_linearB
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord);

  /// Compute position of isosurface vertices using linear interpolation.
  void position_iso_vertices_linearB
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord);

   /// Compute position of isosurface vertices using alpha uncertainty.
  void position_iso_vertices_linearB_alpha
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color);

  /// Compute position of isosurface vertices using alpha_uncertainty.
  void position_iso_vertices_linearB_alpha
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color);


  /// Compute position of isosurface vertices using multilinear interpolation.
  void position_iso_vertices_multilinear
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord);

  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const std::vector<VERTEX_INDEX> & endpoint,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);

  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & endpoint,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info);

  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & endpoint,
     const ISOTABLE_TYPE isotable_type,
     const INTERPOLATION_TYPE interpolation_type,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info);

  // newly added!!
  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & endpoint,
     const ISOTABLE_TYPE isotable_type,
     const INTERPOLATION_TYPE interpolation_type,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color); 


  /// Convert isosurface vertex indices to pairs of endpoints
  /// representing grid edges containing isosurface vertices.
  void convert_indices_to_endpoint_pairs
    (const MC_SCALAR_GRID_BASE & scalar_grid,
     const std::vector<ISO_VERTEX_INDEX> iso_vlist,
     const MERGE_DATA & merge_data,
     std::vector<VERTEX_INDEX> & endpoint);

// **************************************************
// NEP MARCHING CUBES SUBROUTINES
// **************************************************

  /// Identify isosurface patches which lie in a facet.
  /// Set is_in_facet flag for each entry in isosurface lookup table.
  void set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable);

// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

  /// Compute position of interval volume vertices using linear interpolation.
  void position_ivol_vertices_linear
    (const MC_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

};

#endif
