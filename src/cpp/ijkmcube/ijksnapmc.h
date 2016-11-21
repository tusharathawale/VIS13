/// \file ijksnapmc.h
/// Generate snapMC isosurface in arbitrary dimensions

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008 Rephael Wenger

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

#ifndef _IJKSNAPMC_
#define _IJKSNAPMC_

#include <string>

#include "ijk.txx"
#include "ijkgrid.txx"
#include "ijkmerge.txx"

#include "ijkmcube.h"
#include "ijktable.h"
#include "ijkoctree.h"

/// SnapMC classes and routines.
namespace IJKSNAPMC {

// **************************************************
// USING TYPES
// **************************************************

  using IJKTABLE::ISOSURFACE_TABLE;
  using IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON;

  using IJKMCUBE::SCALAR_TYPE;
  using IJKMCUBE::COORD_TYPE;
  using IJKMCUBE::GRID_COORD_TYPE;
  using IJKMCUBE::TABLE_INDEX;
  using IJKMCUBE::AXIS_SIZE_TYPE;
  using IJKMCUBE::VERTEX_INDEX;
  using IJKMCUBE::EDGE_INDEX;
  using IJKMCUBE::MERGE_INDEX;
  using IJKMCUBE::ISO_VERTEX_INDEX;
  using IJKMCUBE::ISOTABLE_TYPE;
  using IJKMCUBE::MC_SCALAR_GRID_BASE;
  using IJKMCUBE::NEP_ISOSURFACE_TABLE;
  using IJKMCUBE::MERGE_DATA;
  using IJKMCUBE::MINMAX_REGIONS;
  using IJKMCUBE::MC_DATA;
  using IJKMCUBE::MC_ISOSURFACE;
  using IJKMCUBE::SNAP_TYPE;
  using IJKMCUBE::SNAP_OCTREE;
  using IJKMCUBE::SNAP_MINMAX;

// **************************************************
// DATA STRUCTURES
// **************************************************

  class SNAP_INFO;
  class SNAP_GRID_BASE;
  class SNAP_GRID;

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

  /// SnapMC Algorithm.
  /// Improves quality of isosurface simplices by "snapping" isosurface vertices to grid vertices for isosurface triangulation construction.
  /// Isosurface vertices constructed by SnapMC are a subset of the isosurface vertices constructed by regular Marching Cubes.
  void snapMC
    (const MC_DATA & mc_data, const SCALAR_TYPE isovalue, 
     MC_ISOSURFACE & mc_isosurface, SNAP_INFO & snap_info);

  /// SnapMC Algorithm.
  /// Improves quality of isosurface simplices by "snapping" isosurface vertices to grid vertices for isosurface triangulation construction.
  /// Isosurface vertices constructed by SnapMC are a subset of the isosurface vertices constructed by regular Marching Cubes.
  /// @param cube_list = list of cubes intersecting the isosurface.
  void snapMC
    (const MC_DATA & mc_data, const SCALAR_TYPE isovalue,
     const std::vector<VERTEX_INDEX> & cube_list,
     MC_ISOSURFACE & mc_isosurface, SNAP_INFO & snap_info);

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

  /// SnapMC Algorithm.
  /// Improves quality of isosurface simplices by "snapping" isosurface vertices to grid vertices for isosurface triangulation construction.
  /// Isosurface vertices constructed by SnapMC are a subset of the isosurface vertices constructed by regular Marching Cubes.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  /// SnapMC Algorithm.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  /// SnapMC Algorithm.
  void snapMC_from_list
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     const std::vector<VERTEX_INDEX> & snap_vlist,
     const bool extract_from_boundary,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

  /// SnapMC Algorithm using partition into uniform regions.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
     SNAP_INFO & snap_info);

  /// SnapMC Algorithm using partition into uniform regions.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
     SNAP_INFO & snap_info);

  /// SnapMC Algorithm using octree.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const SNAP_TYPE snap_value,
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
     SNAP_INFO & snap_info);

  /// SnapMC Algorithm using octree.
  void snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
     const int nep_num_dup,
     std::vector<VERTEX_INDEX> & simplex_vert, 
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
     SNAP_INFO & snap_info);

// **************************************************
// SNAPMC SUBROUTINES
// **************************************************

  /// Snap scalar values which are close to isovalue.
  void snap_scalar_values
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value0, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snaptoward);

  /// Snap scalar values which are close to isovalue.
  void snap_scalar_values
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snap_back, const MINMAX_REGIONS & minmax);

  /// Snap scalar values which are close to isovalue.
  void snap_scalar_values
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
     VERTEX_INDEX * snap_back, const SNAP_OCTREE & octree);

  /// Compute position of isosurface vertices using linear interpolation.
  void position_snap_vertices_linear
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE * snap_grid, 
     const VERTEX_INDEX * snap_toward, const SCALAR_TYPE isovalue,
     const std::vector<ISO_VERTEX_INDEX> & vlist, 
     COORD_TYPE * coord, VERTEX_INDEX & num_snapped);

  /// Merge identical isosurface vertices and compute their position.
  void merge_and_position_vertices_snapMC
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
     const SCALAR_TYPE isovalue, 
     const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, SNAP_INFO & snap_info);

// **************************************************
// UTILITY FUNCTIONS
// **************************************************

  /// Get a list of cubes which intersect the isosurface.
  void get_nonempty_snap_cubes
    (const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & snap_vlist,
     float & procedure_time);

// **************************************************
// SNAP_GRID
// **************************************************

  /// Class storing snap grid information.
  /// scalar[] (inherited from MC_SCALAR_GRID_BASE) stores snapped scalar values.
  /// Note: SNAP_GRID_BASE does not provide any way to allocate memory
  ///    for scalar[] or snap_back[].
  class SNAP_GRID_BASE:public MC_SCALAR_GRID_BASE {
    
  protected:
    VERTEX_INDEX * snap_back;
    // grid vertex (x0,x1,x2,...) snaps back to isosurface vertex
    //     snap_toward[x0 + x1*axis_size[0] + 
    //                   x2*axis_size[0]*axis_size[1] + ...]

  public:
    // *** SHOULD USE MC_SCALAR_GRID_BASE::DIMENSION_TYPE ***
    SNAP_GRID_BASE(const int dimension, const AXIS_SIZE_TYPE * axis_size):
    MC_SCALAR_GRID_BASE(dimension, axis_size) 
      { this->snap_back = NULL; };  
    ~SNAP_GRID_BASE(){ this->snap_back = NULL; };
    // Note: constructor and destructor do not allocate or free any memory

    const VERTEX_INDEX * SnapBack() const { return(snap_back); };

    bool MatchesSize(const MC_SCALAR_GRID_BASE & scalar_grid,
		     IJK::ERROR & error_msg) const;
  };

  /// Class for creating snap grid information
  /// Note: In contrast with MC_SCALAR_GRID_BASE and SNAP_GRID_BASE,
  ///       this class creates and destroys memory to store the grid information
  class SNAP_GRID:public SNAP_GRID_BASE {

  protected:
    void Create(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
		const SNAP_TYPE snap_value);
    void Create
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const MINMAX_REGIONS & minmax);
    void Create
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const SNAP_OCTREE & octree);

  public:
    SNAP_GRID
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value);
    SNAP_GRID
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, float & creation_time);
    SNAP_GRID
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const MINMAX_REGIONS & minmax,
       float & creation_time);
    SNAP_GRID
      (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
       const SNAP_TYPE snap_value, const SNAP_OCTREE & octree,
       float & creation_time);
    ~SNAP_GRID();
    // Note: constructor and destructor ALLOCATE and FREE memory
  };


// **************************************************
// SNAP INFO
// **************************************************

  /// Snap information.
  class SNAP_INFO:public IJKMCUBE::MCUBE_INFO {
  public:
    VERTEX_INDEX num_snapped_iso_vertices;
    // number of snapped isosurface vertices

    SNAP_INFO() { Clear(); };
    SNAP_INFO(const int dimension):MCUBE_INFO(dimension) 
      { Clear(); };

    void Clear();     // clear all data
  };

};

#endif
