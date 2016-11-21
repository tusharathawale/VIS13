/// \file ijkmcube_util.h
/// ijkmcube utility functions.

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009 Rephael Wenger

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

#ifndef _IJKMCUBE_UTIL_
#define _IJKMCUBE_UTIL_

#include <string>

#include "ijk.txx"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"
#include "ijktable.h"

namespace IJKMCUBE {


// **************************************************
// UTILITY FUNCTIONS
// **************************************************

  /// Get grid cubes intersected by the isosurface
  void get_mixed_cubes
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const MINMAX_REGIONS & minmax,  const SCALAR_TYPE isovalue, 
     VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length);
  /// Get grid cubes intersected by the isosurface
  void get_mixed_cubes
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const IJKOCTREE::OCTREE & octree,  const SCALAR_TYPE isovalue, 
     VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length);

  /// Return number of isosurface vertices per each grid vertex.
  ///   Each isosurface vertex lies on a unique grid edge and 
  ///   each grid edge is associated with its "lower" endpoint.
  inline int get_num_iso_vertices_per_grid_vertex(const int dimension)
    { return(dimension); };
  /// Return number of isosurface vertices per each grid vertex
  ///   with NEP (negative-equals-positive) isosurface lookup table.
  ///   Each isosurface vertex lies on a unique grid edge 
  ///   or lies on a grid vertex.
  inline int get_num_nep_iso_vertices_per_grid_vertex(const int dimension)
    { return(dimension+1); };
  /// Return number of interval volume vertices per each grid vertex.
  ///   Each grid edge contains at most one "upper" interval volume vertex
  ///   and at most one "lower" interval volume vertex.
  ///   Each grid edge is associated with its "lower" endpoint.
  inline int get_num_ivol_vertices_per_grid_vertex(const int dimension)
    { return(2*dimension+1); };

  /// Compute isosurface vertex increment.
  ///   isov0+increment[k] is the isosurface vertex index of the (k-1)'th
  ///   isosurface vertex in the cube where isov0 is the vertex index
  ///   of the first isosurface vertex in the cube.
  void compute_iso_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  /// Compute isosurface endpoint increment.
  ///   (iv0+increment[2*i],iv0+increment[2*i+1]) are the vertex indices 
  ///   of the (i-1)'th edge of the cube where iv0 is the primary cube vertex.
  void compute_iso_endpoint_increment
    (const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
     VERTEX_INDEX * increment);
  /// Compute NEP isosurface vertex increment.
  ///   isov0+increment[k] is the NEP isosurface vertex index of the (k-1)'th
  ///   isosurface vertex in the cube where isov0 is the vertex index
  ///   of the first isosurface vertex in the cube.
  void compute_nep_iso_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  /// Compute interval volume vertex increment.
  ///   ivolv0+increment[k] is the interval volume vertex index of the (k-1)'th
  ///   interval volume vertex in the cube where ivolv0 is the vertex index
  ///   of the first interval volume vertex in the cube.
  void compute_ivol_vertex_increment
    (const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * vertex_increment, ISO_VERTEX_INDEX * increment);
  /// Compute hypercube edge increment.
  void compute_hypercube_edge_increment
    (const int dimension, const ISOSURFACE_TABLE & isotable, 
     const VERTEX_INDEX * hcube_vertex_increment, EDGE_INDEX * increment);
  /// Compute facet increment.
  void compute_facet_increment
    (const ISOSURFACE_TABLE_POLYHEDRON & poly, 
     const AXIS_SIZE_TYPE * axis_size, VERTEX_INDEX * increment);
  /// Compute facet vertex increment.
  void compute_facet_vertex_increment
    (const ISOSURFACE_TABLE_POLYHEDRON & poly, 
     const int orth_dir, const int side,
     const AXIS_SIZE_TYPE * axis_size, VERTEX_INDEX * facet_vertex_increment,
     VERTEX_INDEX * cube_vertex_index);

  /// Return file name of cube isosurface lookup table for given type and dimension.
  std::string get_isotable_filename
    (const ISOTABLE_TYPE type, const int dimension);
  /// Return file name of isosurface lookup table for given type, dimension and polytope.
  std::string get_isotable_filename
    (const ISOTABLE_TYPE type, const int dimension, 
     const char * poly_name);

// **************************************************
// TIMING FUNCTIONS
// **************************************************

  /// Convert elapsed CPU clock time to seconds
  inline float clock2seconds(clock_t t)
    { return(float(t)/CLOCKS_PER_SEC); }

// **************************************************
// STRINGS
// **************************************************

  /// Return string for given isosurface topology.
  std::string get_topology_string
    (const ISOSURFACE_TOPOLOGY isosurface_topology);

// **************************************************
// CHECK SUBROUTINES
// **************************************************

  /// Check NEP variable \a num_dup has valid values.
  /// If \a num_dup has invalid values, return false and set error message in \a error.
  /// @param num_dup = Number of simplices representing duplicate isosurface simplices.
  ///   If num_dup = 0, then duplicate simplices "cancel" and both are deleted.
  ///   If num_dup = 1, then duplicate simplices are replaced by a single simplex.
  ///   If num_dup = 2, then duplicate simplices are left in the isosurface.
  bool check_nep_num_dup(const int num_dup, IJK::ERROR & error);

  /// Check order of endpoints in list of endpoints.
  /// If endpoint[2*i] > endpoint[2*i+1], return false and set error message in \a error.
  /// @param endpoint[] = Array of endoints.  
  ///   Edge i has endpoints (endpoint[2*i], endpoint[2*i+1]).
  bool check_order(const std::vector<VERTEX_INDEX> & endpoint,
		   IJK::ERROR & error);

  /// Check isosurface table encoding.
  /// Returns false if isotable encoding does not match encoding.
  bool check_isotable_encoding
    (const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE::ENCODING encoding, IJK::ERROR & error);

  /// Check isosurface table encodings of all isosurface lookup tables in \a poly_isotable.
  /// Returns false if some isosurface table encoding does not match encoding.
  bool check_isotable_encoding
    (const POLY_ISOTABLE & poly_isotable, 
     const ISOSURFACE_TABLE::ENCODING encoding, IJK::ERROR & error);

  /// Check that ambiguity information table matches isosurface lookup table.
  /// Returns false if tables do not match.
  bool check_isotable_ambig_info
    (const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, IJK::ERROR & error);

  /// Check cube isosurface lookup table matches topology.
  ///   If isosurface topology is ASYMPTOTIC_DECIDER_TOPOLOGY or variants,
  ///   then check that configuration with two '-' vertices at opposing corners
  ///   should have exactly two simplices.
  bool check_cube_isotable_fits_topology
    (const ISOSURFACE_TABLE & cube_isotable, 
     ISOSURFACE_TOPOLOGY isosurface_topology, IJK::ERROR & error);

  /// Check isosurface table dimension matches scalar field dimension
  bool check_dimension
    (const ISOSURFACE_TABLE & isotable, 
     const MC_SCALAR_GRID_BASE & scalar_grid,
     IJK::ERROR & error);
};

#endif
