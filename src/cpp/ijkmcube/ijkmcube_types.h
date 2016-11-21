/// \file ijkmcube_types.h
/// Type definitions for ijkmcube.

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

#ifndef _IJKMCUBE_TYPES_
#define _IJKMCUBE_TYPES_

#include <string>

#include "ijk.txx"
#include "ijkscalar_grid.txx"
#include "ijkmerge.txx"

#include "ijktable.h"
#include "ijkoctree.h"


namespace IJKMCUBE {

// **************************************************
// IJKTABLE TYPES
// **************************************************

  using IJKTABLE::ISOSURFACE_TABLE;
  using IJKTABLE::ISOSURFACE_TABLE_AMBIG_INFO;
  using IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON;
  using IJKTABLE::TABLE_INDEX;

// **************************************************
// SCALAR TYPES
// **************************************************

  typedef float SCALAR_TYPE;     ///< Scalar value type.
  typedef float COORD_TYPE;      ///< Isosurface vertex coordinate type.
  typedef int VERTEX_INDEX;      ///< Grid vertex index type.
  typedef float GRADIENT_TYPE;   ///< Gradient coordinate type.
  typedef int GRID_COORD_TYPE;   ///< Grid vertex coordinate type.
  typedef int AXIS_SIZE_TYPE;    ///< Axis size type.
  typedef int ISO_VERTEX_INDEX;  ///< Isosurface vertex index type.
  typedef float SNAP_TYPE;       ///< Snap parameter type.
  typedef int MERGE_INDEX;       ///< Merge index type.

  /// Edge index type.
  /// Vertex and edge indices must have the same type.
  typedef VERTEX_INDEX EDGE_INDEX;

// **************************************************
// ARRAY TYPES
// **************************************************

  typedef std::vector<COORD_TYPE> COORD_ARRAY;   ///< Grid coordinate array.
  typedef std::vector<VERTEX_INDEX>              /// Vertex index array.
    VERTEX_INDEX_ARRAY; 
  typedef std::vector<SCALAR_TYPE> SCALAR_ARRAY; ///< Scalar array.

// **************************************************
// ENUMERATED TYPES
// **************************************************

  /// Interval volume vertex types.
  typedef enum { LOWER, MIDDLE, UPPER} INTERVAL_VOLUME_VERTEX_TYPE;

  /// Type of mesh edge representation
  /// - EDGE_ID = Represent edge by single integer identifier.
  /// - EDGE_ENDPOINT_PAIR = Represent edge by a pair of integer indentifiers of the endpoints.
  typedef enum { EDGE_ID, EDGE_ENDPOINT_PAIR } MESH_EDGE_REPRESENTATION;

  /// Isotable type.
  typedef enum { BINARY, NEP, IVOL } ISOTABLE_TYPE;

  /// Isosurface topology type.
  /// - ISOTABLE_TOPOLOGY: Isosurface lookup table specifies isosurface patch.
  /// - CUBE_DECIDER_TOPOLOGY: In ambiguous cubes with no ambiguous facets,
  ///   resolve ambiguity based on average scalar values of cube vertices.  
  /// - ASYMPTOTIC_DECIDER_TOPOLOGY: Break cubes with ambiguous facets 
  ///    into pyramids to resolve facet ambiguities.
  /// - SADDLE_TOPOLOGY: Break cubes at saddle points to APPROXIMATE trilinear topology.
  /// - LINEAR_TOPOLOGY: Generate topology matching the multilinear interpolant.
  /// - *** SADDLE_TOPOLOGY is experimental and LINEAR_TOPOLOGY is not yet implemented.
  typedef enum { ISOTABLE_TOPOLOGY, CUBE_DECIDER_TOPOLOGY,
		 ASYMPTOTIC_DECIDER_TOPOLOGY,
		 SADDLE_TOPOLOGY,
		 LINEAR_TOPOLOGY } ISOSURFACE_TOPOLOGY;

  /// Interpolation type.
  /// - LINEAR_INTERPOLATION: Determine the location of an isosurface vertex
  ///   using linear interpolation on the endpoints of the cube, pyramid
  ///   or simplex containing the isosurface vertex.
  /// - MULTILINEAR_INTERPOLATION: Determine the location of an 
  ///   isosurface vertex using multilinear interpolation on the cube vertices.
  /// change!

  typedef enum { LINEAR_INTERPOLATION, MULTILINEAR_INTERPOLATION, ALPHA_UNCERTAINTY}
  INTERPOLATION_TYPE;

  /// Type of vertices of multiresolution grid
  typedef enum { NOT_MULTIRES, MULTIRES_A, MULTIRES_B } 
  MULTIRES_VERTEX_TYPE;

// **************************************************
// CLASSES
// **************************************************

  typedef IJK::BOX<VERTEX_INDEX> GRID_BOX;  ///< Grid box type.

  typedef IJK::MERGE_PAIRS_PARAMETERS<int> 
    MERGE_EDGES_PARAMETERS;                 ///< Parameters for merging edges.
};

#endif
