/// \file ijksnapmc.cxx
/// Snap MC isosurface generation

/*
  IJK: Isosurface Jeneration Kode
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

#include <assert.h>
#include <sstream>
#include <string>
#include <vector>

#include "ijksnapmc.h"
#include "ijkmcube_util.h"
#include "ijkmcube_extract.h"
#include "ijkoctree.h"
#include "ijktable.h"

#include "ijkinterpolate.txx"
#include "ijkmerge.txx"

using namespace IJK;
using namespace IJKSNAPMC;
using namespace IJKTABLE;

using IJK::linear_interpolate_coord;
using IJKMCUBE::compute_iso_vertex_increment;
using IJKMCUBE::compute_nep_iso_vertex_increment;
using IJKMCUBE::compute_ivol_vertex_increment;
using IJKMCUBE::get_num_iso_vertices_per_grid_vertex;
using IJKMCUBE::get_num_ivol_vertices_per_grid_vertex;
using IJKMCUBE::get_num_nep_iso_vertices_per_grid_vertex;
using IJKMCUBE::check_nep;
using IJKMCUBE::check_isotable_encoding;
using IJKMCUBE::clock2seconds;

using IJKMCUBE::NEP;
using IJKMCUBE::BINARY;
using IJKMCUBE::IVOL;

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

/// SnapMC Algorithm.
/// Improves quality of isosurface simplices by "snapping" isosurface vertices to grid vertices for isosurface triangulation construction.
/// Isosurface vertices constructed by SnapMC are a subset of the isosurface vertices constructed by regular Marching Cubes.
void IJKSNAPMC::snapMC
(const MC_DATA & mc_data, const SCALAR_TYPE isovalue, 
 MC_ISOSURFACE & mc_isosurface, SNAP_INFO & snap_info)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = mc_data.ScalarGrid().AxisSize();
  float merge_time = 0.0;
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!mc_data.Check(error)) { throw error; };

  mc_isosurface.Clear();
  snap_info.time.Clear();

  IJKMCUBE::NEP_ISO_MERGE_DATA merge_data(dimension, axis_size);
  clock_t t1 = clock();
  merge_time = clock2seconds(t1-t_start);

  const bool snap = mc_data.Snap();
  const int nep_num_dup = mc_data.NEPNumDup();
  const SNAP_TYPE snap_value = mc_data.SnapValue();

  if (mc_data.UseOctree()) {
    snapMC(mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
	   isovalue, snap_value, nep_num_dup,
	   mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
	   merge_data, *(mc_data.SnapOctree()), snap_info);
  }
  else if (mc_data.UseMinmaxRegions()) {
    snapMC(mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
	   isovalue, snap_value, nep_num_dup,
	   mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
	   merge_data, *(mc_data.SnapMinmaxRegions()), snap_info);
  }
  else {
    snapMC(mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
	   isovalue, snap_value, nep_num_dup,
	   mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
	   merge_data, snap_info);
  }

  snap_info.time.extract += merge_time;

  // store times
  clock_t t_end = clock();
  snap_info.time.total = clock2seconds(t_end-t_start);
}

/// SnapMC Algorithm.
/// Improves quality of isosurface simplices by "snapping" isosurface vertices to grid vertices for isosurface triangulation construction.
/// Isosurface vertices constructed by SnapMC are a subset of the isosurface vertices constructed by regular Marching Cubes.
void IJKSNAPMC::snapMC
(const MC_DATA & mc_data, const SCALAR_TYPE isovalue, 
 const std::vector<VERTEX_INDEX> & cube_list,
 MC_ISOSURFACE & mc_isosurface, SNAP_INFO & snap_info)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = mc_data.ScalarGrid().AxisSize();
  float merge_time = 0.0;
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!mc_data.Check(error)) { throw error; };

  mc_isosurface.Clear();
  snap_info.time.Clear();

  IJKMCUBE::NEP_ISO_MERGE_DATA merge_data(dimension, axis_size);
  clock_t t1 = clock();
  merge_time = clock2seconds(t1-t_start);

  const bool snap = mc_data.Snap();
  const int nep_num_dup = mc_data.NEPNumDup();
  const SNAP_TYPE snap_value = mc_data.SnapValue();

  if (mc_data.UseList()) {

    SNAP_GRID snap_grid
      (mc_data.ScalarGrid(), isovalue, snap_value, snap_info.time.snap);

    snapMC_from_list
      (mc_data.ScalarGrid(), snap_grid, mc_data.isotable.cube_nep,
       isovalue, nep_num_dup, cube_list, true, 
       mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
       merge_data, snap_info);
  }
  else {
    snapMC(mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
	   isovalue, snap_value, nep_num_dup,
	   mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
	   merge_data, snap_info);
  }

  snap_info.time.extract += merge_time;

  // store times
  clock_t t_end = clock();
  snap_info.time.total = clock2seconds(t_end-t_start);
}

// **************************************************
// SNAPMC: MARCHING CUBES WITH QUALITY TRIANGLES
// **************************************************

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, 
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, 
		      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  if (scalar_grid.NumVertices() != snap_grid.NumVertices()) {
    error.AddMessage
      ("Programming error. Snap grid does not match scalar grid.");
    error.AddMessage
      ("  Snap grid has ", snap_grid.NumVertices(), " vertices.");
    error.AddMessage
      ("  Scalar grid has ", scalar_grid.NumVertices(), " vertices.");
    throw error;
  }

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC_from_list
(const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 const std::vector<VERTEX_INDEX> & snap_vlist,
 const bool extract_from_boundary,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("snapMC_from_list");

  if (!snap_grid.MatchesSize(scalar_grid, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };
  if (!merge_data.Check(NEP, error)) { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  clock_t t0 = clock();

  extract_iso_simplices_from_list
    (snap_grid, isotable, isovalue, &snap_vlist[0], snap_vlist.size(), 
     iso_simplices, snap_info);

  GRID_SIZE_TYPE num_cube_facet_vertices = 
    compute_num_cube_facet_vertices(dimension);

  snap_info.nep.num_non_empty_boundary_facets = 0;
  if (extract_from_boundary) {
    VERTEX_INDEX num_mixed_cubes;
    extract_iso_simplices_nep_boundary
      (snap_grid, isotable, isovalue, num_cube_facet_vertices,
       iso_simplices, num_mixed_cubes);
    snap_info.nep.num_non_empty_boundary_facets = num_mixed_cubes;
  };

  clock_t t2 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size,
			   iso_simplices, iso_vlist, simplex_vert, merge_data);
  clock_t t3 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  VERTEX_INDEX num_snapped;
  position_snap_vertices_linear
    (scalar_grid, snap_grid.ScalarPtrConst(), snap_grid.SnapBack(),
     isovalue, iso_vlist, &(vertex_coord[0]), num_snapped);
  snap_info.num_snapped_iso_vertices = num_snapped;

  clock_t t4 = clock();

  // store times
  snap_info.time.extract = float(t2-t0)/CLOCKS_PER_SEC;
  snap_info.time.merge = float(t3-t2)/CLOCKS_PER_SEC;
  snap_info.time.position = float(t4-t3)/CLOCKS_PER_SEC;
  snap_info.time.total = float(t4-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_MINMAX & minmax,
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// minmax = snapMC min and max of regions
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, minmax,
		      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, 
	 minmax, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_MINMAX & minmax, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_minmax_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info,
     minmax);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,  const SNAP_TYPE snap_value,
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_OCTREE & octree,
 SNAP_INFO & snap_info)
// extract isosurface using Marching Cubes algorithm
// returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// snap_value = controls snapping. 0 = no snapping. 0.5 = max snapping
// nep_num_dup = controls duplicate isosurface patches.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// octree = snapMC octree
// snap_info = information about running time and grid cubes and edges
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  SNAP_GRID snap_grid(scalar_grid, isovalue, snap_value, octree,
		      snap_info.time.snap);

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  snapMC(scalar_grid, snap_grid, isotable, isovalue,
	 nep_num_dup, simplex_vert, vertex_coord, merge_data, 
	 octree, snap_info);

  clock_t t_end = clock();

  // store time
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const NEP_ISOSURFACE_TABLE & isotable, const SCALAR_TYPE isovalue, 
 const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const SNAP_OCTREE & octree, SNAP_INFO & snap_info)
{
  PROCEDURE_ERROR error("snapMC");

  clock_t t_start = clock();

  if (!check_nep(isotable, merge_data, nep_num_dup, error)) { throw error; }

  simplex_vert.clear();
  vertex_coord.clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_octree_nep
    (snap_grid, isotable, isovalue, nep_num_dup, iso_simplices, snap_info,
     octree);

  merge_and_position_vertices_snapMC
    (scalar_grid, snap_grid, isovalue, iso_simplices,
     simplex_vert, vertex_coord, merge_data, snap_info);

  // store times
  clock_t t_end = clock();
  snap_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// **************************************************
// LINEAR INTERPOLATION
// **************************************************

namespace {

  inline void linear_weights
  (const int dimension, const SCALAR_TYPE s0, const SCALAR_TYPE s1, 
   const SCALAR_TYPE isovalue, COORD_TYPE & w0, COORD_TYPE & w1)
  {
    SCALAR_TYPE s_diff = s1 - s0;
    const double EPSILON = 0.00001;
    if (s_diff > EPSILON || s_diff < -EPSILON) { 
      w0 = (s1 - isovalue) / s_diff;
      w1 = (isovalue - s0) / s_diff;
    }
    else {
      // arbitrarily set weights to 0.5
      w0 = w1 = 0.5;
    };
  }

}


// **************************************************
// SNAPMC SUBROUTINES
// **************************************************

namespace {
  typedef ARRAY<VERTEX_INDEX> VERTEX_ARRAY;

  inline void snap_vertex
  (const VERTEX_INDEX v0, const VERTEX_INDEX v1, const SNAP_TYPE dist,
   const SCALAR_TYPE * scalar, const SCALAR_TYPE isovalue, 
   SCALAR_TYPE * scalar_snap, VERTEX_INDEX * snap_back, 
   SNAP_TYPE * snap_dist)
  {
    if (scalar_snap[v0] != isovalue) {
      scalar_snap[v0] = isovalue;
      snap_back[v0] = v1;
      snap_dist[v0] = dist;
    }
    else if (snap_dist[v0] > dist) {
      scalar_snap[v0] = isovalue;
      snap_back[v0] = v1;
      snap_dist[v0] = dist;
    }
    else if (snap_dist[v0] < dist) {
      return;
    }
    else {
      // snap_dist[v0] == dist
      // break distance tie by using lower vertex index for snap_back[v0]
      if (v1 < snap_back[v0]) {
	snap_back[v0] = v1;
      }
    }
  }

  inline void if_close_then_snap
  (const int dimension,
   const VERTEX_INDEX v0, const SCALAR_TYPE s0,
   const VERTEX_INDEX v1,const SCALAR_TYPE s1,
   const SCALAR_TYPE * scalar, const SCALAR_TYPE isovalue,
   const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap, 
   VERTEX_INDEX * snap_back, SNAP_TYPE * snap_dist)
  {
    COORD_TYPE w0, w1;

    if ((s0 < isovalue && s1 > isovalue) ||
	(s0 > isovalue && s1 < isovalue)) {
      // isosurface intersects edge

      linear_weights(dimension, s0, s1, isovalue, w0, w1);

      if (w1 < snap_value) {

	// snap isosurface vertex to v0.  Distance = w1.
	snap_vertex(v0, v1, w1, scalar, isovalue,
		    scalar_snap, snap_back, snap_dist);
      }
      else if (w0 < snap_value) {

	// snap isosurface vertex to v1.  Distance = w1.
	snap_vertex(v1, v0, w0, scalar, isovalue,
		    scalar_snap, snap_back, snap_dist);
      }
    }
  }

  inline void snap_incident_edges
  (const int dimension, const AXIS_SIZE_TYPE * axis_increment,
   const SCALAR_TYPE * scalar, const VERTEX_INDEX v0,
   const SCALAR_TYPE isovalue,
   const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap, 
   VERTEX_INDEX * snap_back, SNAP_TYPE * snap_dist)
  {
    SCALAR_TYPE s0 = scalar[v0];
    for (int d = 0; d < dimension; d++) {
      VERTEX_INDEX v1 = v0 + axis_increment[d];
      SCALAR_TYPE s1 = scalar[v1];

      if_close_then_snap
	(dimension, v0, s0, v1, s1, scalar, isovalue,
	 snap_value, scalar_snap, snap_back, snap_dist);
    }
  }

  void snap_scalar_values_in_upper_facet
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
   const SNAP_TYPE snap_value, int orth_dir,
   SCALAR_TYPE * scalar_snap,
   VERTEX_INDEX * snap_back, SNAP_TYPE * snap_dist)
  // snap scalar values in upper facet
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    VERTEX_INDEX axis_increment[scalar_grid.Dimension()];

    if (scalar_grid.Dimension() < 2) { return; };
    if (scalar_grid.AxisSize(orth_dir) < 1) { return; };

    compute_increment(scalar_grid, axis_increment);

    const VERTEX_INDEX lower2upper_increment =
      axis_increment[orth_dir]*(scalar_grid.AxisSize(orth_dir)-1);

    VERTEX_INDEX v0, v1;
    SCALAR_TYPE s0, s1;
    COORD_TYPE w0, w1;
    for (int d = 0; d < dimension; d++) {

      if (d != orth_dir) {
	GRID_SIZE_TYPE num_ridge_vertices;
	compute_num_vertices_in_grid_ridge
	  (dimension, axis_size, orth_dir, d, num_ridge_vertices);

	VERTEX_ARRAY ridge_vlist(num_ridge_vertices);
	get_vertices_in_grid_ridge(dimension, axis_size, orth_dir, d, 
				   ridge_vlist.Ptr());


	for (int i = 0; i < num_ridge_vertices; i++) {
	  // convert ridge on lower facet to ridge on upper facet
	  ridge_vlist[i] += lower2upper_increment;
	}

	for (VERTEX_INDEX iv = 0; iv < num_ridge_vertices; iv++) {

	  v0 = ridge_vlist[iv];
	  s0 = scalar[v0];

	  for (int k = 0; k < axis_size[d]-1; k++) {

	    v1 = v0 + axis_increment[d];
	    s1 = scalar[v1];

	    if_close_then_snap(dimension, v0, s0, v1, s1, scalar, isovalue,
			       snap_value, scalar_snap, snap_back, snap_dist);

	    // set v0 to v1
	    v0 = v1;
	    s0 = s1;
	  }
	}

      }
    }
  }

}


void IJKSNAPMC::snap_scalar_values
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  if (dimension < 1) { return; };

  const int d_last = dimension-1;
  GRID_SIZE_TYPE num_cubes_in_facet;
  compute_num_cubes_in_grid_facet
    (dimension, axis_size, d_last, num_cubes_in_facet);

  VERTEX_ARRAY facet_cube_list(num_cubes_in_facet);
  get_cubes_in_grid_facet(dimension, axis_size, d_last, 0, 
			  facet_cube_list.Ptr());

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  VERTEX_INDEX v0, v1;
  SCALAR_TYPE s0, s1;
  COORD_TYPE w0, w1;

  for (int k = 0; k+1 < axis_size[d_last]; k++) {

    VERTEX_INDEX facet_increment = k*axis_increment[d_last];
    for (VERTEX_INDEX icube = 0; icube < num_cubes_in_facet; icube++) {
      v0 = facet_cube_list[icube] + facet_increment;
      s0 = scalar[v0];
      for (int d = 0; d < dimension; d++) {
	v1 = v0 + axis_increment[d];
	s1 = scalar[v1];

	if_close_then_snap(dimension, v0, s0, v1, s1, scalar, isovalue,
			   snap_value, scalar_snap, snap_back, 
			   snap_dist.Ptr());
      }
    }
  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d,
       scalar_snap, snap_back, snap_dist.Ptr());
  }
}


void IJKSNAPMC::snap_scalar_values
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back, const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  const AXIS_SIZE_TYPE region_edge_length = minmax.RegionEdgeLength();
  const GRID_SIZE_TYPE num_cubes = scalar_grid.ComputeNumCubes();
  VERTEX_INDEX max_num_cubes_in_region;
  compute_num_grid_cubes_in_region
    (dimension, region_edge_length, max_num_cubes_in_region);

  VERTEX_INDEX num_regions;
  compute_num_regions(dimension, axis_size, region_edge_length, num_regions);

  IJK::ARRAY<VERTEX_INDEX> rlist(num_regions);
  IJK::ARRAY<bool> is_full(num_regions);
  IJK::ARRAY<VERTEX_INDEX> region_cube_increment(max_num_cubes_in_region);

  compute_region_cube_increment
    (dimension, axis_size, region_edge_length, region_cube_increment.Ptr());

  get_region_primary_vertices
    (dimension, axis_size, region_edge_length, rlist.Ptr(), is_full.Ptr());

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  IJK::ARRAY<VERTEX_INDEX> vlist(max_num_cubes_in_region);
  VERTEX_INDEX num_cubes_in_region;
  for (VERTEX_INDEX iregion = 0; iregion < num_regions; iregion++) {

    if (minmax.Min(iregion) < isovalue &&
	minmax.Max(iregion) >= isovalue) {

      VERTEX_INDEX v0 = rlist[iregion];
      if (is_full[iregion]) {
	for (VERTEX_INDEX j = 0; j < max_num_cubes_in_region; j++) {
	  VERTEX_INDEX v1 = v0 + region_cube_increment[j];

	  snap_incident_edges
	    (dimension, axis_increment, scalar, v1, isovalue,
	     snap_value, scalar_snap, snap_back, snap_dist.Ptr());
	}
      }
      else {
	get_grid_cubes_in_region
	  (dimension, axis_size, v0, region_edge_length,
	   vlist.Ptr(), num_cubes_in_region);
	for (VERTEX_INDEX j = 0; j < num_cubes_in_region; j++) {
	  VERTEX_INDEX v1 = vlist[j];

	  snap_incident_edges
	    (dimension, axis_increment, scalar, v1, isovalue,
	     snap_value, scalar_snap, snap_back, snap_dist.Ptr());
	}
      }
    }

  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d, scalar_snap, snap_back, 
       snap_dist.Ptr());
  }
}

void IJKSNAPMC::snap_scalar_values
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, SCALAR_TYPE * scalar_snap,
 VERTEX_INDEX * snap_back, const SNAP_OCTREE & octree)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  PROCEDURE_ERROR error("snape_scalar_values");

  assert(scalar != NULL && snap_back != NULL);

  if (snap_value < 0.0 || snap_value > 0.5) {
    error.AddMessage("Snap value must be in range [0.0, 0.5].");
    throw error;
  }

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();

  // copy scalar_grid to scalar_snap
  std::copy(scalar, scalar+num_grid_vertices, scalar_snap);

  const GRID_SIZE_TYPE num_cubes = scalar_grid.ComputeNumCubes();

  ARRAY<SNAP_TYPE> snap_dist(num_grid_vertices);

  int num_levels = octree.NumLevels();
  if (num_levels == 0) return;
  int ileaf_level = num_levels-1;

  // Depth first search through octree
  IJKOCTREE::OCTREE_STACK stack(num_levels);

  IJKOCTREE::CONST_NODE_PTR root = octree.Root();
  stack.PushRoot(root);

  VERTEX_INDEX v0, v1;
  SCALAR_TYPE s0, s1;
  while (!stack.IsEmpty()) {

    IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
    if (node->MinValue() < isovalue && isovalue <= node->MaxValue()) {

      if (stack.TopIsLeaf()) {
	// leaf node
	IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
	for (int i = 0; i < node->NumChildren(); i++) {

	  v0 = octree.ComputeChildVertex0(ileaf_level, node, i);
	  s0 = scalar[v0];
	  for (int d = 0; d < dimension; d++) {
	    v1 = v0 + axis_increment[d];
	    s1 = scalar[v1];

	    if_close_then_snap
	      (dimension, v0, s0, v1, s1, scalar, isovalue,
	       snap_value, scalar_snap, snap_back, snap_dist.Ptr());
	  }

	}

	// pop stack
	stack.Pop();
      }
      else {
	// internal node
	if (stack.NoNextChild()) {
	  // pop stack
	  stack.Pop();
	}
	else {
	  // push child_node onto stack
	  stack.PushChild();
	}
      }
    }
    else {
      // pop stack
      stack.Pop();
    }
  }

  // snap scalar value in upper facets
  for (int d = 0; d < dimension; d++) {
    snap_scalar_values_in_upper_facet
      (scalar_grid, isovalue, snap_value, d,
       scalar_snap, snap_back, snap_dist.Ptr());
  }
}

void IJKSNAPMC::position_snap_vertices_linear
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE * scalar_snap, 
 const VERTEX_INDEX * snap_back, const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & vlist, 
 COORD_TYPE * coord, VERTEX_INDEX & num_snapped)
// compute position of isosurface vertices using linear interpolation
//   vertices can have positive, negative or equals values
// scalar_grid = scalar grid data
// scalar_snap[] = array of snapped values
// snap_back[i] = move i'th isosurface vertex toward snap_back[i]
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
// num_snapped = number of snapped isosurface vertices
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  const int VERTEX_OFFSET = dimension;

  num_snapped = 0;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_isov_per_gridv = 
    get_num_nep_iso_vertices_per_grid_vertex(dimension);

  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_isov_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;

    if (k == VERTEX_OFFSET) {
      if (scalar_snap[v0] == scalar[v0]) {
	// isosurface vertex lies on grid vertex v0

	compute_coord(v0, dimension, axis_size, coord0);
	for (int d = 0; d < dimension; d++)
	  coord[i*dimension+d] = coord0[d];
      }
      else {
	// isosurface vertex was snapped to grid vertex v0
	// move isosurface vertex to interior of edge [v0, v1]

	VERTEX_INDEX v1 = snap_back[v0];

	SCALAR_TYPE s0 = scalar[v0];
	SCALAR_TYPE s1 = scalar[v1];

	compute_coord(v0, dimension, axis_size, coord0);
	compute_coord(v1, dimension, axis_size, coord1);

	linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue, 
			   coord2);
	for (int d = 0; d < dimension; d++)
	  coord[i*dimension+d] = coord2[d];

	num_snapped++;
      }
    }
    else {
      // isosurface vertex lies on grid edge
      VERTEX_INDEX v1 = v0+axis_increment[k];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue, 
			 coord2);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

}

void IJKSNAPMC::merge_and_position_vertices_snapMC
(const MC_SCALAR_GRID_BASE & scalar_grid, const SNAP_GRID_BASE & snap_grid,
 const SCALAR_TYPE isovalue, 
 const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, SNAP_INFO & snap_info)
// call merge_identical_vertices and then position_vertices_linear
// scalar_grid = scalar grid data
// snap_grid = snap grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// snap_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices_snapMC");

  if (!merge_data.Check(NEP, error)) { throw error; };

  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size,
			   iso_simplices, iso_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  VERTEX_INDEX num_snapped;
  position_snap_vertices_linear
    (scalar_grid, snap_grid.ScalarPtrConst(), snap_grid.SnapBack(),
     isovalue, iso_vlist, &(vertex_coord[0]), num_snapped);
  snap_info.num_snapped_iso_vertices = num_snapped;

  clock_t t3 = clock();

  // store times
  snap_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  snap_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}

void IJKSNAPMC::get_nonempty_snap_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, std::vector<VERTEX_INDEX> & snap_vlist,
 float & procedure_time)
// get all possible nonempty cubes under all possible snap values
// snap_vlist may contain some extra cubes
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  VERTEX_INDEX axis_increment[dimension];
  SCALAR_TYPE * scalar_snap = NULL;  // snapped scalar values
  VERTEX_INDEX * snap_back = NULL;   // snap grid vertex iv toward 
                                     //   grid vertex snapto[iv]
  const SNAP_TYPE max_snap = 0.5;

  PROCEDURE_ERROR error("get_nonempty_snap_cubes");

  clock_t t_start = clock();

  snap_vlist.clear();

  // create grid of snapped scalar values
  GRID_SIZE_TYPE num_grid_vertices;
  compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
  scalar_snap = new SCALAR_TYPE[num_grid_vertices];
  snap_back = new VERTEX_INDEX[num_grid_vertices];

  // NOTE: SHOULD REPLACE THIS BY A FAST ROUTINE WHICH SIMPLY
  //  SETS TO SCALAR VALUE ANY VERTEX WITHIN 0.5 OF AN ISOSURFACE VERTEX
  snap_scalar_values
    (scalar_grid, isovalue, max_snap, scalar_snap, snap_back);

  VERTEX_INDEX num_cube_vertices = isotable.Polyhedron().NumVertices();
  VERTEX_ARRAY vertex_increment(num_cube_vertices);
  compute_increment(dimension, axis_size, axis_increment);
  compute_cube_vertex_increment
    (dimension, axis_increment, vertex_increment.Ptr());

  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_INDEX * facet_vlist = new VERTEX_INDEX[num_cubes_in_facet0];
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist);

  for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
    for (VERTEX_INDEX iv0 = facet_vlist[j]; 
	 iv0 < facet_vlist[j] + axis_size[0]-1; iv0++) {

      bool found_equals = false;
      for (int k = 0; k < num_cube_vertices; k++) {
	VERTEX_INDEX iv1 = iv0 + vertex_increment[k];
	if (scalar_snap[iv1] == isovalue) {
	  found_equals = true;
	  break;
	};
      };

      if (found_equals) {
	snap_vlist.push_back(iv0);
      }
      else {
	bool found_negative = false;
	for (int k = 0; k < num_cube_vertices; k++) {
	  VERTEX_INDEX iv1 = iv0 + vertex_increment[k];
	  if (scalar_snap[iv1] < isovalue) {
	    found_negative = true;
	    break;
	  };
	};

	bool found_positive = false;
	for (int k = 0; k < num_cube_vertices; k++) {
	  VERTEX_INDEX iv1 = iv0 + vertex_increment[k];
	  if (scalar_snap[iv1] >= isovalue) {
	    found_positive = true;
	    break;
	  };
	};

	if (found_positive && found_negative)
	  { snap_vlist.push_back(iv0); }
      }
    }
  }

  delete [] facet_vlist;
  delete [] scalar_snap;
  delete [] snap_back;

  clock_t t_end = clock();
  procedure_time = clock2seconds(t_end-t_start);
}

// **************************************************
// SNAP_GRID
// **************************************************

bool SNAP_GRID_BASE::MatchesSize
(const MC_SCALAR_GRID_BASE & scalar_grid, ERROR & error) const
{
  const AXIS_SIZE_TYPE * axis_size = AxisSize();
  const AXIS_SIZE_TYPE * scalar_axis_size = scalar_grid.AxisSize();

  if (Dimension() != scalar_grid.Dimension()) {
    error.AddMessage("Snap grid dimension (", Dimension(), 
		     ") does not equal scalar grid dimension(",
		     scalar_grid.Dimension(), ").");
    return(false);
  }

  for (int d = 0; d < Dimension(); d++) {
    if (axis_size[d] != scalar_axis_size[d]) {
      error.AddMessage
	("Snap grid axis ", d , " size does not equal scalar grid axis ", 
	 d, " size.");
      error.AddMessage
	("Snap grid axis size: ", axis_size[d], 
	 ".  Scalar grid axis size: ", scalar_axis_size[d], ".");
      return(false);
    }
  }

  return(true);
}

// class for creating snap grid information
// Note: In contrast with MC_SCALAR_GRID_BASE and SNAP_GRID_BASE,
//       the constructor creates memory to store the grid information
SNAP_GRID::SNAP_GRID
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  Create(scalar_grid, isovalue, snap_value);
}

SNAP_GRID::SNAP_GRID
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::SNAP_GRID
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value,  const IJKMCUBE::MINMAX_REGIONS & minmax, 
 float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value, minmax);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::SNAP_GRID
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value,  const SNAP_OCTREE & octree, 
 float & creation_time):
  SNAP_GRID_BASE(scalar_grid.Dimension(), scalar_grid.AxisSize())
{
  clock_t t0 = clock();
  Create(scalar_grid, isovalue, snap_value, octree);
  clock_t t1 = clock();
  creation_time = float(t1-t0)/CLOCKS_PER_SEC;
}

SNAP_GRID::~SNAP_GRID()
// destructor
{
  delete [] scalar;
  scalar = NULL;
  delete [] snap_back;
  snap_back = NULL;
}

void SNAP_GRID::Create
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back);
}

void SNAP_GRID::Create
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back,
     minmax);
}

void SNAP_GRID::Create
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const SNAP_TYPE snap_value, const SNAP_OCTREE & octree)
{
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  this->scalar = new SCALAR_TYPE[num_grid_vertices];
  this->snap_back = new VERTEX_INDEX[num_grid_vertices];
  snap_scalar_values
    (scalar_grid, isovalue, snap_value, this->scalar, this->snap_back,
     octree);
}

// **************************************************
// SNAP INFO
// **************************************************

void IJKSNAPMC::SNAP_INFO::Clear()
{
  MCUBE_INFO::Clear();
  num_snapped_iso_vertices = 0;
}


