/// \file ijkmcube_util.cxx
/// Utility functionsfor ijkmcube

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

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ijkmcube_util.h"

#include "ijkcoord.txx"
#include "ijkinterpolate.txx"

using namespace IJK;
using namespace IJKTABLE;

// **************************************************
// UTILITY FUNCTIONS
// **************************************************

void IJKMCUBE::compute_iso_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of isosurface vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_iso_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BINARY);
  assert(isotable.NumIsosurfaceVertices() == isotable.Polyhedron().NumEdges());

  const ISO_VERTEX_INDEX num_iso_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_isov_per_gridv = 
    get_num_iso_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_iso_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() != ISOSURFACE_VERTEX::EDGE) {
      error.AddMessage("Illegal isosurface vertex ", j, ".");
      error.AddMessage("Isosurface vertex type must be EDGE.");
      throw error;
    };

    IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
    int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
    int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

    int d = 0;
    while (d < dimension && 
	   isotable.Polyhedron().VertexCoord(iv0, d) ==
	   isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

    // error checking
    if (d == dimension) {
      error.AddMessage("Illegal cube edge ", ie, ".");
      error.AddMessage("Endpoints ", iv0, " and ", iv1, 
		       " have the exact same coordinates.");
      throw error;
    };

    for (int d2 = 0; d2 < dimension; d2++) {
      if (d2 != d && 
	  isotable.Polyhedron().VertexCoord(iv0, d2) !=
	  isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " differ in more than one coordinate.");
	throw error;
      };
    }

    if (isotable.Polyhedron().VertexCoord(iv0, d) >
	isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

    increment[j] = vertex_increment[iv0]*num_isov_per_gridv + d;
  }
}

void IJKMCUBE::compute_iso_endpoint_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of endpoints of grid edges containing isosurface vertices
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // (increment[2*i], increment[2*i+1) = 
  //   increment for endpoints of grid edges containing iso vertex i
  // Precondition: array increment is allocated with size at least
  //               2*isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_iso_endpoint_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BINARY);
  assert(isotable.NumIsosurfaceVertices() == isotable.Polyhedron().NumEdges());

  const ISO_VERTEX_INDEX num_iso_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_isov_per_gridv = 
    get_num_iso_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_iso_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() != ISOSURFACE_VERTEX::EDGE) {
      error.AddMessage("Illegal isosurface vertex ", j, ".");
      error.AddMessage("Isosurface vertex type must be EDGE.");
      throw error;
    };

    IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
    int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
    int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

    increment[2*j] = vertex_increment[iv0];
    increment[2*j+1] = vertex_increment[iv1];

    if (increment[2*j] > increment[2*j+1])
      { std::swap(increment[2*j], increment[2*j+1]); }
  }
}

void IJKMCUBE::compute_nep_iso_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of isosurface vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_nep_iso_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.NumIsosurfaceVertices() == 
	 isotable.Polyhedron().NumEdges()+isotable.Polyhedron().NumVertices());

  const ISO_VERTEX_INDEX num_iso_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_isov_per_gridv = 
    get_num_nep_iso_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_iso_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() == ISOSURFACE_VERTEX::EDGE) {

      IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
      int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
      int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

      int d = 0;
      while (d < dimension && 
	     isotable.Polyhedron().VertexCoord(iv0, d) ==
	     isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

      // error checking
      if (d == dimension) {
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " have the exact same coordinates.");
	throw error;
      };

      for (int d2 = 0; d2 < dimension; d2++) {
	if (d2 != d && 
	    isotable.Polyhedron().VertexCoord(iv0, d2) !=
	    isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	  error.AddMessage("Illegal cube edge ", ie, ".");
	  error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			   " differ in more than one coordinate.");
	  throw error;
	};
      }

      if (isotable.Polyhedron().VertexCoord(iv0, d) >
	  isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

      increment[j] = vertex_increment[iv0]*num_isov_per_gridv + d;
    }
    else if (isotable.IsosurfaceVertex(j).Type() == 
	     ISOSURFACE_VERTEX::VERTEX) {
      VERTEX_INDEX iv = isotable.IsosurfaceVertex(j).Face();
      increment[j] = vertex_increment[iv]*num_isov_per_gridv + dimension;
    }
    else {
      error.AddMessage("Illegal isosurface vertex ", j, ".");
      error.AddMessage("Isosurface vertex type must be EDGE or VERTEX.");
      throw error;
    }
  }
}

void IJKMCUBE::compute_ivol_vertex_increment
(const ISOSURFACE_TABLE & isotable, const VERTEX_INDEX * vertex_increment,
 ISO_VERTEX_INDEX * increment)
  // compute increment to add to primary vertex to get index 
  //   of interval volume vertex i of hypercube
  // vertex_increment[k] = vertex increment for hypercube vertex k
  // increment[i] = increment for isosurface vertex i
  //              = index of edge containing vertex i
  // Precondition: array increment is allocated with size at least
  //               isotable.NumIsosurfaceVertices()
{
  PROCEDURE_ERROR error("IJKMCUBE::compute_ivol_vertex_increment");

  assert(vertex_increment != NULL);
  assert(increment != NULL);
  assert(isotable.Encoding() == ISOSURFACE_TABLE::BASE3);
  assert(isotable.NumIsosurfaceVertices() == 
	 2*isotable.Polyhedron().NumEdges()+
	 isotable.Polyhedron().NumVertices());

  const ISO_VERTEX_INDEX num_ivol_vertices = 
    isotable.NumIsosurfaceVertices();
  const int dimension = isotable.Dimension();
  const int num_ivolv_per_gridv = 
    get_num_ivol_vertices_per_grid_vertex(dimension);

  for (long j = 0; j < num_ivol_vertices; j++) {

    if (isotable.IsosurfaceVertex(j).Type() == ISOSURFACE_VERTEX::EDGE) {

      IJKTABLE::EDGE_INDEX ie = isotable.IsosurfaceVertex(j).Face();
      int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
      int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

      int d = 0;
      while (d < dimension && 
	     isotable.Polyhedron().VertexCoord(iv0, d) ==
	     isotable.Polyhedron().VertexCoord(iv1, d)) { d++; };

      // error checking
      if (d == dimension) {
	error.AddMessage("Illegal cube edge ", ie, ".");
	error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			 " have the exact same coordinates.");
	throw error;
      };

      for (int d2 = 0; d2 < dimension; d2++) {
	if (d2 != d && 
	    isotable.Polyhedron().VertexCoord(iv0, d2) !=
	    isotable.Polyhedron().VertexCoord(iv1, d2)) { 
	  error.AddMessage("Illegal cube edge ", ie, ".");
	  error.AddMessage("Endpoints ", iv0, " and ", iv1, 
			   " differ in more than one coordinate.");
	  throw error;
	};
      }

      if (isotable.Polyhedron().VertexCoord(iv0, d) >
	  isotable.Polyhedron().VertexCoord(iv1, d)) { std::swap(iv0, iv1); };

      if (isotable.IsosurfaceVertex(j).IsLabelSet()) {
	if (isotable.IsosurfaceVertex(j).Label() == "0") {
	  increment[j] = vertex_increment[iv0]*num_ivolv_per_gridv + 2*d;
	}
	else if (isotable.IsosurfaceVertex(j).Label() == "1") {
	  increment[j] = vertex_increment[iv0]*num_ivolv_per_gridv + 2*d+1;
	}
	else {
	  error.AddMessage("Interval volume vertex ", j, " has improper label.");
	  error.AddMessage("Interval volume vertices on grid edges must have labels 0 or 1.");
	  throw error;
	};
      }
      else {
	error.AddMessage("Interval volume vertex ", j, " has no label.");
	error.AddMessage("Interval volume vertices on grid edges must have labels 0 or 1.");
	throw error;
      };

    }
    else if (isotable.IsosurfaceVertex(j).Type() == 
	     ISOSURFACE_VERTEX::VERTEX) {
      VERTEX_INDEX iv = isotable.IsosurfaceVertex(j).Face();
      increment[j] = vertex_increment[iv]*num_ivolv_per_gridv + 2*dimension;
    } else {
      error.AddMessage("Interval volume vertex ", j, " has improper type.");
      error.AddMessage("Interval volume vertices must lie on grid vertices or grid edges.");
      throw error;
    };

  }

}

void IJKMCUBE::compute_hypercube_edge_increment
(const int dimension, const ISOSURFACE_TABLE & isotable, 
 const VERTEX_INDEX * hcube_vertex_increment, EDGE_INDEX * increment)
// compute increment from first edge in hypercube
// increment[i] = increment of i'th edge
//   Precondition: array increment is allocated with size at least
//                 isotable.Polyhedron().NumEdges()
{
  const int nume = isotable.Polyhedron().NumEdges();
  PROCEDURE_ERROR error("ijkmcube");

  assert(increment != NULL);

  for (int ie = 0; ie < nume; ie++) {
    int iv0 = isotable.Polyhedron().EdgeEndpoint(ie, 0);
    int iv1 = isotable.Polyhedron().EdgeEndpoint(ie, 1);

    if (iv0 > iv1) {
      std::swap(iv0, iv1);
    };

    bool edge_found = false;
    for (int d = 0; d < dimension; d++) {
      int jv1 = iv0 + (1L << d);
      if (jv1 == iv1) {
	increment[ie] = hcube_vertex_increment[iv0]*dimension+d;
	edge_found = true;
	break;
      }
    }

    if (!edge_found) {
      error.AddMessage("Error computing hypercube edge increments.");
      error.AddMessage
	("Hypercube vertices may not be in lexicographic order.");
      throw error;
    };
  }

}

void IJKMCUBE::compute_facet_increment
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const AXIS_SIZE_TYPE * axis_size,
 VERTEX_INDEX * facet_increment)
  // compute increment to add to index of current voxel to get
  //   facet index
  // Precondition: Polyhedron is a cube.
  // Precondition: array increment is allocated with size 
  //               at least poly.NumFacets();
{
  const int dimension = poly.Dimension();
  const int numf = poly.NumFacets();
  const int numv = poly.NumVertices();
  VERTEX_INDEX axis_increment[dimension];
  int orthogonal_axis[numf];
  int axis_intercept[numf];
  int num_orthogonal[dimension];  // number of facets orthogonal to dimension
  int facet0[dimension];
  int facet1[dimension];


  PROCEDURE_ERROR error("compute_facet_increment");
  error.AddMessage("Programming error.");

  assert(numf == 2*dimension);
  assert(numv > 0 && numv == (1L << dimension));

  compute_increment(dimension, axis_size, axis_increment);

  // find orthogonal axis
  for (int jf = 0; jf < poly.NumFacets(); jf++) {
    bool is_orthogonal_axis_set = false;
    for (int d = 0; d < dimension; d++) {
      int coord;
      bool coord_is_set = false;
      bool is_orthogonal = true;
      for (int iv = 0; iv < numv; iv++) {
	if (poly.IsVertexInFacet(jf, iv)) {
	  if (coord_is_set) {
	    if (poly.VertexCoord(iv, d) != coord) {
	      is_orthogonal = false;
	      break;
	    }
	  }
	  else {
	    coord = poly.VertexCoord(iv, d); 
	    coord_is_set = true;
	  }
	}
      }

      if (coord_is_set && is_orthogonal) {
	orthogonal_axis[jf] = d;
	axis_intercept[jf] = coord;
	is_orthogonal_axis_set = true;
	break;
      }

      if (!coord_is_set) {
	error.AddMessage("No vertices lie on polyhedron facet ", jf, ".");
	throw error;
      }
    }

    if (!is_orthogonal_axis_set) {
      error.AddMessage("Unable to find axis orthogonal to polyhedron facet ", 
		       jf, ".");
      throw error;
    }
  }

  // initialize num_orthogonal[d] to 0
  for (int d = 0; d < dimension; d++) {
    num_orthogonal[d] = 0;
  }

  // set facet0 and facet1
  for (int jf = 0; jf < poly.NumFacets(); jf++) {
    int d = orthogonal_axis[jf];
    if (num_orthogonal[d] == 0) {
      facet0[d] = jf;
    }
    else if (num_orthogonal[d] == 1) {
      facet1[d] = jf;
    }
    else {
      error.AddMessage("More than two facets are orthogonal to axis ", 
		       d, ".");
      throw error;
    }

    num_orthogonal[d]++;
  };

  for (int d = 0; d < dimension; d++) {

    if (num_orthogonal[d] != 2) {
      error.AddMessage("Axis ", d, " does not have two orthogonal facets.");
      throw error;
    }

    int jf0 = facet0[d];
    int jf1 = facet1[d];
    if (axis_intercept[jf0] > axis_intercept[jf1])
      { std::swap(jf0, jf1); };

    facet_increment[jf0] = d;
    facet_increment[jf1] = dimension*axis_increment[d] + d;
  }

}

namespace {

  void get_min_max_coord
  (const ISOSURFACE_TABLE_POLYHEDRON & poly, const int d,
   int & minc, int & maxc)
  {
    minc = maxc = 0;
    if (poly.NumVertices() < 1) { return; };

    minc = poly.VertexCoord(0, d);
    maxc = poly.VertexCoord(0, d);
    for (int iv = 1; iv < poly.NumVertices(); iv++)
      {
	int coord = poly.VertexCoord(iv, d);
	if (coord < minc) { minc = coord; };
	if (coord > maxc) { maxc = coord; };
      }
  }

}

void IJKMCUBE::compute_facet_vertex_increment
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const int orth_dir, const int side,
 const AXIS_SIZE_TYPE * axis_size, 
 VERTEX_INDEX * facet_vertex_increment, VERTEX_INDEX * cube_vertex_index)
// compute increment to add to primary facet vertex
//   to get other facet vertices
// poly = isosurface table polyhedron.
//   Precondition: Polyhedron is a cube.
// orth_dir = direction orthogonal to facet. 
// side = 0. Facet orthogonal to orth_dir with min orth_dir coordinates.
// side = 1. Facet orthogonal to orth_dir with max orth_dir coordinates.
// facet_vertex_increment[k]: increment to add to primary vertex
//               to compute k'th facet vertex
// Precondition: facet_vertex increment is allocated with size 
//               num_facet_vertices
// cube_vertex_index[k]: Index of k'th facet vertex in cube
//   Precondition: cube_vertex_index is preallocated with size 
//               num_facet_vertices
{
  const int dimension = poly.Dimension();
  VERTEX_INDEX axis_increment[dimension];
  int minc[dimension];
  int maxc[dimension];
  PROCEDURE_ERROR error("compute_facet_vertex_increment");

  if (dimension == 0) { return; };

  assert(0 <= orth_dir && orth_dir < dimension);

  const GRID_SIZE_TYPE num_cube_facet_vertices = 
    compute_num_cube_facet_vertices(dimension);
  compute_increment(dimension, axis_size, axis_increment);

  for (int d = 0; d < dimension; d++) 
    { get_min_max_coord(poly, d, minc[d], maxc[d]); }

  int facet_coord = minc[orth_dir];
  if (side == 0)
    { facet_coord = maxc[orth_dir]; };

  int k = 0;
  for (int iv = 0; iv < poly.NumVertices(); iv++) {
    int coord0 = poly.VertexCoord(iv, orth_dir);
    if (coord0 == facet_coord) {
      facet_vertex_increment[k] = 0;
      cube_vertex_index[k] = iv;
      for (int d = 0; d < dimension; d++) {
	if (d != orth_dir) {
	  int coord1 = poly.VertexCoord(iv, d);
	  if (coord1 == maxc[d])
	    { facet_vertex_increment[k] += axis_increment[d]; }
	}
      }
      k++;
    }
  }

  if (k != num_cube_facet_vertices) {
    error.AddMessage("Programming error.  Failed to process correct number of facet vertices.");
    error.AddMessage("  Processed ", k, " vertices.  Number of cube facet vertices = ", num_cube_facet_vertices, ".");
    throw error;
  }
}


/// Get mixed cubes (cubes with positive and negative vertices.)
/// Cubes are represented by primary vertices (lowest coord vertices in cubes)
/// @param minmax = minmax data structures
/// @param[out] vlist[] = list of primary cube vertices.
/// @param[out] vlist_length = Number of vertices in vlist.
/// @pre Array vlist[] is preallocated to length at least number of cubes.
void IJKMCUBE::get_mixed_cubes
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MINMAX_REGIONS & minmax,  const SCALAR_TYPE isovalue, 
 VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length)
{
  assert(axis_size != NULL && vlist != NULL);

  const AXIS_SIZE_TYPE region_edge_length = minmax.RegionEdgeLength();
  VERTEX_INDEX max_num_cubes_in_region;
  compute_num_grid_cubes_in_region
    (dimension, region_edge_length, max_num_cubes_in_region);

  VERTEX_INDEX num_regions;
  compute_num_regions(dimension, axis_size, region_edge_length, num_regions);

  GRID_SIZE_TYPE num_cubes;
  compute_num_grid_cubes(dimension, axis_size, num_cubes);

  vlist_length = 0;
  IJK::ARRAY<VERTEX_INDEX> rlist(num_regions);
  IJK::ARRAY<bool> is_full(num_regions);
  IJK::ARRAY<VERTEX_INDEX> region_cube_increment(max_num_cubes_in_region);

  compute_region_cube_increment
    (dimension, axis_size, region_edge_length, region_cube_increment.Ptr());

  get_region_primary_vertices
    (dimension, axis_size, region_edge_length, rlist.Ptr(), is_full.Ptr());

  VERTEX_INDEX num_cubes_in_region;
  for (VERTEX_INDEX iregion = 0; iregion < num_regions; iregion++) {

    if (minmax.Min(iregion) < isovalue &&
	minmax.Max(iregion) >= isovalue) {

      VERTEX_INDEX iv0 = rlist[iregion];
      if (is_full[iregion]) {
	for (VERTEX_INDEX j = 0; j < max_num_cubes_in_region; j++) {
	  vlist[vlist_length] = iv0 + region_cube_increment[j];
	  vlist_length++;
	}
      }
      else {
	get_grid_cubes_in_region
	  (dimension, axis_size, iv0, region_edge_length,
	   vlist+vlist_length, num_cubes_in_region);

	vlist_length += num_cubes_in_region;
      }
    }

  }

  if (num_cubes < vlist_length) {
    PROCEDURE_ERROR error("get_mixed_cubes");
    error.AddMessage("Programming error.");
    error.AddMessage
      ("Added to vlist more primary cube vertices than exist in grid.");
    throw error;
  }
}


void IJKMCUBE::get_mixed_cubes
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const IJKOCTREE::OCTREE & octree,  const SCALAR_TYPE isovalue, 
 VERTEX_INDEX * vlist, VERTEX_INDEX & vlist_length)
// get mixed cubes (cubes with positive and negative vertices)
// cubes are represented by primary vertices (lowest coord vertices in cubes)
// minmax = minmax data structures
// vlist[k] = list of primary cube vertices
//    Precondition: vlist[] is preallocated to length at least 
//                  number of cubes
// vlist_length = number of vertices in vlist
{
  assert(axis_size != NULL && vlist != NULL);

  GRID_SIZE_TYPE num_cubes;
  compute_num_grid_cubes(dimension, axis_size, num_cubes);

  vlist_length = 0;

  int num_levels = octree.NumLevels();
  if (num_levels == 0) return;
  int ileaf_level = num_levels-1;

  // Depth first search through octree
  IJKOCTREE::OCTREE_STACK stack(num_levels);

  IJKOCTREE::CONST_NODE_PTR root = octree.Root();
  stack.PushRoot(root);

  while (!stack.IsEmpty()) {

    IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
    if (node->MinValue() < isovalue && isovalue <= node->MaxValue()) {

      if (stack.TopIsLeaf()) {
	// leaf node
	IJKOCTREE::CONST_NODE_PTR node = stack.TopNode();
	for (int i = 0; i < node->NumChildren(); i++) {
	  VERTEX_INDEX vertex0 = 
	    octree.ComputeChildVertex0(ileaf_level, node, i);
	  vlist[vlist_length] = vertex0;
	  vlist_length++;
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


  if (num_cubes < vlist_length) {
    PROCEDURE_ERROR error("get_mixed_cubes");
    error.AddMessage("Programming error.");
    error.AddMessage
      ("Procedure get_mixed_cubes returns more cubes from minmax than exist in grid.");
  }
}


std::string IJKMCUBE::get_isotable_filename
(const ISOTABLE_TYPE isotable_type, const int dimension,
 const char * poly_name)
// return isotable file name base on isotable_type and dimension
{
  PROCEDURE_ERROR error("get_isotable_filename");
  std::string isotable_filename;

  switch(isotable_type){

  case BINARY:
    isotable_filename = "iso";
    break;

  case NEP:
    isotable_filename = "iso.nep";
    break;

  case IVOL:
    isotable_filename = "ivol";
    break;

  default:
    error.AddMessage("Programming error. Illegal isosurface table type.");
    throw error;
  };

  isotable_filename += "." + std::string(poly_name);

  std::ostringstream dimension_stringstream;
  dimension_stringstream << dimension;

  isotable_filename += "." + dimension_stringstream.str() + "D.xit";

  return(isotable_filename);
}

std::string IJKMCUBE::get_isotable_filename
(const ISOTABLE_TYPE isotable_type, const int dimension)
// return isotable file name base on isotable_type and dimension
{
  return(get_isotable_filename(isotable_type, dimension, "cube"));
}

// **************************************************
// STRINGS
// **************************************************

std::string IJKMCUBE::get_topology_string
(const ISOSURFACE_TOPOLOGY isosurface_topology)
{
  using std::string;

  string topology_str;

  switch (isosurface_topology) {

  case ISOTABLE_TOPOLOGY:
    topology_str = string("isosurface table");
    break;

  case CUBE_DECIDER_TOPOLOGY:
    topology_str = string("cube decider");
    break;

  case ASYMPTOTIC_DECIDER_TOPOLOGY:
    topology_str = string("asymptotic decider");
    break;

  case SADDLE_TOPOLOGY:
    topology_str = string("saddle topology");
    break;

  case LINEAR_TOPOLOGY:
    topology_str = string("linear topology");
    break;

  default:
    topology_str = string("unknown topology");
    break;
  }

  return(topology_str);
}

// **************************************************
// CHECK SUBROUTINES
// **************************************************

bool IJKMCUBE::check_nep_num_dup(const int nep_num_dup, ERROR & error)
{
  if (nep_num_dup < 0 || nep_num_dup > 2) {
    error.AddMessage("Programming error.  Illegal value ", nep_num_dup,
		     " for nep_num_dup.");
    error.AddMessage("   Value should be 0, 1 or 2.");

    return(false);
  }

  return(true);
}

bool IJKMCUBE::check_order(const std::vector<VERTEX_INDEX> & endpoint,
			   ERROR & error)
{
  const int nume = endpoint.size()/2;
  for (int i = 0; i < nume; i++) {
    if (endpoint[2*i] > endpoint[2*i+1]) {
      error.AddMessage("Programming error.  Illegal endpoint order.");
      error.AddMessage("Edge ", i, " has endpoints (",
		       endpoint[2*i], ",", endpoint[2*i+1], ")");
      error.AddMessage("Endpoints should be listed in reverse order.");
      return(false);
    }
  }

  return(true);
}

bool IJKMCUBE::check_isotable_encoding
(const ISOSURFACE_TABLE & isotable, 
 const ISOSURFACE_TABLE::ENCODING encoding, ERROR & error)
{
  if (isotable.Encoding() != encoding) {
    error.AddMessage("Illegal isosurface table encoding.");
    if (encoding == ISOSURFACE_TABLE::BINARY)
      { error.AddMessage("  BINARY encoding required."); }
    else if (encoding == ISOSURFACE_TABLE::BASE3)
      { error.AddMessage("  BASE3 encoding required."); }

    return(false);
  }

  return(true);
}

bool IJKMCUBE::check_isotable_encoding
(const POLY_ISOTABLE & poly_isotable, 
 const ISOSURFACE_TABLE::ENCODING encoding, ERROR & error)
{
  if (poly_isotable.cube.IsTableAllocated()) {
    if (!check_isotable_encoding(poly_isotable.cube, encoding, error))
      { return(false); };
  }

  if (poly_isotable.pyramid.IsTableAllocated()) {
    if (!check_isotable_encoding(poly_isotable.pyramid, encoding, error))
      { return(false); };
  }

  if (poly_isotable.simplex.IsTableAllocated()) {
    if (!check_isotable_encoding(poly_isotable.simplex, encoding, error))
      { return(false); };
  }

  return(true);
}

bool IJKMCUBE::check_isotable_ambig_info
(const ISOSURFACE_TABLE & isotable, 
 const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, ERROR & error)
{
  if (isotable.NumTableEntries() != ambig_info.NumTableEntries()) {
    if (ambig_info.NumTableEntries() == 0) {
      error.AddMessage("Isosurface table ambiguity information has not been computed.");
    }
    else {
      error.AddMessage("Isosurface table ambiguity information does not match isosurface table.");
      error.AddMessage("  Isosurface table ambiguity information has ",
		       ambig_info.NumTableEntries(), " table entries.");
      error.AddMessage("  Isosurface table has ",
		       isotable.NumTableEntries(), " table entries.");
    }
    return(false);
  }

  return(true);
}

/// Check cube isosurface lookup table fits topology
bool IJKMCUBE::check_cube_isotable_fits_topology
(const ISOSURFACE_TABLE & cube_isotable, 
 ISOSURFACE_TOPOLOGY isosurface_topology, IJK::ERROR & error)
{
  const int dimension = cube_isotable.Dimension();
  const VERTEX_INDEX numv = 
    compute_num_cube_vertices<VERTEX_INDEX>(dimension);

  if (numv != cube_isotable.Polyhedron().NumVertices()) {
    error.AddMessage("Isotable polyhedron is not a cube.");
    return(false);
  }

  if (isosurface_topology == ISOTABLE_TOPOLOGY) {
    // Nothing to check.
    return(true);
  };

  std::string topology_str = get_topology_string(isosurface_topology);

  for (TABLE_INDEX ifirst = 0; ifirst < numv/2; ifirst++) {
    TABLE_INDEX ilast = (numv-ifirst-1);
    TABLE_INDEX i0 = (TABLE_INDEX(1) << ifirst);
    TABLE_INDEX i1 = (TABLE_INDEX(1) << ilast);

    TABLE_INDEX it = (i0 | i1);
    it = cube_isotable.NumTableEntries()-1-it;

    if (cube_isotable.NumSimplices(it) != 2) {
      std::string msg = "Isosurface table may not be appropriate for " + topology_str + ".";
      error.AddMessage(msg);
      error.AddMessage("  Table entry ", it, " has ",
		       cube_isotable.NumSimplices(it), " simplices.");
      error.AddMessage("  Expected two simplices.");
      return(false);
    }
  }

  return(true);
}

/// Check isosurface table dimension matches scalar field dimension
bool IJKMCUBE::check_dimension
(const ISOSURFACE_TABLE & isotable, const MC_SCALAR_GRID_BASE & scalar_grid,
 IJK::ERROR & error)
{
  if (isotable.Dimension() != scalar_grid.Dimension()) {
    error.AddMessage("Isosurface table dimension does not match scalar field dimension.");
    error.AddMessage("Isosurface table dimension = ", isotable.Dimension(), 
		     ".");
    error.AddMessage("Scalar grid dimension = ", scalar_grid.Dimension(), 
		     ".");
    return(false);
  }

  return(true);
}
