/// \file ijkmcube_sub.cxx
/// Subroutines for ijkmcube

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006,2007,2008,2009 Rephael Wenger

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

#include "ijkmcube_sub.h"
#include "ijkmcube_util.h"
#include "ijkoctree.h"
#include "ijktable.h"
//#include <Algorithm>

#include "ijkcoord.txx"
#include "ijkinterpolate.txx"
#include "ijkpoly.txx"
#include "ijkmerge.txx"
#include "ijkisopoly.txx"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKTABLE;
using namespace std;

// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
// call merge_identical_vertices and then position_vertices_linear
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// simplex vert[] = list of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of simplex is
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");
  

  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical_vertices(dimension, axis_size, iso_simplices, 
			   iso_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linear
    (scalar_grid, isovalue, isotable_type, iso_vlist, &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}

void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const std::vector<VERTEX_INDEX> & endpoint,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_edges_parameter, 
 MCUBE_INFO & mcube_info)
// call merge_identical_edges and then position_iso_vertices_linearB
// Input includes an array of edge endpoints.  Isosurface vertices lie on these edges.
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// endpoint[] = list of isosurface simplex vertices
//   (endpoint[2*i], endpoint[2*i+1]) = endoints
//    of edge containing isosurface vertex i.
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;

  merge_identical_edges
    (endpoint, merged_endpoint, iso_simplices, merge_edges_parameter);

  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linearB
    (scalar_grid, isovalue, isotable_type, merged_endpoint, 
     &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & endpoint,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
// call merge_identical_edges and then position_iso_vertices_linearB
// Input includes an array of edge endpoints.  Isosurface vertices lie on these edges.
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// endpoint[] = list of isosurface simplex vertices
//   (endpoint[2*i], endpoint[2*i+1]) = endoints
//    of edge containing isosurface vertex i.
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;
  merge_identical_edges
    (endpoint, merged_endpoint, iso_simplices, merge_parameters);
			
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linearB
    (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
     merged_endpoint, &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}

void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & endpoint,
 const ISOTABLE_TYPE isotable_type,
 const INTERPOLATION_TYPE interpolation_type,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
// call merge_identical_edges and then position_iso_vertices_linearB
// Input includes an array of edge endpoints.  Isosurface vertices lie on these edges.
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// endpoint[] = list of isosurface simplex vertices
//   (endpoint[2*i], endpoint[2*i+1]) = endoints
//    of edge containing isosurface vertex i.
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;
  merge_identical_edges
    (endpoint, merged_endpoint, iso_simplices, merge_parameters);
			
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  switch(interpolation_type) {

  case LINEAR_INTERPOLATION:
    position_iso_vertices_linearB
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  case MULTILINEAR_INTERPOLATION:
    position_iso_vertices_multilinear
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  case ALPHA_UNCERTAINTY:
    cout<<"Here in alpha uncertainty!"; 

    /* This was called before
      position_iso_vertices_linearB_alpha
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));    */
    break;	


  default:
    error.AddMessage("Illegal interpolation type.");
    throw error;
  }

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


// newly added!!!
// Take mu (scalar_grid) and delta (scalar_grid_1) grids as input and compute geometry and its variance
void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & endpoint,
 const ISOTABLE_TYPE isotable_type,
 const INTERPOLATION_TYPE interpolation_type,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info,const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color)
// call merge_identical_edges and then position_iso_vertices_linearB
// Input includes an array of edge endpoints.  Isosurface vertices lie on these edges.
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// endpoint[] = list of isosurface simplex vertices
//   (endpoint[2*i], endpoint[2*i+1]) = endoints
//    of edge containing isosurface vertex i.
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
{
   
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;
  merge_identical_edges
    (endpoint, merged_endpoint, iso_simplices, merge_parameters);
 
			
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  switch(interpolation_type) {

  case LINEAR_INTERPOLATION: 
    position_iso_vertices_linearB
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  case MULTILINEAR_INTERPOLATION:
    position_iso_vertices_multilinear
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  case ALPHA_UNCERTAINTY:
    position_iso_vertices_linearB_alpha
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]),scalar_grid_1, f_color);    
    break;	


  default:
    error.AddMessage("Illegal interpolation type.");
    throw error;
  }

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}






void IJKMCUBE::merge_identical_vertices
(const int dimension, const AXIS_SIZE_TYPE * grid_length,
 const std::vector<ISO_VERTEX_INDEX> & vlist0,
 std::vector<ISO_VERTEX_INDEX> & vlist1, 
 std::vector<ISO_VERTEX_INDEX> & vlist0_map, MERGE_DATA & merge_data)
// merge identical vertices in list of isosurface vertices
// return vlist1 and vlist0_map
// dimension = volume dimension
// vlist0[] = input list of isosurface vertices
//   Precondition: Vertex indices are in range [0,merged_data.NumEdges()-1]
// vlist1[] = output list of isosurface vertices
// vlist0_map[] = vlist0[iv] is mapped to vlist1[vlist0_map[iv]]
// merge_data = internal merge data structure
{
  PROCEDURE_ERROR error("merge_identical_vertices");

  if (!merge_data.Check(error)) { throw error; }

  // initialize
  merge_data.ClearList();
  vlist1.clear();
  vlist0_map.clear();

  ISO_VERTEX_INDEX v;
  for (int i = 0; i < vlist0.size(); i++) {
    v = vlist0[i];

    if (!merge_data.InList(v)) {
      vlist1.push_back(v);
      merge_data.Insert(v);
    };
    vlist0_map.push_back(merge_data.ListLoc(v));
  };
}

void IJKMCUBE::merge_identical_edges
(const std::vector<VERTEX_INDEX> & elist0,
 std::vector<VERTEX_INDEX> & elist1, 
 std::vector<ISO_VERTEX_INDEX> & eindex,
 const MERGE_EDGES_PARAMETERS & merge_edges_parameter)
// merge identical vertices in list of isosurface vertices
// return elist1, modifies eindex
{
  const int nume = elist0.size()/2;
  PROCEDURE_ERROR error("merge_identical_edges");

  if (elist0.size() == 0 && eindex.size() > 0) {
    error.AddMessage("Programming error.  elist0 has zero edges but eindex is not empty.");
    throw error;
  }

  if (!check_order(elist0, error)) { throw error; };

  // initialize
  elist1.clear();

  ARRAY<ISO_VERTEX_INDEX> elist0_map(nume);

  merge_pairs(elist0, elist1, elist0_map.Ptr(), 
	      merge_edges_parameter);

  if (!check_pair_list
      (elist0, elist1, elist0_map.PtrConst(), error)) 
    { throw error; }

  for (int i = 0; i < eindex.size(); i++) 
    { eindex[i] = elist0_map[eindex[i]]; }
}

void IJKMCUBE::set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable)
{
  const int numf = isotable.Polyhedron().NumFacets();
  const int numv_per_simplex = isotable.NumVerticesPerSimplex();

  isotable.SetIsInFacet(false);

  for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
    int nums = isotable.NumSimplices(it);
    if (nums == 0)   // no isosurface patch for table entry it
      { continue; }

    for (int jf = 0; jf < numf; jf++) {
      bool in_facet = true;
      for (int is = 0; is < nums && in_facet; is++) {
	for (int k = 0; k < numv_per_simplex && in_facet; k++) {
	  ISOSURFACE_VERTEX_INDEX isov = isotable.SimplexVertex(it, is, k);
	  if (isotable.IsosurfaceVertex(isov).Type() ==
	      ISOSURFACE_VERTEX::VERTEX) {
	    int iv = isotable.IsosurfaceVertex(isov).Face();
	    if (!isotable.Polyhedron().IsVertexInFacet(jf, iv)) 
	      { in_facet = false; };
	  }
	  else {
	    in_facet = false;
	  }
	}
      }

      if (in_facet) {
	isotable.SetContainingFacet(it, jf);
	break;
      }
    }
  }

}


// **************************************************
// POSITION ISO VERTICES
// **************************************************

namespace {

  /// Get scalar values of polyhedron vertices.
  /// *** MAKE THIS A TEMPLATE ***
  inline void get_poly_scalar
  (const SCALAR_TYPE * scalar, const VERTEX_INDEX iv0,
   const VERTEX_INDEX * increment, const int num_vertices, 
   SCALAR_TYPE * poly_scalar)
  {
    for (int i = 0; i < num_vertices; i++) {
      VERTEX_INDEX iv1 = iv0 + increment[i];
      poly_scalar[i] = scalar[iv1];
    }
  }

  void interpolate_using_multilinear
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const VERTEX_INDEX * increment, const int num_cube_vertices,
   const COORD_TYPE * coord0, const COORD_TYPE * coord1,
   COORD_TYPE * coord2)
  {
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv_base;
    GRID_COORD_TYPE base_coord[dimension];
    COORD_TYPE coordA[dimension];
    COORD_TYPE coordB[dimension];
    COORD_TYPE coordC[dimension];
    SCALAR_TYPE cube_scalar[num_cube_vertices];

    for (int d = 0; d < dimension; d++) {
      COORD_TYPE x = std::min(coord0[d], coord1[d]);
      base_coord[d] = GRID_COORD_TYPE(floor(x));
    }

    iv_base = scalar_grid.ComputeVertexIndex(base_coord);

    get_poly_scalar(scalar_grid.ScalarPtrConst(), iv_base, increment,
		    num_cube_vertices, cube_scalar);


    subtract_coord(dimension, coord0, base_coord, coordA);
    subtract_coord(dimension, coord1, base_coord, coordB);

    const int NUM_ITER = 5;
    multilinear_interpolate_coord
    (dimension, coordA, coordB, isovalue, num_cube_vertices,
     cube_scalar, NUM_ITER, coordC);

    add_coord(dimension, coordC, base_coord, coord2);
  }

  void position_iso_vertices_linear_binary
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
  // compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // vlist[] = list of isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    AXIS_SIZE_TYPE axis_increment[dimension];
    const int numv = vlist.size();
    GRID_COORD_TYPE coord0[dimension];
    GRID_COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];

    // compute_axis_increment
    compute_increment(scalar_grid, axis_increment);

    const int num_isov_per_gridv = 
      get_num_iso_vertices_per_grid_vertex(dimension);


    for (int i = 0; i < numv; i++) {
      int d = vlist[i]%num_isov_per_gridv;
      VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;
      VERTEX_INDEX v1 = v0+axis_increment[d];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue, coord2);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

  void position_iso_vertices_linear_nep
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
  // compute position of isosurface vertices using linear interpolation
  //   vertices can have positive, negative or equals values
  // dimension = volume dimension
  // scalar_grid[iv] = scalar value at grid vertex iv
  // grid_length[i] = # grid vertices along axis i
  // vlist[] = list of isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
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

    // compute_axis_increment
    compute_increment(scalar_grid, axis_increment);

    const int num_isov_per_gridv = 
      get_num_nep_iso_vertices_per_grid_vertex(dimension);

    for (int i = 0; i < numv; i++) {
      int k = vlist[i]%num_isov_per_gridv;
      VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;

      if (k == VERTEX_OFFSET) {
	// isosurface vertex lies on grid vertex v0
	compute_coord(v0, dimension, axis_size, coord0);
	for (int d = 0; d < dimension; d++)
	  coord[i*dimension+d] = coord0[d];
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

  void position_iso_vertices_linear_binaryB
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
  // compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // elist[] = list of edges containing isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist.size()/2;
    GRID_COORD_TYPE coord0[dimension];
    GRID_COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist[2*i];
      VERTEX_INDEX v1 = elist[2*i+1];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      scalar_grid.ComputeCoord(v0, coord0);
      scalar_grid.ComputeCoord(v1, coord1);

      linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue, coord2);
      for (int d = 0; d < dimension; d++)
	{ coord[i*dimension+d] = coord2[d]; }
    }
  }


  void position_iso_vertices_linear_binaryB
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const MC_MESH_VERTEX_LIST & new_mesh_vertices,
   const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
  // compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // elist[] = list of edges containing isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist.size()/2;
    COORD_TYPE coord0[dimension];
    COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];
    VERTEX_INDEX axis_increment[dimension];

    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VERTEX_INDEX> increment(num_cube_vertices);
    compute_increment(scalar_grid, axis_increment);
    compute_cube_vertex_increment
      (dimension, axis_increment, increment.Ptr());

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist[2*i];
      VERTEX_INDEX v1 = elist[2*i+1];

      SCALAR_TYPE s0, s1;

      if (v0 < scalar_grid.NumVertices()) {
	s0 = scalar[v0];
	scalar_grid.ComputeCoord(v0, coord0);
      }
      else {
	VERTEX_INDEX k = v0 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord0);
	s0 = new_mesh_vertices.Scalar(k);
      }

      if (v1 < scalar_grid.NumVertices()) {
	s1 = scalar[v1];
	scalar_grid.ComputeCoord(v1, coord1);
      }
      else {
	VERTEX_INDEX k = v1 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord1);
	s1 = new_mesh_vertices.Scalar(k);
      }

      linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue, coord2);
      for (int d = 0; d < dimension; d++)
	{ coord[i*dimension+d] = coord2[d]; }
    }
  }


   void position_iso_vertices_linear_binaryB_alpha
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color)
  // compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // elist[] = list of edges containing isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist.size()/2;
    GRID_COORD_TYPE coord0[dimension];
    GRID_COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist[2*i];
      VERTEX_INDEX v1 = elist[2*i+1];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      scalar_grid.ComputeCoord(v0, coord0);
      scalar_grid.ComputeCoord(v1, coord1);

    // This was called!  alpha_uncertainty_coord(dimension, s0, coord0, s1, coord1, isovalue, coord2);
      for (int d = 0; d < dimension; d++)
	{ coord[i*dimension+d] = coord2[d]; }
    }
  }

  void getHistogram(float* area_hist, float* input, int input_length, int bin_count, int* hist_flag, float* mnarea, float* mxarea, float* expected_value, float* mode, float* mode_low, float* mode_high)
{		
	//find minimum of input
	float min_input = input[0];
	for(int i=1; i<input_length; i++)
	{	
		if(input[i] < min_input)
			min_input = input[i];				
        }		

	//find maximum cell area
	float max_input = input[0];
	for(int i=1; i<input_length; i++)
	{	
		if(input[i] > max_input)
			max_input = input[i];				
        }	
		
	if(max_input > min_input)
	{
		*mnarea = min_input;
		*mxarea = max_input;
		*hist_flag = 1;		
		//Create array of size of number of bins to store frequencies			
		for(int i=0; i<bin_count; i++)
		{
			area_hist[i] = 0;
		}
		int bin_idx;
		float bin_width = (float)((max_input - min_input)/bin_count);			
		for(int i=0; i<input_length; i++)
		{						
			bin_idx = (int)((input[i] - min_input)/bin_width);				
			if(bin_idx == bin_count)
				bin_idx = bin_idx - 1;
			area_hist[bin_idx] = area_hist[bin_idx] + 1;				
		}
		int SumCount = 0;
		int md = 0;		
		int max_freq = 0;	
		for(int i=0; i<bin_count; i++)
		{
			if (area_hist[i] > max_freq)
			{
				md = i;
				max_freq = area_hist[i];
			}
			SumCount += area_hist[i];				
		}		
		*mode_low = bin_width*md;
		*mode_high = bin_width*(md+1);  
		
		for(int i=0; i<bin_count; i++)
		{
			area_hist[i] = (float)(area_hist[i]/SumCount);				
		}
		*expected_value = 0;			
		for(int i=0; i<bin_count; i++)
		{		
			float ar = min_input +  i*bin_width + (float)(bin_width/2);
			*expected_value += ar*area_hist[i];
		}

	}
	else
	{
		*hist_flag = 0;	
		*mnarea = min_input;
		*expected_value = min_input;
	}
}

void equalizeHistogram(float* intensity_levels,int* count_intensity,int track_level,float* equalized_intensity)
{


}

  
  // Find expeced level-crossing location and its variance (with proportional colors)
  // perform histogram equalization of colors

  void position_iso_vertices_linear_binaryB_alpha
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const MC_MESH_VERTEX_LIST & new_mesh_vertices,
   const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color)
  // compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // elist[] = list of edges containing isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    // mean grid
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist.size()/2;
    COORD_TYPE coord0[dimension];
    COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];
    VERTEX_INDEX axis_increment[dimension];

    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VERTEX_INDEX> increment(num_cube_vertices);
    compute_increment(scalar_grid, axis_increment);
    compute_cube_vertex_increment
      (dimension, axis_increment, increment.Ptr());

    // newly added!!
    // delta grid
    const int dimension_1 = scalar_grid_1.Dimension();
    const SCALAR_TYPE * scalar_1 = scalar_grid_1.ScalarPtrConst();
    f_color = new float[4*nume];
    float normalized_variance[nume];
    float copy_normalized_variance[nume];	
    float intensity_levels[nume]; 	
    int count_intensity[nume];	

    float* max_variance = new float; 
    float* color = new float;

    *max_variance = 0;

    // Loop through all edges
    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist[2*i];
      VERTEX_INDEX v1 = elist[2*i+1];

      // mean
      SCALAR_TYPE s0, s1;
 
      // delta
      SCALAR_TYPE s2, s3;  	

      // Read mean values at edge vertices
      if (v0 < scalar_grid.NumVertices()) {
	s0 = scalar[v0];
	scalar_grid.ComputeCoord(v0, coord0);
      }
      else {
	VERTEX_INDEX k = v0 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord0);
	s0 = new_mesh_vertices.Scalar(k);
      }

      if (v1 < scalar_grid.NumVertices()) {
	s1 = scalar[v1];
	scalar_grid.ComputeCoord(v1, coord1);
      }
      else {
	VERTEX_INDEX k = v1 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord1);
	s1 = new_mesh_vertices.Scalar(k);
      }

      // Read delta values at edge vertices
      if (v0 < scalar_grid_1.NumVertices()) {
	s2 = scalar_1[v0];
	scalar_grid_1.ComputeCoord(v0, coord0);
      }
      else {
	VERTEX_INDEX k = v0 - scalar_grid_1.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord0);
	s2 = new_mesh_vertices.Scalar(k);
      }

      if (v1 < scalar_grid_1.NumVertices()) {
	s3 = scalar_1[v1];
	scalar_grid_1.ComputeCoord(v1, coord1);
      }
      else {
	VERTEX_INDEX k = v1 - scalar_grid_1.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord1);
	s3 = new_mesh_vertices.Scalar(k);
      }	

     // cout<<"in second alpha uncertainty!!\n";	

       //s2 = .01;
       //s3 = .01;		

      // Compute expected level-crossing location for the edge and varince in level-crossing
      alpha_uncertainty_coord(dimension, s0, coord0, s1, coord1, isovalue, coord2, s2, s3, max_variance, color);
      normalized_variance[i] = *color;

      for (int d = 0; d < dimension; d++)
	{ coord[i*dimension+d] = coord2[d]; }
    }

    int color_scheme = 1;	
    float l1;
    float l2;
    float l3;	

    if(color_scheme == 0)
    {		
	// normalize the variance and copy into another array
     	for (int i = 0; i < nume; i++) 
	{
	  normalized_variance[i] = (double)(normalized_variance[i]/(*max_variance));	
    	}

        // Get the histogram for the normalized variances
	int bin_count = 15;
	float* var_hist = new float[bin_count];
	int* hist_flag = new int;
	float* mnvar = new float; 
	float* mxvar = new float; 
 	float* expected_value = new float;
	float* mode = new float;
    	float* mode_low = new float;	
	float* mode_high = new float;		
	 
	getHistogram(var_hist, normalized_variance, nume, bin_count,  hist_flag, mnvar, mxvar, expected_value,mode,mode_low,mode_high);

	cout<<"Expected normalized variance:"<<*expected_value<<"\n";
	cout<<"Minimum variance:"<<*mnvar<<"\n";
	cout<<"mode_low:"<<*mode_low<<"\n";
	cout<<"mode_high:"<<*mode_high<<"\n";
	cout<<"variance histogram:"<<"\n";
	double sm = 0;
	for (int i = 0; i < bin_count; i++)
	{ 
		cout<<var_hist[i]<<" ";
		sm+=var_hist[i];
	}
	cout<<"\n";
	cout<<"sum:"<<sm<<"\n";


	// set limits
	l1 = (double)((*mode_low - *mnvar)/2);
	l2 = (double)((*mxvar - *mode_high)/2);

	int lin_interp = 1;
	
	for (int i = 0; i < nume; i++) 
	{
		if(lin_interp == 0)
		{	

				if(normalized_variance[i] >= 0 && normalized_variance[i]<=l1 )
   				{	
					  //interpolate between green(0,1,0) & yellow(1,1,0)
					  f_color[4*i + 0] = 0;
				 	  f_color[4*i + 1] = 0.5;	
					  f_color[4*i + 2] = 0;	
					  f_color[4*i + 3] = 1;			
			        }     

				else if(normalized_variance[i] > l1 && normalized_variance[i]<=l2 )
    				{	
				          f_color[4*i + 0] = 0;
				 	  f_color[4*i + 1] = 0;	
					  f_color[4*i + 2] = 0.5;	
					  f_color[4*i + 3] = 1;			
			        }     

				else if(normalized_variance[i] > l2 && normalized_variance[i]<=1)	
    				{
				          f_color[4*i + 0] = 0.5;
				 	  f_color[4*i + 1] = 0;	
					  f_color[4*i + 2] = 0;	
					  f_color[4*i + 3] = 1;			
			        }     

		}

		else if(lin_interp == 1)
		{

			if((normalized_variance[i] >= 0) && (normalized_variance[i]<=1) )
   			{
	  			double t = normalized_variance[i];
	  			//interpolate between green(0,1,0) & red(1,0,0)
	  			f_color[4*i + 0] = t*1 + (1-t)*0;
	 	  		f_color[4*i + 1] = t*0 + (1-t)*1;	
		  		f_color[4*i + 2] = t*0 + (1-t)*0;	
		  		f_color[4*i + 3] = 1;			
    			}	
			
		}     

    	 }

     }

     // perform histogram equalization
     else if (color_scheme == 1)
     {

     		// normalize the variance and copy into another array
     		for (int i = 0; i < nume; i++)
	        {
	 		 normalized_variance[i] = (double)(normalized_variance[i]/(*max_variance));	
	  		 copy_normalized_variance[i] = normalized_variance[i];	  	 
    		}

		// sort the copied normalized variance
		std::sort(copy_normalized_variance, copy_normalized_variance+nume);

		// Count frequency for each intensity
	  	int track_level = 0;
	  	for (int i = 0; i < nume; i++) 
		{
		
			if( i == 0)
			{
				//copy the intensity level
				intensity_levels[track_level] = copy_normalized_variance[i];
				count_intensity[track_level] = 1;  
			}
			else
			{
				if( copy_normalized_variance[i] != copy_normalized_variance[i-1] )
				{
					track_level++;
					intensity_levels[track_level] = copy_normalized_variance[i];
					count_intensity[track_level] = 1; 
				}
				else
				{
					count_intensity[track_level]++;
				}
			}

		}

		// Find cdf
		int  cdf[track_level];
		float cdf_min;
		for (int i = 0; i <= track_level; i++) 
		{
			if ( i == 0)
			{
				cdf[i] = count_intensity[i];
				cdf_min = cdf[i];
			}
			else
				cdf[i] = count_intensity[i] + cdf[i-1];		
		}

		float cdf_i;	

		// create new variance array according to histogram equalization
		for (int i = 0; i < nume; i++) 
		{

			for(int j=0; j<=track_level; j++)
			{
				if (normalized_variance[i]  == intensity_levels[j])
					cdf_i = cdf[j];
			}
			
			normalized_variance[i] = (double)(((cdf_i - cdf_min)*track_level)/(nume - cdf_min));

		}	

		// transform the intensity values in equalized histogram space
		for(int j=0; j<=track_level; j++)
		{
			cdf_i = cdf[j];
			intensity_levels[j] = (double)(((cdf_i - cdf_min)*track_level)/(nume - cdf_min));	
		}	
	
		
		// set limits
		int lim = (intensity_levels[track_level] - intensity_levels[0])/3;
		l1 = lim;
		l2 = 2*lim;
		
		// For small intensity magnitudes, for appropriate color mapping change the denominator to the higher values
		l3 = (intensity_levels[track_level] - intensity_levels[0])/2;

		//cout<<"\nl1: "<<l1<<"\n";
		//cout<<"\nl2: "<<l2<<"\n";

		int lin_interp = 1;
		
    		//linearly interpolate from green to yellow to red
    		for (int i = 0; i < nume; i++) 
		{

			if(lin_interp == 0)
			{	

				if(normalized_variance[i] >= intensity_levels[0] && normalized_variance[i]<=l1 )
   				{	
					  //interpolate between green(0,1,0) & yellow(1,1,0)
					  f_color[4*i + 0] = 0;
				 	  f_color[4*i + 1] = 0.5;	
					  f_color[4*i + 2] = 0;	
					  f_color[4*i + 3] = 1;			
			        }     

				else if(normalized_variance[i] > l1 && normalized_variance[i]<=l2 )
    				{	
				          f_color[4*i + 0] = 0;
				 	  f_color[4*i + 1] = 0;	
					  f_color[4*i + 2] = 0.5;	
					  f_color[4*i + 3] = 1;			
			        }     

				else if(normalized_variance[i] > l2 && normalized_variance[i]<=intensity_levels[track_level])	
    				{
				          f_color[4*i + 0] = 0.5;
				 	  f_color[4*i + 1] = 0;	
					  f_color[4*i + 2] = 0;	
					  f_color[4*i + 3] = 1;			
			        }     

			}
			
			else if(lin_interp == 1)
			{

				if((normalized_variance[i] >= intensity_levels[0]) && (normalized_variance[i]<=l3) )
   				{
	  				double t = (double)((normalized_variance[i] - intensity_levels[0])/l3);
	  				//interpolate between green(0,1,0) & blue(0,0,1)
					f_color[4*i + 0] = t*0 + (1-t)*0;
			 		f_color[4*i + 1] = t*0 + (1-t)*1;	
					f_color[4*i + 2] = t*1 + (1-t)*0;	
					f_color[4*i + 3] = 1;			
    				}     

				//else if((normalized_variance[i] > l3) && (normalized_variance[i]<=intensity_levels[track_level]))
					else if((normalized_variance[i] > l3) && (normalized_variance[i]<=2*l3))
   			 	{	
	 				// double t = (double)((normalized_variance[i]-l3)/(intensity_levels[track_level]-l3));	
	 				double t = (double)((normalized_variance[i]-l3)/l3);	
          				//interpolate between blue(0,0,1) & red(1,0,0)
         				 f_color[4*i + 0] =  t*1 + (1-t)*0;
 	 				 f_color[4*i + 1] =  t*0 + (1-t)*0;	
	 				 f_color[4*i + 2] =  t*0 + (1-t)*1;	
					 f_color[4*i + 3] =  1;			
    				}   	
    				
    			else if(normalized_variance[i] > 2*l3)	
    			{
    					 f_color[4*i + 0] =  1;
 	 				 		 f_color[4*i + 1] =  0;	
	 					 	 f_color[4*i + 2] =  0;	
							 f_color[4*i + 3] =  1;			    				
    			}
    		 
			}      
  		} 	
	}
}


  void position_iso_vertices_multilinear_binary
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const MC_MESH_VERTEX_LIST & new_mesh_vertices,
   const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
  // compute position of isosurface vertices using multilinear interpolation
  // scalar_grid = scalar grid data
  // elist[] = list of edges containing isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist.size()/2;
    COORD_TYPE coord0[dimension];
    COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];
    VERTEX_INDEX axis_increment[dimension];

    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VERTEX_INDEX> increment(num_cube_vertices);
    compute_increment(scalar_grid, axis_increment);
    compute_cube_vertex_increment
      (dimension, axis_increment, increment.Ptr());

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist[2*i];
      VERTEX_INDEX v1 = elist[2*i+1];

      SCALAR_TYPE s0, s1;
      bool is_new0 = false;
      bool is_new1 = false;

      if (v0 < scalar_grid.NumVertices()) {
	s0 = scalar[v0];
	scalar_grid.ComputeCoord(v0, coord0);
      }
      else {
	VERTEX_INDEX k = v0 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord0);
	s0 = new_mesh_vertices.Scalar(k);
	is_new0 = true;
      }

      if (v1 < scalar_grid.NumVertices()) {
	s1 = scalar[v1];
	scalar_grid.ComputeCoord(v1, coord1);
      }
      else {
	VERTEX_INDEX k = v1 - scalar_grid.NumVertices();
	new_mesh_vertices.CopyCoord(k, coord1);
	s1 = new_mesh_vertices.Scalar(k);
	is_new1 = true;
      }

      if (is_new0 == is_new1) {
	linear_interpolate_coord(dimension, s0, coord0, s1, coord1, 
			   isovalue, coord2);
      }
      else {
	interpolate_using_multilinear
	  (scalar_grid, isovalue, increment.Ptr(), num_cube_vertices,
	   coord0, coord1, coord2);
      }

      for (int d = 0; d < dimension; d++)
	{ coord[i*dimension+d] = coord2[d]; }
    }

  }

}

void IJKMCUBE::position_iso_vertices_linear
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of isosurface vertices using linear interpolation
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binary(scalar_grid, isovalue, vlist, coord);
    break;

  case NEP:
    position_iso_vertices_linear_nep(scalar_grid, isovalue, vlist, coord);
    break;

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }
}


void IJKMCUBE::position_iso_vertices_linearB
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
// Compute position of isosurface vertices using linear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binaryB(scalar_grid, isovalue, elist, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}


void IJKMCUBE::position_iso_vertices_linearB
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
// Compute position of isosurface vertices using linear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binaryB
      (scalar_grid, isovalue, new_mesh_vertices, elist, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}

// Take mean and delta grids and compute geometry and its variance
void IJKMCUBE::position_iso_vertices_linearB_alpha
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color)
// Compute position of isosurface vertices using linear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear_alpha");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binaryB_alpha(scalar_grid, isovalue, elist, coord,scalar_grid_1,f_color);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}


void IJKMCUBE::position_iso_vertices_linearB_alpha
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord, const MC_SCALAR_GRID_BASE & scalar_grid_1, float* &f_color)
// Compute position of isosurface vertices using linear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_linear_alpha");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binaryB_alpha
      (scalar_grid, isovalue, new_mesh_vertices, elist, coord, scalar_grid_1, f_color);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}



void IJKMCUBE::position_iso_vertices_multilinear
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
// Compute position of isosurface vertices using multilinear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_multilinear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_multilinear_binary
      (scalar_grid, isovalue, new_mesh_vertices, elist, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}

/// Convert isosurface vertex indices to pairs of endpoints
/// representing grid edges containing isosurface vertices.
void IJKMCUBE::convert_indices_to_endpoint_pairs
(const MC_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<ISO_VERTEX_INDEX> iso_vlist,
 const MERGE_DATA & merge_data,
 std::vector<VERTEX_INDEX> & endpoint)
{
  const int dimension = scalar_grid.Dimension();
  AXIS_SIZE_TYPE axis_increment[dimension];

  endpoint.clear();

  compute_increment(scalar_grid, axis_increment);

  for (int i = 0; i < iso_vlist.size(); i++) {
    MERGE_INDEX isov = iso_vlist[i];
    VERTEX_INDEX v0 = merge_data.GetFirstEndpoint(isov);
    MERGE_INDEX dir = merge_data.GetEdgeDir(isov);
    VERTEX_INDEX v1 = v0 + axis_increment[dir];

    endpoint.push_back(v0);
    endpoint.push_back(v1);
  }

}


// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

void IJKMCUBE::position_ivol_vertices_linear
(const MC_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of interval volume vertices using linear interpolation
// scalar_grid = scalar grid data
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE axis_increment[dimension];
  const int numv = vlist.size();
  GRID_COORD_TYPE coord0[dimension];
  GRID_COORD_TYPE coord1[dimension];
  COORD_TYPE coord2[dimension];
  const int VERTEX_OFFSET = 2*dimension;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment);

  const int num_ivolv_per_gridv = 
    get_num_ivol_vertices_per_grid_vertex(dimension);


  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_ivolv_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_ivolv_per_gridv;

    if (k == VERTEX_OFFSET) {
      // isosurface vertex lies on grid vertex v0
      compute_coord(v0, dimension, axis_size, coord0);
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord0[d];
    }
    else {
      // isosurface vertex lies on grid edge
      int d = k/2;
      VERTEX_INDEX v1 = v0+axis_increment[d];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0);
      compute_coord(v1, dimension, axis_size, coord1);

      if (k % 2 == 0) {
	// linear interpolate using isovalue0
	linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue0, 
			   coord2);
      }
      else {
	// linear interpolate using isovalue1
	linear_interpolate_coord(dimension, s0, coord0, s1, coord1, isovalue1, 
			   coord2);
      }
      for (int d = 0; d < dimension; d++)
	coord[i*dimension+d] = coord2[d];
    }
  }

}


