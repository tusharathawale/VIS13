/// \file ijkgrid.txx
/// ijk templates defining regular grid classes and functions.
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008,2009,2010 Rephael Wenger

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

#ifndef _IJKGRID_
#define _IJKGRID_

#include <algorithm>

#include "ijk.txx"

namespace IJK {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  /// Default type for grid size
  typedef long GRID_SIZE_TYPE;

  // **************************************************
  // TEMPLATE CLASS GRID
  // **************************************************

  /// Base grid class
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  class GRID {

  protected:
    DTYPE dimension;     ///< grid dimension
    ATYPE * axis_size;   ///< axis_size[i] = # grid points along axis i
    NTYPE num_vertices;  ///< number of grid vertices

    void Init            ///< Initialize grid.
    (const DTYPE dimension, const ATYPE * axis_size);
    void FreeAll();      ///< Free all allocated memory.


  public:
    typedef DTYPE DIMENSION_TYPE;
    typedef ATYPE AXIS_SIZE_TYPE;
    typedef VTYPE VERTEX_INDEX_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    GRID(const DTYPE dimension, const ATYPE * axis_size);
    GRID() { Init(0, NULL); };
    ~GRID();
    GRID(const GRID & grid);                      // copy constructor
    const GRID & operator = (const GRID & right); // copy assignment

    // set functions
    template <class DTYPE2, class ATYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    DTYPE Dimension() const { return(dimension); }
    const ATYPE * AxisSize() const { return(axis_size); }
    ATYPE AxisSize(const DTYPE i) const { return(axis_size[i]); }
    NTYPE NumVertices() const { return(num_vertices); }

    // compute functions
    NTYPE ComputeNumCubes() const;
    template <class GTYPE>
    VTYPE ComputeVertexIndex(const GTYPE * coord) const;
    template <class GTYPE>
    void ComputeCoord(const VTYPE iv, GTYPE * coord) const;

    // compare
    template <class DTYPE2, class ATYPE2>
      bool CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const;
    template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
      bool CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid) const;

    // check function
    bool Check(const DTYPE dimension, const ATYPE * axis_size,
	       IJK::ERROR & error) const;
  };

  // **************************************************
  // TEMPLATE CLASS GRID_PLUS
  // **************************************************

  /// GRID class plus other indexes and operators for fast accessing of grid vertices.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  class GRID_PLUS:public GRID<DTYPE,ATYPE,VTYPE,NTYPE> {

  protected:

    /// iv+axis_increment[d] is vertex next to iv
    VTYPE * axis_increment;    

    /// iv0+cube_vertex_increment[k] = k'th vertex of cube with primary vertex iv0
    VTYPE * cube_vertex_increment;

    /// unit_cube_coord[dimension*k+j] = j'th coordinate of k'th vertex of unit cube
    NTYPE * unit_cube_coord;

    NTYPE num_cube_vertices;  ///< Number of cube vertices.

    void Init();
    void FreeAll();
    void Create();            ///< Allocate and set data in GRID_PLUS.

  public:
    GRID_PLUS(const DTYPE dimension, const ATYPE * axis_size):
      GRID<DTYPE,ATYPE,VTYPE,NTYPE> (dimension,axis_size) { Init(); };
    GRID_PLUS() { Init(); };
    ~GRID_PLUS();

    // set functions
    template <class DTYPE2, class ATYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    NTYPE NumCubeVertices() const       ///< Return number of cube vertices.
    { return(num_cube_vertices); };
    const VTYPE * AxisIncrement() const ///< Return array axis_increment[]
    { return(axis_increment); }
    const VTYPE AxisIncrement(const DTYPE d) const ///< Return axis_increment[d]
    { return(axis_increment[d]); }
    const VTYPE * CubeVertexIncrement() const ///< Return cube_vertex_increment[]
    { return(cube_vertex_increment); }
    const VTYPE CubeVertexIncrement   ///< Return cube_vertex_increment[k]
      (const DTYPE k) const              
    { return(cube_vertex_increment[k]); }

    /// Return pointer to coordinates of k'th cube vertex
    const NTYPE * UnitCubeCoord(const NTYPE k) 
    { return(unit_cube_coord+this->Dimension()*k); }
    const NTYPE UnitCubeCoord         ///< Return j'th coordinate of k'th vertex
    (const NTYPE k, const NTYPE j) 
    { return(unit_cube_coord[this->Dimension()*k+j]); }

    /// Return next vertex in direction d.
    /// @pre iv is not the last vertex in direction d.
    VTYPE NextVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv+axis_increment[d]); }

    /// Return previous vertex in direction d.
    /// @pre iv is not the first vertex in direction d.
    VTYPE PrevVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv-axis_increment[d]); }

    /// Return next vertex in direction d.
    /// @pre iv is a primary cube vertex.
    /// @pre k is less than the number of unit cube vertices.
    VTYPE CubeVertex(const VTYPE iv0, const int k) const  
    { return(iv0+cube_vertex_increment[k]); }

  };

  // **************************************************
  // TEMPLATE CLASS FACET_LIST
  // **************************************************

  /// Base class for list of facets.
  template <class DTYPE>
  class FACET_LIST {

  protected:
    FACET_LIST(){};

  public:
    inline bool Contains(const DTYPE d) const;
  };

  /// list containing zero facets
  template <class DTYPE>
  class FACET_LIST0:public FACET_LIST<DTYPE> {

  public:
    FACET_LIST0() {};
    inline bool Contains(const DTYPE d) const { return(false); };
  };

  /// List containing one facet.
  template <class DTYPE>
  class FACET_LIST1:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;

  public:
    inline void Set(const DTYPE facet0) { this->facet0 = facet0; };
    FACET_LIST1(const DTYPE facet0) { Set(facet0); };
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0) { return(true); };
      return(false);
    };
  };

  /// List containing two facets.
  template <class DTYPE>
  class FACET_LIST2:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;
    DTYPE facet1;

  public:
    inline void Set(const DTYPE facet0, const DTYPE facet1) 
    { this->facet0 = facet0; this->facet1 = facet1; };
    FACET_LIST2(const DTYPE facet0, const DTYPE facet1) 
    { Set(facet0, facet1); };
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0 || d == facet1) { return(true); };
      return(false);
    };
  };

  // **************************************************
  // inline UTILITY FUNCTIONS
  // **************************************************

  /// Integer divide.
  template <class ATYPE, class BTYPE>
  inline long integer_divide(const ATYPE a, const BTYPE b)
  { return(long(a)/long(b)); }

  /// Integer divide.
  inline long integer_divide(const long a, const long b)
  { return(a/b); }

  /// Integer divide.
  inline int integer_divide(const int a, const int b)
  { return(a/b); }
  
  /// Integer divide.
  inline unsigned long 
  integer_divide(const unsigned long a, const unsigned long b)
  { return(a/b); }

  /// Integer divide.
  inline unsigned int
  integer_divide(const unsigned int a, const unsigned int b)
  { return(a/b); }

  // **************************************************
  // TEMPLATE FUNCTIONS: COUNTING AND INDEXING
  // **************************************************

  ///
  /// \defgroup counting Counting and indexing
  /* \ingroup counting */
  /* \{ */

  /// Return coordinate of vertex \a iv.
  /// @param iv = Vertex index.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param[out] coord[k] = \a k'th coordinate of vertex \a iv.
  /// @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
  /// @pre \li Array coord[] is pre-allocated to size at least \a dimension.
  template <class VTYPE, class DTYPE, class ATYPE, class GTYPE>
  void compute_coord(const VTYPE iv, const DTYPE dimension,
		     const ATYPE * axis_size, GTYPE * coord)
  // compute coordinates of grid vertex iv
  {
    VTYPE k = iv;
    for (DTYPE d = 0; d < dimension; d++) {
      coord[d] = k % axis_size[d];
      k = k / axis_size[d];
    };
  }

  /// Return index of vertex with specified coord
  template <class VTYPE, class GTYPE, class DTYPE, class ATYPE>
  VTYPE compute_vertex_index
  (const GTYPE * coord, const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE iv = 0;
    VTYPE inc = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      iv += inc*coord[d];
      inc = inc*axis_size[d];
    }

    return(iv);
  }

  /// Return number of vertices in grid or subgrid.
  /// @param facet_list = List of facets determining subgrid.
  ///   If empty, return number of vertices in entire grid.
  ///   Otherwise, return number of vertices contained in the intersection of all facets in facet_list.
  template <class DTYPE, class ATYPE, class FTYPE, class NTYPE>
  void compute_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const FTYPE & facet_list,
   NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) 
	{ num_vertices = num_vertices * axis_size[d]; }
    }
  }

  /// Return number of subsampled vertices along axis
  /// @param axis_size = number of vertices along axis
  /// @param subsample_period = only count every k'th vertex along axis where k = subsample_period.
  /// @pre subsample_period is a positive integer.
  template <class ATYPE, class PTYPE>
  inline ATYPE compute_num_subsampled_vertices_along_axis
  (const ATYPE axis_size, const PTYPE subsample_period)
  { 
    return(integer_divide(axis_size+subsample_period-1, subsample_period)); 
  }

  /// Return number of vertices in subsampled grid or grid subspace.
  template <class DTYPE, class ATYPE, class PTYPE, class FTYPE, class NTYPE>
  void compute_num_subsampled_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const PTYPE subsample_period, const FTYPE & facet_list,
   NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
    // subsample_period[d] = only count every k'th vertex along axis d
    //   where k = subsample_period[d]
    // facet_list = ignore facets in facet_list
    // num_vertices = number of vertices.
  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {

	if (subsample_period[d] < 1) 
	  { throw_subsample_period_error
	      ("compute_num_subsampled_vertices", subsample_period[d]); }

	ATYPE num_subsampled_vertices_along_axis =
	  compute_num_subsampled_vertices_along_axis
	  (axis_size[d], subsample_period[d]);
	num_vertices = num_vertices * num_subsampled_vertices_along_axis;
      }
    }
  }

  /// Return number of grid vertices
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  { 
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_vertices(dimension, axis_size, facet_list0, num_vertices);
  }

  /// Return number of grid vertices between two vertices.
  /// @param iv0 = Vertex index.
  /// @param iv1 = Vertex index.
  /// @param[out] num_vertices = Number of grid vertices between \a iv0 and \a iv1.
  /// @pre 0 <= iv0 <= iv1 < (total number of grid vertices)
  template <class VTYPE, class DTYPE, class ATYPE, class NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, NTYPE & num_vertices)
  {
    VTYPE coord0[dimension];
    VTYPE coord1[dimension];
    IJK::PROCEDURE_ERROR error("compute_num_grid_vertices");

    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] > coord1[d]) {
	error.AddMessage("Programming error in calculating ", d,
			 "'th coordinate.");
	error.AddMessage("  coord0 (", coord0[d], 
			 ") > coord1 (", coord1[d], ").");
	throw error;
      }

      num_vertices = num_vertices * (coord1[d]-coord0[d]+1);
    }
  }

  /// Return number of vertices in grid interior.
  template <class DTYPE, class ATYPE, class WTYPE, class NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    ATYPE interior_axis_size[dimension];

    num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 2*boundary_width) { return; };
      interior_axis_size[d] = axis_size[d]-2*boundary_width;
    }
    compute_num_grid_vertices(dimension, interior_axis_size, num_vertices);
  }

  /// Return number of vertices in grid interior for boundary with 1.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  {
    compute_num_interior_grid_vertices
      (dimension, axis_size, 1, num_vertices);
  }

  /// Return number of vertices in grid boundary.
  template <class DTYPE, class ATYPE, class WTYPE, class NTYPE>
  void compute_num_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, NTYPE & num_boundary_vertices)
  {
    NTYPE num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    NTYPE num_interior_vertices;
    compute_num_interior_grid_vertices
      (dimension, axis_size, boundary_width, num_interior_vertices);
    num_boundary_vertices = num_grid_vertices - num_interior_vertices;
  }

  /// Return number of cubes in grid or grid subspace
  template <class DTYPE, class ATYPE, class FTYPE, class NTYPE>
  void compute_num_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   const FTYPE & facet_list, NTYPE & num_cubes)
    // facet_list = ignore facets in facet_list
  {
    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
	if (axis_size[d] < 2) { 
	  num_cubes = 0;
	  return; 
	};
	num_cubes = num_cubes * (axis_size[d]-1);
      }
    }
  }

  /// Return number of grid cubes
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_cubes(dimension, axis_size, facet_list0, num_cubes);
  }

  /// Return number of inner vertices in grid or grid subspace.
  /// Inner vertices are vertices which are not on an outer facet of the grid.
  /// Outer facets are grid facets which do not contain the origin.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_inner_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    if (dimension < 1) { 
      num_vertices = 0;
      return; 
    };

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE asize = axis_size[d];
      if (asize < 2) { 
	num_vertices = 0;
	return; 
      }
      else { num_vertices = num_vertices * (asize-1);  }
    }
  }

  /// Return number of outer vertices in grid or grid subspace.
  /// Outer vertices are vertices which are are on an outer facet of the grid.
  /// Outer facets are grid facets which do not contain the origin.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_outer_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_outer_vertices)
    // dimension = grid dimension
    // axis_size[d] = number of vertices along grid axis d
  {
    NTYPE num_grid_vertices;
    NTYPE num_inner_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    compute_num_inner_vertices(dimension, axis_size, num_inner_vertices);
    num_outer_vertices = num_grid_vertices-num_inner_vertices;
  }

  /// Return number of cube vertices
  template <class DTYPE> 
  long compute_num_cube_vertices(const DTYPE dimension)
  { return(1L << dimension); };

  /// Return number of cube facet vertices
  /// @pre dimension > 0
  template <class DTYPE> 
  long compute_num_cube_facet_vertices(const DTYPE dimension)
  { return(1L << dimension-1); };

  /// Return number of grid vertices in a region all of whose edges have length region_edge_length.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_vertices)
  {
    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_vertices = num_region_vertices * (region_edge_length+1);
    }
  }

  /// Return number of grid cubes in a region all of whose edges have length region_edge_length.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_cubes)
  {
    if (region_edge_length < 1) { 
      num_region_cubes = 0; 
      return;
    };

    num_region_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_cubes = num_region_cubes * region_edge_length;
    }
  }

  /// Return number of vertices in a region.
  /// @param max_region_edge_length = Maximum number of grid edges contained in each region edge.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const ATYPE max_region_edge_length,
   NTYPE & num_region_vertices)
  {
    ATYPE coord[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE numv_along_axis = max_region_edge_length + 1;
      if (coord[d] + max_region_edge_length >= axis_size[d]) {
	if (coord[d] < axis_size[d]) 
	  { numv_along_axis = axis_size[d] - coord[d]; }
	else
	  { numv_along_axis = 0; };
      }
      num_region_vertices = num_region_vertices * numv_along_axis;
    }
  }

  /// Return number of regions along a single axis.
  /// @param axis_size = Number of grid vertices on the axis.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <class ATYPE>
  ATYPE compute_num_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    ATYPE num_regions_along_axis = 
      long(axis_size+region_edge_length-2)/long(region_edge_length);
    return(num_regions_along_axis);
  }

  /// Return total number of regions in grid or subgrid.
  /// @param facet_list = List of facets determining subgrid.
  ///   If empty, return number of regions in entire grid.
  ///   Otherwise, return number of regions in subgrid formed by intersection of all facets in facet_list.
  template <class DTYPE, class ATYPE, class FTYPE, class NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (!facet_list.Contains(d)) {
	if (axis_size[d] <= 1) { 
	  num_regions = 0;
	  return; 
	};

	ATYPE num_regions_along_axis = 
	  compute_num_regions_along_axis(axis_size[d], region_edge_length);
	num_regions = num_regions * num_regions_along_axis; 
      }
    }
  }

  /// Return total number of regions.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_regions
      (dimension, axis_size, region_edge_length, facet_list0, num_regions);
  }

  /// Return number of full regions along a single axis.
  /// @param axis_size = Number of grid vertices on the axis.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <class ATYPE>
  ATYPE compute_num_full_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    ATYPE num_full_regions_along_axis = 
      long(axis_size-1)/long(region_edge_length);
    return(num_full_regions_along_axis);
  }

  /// Return number of full regions in grid or subgrid.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param facet_list = List of facets determining subgrid.
  ///   If empty, return number of vertices in entire grid.
  ///   Otherwise, return number of vertices in subgrid. formed by intersection of all facets in facet_list.
  /// @pre \a region_edge_length > 0.
  template <class DTYPE, class ATYPE, class FTYPE, class NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_full_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_full_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_full_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
	NTYPE k = compute_num_full_regions_along_axis
	  (axis_size[d], region_edge_length);
	num_full_regions = num_full_regions*k;
      }
    }
  }

  /// Return number of full regions in grid.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_full_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_full_regions);
  }

  /// Return number of partial regions along a single axis (0 or 1.)
  /// @param axis_size = Number of grid vertices on the axis.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @pre \a region_edge_length > 0.
  template <class ATYPE>
  ATYPE compute_num_partial_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    if (long(axis_size-1)%long(region_edge_length) == 0) { return(0); };
    return(1);
  }

  /// Return number of partial regions in grid or subgrid.
  template <class DTYPE, class ATYPE, 
	    class FTYPE, class NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_partial_regions)
  {
    NTYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_regions);
    NTYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_full_regions);
    num_partial_regions = num_regions-num_full_regions;
  }

  /// Return number of partial regions in grid.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, 
   NTYPE & num_partial_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_partial_regions);
  }

  /// Return total number of regions in grid facet.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param orth_dir = Direction orthogonal to facet.
  template <class NTYPE, class DTYPE, class ATYPE>
  NTYPE compute_num_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE * region_edge_length, const DTYPE orth_dir)
  {
    NTYPE num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (d != orth_dir) {
	if (axis_size[d] <= 1) { return(0); };

	ATYPE num_regions_along_axis = 
	  compute_num_regions_along_axis(axis_size[d], region_edge_length[d]);
	num_regions = num_regions * num_regions_along_axis; 
      }
    }
    return(num_regions);
  }

  /// Return total number of full regions in grid facet.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param orth_dir = Direction orthogonal to facet.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_full_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_full_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_full_regions);
  }

  /// Return total number of partial regions in grid facet.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_partial_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_partial_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_partial_regions);
  }

  /// Return number of vertices along subsampled axis.
  /// @param axis_size = Number of vertices along axis.
  /// @param subsample_period = Only count every k'th vertex along axis where k = \a subsample period.
  /// @pre \a subsample_period > 0.
  template <class ATYPE, class PTYPE>
  ATYPE compute_subsample_size
  (const ATYPE axis_size, const PTYPE subsample_period)
  {
    ATYPE subsample_size =
      integer_divide(axis_size+subsample_period-1, subsample_period);
    return(subsample_size);
  }

  /// Return number of vertices in subsampled grid.
  /// @param subsample_period[d] = subsample period along axis d
  template <class DTYPE, class ATYPE, class PTYPE, class NTYPE>
  void compute_subsample_size
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, NTYPE & num_vertices)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_subsampled_vertices
      (dimension, axis_size, subsample_period, facet_list0, num_vertices);
  };

  /// Return number of vertices along supersampled axis.
  /// @param axis_size = Number of vertices along axis.
  /// @param supersample_period = Only count every k'th vertex along axis where k = \a supersample period.
  /// @pre \a supersample_period > 0.
  template <class ATYPE, class PTYPE>
  ATYPE compute_supersample_size
  (const ATYPE axis_size, const PTYPE supersample_period)
  {
    if (axis_size < 1) { return(0); };
    ATYPE supersample_size = (axis_size-1)*supersample_period+1;
    return(supersample_size);
  }

  /* \} */

  // **************************************************
  // TEMPLATES TO CHECK VALUES AND ARRAYS.
  // **************************************************

  /// Check dimension
  /// Return true if dimension is non-negative
  template <class DTYPE>
  bool check_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension < 0) {
      error.AddMessage("Illegal dimension ", dimension, ".");
      error.AddMessage("Dimension must be non-negative.");
      return(false);
    }

    return(true);
  }

  /// Check range of vertices
  /// Return true if 0 <= iv0 <= iv1 < total_num_grid_vertices
  template <class VTYPE, class DTYPE, class ATYPE>
  bool check_range
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, IJK::ERROR & error)
  {
    GRID_SIZE_TYPE total_num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, total_num_grid_vertices);

    if (iv0 > iv1) {
      error.AddMessage("Illegal vertex range. Vertex index ", iv0, 
		       " is greater than vertex index ", iv1, ".");
      return(false);
    }

    if (0 > iv0 || iv1 >= total_num_grid_vertices) {
      error.AddMessage("Illegal vertex indices: ", iv0, ", ", iv1, ".");
      error.AddMessage("Vertex indices should be in range: [0,",
		       total_num_grid_vertices, "].");
      return(false);
    }

    return(true);
  }

  /// Check region coordinates.
  template <class DTYPE, class VTYPE, class CTYPE>
  bool check_region_coordinates
  (const DTYPE dimension, const VTYPE iv0, const CTYPE * coord0, 
   const VTYPE iv1, const CTYPE * coord1, IJK::ERROR & error)
  // return true if (coord0[d] <= coord1[d]) for all d
  {
    for (DTYPE d = 0; d < dimension; ++d) {
      if (coord0[d] > coord1[d]) {
	error.AddMessage("Illegal coordinates.  Coordinates of vertex 0 should be less than or equal to coordinates of vertex 1.");
	error.AddMessage(" Vertex 0 = ", iv0, 
			 ".  Coordinate ", d , " = ", coord0[d], ".");
	error.AddMessage(" Vertex 1 = ", iv1, 
			 ".  Coordinate ", d, " = ", coord1[d], ".");
	return(false);
      }
    }

    return(true);
  }

  /// Check that axis size is positive.
  template <class DTYPE, class ATYPE>
  bool check_positive_axis_size
  (const DTYPE dimension, ATYPE * axis_size, IJK::ERROR & error)
  {
    for (DTYPE d = 0; d < dimension; d++) {
      error.AddMessage("Illegal axis size.  Axis size must be positive.");
      error.AddMessage("  Axis ", d, " has size ", axis_size[d], ".");
      return(false);
    }
    return(true);
  }

  /// Check that array vertex_list[] is not NULL.
  template <class VTYPE>
  bool check_vertex_list(const VTYPE * vertex_list, IJK::ERROR & error)
  {
    if (vertex_list == NULL) {
      error.AddMessage("Vertex list is NULL");
      return(false);
    }
    else { return(true); }
  }

  /// Check that region edge length is positive.
  template <class LTYPE>
  bool check_region_edge_length(const LTYPE length, IJK::ERROR & error)
  {
    if (length <= 0) {
      error.AddMessage("Region edge length must be positive.");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices added to vertex list.
  template <class NTYPE0, class NTYPE1>
  bool check_num_vertices_added
  (const NTYPE0 num_added, const NTYPE1 numv, IJK::ERROR & error)
  {
    if (num_added != numv) {
      error.AddMessage("Added ", num_added, " vertices to vertex list.");
      error.AddMessage("Number of vertices i list should be ", numv, ".");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices equals number of grid vertices
  template <class DTYPE, class ATYPE, class NTYPE>
  bool check_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const NTYPE numv,
   IJK::ERROR & error)
  {
    NTYPE numv2;
    compute_num_grid_vertices(dimension, axis_size, numv2);
    if (numv != numv2) {
      error.AddMessage("Programming error. Incorrect number of vertices.");
      return(false);
    }
    return(true);
  }

  /// Throw subsample period error.
  template <class STRING_TYPE, class PTYPE>
  void throw_subsample_period_error
  (const STRING_TYPE proc_name, const PTYPE subsample_period)
  {
    IJK::PROCEDURE_ERROR error(proc_name);
    error.AddMessage("Subsample period must be a positive integer.");
    error.AddMessage("  Subsampling period = ", subsample_period, ".");
    throw error;
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTING INCREMENTS
  // **************************************************

  /// Compute increment to add to index of a vertex to get next vertices along the axes.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[] = Axis size. axis_size[d] is the number of vertices along axis \a d.
  /// @param[out] increment[] = Axis increment. iv+increment[d] is the index of the vertex after iv along axis \a d.
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  // NOTE: *** SHOULD RENAME THIS: compute_axis_increment.
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_increment");

    if (dimension <= 0) { return; };

    if (axis_size == NULL || increment == NULL) {
      error.AddMessage("Programming error. axis_size == NULL or increment == NULL.");
      throw error;
    }

    increment[0] = 1;
    for (DTYPE d = 1; d < dimension; d++)
      { increment[d] = increment[d-1]*axis_size[d-1]; }
  }

  /// Compute increment to add to index of current vertex to get
  ///   next vertex along each axis
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE, class ITYPE>
  void compute_increment
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, ITYPE * increment)
  {
    compute_increment(grid.Dimension(), grid.AxisSize(), increment);
  }

  /// Compute increment to add to index of vertex to get next subsampled vertex along each axis.
  /// @param subsample_period[d] = Subsample period along axis \a d.
  /// @param[out] increment[d] = Increment to add to index of a vertex to get next subsampled vertex along axis \a d.
  /// @pre Array increment[] is pre-allocated to size at least \a dimension.
  template <class DTYPE, class ATYPE, class PTYPE, class ITYPE>
  void compute_subsample_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, ITYPE * increment)
  {
    compute_increment(dimension, axis_size, increment);

    for (DTYPE d = 0; d < dimension; ++d)
      { increment[d] *= subsample_period[d]; };
  }

  /// Compute increment to add to index of primary vertex to get
  ///   k'th corner vertex of region.
  /// @param[out] increment[] = Region corner increment. iv0+increment[k] is k'th corner vertex of region.
  /// @pre Array increment[] is allocated with size at least number of corner regions.
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE * region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_region_corner_increment");
    ITYPE axis_increment[dimension];

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (region_axis_size == NULL) {
      error.AddMessage("Programming error. region_axis_size == NULL.");
      throw error;
    }    
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size[d] < 1) {
	error.AddMessage("Programming error.  Region axis size must be at least 1.");
	throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment);
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
	if ((j0 % 2) == 1) {
	  increment[j] = increment[j] + (region_axis_size[d]-1)*axis_increment[d];
	};
	j0 = j0/2;
      };
    }
  }

  /// Compute increment to add to index of primary vertex to get
  ///   k'th corner vertex of cubic region (all axes have same size.)
  /// @param dimension = Dimension of grid.
  /// @param grid_axis_size[] = Grid axis size. grid_axis_size[i] is the number of grid vertices along axis i.
  /// @param region_axis_size = Region axis size. Number of region vertices along any axis.
  /// @param[out] increment[] = Region corner increment. iv0+increment[k] is the k'th corner vertex of the region.
  /// @pre Array increment[] is allocated with size at least number of corner regions.
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_cubic_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_cubic_region_corner_increment");
    ITYPE axis_increment[dimension];

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size < 1) {
	error.AddMessage("Programming error.  Region axis size must be at least 1.");
	throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment);
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
	if ((j0 % 2) == 1) {
	  increment[j] = increment[j] + (region_axis_size-1)*axis_increment[d];
	};
	j0 = j0/2;
      };
    }
  }


  /// Compute increment to add to index of primary vertex to get k'th vertex in region.
  /// @pre Array increment is allocated with size at least number of region vertices
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_region_vertex_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const ATYPE scale, ITYPE * increment)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(scale);

    for (DTYPE d = 0; d < dimension; d++) {
      subgrid_axis_size[d] = region_edge_length*scale+1; 

      if (subgrid_axis_size[d] > axis_size[d]) {
	IJK::PROCEDURE_ERROR error("compute_region_vertex_increment");
	error.AddMessage("Region is larger than grid.");
	error.AddMessage("  Grid axis size[", d, "] = ", axis_size[d], ".");
	error.AddMessage("  Region size[", d, "] = ", subgrid_axis_size[d],
			 ".");
	throw error;
      }
    }

    subsample_subgrid_vertices
      (dimension, axis_size, 0, subgrid_axis_size, subsample_period, 
       increment);
  }

  /// Compute increment to add to index of primary vertex of region to get
  ///   primary vertex of k'th grid cube in cubic region.
  /// @param region_edge_length = Number of grid edges along every edge of the cubic region.
  /// @param increment = Region increment. iv0+increment[k] is the k'th cube of the region.
  /// @pre Array increment[] is allocated with size at least number of grid cubes contained in the region.
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_region_cube_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, ITYPE * increment)
  {
    ATYPE subgrid_size[dimension];

    if (region_edge_length < 1)  // No cubes.
      { return; }

    for (DTYPE d = 0; d < dimension; d++) 
      { subgrid_size[d] = region_edge_length; }

    get_subgrid_vertices(dimension, axis_size, 0, subgrid_size, increment);
  }

  /// Compute increment to add to vertex 0 to compute vertex i of hypercube.
  /// @param dimension = Dimension of grid.
  /// @param axis_increment[] = Axis increment. iv+axis_increment[i] is the next vertex after vertex iv along axis i.
  /// @param[out] cube_vertex_increment[] = Cube vertex increment. iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube with primary vertex iv0.
  /// @pre Array cube_vertex_increment[] is allocated with size at least number of cube vertices
  template <class DTYPE, class ATYPE, class ITYPE>
  void compute_cube_vertex_increment
  (const DTYPE dimension, const ATYPE * axis_increment, 
   ITYPE * cube_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_vertex_increment");

    if (dimension <= 0) { return; };

    if (axis_increment == NULL) {
      error.AddMessage("Programming error. Array axis_increment[] must be allocated and set before calling compute_cube_vertex_increment.");
      throw error;
    }

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated before calling compute_cube_vertex_increment.");
      throw error;
    }
    
    ITYPE num_cube_vertices = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_cube_vertices; j++) {
      cube_vertex_increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
	if ((j0 % 2) == 1) {
	  cube_vertex_increment[j] += axis_increment[d];
	};
	j0 = j0/2;
      }
    }
  }

  /// Compute unit cube vertex coordinates (0,0,...,0) to (1,1,...,1)
  /// @param dimension = Dimension of grid.
  /// @param[out] coord[] = Unit cube vertex coordinates.
  /// @pre Array coord[] is allocated with size at least (number of cube vertices)*dimension
  template <class DTYPE, class CTYPE>
  void compute_unit_cube_coord(const DTYPE dimension, CTYPE * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_unit_cube_vertex_coord");

    if (dimension <= 0) { return; };

    if (coord == NULL) {
      error.AddMessage("Programming error. Array coord[] must be allocated before calling compute_unit_cube_coord.");
      throw error;
    }
    
    long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
	coord[j*dimension+d] = j0 % 2;
	j0 = j0/2;
      }
    }
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: GETTING VERTICES
  // **************************************************

  /// Subsample vertices in subgrid.
  /// @param dimension = Grid dimension.
  /// @param axis_size[d] = Number of vertices along grid axis d.
  /// @param  subgrid_origin = Subgrid origin.
  ///   Note: subgrid_origin is always reported (unless num_vertices == 0.)
  /// @param subgrid_axis_size[d] = Number of vertices along subgrid axis d.
  /// @param subsample_period[d] = Only report every k'th vertex along subgrid axis \a d where k = \a subsample_period[d].
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li Subgrid is contained in grid, i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li \a subsample_period[d] is a positive integer.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <class DTYPE, class ATYPE, class PTYPE, class VTYPE>
  void subsample_subgrid_vertices
    (const DTYPE dimension, const ATYPE * axis_size, 
     const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
     const PTYPE subsample_period, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("subsample_subgrid_vertices");

    GRID_SIZE_TYPE num_vertices;
    compute_subsample_size
      (dimension, subgrid_axis_size, subsample_period, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: subgrid_axis_size[d] >= 1 for all d

    if (!check_vertex_list(vlist, error)) { throw error; };

    VTYPE subsample_increment[dimension];
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment);

    // initialize prev_num_vertices
    vlist[0] = subgrid_origin;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_subsampled_along_axis =
	compute_num_subsampled_vertices_along_axis
	(subgrid_axis_size[d], subsample_period[d]);

      VTYPE iv0 = subsample_increment[d];

      for (VTYPE i = 1; i < num_subsampled_along_axis; i++) {
	for (VTYPE * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	  *(vcur_ptr) = iv0 + *(vprev_ptr);
	  vcur_ptr++;
	};
	iv0 = iv0 + subsample_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_subsampled_along_axis;
    }

    if (!check_num_vertices_added(prev_num_vertices, num_vertices, error))
      throw error;

  }

  /// Subsample vertices in subgrid.
  template <class DTYPE, class ATYPE, class PTYPE, class VTYPE, class NTYPE>
  void subsample_subgrid_vertices
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices
      (grid.Dimension(), grid.AxisSize(), subgrid_origin, subgrid_axis_size,
       subsample_period, vlist);
  }

  /// Subsample vertices in grid.
  /// @param dimension = Grid dimension.
  /// @param axis_size[d] = Number of vertices along grid axis d.
  /// @param subsample_period[d] = Only report every k'th vertex along subgrid axis \a d where k = \a subsample_period[d].
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li \a subsample_period is a positive integer.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <class DTYPE, class ATYPE, class PTYPE, class VTYPE>
  void subsample_grid_vertices
    (const DTYPE dimension, const ATYPE * axis_size, 
     const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices(dimension, axis_size, VTYPE(0), axis_size, 
			       subsample_period, vlist);
  }

  /// Get vertices in subgrid.
  /// @param dimension = Grid dimension.
  /// @param axis_size[d] = Number of vertices along grid axis d.
  /// @param  subgrid_origin = Subgrid origin.
  ///   Note: subgrid_origin is always reported (unless num_vertices == 0.)
  /// @param subgrid_axis_size[d] = Number of vertices along subgrid axis d.
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li Subgrid is contained in grid, i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
  /// @pre \li \a subsample_period is a positive integer.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_subgrid_vertices
    (const DTYPE dimension, const ATYPE * axis_size, 
     const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
     VTYPE * vlist)
  {
    IJK::CONSTANT<ATYPE,ATYPE> ONE(1);

    subsample_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size, 
       ONE, vlist);
  }

  /// Get vertices in grid.
  /// @param dimension = Grid dimension.
  /// @param axis_size[d] = Number of vertices along grid axis d.
  /// @param[out] vlist[] = List of vertices.
  /// @pre \li Array vlist[] is preallocated to length at least number of vertices in grid or subgrid.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_grid_vertices
    (const DTYPE dimension, const ATYPE * axis_size, 
     VTYPE * vlist)
  {
    subsample_grid_vertices(dimension, axis_size, 1, vlist);
  }

  /// Get grid vertices in region between two grid vertices (inclusive).
  /// @param iv0 = Lower grid vertex.
  /// @param iv1 = Upper grid vertex.
  /// @param[out] vlist = List of vertices between \a iv0 and \a iv1.
  /// @pre \li 0 <= iv0 <= iv1 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_grid_vertices_between
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const VTYPE iv1, VTYPE * vlist)
  {
    ATYPE coord0[dimension];
    ATYPE coord1[dimension];
    ATYPE region_size[dimension];
    IJK::PROCEDURE_ERROR error("get_grid_vertices_between");

    if (dimension < 0) { return; };
    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    for (DTYPE d = 0; d < dimension; ++d) {
      if (!check_region_coordinates
	  (dimension, iv0, coord0, iv1, coord1, error)) { throw error; };

      region_size[d] = coord1[d]-coord0[d]+1;
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }

  /// Get grid vertices in region.
  /// @param iv0 = Primary vertex of region.
  /// @param[out] vlist = List of vertices in region.
  /// @pre \li 0 <= iv0 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist)
  {
    ATYPE coord[dimension];
    ATYPE region_size[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
	{ region_size[d] = max_region_edge_length + 1; }
      else if (coord[d] < axis_size[d])
	{ region_size[d] = axis_size[d] - coord[d]; }
      else {
	// Vertex iv0 does not lie inside grid.
	return;
      }
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }

  /// Get grid cubes in region.
  /// @param iv0 = Primary vertex of region.
  /// @param[out] vlist = List of vertices in region.
  /// @pre \li 0 <= iv0 < total_num_grid_vertices.
  /// @pre \li Array vlist[] is preallocated to size at least number of region vertices.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  void get_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist,
   NTYPE & num_cubes)
  {
    ATYPE coord[dimension];
    ATYPE subgrid_size[dimension];

    num_cubes = 0;
    if (max_region_edge_length < 1) { return; };

    compute_coord(iv0, dimension, axis_size, coord);

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
	{ subgrid_size[d] = max_region_edge_length; }
      else if (coord[d]+1 < axis_size[d])
	{ subgrid_size[d] = axis_size[d] - coord[d] - 1; }
      else {
	// Region contains no cubes
	return;
      }
    }

    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) 
      { num_cubes *= subgrid_size[d]; }

    get_subgrid_vertices(dimension, axis_size, iv0, subgrid_size, vlist);
  }

  /// Get primary cube vertices in subgrid.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_primary_cube_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size,
   VTYPE * vlist)
  {
    ATYPE subgrid_axis_size2[dimension];

    for (DTYPE d = 0; d < dimension; d++) {
      if (subgrid_axis_size[d] < 2) 
	  { return; }                      // zero cubes
      subgrid_axis_size2[d] = subgrid_axis_size[d]-1;
    }
    get_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size2, vlist);
  }

  // ********************************************************
  // TEMPLATE FUNCTIONS: FACET VERTICES, CUBES AND REGIONS
  // ********************************************************

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  /// @param boundary_width = Width of boundary, (Number of vertices. Must be non-negative.)
  template <class DTYPE, class ATYPE, class WTYPE, class NTYPE>
  void compute_num_vertices_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir, const WTYPE boundary_width,
     NTYPE & num_vertices)
  {
    if (dimension < 1) { num_vertices = 0; }
    else { num_vertices = 1; }

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
	if (axis_size[d] > 2*boundary_width) 
	  { num_vertices *= (axis_size[d]-2*boundary_width); }
	else
	  { num_vertices = 0; };
      }
    }
  }

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, 0, num_vertices);
  }

  /// Return number of vertices in specified grid facet
  /// Specify grid facet by the direction orthogonal to the facet.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  void compute_num_vertices_in_grid_facet
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, const DTYPE orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, num_vertices);
  }

  /// Return number of vertices in specified grid ridge
  /// Specify grid facet by the directions orthogonal to ridge
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_vertices_in_grid_ridge
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir0, const DTYPE orth_dir1,
     NTYPE & num_vertices)
  {
    FACET_LIST2<DTYPE> facet_list2(orth_dir0, orth_dir1);
    compute_num_vertices(dimension, axis_size, facet_list2, num_vertices);
  }

  /// Return maximum number of vertices over all grid facets
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_max_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   NTYPE & max_num_vertices)
  {
    max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      NTYPE num_face_vertices; 
      compute_num_vertices_in_grid_facet
	(dimension, axis_size, d, num_face_vertices);
      if (num_face_vertices > max_num_vertices)
	max_num_vertices = num_face_vertices;
    };
  }

  /// Return number of cubes in specified grid facet.
  /// Facet is lower facet orthogonal to the specified direction
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_cubes_in_grid_facet
    (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir,
     NTYPE & num_cubes)
  {
    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) { 
	num_cubes = 0;
	return; 
      };
      if (d != orth_dir) {
	num_cubes = num_cubes*(axis_size[d]-1);
      };
    }
  }

  /// Return number of cubes in grid facet orthogonal to axis 0
  template <class DTYPE, class ATYPE, class NTYPE>
  void compute_num_cubes_in_grid_facet0
    (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    compute_num_cubes_in_grid_facet(dimension, axis_size, DTYPE(0), num_cubes);
  }

  /// Get vertices in specified grid facet.
  /// Does not return any vertices if \a axis_size[orth_dir] == 0.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param orth_dir = Direction orthogonal to facet.
  /// @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
  /// @param boundary_width = Width of boundary, (Number of vertices. Must be non-negative.)
  /// @param[out] vlist[] = List of primary vertices.
  /// @pre Array vlist[] must be pre-allocated to size at least number of vertices in facet.
  template <class DTYPE, class ATYPE, class BTYPE, class VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir, 
   const bool side, const BTYPE boundary_width, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    ATYPE coord[dimension];

    if (axis_size[orth_dir] < 1) { return; };
    VTYPE subgrid_origin = 0;

    if (boundary_width == 0) {
      std::copy(axis_size, axis_size+dimension, subgrid_axis_size);
    
      if (side) {
	for (DTYPE d = 0; d < dimension; d++) { coord[d] = 0; };
	coord[orth_dir] = axis_size[orth_dir]-1;
	subgrid_origin = 
	  compute_vertex_index<VTYPE>(coord, dimension, axis_size);
      }
    }
    else {
      for (DTYPE d = 0; d < dimension; d++) {
	if (d != orth_dir) {
	  if (axis_size[d] < 2*boundary_width) { return; };

	  coord[d] = boundary_width; 
	  subgrid_axis_size[d] = axis_size[d]-2*boundary_width;
	}
      }

      if (side) { coord[orth_dir] = axis_size[orth_dir]-1; }
      else { coord[orth_dir] = 0; }

      subgrid_origin = 
	compute_vertex_index<VTYPE>(coord, dimension, axis_size);
    }
    subgrid_axis_size[orth_dir] = 1;

    get_subgrid_vertices(dimension, axis_size, subgrid_origin, 
			 subgrid_axis_size, vlist);
  }

  /// Get vertices in specified grid facet.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir, 
   const bool side, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, side, 0, vlist);
  }

  /// Get vertices in grid facet 0.
  template <class DTYPE, class ATYPE, class BTYPE, class VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, 
   const BTYPE boundary_width, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, DTYPE(0), false, boundary_width, vlist);
  }

  /// Get vertices in grid facet 0.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0(dimension, axis_size, 0, vlist);
  }

  /// Get vertices in grid facet 0.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  void get_vertices_in_grid_facet0
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0
      (grid.Dimension(), grid.AxisSize(), vlist);
  }

  /// Get vertices in specified grid ridge
  /// Ridge is lower ridge orthogonal to the specified directions
  /// @param[out] vlist[] = List of vertices in grid ridge.
  /// @pre Array vlist[] is preallocated to length at least number of vertices in grid ridge.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_vertices_in_grid_ridge
    (const DTYPE dimension, const ATYPE * axis_size,
     const DTYPE orth_dir0, const DTYPE orth_dir1, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];

    if (axis_size[orth_dir0] < 1 || axis_size[orth_dir1] < 1) { return; }
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size);
    subgrid_axis_size[orth_dir0] = 1;
    subgrid_axis_size[orth_dir1] = 1;
    
    get_subgrid_vertices(dimension, axis_size, 0, subgrid_axis_size, vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in specified grid facet where d = \a dimension.
  /// Does not return any vertices if \a axis_size[orth_dir] == 0.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[d] = Number of vertices along axis \a d.
  /// @param orth_dir = Direction orthogonal to facet.
  /// @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
  /// @param[out] vlist[] = List of primary vertices.
  /// @pre Array vlist[] must be pre-allocated to size at least number of primary vertices in facet.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, const bool side, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];

    if (axis_size[orth_dir] == 0) { return; };
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size);
    subgrid_axis_size[orth_dir] = 2;

    VTYPE subgrid_origin = 0;

    if (side) {
      VTYPE axis_increment[dimension];
      compute_increment(dimension, axis_size, axis_increment);

      subgrid_origin = axis_increment[orth_dir]*(axis_size[orth_dir]-1);
    }

    get_primary_cube_vertices
    (dimension, axis_size, subgrid_origin, subgrid_axis_size, vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in specified grid facet where d = \a dimension.
  /// @param orth_dir = Direction orthogonal to facet.  Facet contains (0,0,...,0).
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, VTYPE * vlist)
    // orth_dir = direction orthogonal to facet
  {
    get_cubes_in_grid_facet(dimension, axis_size, orth_dir, false, vlist);
  }

  /// Get primary vertices of (d-1)-dimensional cubes in grid facet 0 where d = \a dimension.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, DTYPE(0), vlist);
  }

  /// Get outer grid vertices
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_outer_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
    // Precondition: vlist[] is preallocated to size 
    //   at least num_outer_vertices
  {
    VTYPE axis_increment[dimension];

    if (dimension < 1) { return; };

    DTYPE d_last = dimension - 1;

    if (dimension == 1) {
      if (axis_size[0] < 2) { return; }
      else {
	vlist[0] = axis_size[0]-1;
	return;
      }
    }
    else {
      if (axis_size[d_last] < 2) { return; }

      get_outer_grid_vertices(d_last, axis_size, vlist);
      VTYPE num_vertices; 
      compute_num_outer_vertices(d_last, axis_size, num_vertices);

      compute_increment(dimension, axis_size, axis_increment);

      for (VTYPE i = 1; i < axis_size[d_last]-1; i++) {
	VTYPE k = i*num_vertices;
	VTYPE k_increment = i*axis_increment[d_last];
	for (VTYPE j = 0; j < num_vertices; j++) {
	  vlist[k+j] = vlist[j]+k_increment;
	} 
      }

      num_vertices = num_vertices * (axis_size[d_last]-1);

      VTYPE * vlist2 = vlist + num_vertices;
      FACET_LIST1<DTYPE> facet_list1(d_last);

      get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist2);
    }
  }

  /// Get primary vertices of regions in grid or subgrid.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of full regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of full regions.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);
    
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < 2) { return; };

      subgrid_axis_size[d] = axis_size[d]-1;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }

  /// Get primary vertices of regions in grid or subgrid.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of full regions.
  /// @param[out] is_full[k] = True if region \a k is a full \a LxLx...xL region where \a L = \a region_edge_length.
  /// @pre Array vlist[] must be pre-allocated to size at least number of regions.
  /// @pre Array is_full[] must be pre-allocated to size at least number of regions.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist, bool * is_full)
  {
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    GRID_SIZE_TYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, num_regions);

    if (num_regions < 1) { return; };
    // Note: axis_size[d] >= 2 for all d

    GRID_SIZE_TYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, num_full_regions);

    if (!check_vertex_list(vlist, error)) { throw error; };

    VTYPE subsample_increment[dimension];
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment);

    // set vlist[0], is_full[0] and initialize prev_num_vertices
    vlist[0] = 0;
    if (num_full_regions > 0) { is_full[0] = true; }
    else { is_full[0] = false; };
    VTYPE prev_num_regions = 1;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_regions_along_axis =
	compute_num_regions_along_axis(axis_size[d], region_edge_length);

      ATYPE num_full_along_axis =
	compute_num_full_regions_along_axis(axis_size[d], region_edge_length);

      VTYPE iv0 = subsample_increment[d];
      for (VTYPE i = 1; i < num_full_along_axis; i++) {
	VTYPE i2 = i * prev_num_regions;
	for (VTYPE j = 0; j < prev_num_regions; j++) {
	  VTYPE k = j + i2;
	  vlist[k] = vlist[j] + iv0;
	  is_full[k] = is_full[j];
	}
	iv0 = iv0 + subsample_increment[d];
      }

      if (num_regions_along_axis != num_full_along_axis &&
	  num_regions_along_axis > 1) {
	VTYPE i2 = num_full_along_axis * prev_num_regions;
	for (VTYPE j = 0; j < prev_num_regions; j++) {
	  VTYPE k = j + i2;
	  vlist[k] = vlist[j] + iv0;
	  is_full[k] = false;
	}
      }

      prev_num_regions = prev_num_regions*num_regions_along_axis;
    }

    if (!check_num_vertices_added(prev_num_regions, num_regions, error))
      throw error;

  }

  /// Get primary vertices of full regions in grid.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of full regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of full regions.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_full_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    IJK::PROCEDURE_ERROR error("get_full_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= region_edge_length) { return; };

      subgrid_axis_size[d] = axis_size[d]-region_edge_length;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }

  /// Get primary vertices of partial regions in grid.
  /// @param region_edge_length = Number of grid edges contained in each region edge.
  /// @param[out] vlist[] = Array of primary vertices of partial regions.
  /// @pre Array vlist[] must be pre-allocated to size at least number of partial regions.
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_partial_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("get_partial_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    GRID_SIZE_TYPE num_vertices;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: axis_size[d] >= 1 for all d not in facet_list

    if (!check_vertex_list(vlist, error)) { throw error; };

    VTYPE region_axis_increment[dimension];

    compute_subsample_increment(dimension, axis_size, subsample_period,
				region_axis_increment);

    // initialize prev_num_vertices
    VTYPE prev_num_vertices = 0;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {
      VTYPE * vcur_ptr = vlist + prev_num_vertices;

      ATYPE num_full_regions_along_axis = 
	compute_num_full_regions_along_axis
	(axis_size[d], region_edge_length);

      VTYPE iv0 = region_axis_increment[d];
      for (VTYPE i = 1; i < num_full_regions_along_axis; i++) {
	for (VTYPE * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
	  *(vcur_ptr) = iv0 + *(vprev_ptr);
	  vcur_ptr++;
	};
	iv0 = iv0 + region_axis_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_full_regions_along_axis;

      ATYPE num_partial_regions_along_axis =
	compute_num_partial_regions_along_axis
	(axis_size[d], region_edge_length);

      if (num_partial_regions_along_axis > 0) {

	VTYPE inc = 
	  region_axis_increment[d]*num_full_regions_along_axis;

	if (d == 0) {
	  vlist[prev_num_vertices] = inc;
	  prev_num_vertices++;
	}
	else {

	  get_region_primary_vertices(d, axis_size, region_edge_length, 
				      vlist + prev_num_vertices);

	  ATYPE k;
	  compute_num_regions(d, axis_size, region_edge_length, k);
	  for (VTYPE i = prev_num_vertices; i < prev_num_vertices+k; i++)
	    { vlist[i] += inc; }

	  prev_num_vertices += k;
	}
      }
    }

    if (prev_num_vertices != num_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
		       " vertices to vertex list.");
      error.AddMessage("Number of vertices in list should be ", 
		       num_vertices, ".");
      throw error;
    }

  }

  /// Get boundary grid vertices
  template <class DTYPE, class ATYPE, class VTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    if (dimension < 1) { return; }
    if (dimension == 1) {
      if (axis_size[0] > 0) { vlist[0] = 0; }
      if (axis_size[0] > 1) { vlist[1] = axis_size[0]-1; };
      return;
    }

    DTYPE d_last = dimension - 1;
    if (axis_size[d_last] < 1) { return; };

    // get vertices in lower facet
    VTYPE num_vertices_in_grid_facet =
      compute_num_vertices_in_grid_facet(dimension, axis_size, d_last);
    get_vertices_in_grid_facet(dimension, axis_size, d_last, false, vlist);

    VTYPE * vlist2 = vlist+num_vertices_in_grid_facet;
    VTYPE * vlist3 = vlist2;
    if (axis_size[d_last] > 2) {
      ATYPE axis_increment[dimension];

      compute_increment(dimension, axis_size, axis_increment);
      get_boundary_grid_vertices(dimension-1, axis_size, vlist2);

      VTYPE n = compute_num_boundary_grid_vertices(dimension-1, axis_size);
      for (VTYPE * vcur_ptr = vlist2; vcur_ptr != vlist2+n; vcur_ptr++)
	{ *vcur_ptr += axis_increment[d_last]; }

      vlist3 = vlist2 + n;
      for (ATYPE j = 2; j < axis_size[d_last]-1; j++) {
	VTYPE inc = axis_increment[d_last]*(j-1);

	for (VTYPE i = 0; i < n; i++)  { vlist3[i] = vlist2[i] + inc; }

	vlist3 += n;
      }
    }


    get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist3);
  }


  /// Get boundary grid vertices
  template <class DTYPE, class ATYPE, class VTYPE, class WTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, VTYPE * vlist)
  {
    if (dimension < 1) { return; }
    if (boundary_width < 1) { return; }

    const DTYPE d_last = dimension - 1;

    if (axis_size[d_last] <= 2*boundary_width) {
      // all vertices are on the boundary
      VTYPE num_grid_vertices;
      compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
      for (VTYPE j = 0; j < num_grid_vertices; j++) 
	{ vlist[j] = j; }
      return;
    }

    if (dimension == 1) {
      for (VTYPE j = 0; j < boundary_width; j++) 
	{ vlist[j] = j; }

      for (VTYPE j = 0; j < boundary_width; j++) {
	VTYPE iv = axis_size[0]-boundary_width + j;
	vlist[j+boundary_width] = iv;
      };
      return;
    }

    ATYPE axis_increment[dimension];
    compute_increment(dimension, axis_size, axis_increment);

    // get vertices in lower facet
    VTYPE num_vertices_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d_last, num_vertices_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d_last, false, vlist);

    // get remaining vertices in lower boundary
    for (VTYPE i = 1; i < boundary_width; i++) {
      VTYPE inc = i*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
	vlist[i*num_vertices_in_grid_facet + j] = vlist[j] + inc;
      }
    }

    VTYPE * vlist2 = vlist+boundary_width*num_vertices_in_grid_facet;
    get_boundary_grid_vertices(dimension-1, axis_size, 
			       boundary_width, vlist2);

    VTYPE num_boundary_grid_vertices;
    compute_num_boundary_grid_vertices
      (dimension-1, axis_size, boundary_width, num_boundary_grid_vertices);
    for (VTYPE * vcur_ptr = vlist2; 
	 vcur_ptr != vlist2+num_boundary_grid_vertices; vcur_ptr++)
      { *vcur_ptr += boundary_width*axis_increment[d_last]; }

    VTYPE * vlist3 = vlist2+num_boundary_grid_vertices;
    for (ATYPE j = boundary_width+1; j+boundary_width < axis_size[d_last]; 
	 j++) {
      VTYPE inc = axis_increment[d_last]*(j-boundary_width);

      for (VTYPE i = 0; i < num_boundary_grid_vertices; i++)  
	{ vlist3[i] = vlist2[i] + inc; }

      vlist3 += num_boundary_grid_vertices;
    }

    VTYPE * vlist4 = vlist3 + (boundary_width-1)*num_vertices_in_grid_facet;

    // get vertices in upper facet
    get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist4);

    // get remaining vertices in upper boundary
    for (VTYPE i = 0; i+1 < boundary_width; i++) {
      VTYPE inc = (boundary_width-i-1)*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
	vlist3[i*num_vertices_in_grid_facet + j] = vlist4[j] - inc;
      }
    }
  
  }

  // **************************************************
  // TEMPLATE CLASS GRID MEMBER FUNCTIONS
  // **************************************************

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::GRID
  (const DTYPE dimension, const ATYPE * axis_size)
  // constructor
  {
    Init(dimension, axis_size);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID<DTYPE,ATYPE,VTYPE,NTYPE>()
  // destructor
  {
    FreeAll();
  }

  /// Initialize grid.
  /// @param dimension = Dimension of grid.
  /// @param axis_size[i] = Number of grid vertices along axis i.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Init
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    this->axis_size = NULL;
    this->dimension = 0;
    this->num_vertices = 1;
    if (dimension > 0) 
      { SetSize(dimension, axis_size); };
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  template <class DTYPE2, class ATYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeAll();
    this->dimension = dimension;
    this->axis_size = new ATYPE[dimension];
    for (DTYPE d = 0; d < dimension; d++)
      { this->axis_size[d] = axis_size[d]; }

    compute_num_grid_vertices(dimension, axis_size, this->num_vertices);

    if (dimension < 0) {
      IJK::PROCEDURE_ERROR error("Grid::SetSize");
      error.AddMessage("Programming error.  Illegal dimension ",
		       dimension, ".");
      error.AddMessage("Dimension should be non-negative.");
      throw error;
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }


  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::FreeAll()
  {
    if (axis_size != NULL) { delete [] axis_size; };
    axis_size = NULL;
    dimension = 0;
    num_vertices = 0;
  }

  /// Copy constructor.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID(const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid)
  {
    Init(grid.Dimension(), grid.AxisSize());
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::operator = (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & right)
  // copy assignment	
  {
    if (&right != this) {         // avoid self-assignment
      SetSize(right.Dimension(), right.AxisSize());
    }
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumCubes() const
  {
    NTYPE num_grid_cubes;
    compute_num_grid_cubes(Dimension(), AxisSize(), num_grid_cubes);
    return(num_grid_cubes);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  template <class GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeVertexIndex(const GTYPE * coord) const
  {
    return(compute_vertex_index<VTYPE>(coord, Dimension(), AxisSize()));
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  template <class GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCoord(const VTYPE iv, GTYPE * coord) const
  // Precondition: coord[] is preallocated to length at least dimension
  {
    compute_coord(iv, Dimension(), AxisSize(), coord);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  template <class DTYPE2, class ATYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const
  // return true if grid dimension and axis size match parameters
  {
    if (dimension != this->Dimension()) { return(false); };
    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->AxisSize(d)) { return(false); };
    }

    return(true);
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE>
  template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid) const
  // return true if grid dimension and axis size match parameters
  {
    return(CompareSize(grid.Dimension(), grid.AxisSize()));
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Check
  (const DTYPE dimension, const ATYPE * axis_size, IJK::ERROR & error) const
  // return true if grid dimension and axis size match parameter dimension
  //      and parameter axis size, respectively
  {
    if (dimension != this->dimension) {
      error.AddMessage("Incorrect grid dimension ", this->dimension, ".");
      return(false);
    }

    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->axis_size[d]) {
	error.AddMessage("Illegal axis size[", d, "] = ", 
			 this->axis_size[d], ".");
	return(false);
      }
    }

    NTYPE num_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_vertices);

    if (num_vertices != this->num_vertices) {
      error.AddMessage("Incorrect number of grid vertices ", 
		       this->num_vertices, ".");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS GRID_PLUS MEMBER FUNCTIONS
  // **************************************************

  /// Initialize data structures in GRID_PLUS
  /// Precondition: dimension and axis_size are already set
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::Init()
  {
    this->axis_increment = NULL;
    this->cube_vertex_increment = NULL;
    this->unit_cube_coord = NULL;
    this->num_cube_vertices = 0;

    if (this->Dimension() > 0) { Create(); };
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>()
  // destructor
  {
    FreeAll();
  }

  /// Allocate arrays and compute data in GRID_PLUS.
  /// Precondition: dimension and axis_size[] are already set
  /// Precondition: axis_increment, cube_vertex_increment and unit_cube_coord
  ///   contain NULL pointers.
  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::Create()
  {
    IJK::PROCEDURE_ERROR error("GRID_PLUS::Create");

    if (!check_dimension(this->Dimension(), error)) { throw error; };
    if (axis_increment != NULL || cube_vertex_increment != NULL ||
	unit_cube_coord != NULL) {
      error.AddMessage("Programming error. Previously allocated memory not released.");
      throw error;
    }

    this->num_cube_vertices = compute_num_cube_vertices(this->Dimension());
    this->axis_increment = new VTYPE[this->Dimension()];
    this->cube_vertex_increment = new VTYPE[num_cube_vertices];
    this->unit_cube_coord = new NTYPE[num_cube_vertices*this->Dimension()];

    compute_increment
      (this->Dimension(), this->AxisSize(), this->axis_increment);
    compute_cube_vertex_increment
      (this->Dimension(), this->AxisIncrement(), this->cube_vertex_increment);

    compute_unit_cube_coord(this->Dimension(), this->unit_cube_coord);
  }


  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::FreeAll()
  {
    if (axis_increment != NULL) { delete [] axis_increment; }
    if (cube_vertex_increment != NULL) 
      { delete [] cube_vertex_increment; };
    if (unit_cube_coord != NULL) { delete [] unit_cube_coord; };
    axis_increment = NULL;
    cube_vertex_increment = NULL;
    unit_cube_coord = NULL;
    num_cube_vertices = 0;
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  template <class DTYPE2, class ATYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    IJK::PROCEDURE_ERROR error("GRID_PLUS::SetSize");

    FreeAll();

    GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize(dimension, axis_size);
    Create();
  }

  template <class DTYPE, class ATYPE, class VTYPE, class NTYPE> 
  template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

}

#endif
