/// \file ijkscalar_grid.txx
/// ijk templates defining scalar grid classes and functions.
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

#ifndef _IJKSCALAR_GRID_
#define _IJKSCALAR_GRID_

#include <algorithm>

#include "ijk.txx"
#include "ijkgrid.txx"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_BASE
  // **************************************************

  /// Class storing scalar grid information.
  /// Does not allocate or destroy any scalar information.
  /// Note: SCALAR_GRID_BASE does not provide any way to allocate memory
  ///    for scalar or set scalar to point to existing memory.
  template <class GRID_CLASS, class STYPE>
  class SCALAR_GRID_BASE:public GRID_CLASS {

  protected:
    STYPE * scalar;
    // point (x0,x1,x2,...) has scalar value
    //     scalar_grid[x0 + x1*axis_size[0] + 
    //                   x2*axis_size[0]*axis_size[1] + ...]

    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;
    
  public:
    SCALAR_GRID_BASE(const DTYPE dimension, const ATYPE * axis_size):
      GRID_CLASS(dimension,axis_size){ scalar = NULL; };
    // Note: constructor and destructor do not allocate or free scalar memory
    ~SCALAR_GRID_BASE() { scalar = NULL; };

    // copy constructor and assignment: NOT IMPLEMENTED
    SCALAR_GRID_BASE(const SCALAR_GRID_BASE & scalar_grid); 
    const SCALAR_GRID_BASE & operator = (const SCALAR_GRID_BASE & right);

    // set functions
    void Set               ///< Set scalar value of vertex \a iv to \a s.
    (const VTYPE iv, const STYPE s) { scalar[iv] = s; };
    void SetAll(const STYPE s);  ///< Set all scalar values to \a s.
    void CopyScalar        ///< Copy scalar values from \a scalar_grid.
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid);
    void SetCorners        ///< Set scalar values at grid corners to \a s.
    (const STYPE s); 

    /// Set scalar values of vertices in subsampled grid to s.
    void SetSubsample(const ATYPE scale, const STYPE s);


    /// Set scalar values of vertices in region to \a s.
    template <class GTYPE>
    void SetRegion(const IJK::BOX<GTYPE> & box, 
		   const ATYPE scale, const STYPE s);

    /// Set scalar values of vertices in list of regions to \a s.
    template <class GTYPE>
    void SetRegion(const std::vector< IJK::BOX<GTYPE> > & box_list, 
		   const ATYPE scale, const STYPE s);

    /// Replace region from \a iv0 to \a iv1 with values from \a scalar_grid2.
    /// @pre scalar_grid2 has exactly the same dimension and axis_size as the current grid.
    void Replace
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid2,
     const VTYPE iv0, const VTYPE iv1);

    // get functions
    STYPE * ScalarPtr() { return(scalar); };
    const STYPE * ScalarPtrConst() const { return(scalar); };
    const STYPE * End() const { return(scalar+this->NumVertices()); };
    STYPE Scalar(const VTYPE iv) const { return(scalar[iv]); };

    // compute/find functions
    STYPE FindMinScalar() const;
    STYPE FindMaxScalar() const;
    NTYPE CountScalar(const STYPE s) const; 
  };

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID
  // **************************************************

  /// Class storing scalar grid information
  /// Allocates and destroys array scalar[]
  template <class GRID_CLASS, class STYPE>
  class SCALAR_GRID:public SCALAR_GRID_BASE<GRID_CLASS,STYPE> {

  protected:
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

    void Allocate(const DTYPE dimension, const ATYPE * axis_size);
    void FreeAll();

    /// Copy (x0,x1,...,xd) in scalar_grid to p*(x0,x1,...,xd) in current grid where p is supersample_period. 
    template <class STYPE2, class PTYPE>
    void SupersampleCopy
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & scalar_grid,
     const PTYPE supersample_period);

    /// Linear interpolate between supersampled vertices.
    template <class PTYPE>
    void LinearInterpolate(const PTYPE supersample_period);


  public:
    SCALAR_GRID();
    SCALAR_GRID(const DTYPE dimension, const ATYPE * axis_size);
    ~SCALAR_GRID() { FreeAll(); };

    /// Copy scalar grid
    template <class STYPE2>
    void Copy(const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & 
	      scalar_grid2); 
    SCALAR_GRID(const SCALAR_GRID & scalar_grid2) ///< Copy SCALAR_GRID
    { Copy(scalar_grid2); }
    const SCALAR_GRID & operator =                ///< Copy SCALAR_GRID
    (const SCALAR_GRID & right);

    // set functions
    template <class DTYPE2, class ATYPE2>
      void SetSize     ///< Set dimensions and axis sizes.
      (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);


    /// Subsample \a scalar_grid2. Resizes current grid.
    void Subsample   
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid2,
     const ATYPE subsample_period);

    /// Supersample \a scalar_grid2.  Resizes current grid.
    template <class STYPE2, class PTYPE>
    void Supersample
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & scalar_grid2,
     const PTYPE supersample_period);

  };

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_WRAPPER
  // **************************************************

  /// scalar grid wrapper for scalar array
  /// Does not allocate or destroy any scalar information
  template <class GRID_CLASS, class STYPE>
  class SCALAR_GRID_WRAPPER:
    public SCALAR_GRID_BASE<GRID_CLASS,STYPE> {

  protected:
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    SCALAR_GRID_WRAPPER(const DTYPE dimension, const ATYPE * axis_size,
			STYPE * scalar);
    ~SCALAR_GRID_WRAPPER() { this->scalar = NULL; };
    // Note: constructor and destructor do not allocate or free scalar memory
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_BASE
  // **************************************************

  /// Base class for representing min & max of scalar grid regions.
  template <class GRID_CLASS, class STYPE>
  class MINMAX_BASE:public GRID_CLASS {

  protected:
    STYPE * scalar_min;
    STYPE * scalar_max;

    void FreeAll();

    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    MINMAX_BASE();                      ///< Constructor.
    MINMAX_BASE                         /// Constructor.
    (const DTYPE dimension, const ATYPE * axis_size);
    ~MINMAX_BASE(){ FreeAll(); };       ///< Destructor.

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_BASE(const MINMAX_BASE & minmax_grid); 
    const MINMAX_BASE & operator = (const MINMAX_BASE & right);

    // set functions
    template <class DTYPE2, class ATYPE2>
      void SetSize     ///< Set dimensions and axis sizes.
      (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
      void SetSize       ///< Set dimensions and axis size.
      (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // copy min & max values
    void Copy(const STYPE * scalar_min, const STYPE * scalar_max);

    // get functions
    const STYPE * Min() const { return(scalar_min); };
    const STYPE * Max() const { return(scalar_max); };
    STYPE Min(const VTYPE iv) const { return(scalar_min[iv]); };
    STYPE Max(const VTYPE iv) const { return(scalar_max[iv]); };
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_GRID
  // **************************************************

  /// class for computing min & max of scalar grid
  template <class GRID_CLASS, class STYPE>
  class MINMAX_GRID:public MINMAX_BASE<GRID_CLASS,STYPE> {

  protected:
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

  public:
    MINMAX_GRID(){};
    MINMAX_GRID(const DTYPE dimension, const ATYPE * axis_size):
      MINMAX_BASE<GRID_CLASS,STYPE>(dimension, axis_size) {};
    ~MINMAX_GRID(){};

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_GRID(const MINMAX_GRID & minmax_grid); 
    const MINMAX_GRID & operator = (const MINMAX_GRID & right);

    /// Compute min and max of each region and store in primary vertices
    void OverwriteWithMinMax(const ATYPE region_edge_length);

    void OverwriteWithMinMax(const ATYPE region_edge_length,
			     const ATYPE offset_edge_length);

    /// Project scalar grid onto minmax grid
    /// Compute min and max of projection
    void ProjectMinMax
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const ATYPE num_region_edges, const ATYPE facet_index, 
     const ATYPE axis_increment);

    void ProjectMinMax
    (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
     const ATYPE istart, const ATYPE iend, const ATYPE ionto,
     const ATYPE axis_increment);

    void ProjectMinMax
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid,
     const ATYPE num_region_edges, const ATYPE facet_index, 
     const ATYPE axis_increment);
  };

  // **************************************************
  // TEMPLATE CLASS MINMAX_REGIONS
  // **************************************************

  /// class for computing min & max of scalar grid regions
  template <class GRID_CLASS, class STYPE>
  class MINMAX_REGIONS:
    public MINMAX_BASE<GRID_CLASS,STYPE> {

  protected:
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

    ATYPE region_edge_length;

  public:
    MINMAX_REGIONS();
    MINMAX_REGIONS(const DTYPE dimension, const ATYPE * axis_size,
		   const ATYPE region_edge_length);
    ~MINMAX_REGIONS();

    // copy constructor and assignment: NOT IMPLEMENTED
    MINMAX_REGIONS(const MINMAX_REGIONS & minmax_regions); 
    const MINMAX_REGIONS & operator = (const MINMAX_REGIONS & right);

    // get functions
    ATYPE RegionEdgeLength() const { return(region_edge_length); };
    NTYPE NumRegions() const 
    { return(MINMAX_BASE<GRID_CLASS,STYPE>::NumVertices()); };

    NTYPE NumVertices() const;

    // Compute min and max of each region
    // Creates a grid of regions and stores min and max for each region
    void ComputeMinMax
    (const DTYPE dimension, const ATYPE * axis_size,
     const STYPE * scalar, const ATYPE region_edge_length);

    void ComputeMinMax
    (const DTYPE dimension, const ATYPE * axis_size,
     const STYPE * scalar, const ATYPE region_edge_length,
     const ATYPE offset_edge_length);

    void ComputeMinMax
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid, 
     const ATYPE region_edge_length);

    void ComputeMinMax
    (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid, 
     const ATYPE region_edge_length, const ATYPE offset_edge_length);
  };

  // **************************************************
  // LINEAR INTERPOLATION
  // **************************************************
  
  /// Linear interpolate.
  template <class STYPE, class CTYPE0, class CTYPE1, class CTYPE2>
  inline STYPE linear_interpolate
  (const STYPE s0, const CTYPE0 x0, const STYPE s1, const CTYPE1 x1,
   const CTYPE2 x2)
  {
    CTYPE2 w0, w1;
    CTYPE2 x_diff = STYPE(x1 - x0);
    const double EPSILON = 0.00001;
    if (x_diff > EPSILON || x_diff < -EPSILON) { 
      w0 = (x1-x2) / x_diff;
      w1 = (x2-x0) / x_diff;
    }
    else {
      // arbitrarily set normalize weights to 0.5
      w0 = w1 = 0.5;
    };

    STYPE s = s0 * w0 + s1 * w1;
    return(s);
  };

  /// Linear interpolate.  Force use of float in calculations.
  template <class STYPE>
  inline STYPE linear_interpolate
  (const STYPE s0, const int x0, const STYPE s1, const int x1,
   const int x2)
  {
    return(linear_interpolate(s0, float(x0), s1, float(x1), float(x2)));
  }

  /// Linear interpolate.  Force use of double in calculations.
  template <class STYPE>
  inline STYPE linear_interpolate
  (const STYPE s0, const long x0, const STYPE s1, const long x1,
   const long x2)
  {
    return(linear_interpolate(s0, double(x0), s1, double(x1), double(x2)));
  }

  // **************************************************
  // TEMPLATE FUNCTIONS: SORTING GRID VERTICES
  // **************************************************

  /// Comparison function template.
  template <class STYPE> class SCALAR_LESS_THAN {
  protected:
    const STYPE * scalar;

  public:
    SCALAR_LESS_THAN(const STYPE * s):
      scalar(s) {};

    bool operator()(const int iv0, const int iv1) const
    { return(scalar[iv0] < scalar[iv1]); }
  };

  /// Sort grid vertices by increasing scalar value.
  template <class STYPE, class NTYPE, class ITYPE>
  void sort_grid_vertices
  (const STYPE * scalar, const NTYPE num_vertices, ITYPE * index_sorted)
  {
    SCALAR_LESS_THAN<STYPE> scalar_less_than(scalar);

    for (NTYPE i = 0; i < num_vertices; i++)
      { index_sorted[i] = i; };

    std::sort(index_sorted, index_sorted+num_vertices, scalar_less_than);
  }

  /// Find and store locations of elements in list
  template <class ETYPE, class NTYPE, class LTYPE>
  void list_locate(const ETYPE * list, const NTYPE num_elements,
		   LTYPE * list_loc)
  {
    for (NTYPE i = 0; i < num_elements; i++) {
      ETYPE e = list[i];
      list_loc[e] = i;
    }
  }

  // ********************************************************
  // TEMPLATE FUNCTIONS: COMPUTE CUBE/REGION MIN AND MAX
  // ********************************************************

  /// Compute minimum and maximum scalar value of each region
  template <class DTYPE, class ATYPE, class STYPE>
  void compute_region_minmax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length, STYPE * region_min, STYPE * region_max)
  // region_edge_length = number of grid edges per region edge
  // region_min[j] = min scalar value of region j
  // region_max[j] = max scalar value of region j
  //   Precondition: region_min[] and region_max[] are preallocated to size 
  //                 at least number of regions
  {
    IJK::PROCEDURE_ERROR error("compute_region_minmax");
    GRID_SIZE_TYPE axis_increment[dimension];
    GRID_SIZE_TYPE vprimary[dimension];
    std::vector<STYPE> facet_min;
    std::vector<STYPE> facet_max;

    GRID_SIZE_TYPE num_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_vertices);

    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    GRID_SIZE_TYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, num_regions);

    if (num_regions < 1 || num_vertices < 1) { return; };

    region_min[0] = scalar[0];
    region_max[0] = scalar[0];
    
    if (dimension < 1) { return;  };

    MINMAX_GRID< GRID<DTYPE,ATYPE,GRID_SIZE_TYPE,GRID_SIZE_TYPE>,STYPE> 
      minmax_grid(dimension-1, axis_size);

    compute_increment(dimension, axis_size, axis_increment);

    const DTYPE d_last = dimension-1;
    GRID_SIZE_TYPE num_facet_vertices;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d_last, num_facet_vertices);

    ATYPE num_regions_along_axis = 
      compute_num_regions_along_axis(axis_size[d_last], region_edge_length);
    GRID_SIZE_TYPE num_regions_in_facet;
    compute_num_regions
      (dimension-1, axis_size, region_edge_length, num_regions_in_facet);

    GRID_SIZE_TYPE iregion = 0;
    if (dimension > 1) {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, region_edge_length, j*region_edge_length,
	   axis_increment[d_last]);

	minmax_grid.OverwriteWithMinMax(region_edge_length);

	GRID_SIZE_TYPE facet_increment = (j*region_edge_length)*axis_increment[d_last];
	for (DTYPE d = 0; d < dimension; d++) 
	  { vprimary[d] = facet_increment; }

	for (int jregion = 0; jregion < num_regions_in_facet; jregion++) {

	  GRID_SIZE_TYPE iv = vprimary[0]-facet_increment;
	  region_min[iregion] = minmax_grid.Min(iv);
	  region_max[iregion] = minmax_grid.Max(iv);

	  iregion++;

	  DTYPE d = 0;
	  while (d+1 < dimension && 
		 (vprimary[d] + (region_edge_length+1)*axis_increment[d] >=
		  vprimary[d+1] + axis_increment[d+1])) {
	    d++;
	  }

	  // go to next subspace of dimension d
	  vprimary[d] = vprimary[d] + region_edge_length*axis_increment[d];

	  for (DTYPE d2 = 0; d2 < d; d2++) 
	    { vprimary[d2] = vprimary[d]; };
	}
      }
    }
    else {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, region_edge_length, j*region_edge_length,
	   axis_increment[d_last]);

	region_min[iregion] = minmax_grid.Min(0);
	region_max[iregion] = minmax_grid.Max(0);
	iregion++;
      }
    }

    if (iregion != num_regions) {
      error.AddMessage("Programming error.  Computed min value for ", iregion,
		       " regions.");
      error.AddMessage("Number of regions is ", num_regions, ".");
      throw error;
    }
  }

  /// Compute minimum and maximum scalar value of each region
  template <class DTYPE, class ATYPE, class STYPE>
  void compute_region_minmax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length, const ATYPE offset_edge_length,
   STYPE * region_min, STYPE * region_max)
  // region_edge_length = number of grid edges per region edge
  // offset_edge_length = offset edge length
  // region_min[j] = min scalar value of region j
  // region_max[j] = max scalar value of region j
  //   Precondition: region_min[] and region_max[] are preallocated to size 
  //                 at least number of regions
  {
    IJK::PROCEDURE_ERROR error("compute_region_minmax");
    GRID_SIZE_TYPE axis_increment[dimension];
    GRID_SIZE_TYPE vprimary[dimension];
    std::vector<STYPE> facet_min;
    std::vector<STYPE> facet_max;

    GRID_SIZE_TYPE num_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_vertices);

    if (region_edge_length <= 0) {
      error("Programming error.  Number of region edges must be positive.");
      throw error;
    }

    GRID_SIZE_TYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, num_regions);

    if (num_regions < 1 || num_vertices < 1) { return; };

    region_min[0] = scalar[0];
    region_max[0] = scalar[0];
    
    if (dimension < 1) { return;  };

    MINMAX_GRID< GRID<DTYPE,ATYPE,GRID_SIZE_TYPE,GRID_SIZE_TYPE>,STYPE> 
      minmax_grid(dimension-1, axis_size);

    compute_increment(dimension, axis_size, axis_increment);

    const DTYPE d_last = dimension-1;
    GRID_SIZE_TYPE num_facet_vertices;
    compute_num_vertices_in_grid_facet(dimension, axis_size, d_last, num_facet_vertices);

    ATYPE num_regions_along_axis = 
      compute_num_regions_along_axis(axis_size[d_last], region_edge_length);
    GRID_SIZE_TYPE num_regions_in_facet;
    compute_num_regions
      (dimension-1, axis_size, region_edge_length, num_regions_in_facet);

    GRID_SIZE_TYPE iregion = 0;
    if (dimension > 1) {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	ATYPE istart = 0;
	if (j > 0) { istart = j*region_edge_length-offset_edge_length; };
	ATYPE iend = (j+1)*region_edge_length+offset_edge_length+1;
	if (iend > axis_size[d_last]) { iend = axis_size[d_last]; };

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, istart, iend, j*region_edge_length,
	   axis_increment[d_last]);

	minmax_grid.OverwriteWithMinMax(region_edge_length, offset_edge_length);

	GRID_SIZE_TYPE facet_increment = (j*region_edge_length)*axis_increment[d_last];
	for (DTYPE d = 0; d < dimension; d++) 
	  { vprimary[d] = facet_increment; }

	for (int jregion = 0; jregion < num_regions_in_facet; jregion++) {

	  GRID_SIZE_TYPE iv = vprimary[0]-facet_increment;
	  region_min[iregion] = minmax_grid.Min(iv);
	  region_max[iregion] = minmax_grid.Max(iv);

	  iregion++;

	  DTYPE d = 0;
	  while (d+1 < dimension && 
		 (vprimary[d] + (region_edge_length+1)*axis_increment[d] >=
		  vprimary[d+1] + axis_increment[d+1])) {
	    d++;
	  }

	  // go to next subspace of dimension d
	  vprimary[d] = vprimary[d] + region_edge_length*axis_increment[d];

	  for (DTYPE d2 = 0; d2 < d; d2++) 
	    { vprimary[d2] = vprimary[d]; };
	}
      }
    }
    else {
      for (ATYPE j = 0; j < num_regions_along_axis; j++) {

	ATYPE istart = 0;
	if (j > 0) { istart = j*region_edge_length-offset_edge_length; };
	ATYPE iend = (j+1)*region_edge_length+offset_edge_length+1;
	if (iend > axis_size[d_last]) { iend = axis_size[d_last]; };

	minmax_grid.ProjectMinMax
	  (dimension, axis_size, scalar, istart, iend, j*region_edge_length,
	   axis_increment[d_last]);

	region_min[iregion] = minmax_grid.Min(0);
	region_max[iregion] = minmax_grid.Max(0);
	iregion++;
      }
    }

    if (iregion != num_regions) {
      error.AddMessage("Programming error.  Computed min value for ", iregion,
		       " regions.");
      error.AddMessage("Number of regions is ", num_regions, ".");
      throw error;
    }
  }

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID_BASE MEMBER FUNCTIONS
  // **************************************************

  /// Set all scalar values to s.
  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  SetAll(const STYPE s)
  {
    for (VTYPE i = 0; i < this->NumVertices(); i++)
      { scalar[i] = s; }
  }

  /// Copy scalar values of scalar_grid to current grid.
  /// Precondition: scalar_grid has same axis_size as current grid.
  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  CopyScalar
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid)
  {
    std::copy(scalar_grid.ScalarPtrConst(),
	      scalar_grid.ScalarPtrConst() + this->NumVertices(),
	      this->scalar);
  }

  /// Set scalar value at grid corners to s.
  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  SetCorners(const STYPE s)
  {
    const DTYPE dimension = this->Dimension();
    const ATYPE * axis_size = this->AxisSize();

    VTYPE axis_increment[dimension];

    compute_increment(*this, axis_increment);
    NTYPE num_corners = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VTYPE> grid_corner_increment(num_corners);
    compute_region_corner_increment
      (dimension, axis_size, axis_size, grid_corner_increment.Ptr());

    for (VTYPE i = 0; i < num_corners; i++) {
      VTYPE iv = grid_corner_increment[i];
      this->Set(iv, s);
    }
  }

  /// Set scalar values of vertices in subsampled grid to s.
  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  SetSubsample(const ATYPE scale, const STYPE s)
  {
    const DTYPE dimension = this->Dimension();
    const ATYPE * axis_size = this->AxisSize();
    ATYPE coord[dimension];
    ATYPE axis_increment[dimension];

    if (scale <= 1) {
    // set scalar value of all vertices to s
      SetAll(s);
      return;
    }

    if (dimension < 1) { return; };

    compute_increment(dimension, axis_size, axis_increment);

    for (int d = 0; d < dimension; d++) { coord[d] = 0; };

    VTYPE iv = 0;
    bool done = false;
    while (!done) {
      this->Set(iv, s);

      // increment coord and iv
      coord[0] = coord[0] + scale;
      iv = iv + scale*axis_increment[0];

      DTYPE d = 0;
      while (d < dimension && coord[d] >= axis_size[d]) {
	iv = iv - coord[d] * axis_increment[d];
	coord[d] = 0;
	d++;
	if (d < dimension) { 
	  coord[d] = coord[d] + scale; 
	  iv = iv + scale*axis_increment[d];
	}
	else {
	  done = true;
	}
      }
    }

  }

  /// Set scalar values of vertices in region to s.
  template <class GRID_CLASS, class STYPE>
  template <class GTYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  SetRegion(const IJK::BOX<GTYPE> & box,
	    const ATYPE scale, const STYPE s)
  // scale = subsample region vertices at given scale
  {
    const DTYPE dimension = this->Dimension();
    const ATYPE * axis_size = this->AxisSize();
    ATYPE axis_increment[dimension];
    GTYPE coord[dimension];
    GTYPE start_coord[dimension];
    GTYPE end_coord[dimension];
    IJK::PROCEDURE_ERROR error("SCALAR_GRID_BASE::SetRegion");

    NTYPE num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);

    if (num_grid_vertices == 0) { return; };

    compute_increment(dimension, axis_size, axis_increment);

    ATYPE scale2 = scale;
    if (scale2 < 1) { scale2 = 1; };

    if (box.Dimension() != dimension) {
      error.AddMessage("Box dimension ", box.Dimension(),
		       " does not match grid dimension ", dimension, ".");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      GTYPE c = box.MinCoord(d);
      c = (c+scale2-1)/scale2;
      start_coord[d] = c*scale2;
    }

    VTYPE iv = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      coord[d] = start_coord[d];
      end_coord[d] = std::min(box.MaxCoord(d)+1, GTYPE(axis_size[d]));
      iv = iv + coord[d]*axis_increment[d];
    };

    bool done = false;
    for (int d = 0; d < dimension; d++) {
      if (coord[d] >= end_coord[d]) {
	done = true;
      }
    }

    while (!done) {
      this->Set(iv, s);

      // increment coord and iv
      coord[0] = coord[0]+scale2;
      iv = iv + scale2*axis_increment[0];

      DTYPE d = 0;
      while (d < dimension && coord[d] >= end_coord[d]) {
	iv = iv - (coord[d]-start_coord[d]) * axis_increment[d];
	coord[d] = start_coord[d];
	d++;
	if (d < dimension) { 
	  coord[d] = coord[d]+scale2; 
	  iv = iv + scale2*axis_increment[d];
	}
	else {
	  done = true;
	}
      }
    }

  }

  /// Set scalar values of vertices in region to s.
  template <class GRID_CLASS, class STYPE>
  template <class GTYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  SetRegion(const std::vector< IJK::BOX<GTYPE> > & box_list,
	    const ATYPE scale, const STYPE s)
  // scale = subsample region vertices at given scale
  {
    for (NTYPE i = 0; i < box_list.size(); i++) 
      { SetRegion(box_list[i], scale, s); }
  }

  /// Replace region from iv0 to iv1 with values from scalar_grid2.
  // Precondition: scalar_grid2 has exactly the same dimension
  //   and axis_size as the current scalar grid.
  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID_BASE<GRID_CLASS,STYPE>::Replace
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid2, 
   const VTYPE iv0, const VTYPE iv1)
  {
    const DTYPE dimension = this->dimension;
    const ATYPE * axis_size = this->axis_size;
    const STYPE * scalar2 = scalar_grid2.ScalarPtrConst();
    IJK::PROCEDURE_ERROR error("SCALAR_GRID_BASE::Replace");

    if (!Check(scalar_grid2.Dimension(), scalar_grid2.AxisSize(), error))
      { throw error; };

    if (dimension <= 0) { return; }

    ATYPE c0 = iv0%(this->AxisSize(0));
    ATYPE c1 = iv1%(this->AxisSize(0));

    // iv2 is the projection of iv1 onto the hyperplane orthogonal 
    // to axis 0 through iv0
    VTYPE iv2 = iv1 - c1 + c0;
  
    GRID_SIZE_TYPE num_slice_vertices;
    compute_num_grid_vertices(dimension, axis_size, iv0, iv2,
			      num_slice_vertices);

    // vlist contains a region slice orthogonal to the 0 axis
    VTYPE * vlist = new VTYPE[num_slice_vertices];
    get_vertices_between(dimension, axis_size, iv0, iv2, vlist);

    for (ATYPE c = c0; c <= c1; c++) {
      for (VTYPE j = 0; j < num_slice_vertices; j++) {
	VTYPE jv = vlist[j];
	scalar[jv] = scalar2[jv];
      }

      // go to next slice
      for (VTYPE j = 0; j < num_slice_vertices; j++) 
	{ vlist[j]++; }
    }

    delete [] vlist;
  }

  /// Return minimum scalar value
  template <class GRID_CLASS, class STYPE>
  STYPE SCALAR_GRID_BASE<GRID_CLASS,STYPE>::FindMinScalar() const
  {
    STYPE min_scalar = *(min_element(ScalarPtrConst(), End()));
    return(min_scalar);
  }

  /// Return maximum scalar value
  template <class GRID_CLASS, class STYPE>
  STYPE SCALAR_GRID_BASE<GRID_CLASS,STYPE>::FindMaxScalar() const
  {
    STYPE max_scalar = *(max_element(ScalarPtrConst(), End()));
    return(max_scalar);
  }


  /// Count number of vertices with scalar value equal to s.
  template <class GRID_CLASS, class STYPE>
  typename GRID_CLASS::NUMBER_TYPE SCALAR_GRID_BASE<GRID_CLASS,STYPE>::
  CountScalar(const STYPE s) const
  {
    NTYPE num_scalar = 0;
    for (const STYPE * scalar_ptr = ScalarPtrConst(); scalar_ptr != End();
	 scalar_ptr++) {
      if (*scalar_ptr == s) { ++num_scalar; }
    }

    return(num_scalar);
  }

  // **************************************************
  // TEMPLATE CLASS SCALAR_GRID MEMBER FUNCTIONS
  // **************************************************

  template <class GRID_CLASS, class STYPE>
  SCALAR_GRID<GRID_CLASS,STYPE>::SCALAR_GRID():
    SCALAR_GRID_BASE <GRID_CLASS,STYPE> (0, NULL)
  // constructor
  {
    this->scalar = NULL;
  }

  template <class GRID_CLASS, class STYPE>
  SCALAR_GRID<GRID_CLASS,STYPE>::
  SCALAR_GRID
  (const DTYPE dimension, const ATYPE * axis_size) :
    SCALAR_GRID_BASE <GRID_CLASS,STYPE> (dimension, axis_size)
  // constructor
  {
    this->scalar = new STYPE[this->NumVertices()];
  }

  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID<GRID_CLASS,STYPE>::FreeAll()
  {
    delete [] this->scalar;
    this->scalar = NULL;
    GRID_CLASS::FreeAll();
  }

  template <class GRID_CLASS, class STYPE>
  template <class DTYPE2, class ATYPE2>
  void SCALAR_GRID<GRID_CLASS,STYPE>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    GRID_SIZE_TYPE numv = this->NumVertices();
    if (this->scalar == NULL || !CompareSize(dimension, axis_size)) {
      GRID_CLASS::SetSize(dimension, axis_size);
      if (this->scalar == NULL || numv != this->NumVertices()) {
	if (this->scalar != NULL) { delete [] this->scalar; };
	this->scalar = new STYPE[this->NumVertices()];
      }
    }
  }

  template <class GRID_CLASS, class STYPE> 
  template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
  void SCALAR_GRID<GRID_CLASS,STYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  /// *** NOTE USE SCALAR_GRID_BASE
  /// Copy scalar grid
  template <class GRID_CLASS, class STYPE>
  template <class STYPE2>
  void SCALAR_GRID<GRID_CLASS,STYPE>::Copy
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & scalar_grid2)
  {
    this->SetSize(scalar_grid2);
    CopyScalar(scalar_grid2);
  }

  /// *** NOTE USE SCALAR_GRID_BASE
  template <class GRID_CLASS, class STYPE>
  const SCALAR_GRID<GRID_CLASS,STYPE> & 
  SCALAR_GRID<GRID_CLASS,STYPE>::operator = 
  (const SCALAR_GRID<GRID_CLASS,STYPE> & right)
  // copy assignment	
  {
    if (&right != this) {         // avoid self-assignment
      Copy(right);
    }
  }

  template <class GRID_CLASS, class STYPE>
  void SCALAR_GRID<GRID_CLASS,STYPE>::Subsample
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid,
   const ATYPE subsample_period)
  {
    const DTYPE dimension = scalar_grid.Dimension();
    VTYPE axis_increment[dimension];
    VTYPE vprimary[dimension];
    ATYPE subsampled_axis_size[dimension];
    IJK::PROCEDURE_ERROR error("SCALAR_GRID::Subsample");

    for (DTYPE d = 0; d < dimension; d++) {
      subsampled_axis_size[d] = 
	compute_subsample_size(scalar_grid.AxisSize(d), subsample_period);
    }

    SetSize(dimension, subsampled_axis_size);

    if (this->NumVertices() < 1) { return; };

    compute_increment(dimension, scalar_grid.AxisSize(), axis_increment);

    for (DTYPE d = 0; d < dimension; d++) { vprimary[d] = 0; };

    VTYPE iv = 0;
    while (vprimary[0] < scalar_grid.NumVertices()) {
      VTYPE jv = vprimary[0];
      this->scalar[iv] = scalar_grid.Scalar(jv);
      iv++;

      DTYPE d = 0;
      while (d+1 < dimension && 
	     (vprimary[d] + subsample_period*axis_increment[d] >=
	      vprimary[d+1] + axis_increment[d+1])) {
	d++;
      }

      // go to next subspace of dimension d
      vprimary[d] = vprimary[d] + subsample_period*axis_increment[d];

      for (DTYPE d2 = 0; d2 < d; d2++) 
	{ vprimary[d2] = vprimary[d]; };
    }

    if (iv != this->NumVertices()) {
      error.AddMessage("Programming error.  Subsampled ", iv, 
		       " grid vertices.");
      error.AddMessage("Subsampled grid should have ", this->NumVertices(),
		       " vertices.");
      throw error;
    }

  }


  template <class GRID_CLASS, class STYPE>
  template <class STYPE2, class PTYPE>
  void SCALAR_GRID<GRID_CLASS,STYPE>::Supersample
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & scalar_grid2,
   const PTYPE supersample_period)
  {
    const DTYPE dimension = scalar_grid2.Dimension();
    ATYPE supersampled_axis_size[dimension];
    IJK::PROCEDURE_ERROR error("SUPERSAMPLE_GRID::Supersample");

    for (DTYPE d = 0; d < dimension; d++) {
      supersampled_axis_size[d] = 
	compute_supersample_size(scalar_grid2.AxisSize(d), supersample_period);
    }

    SetSize(dimension, supersampled_axis_size);

    if (this->NumVertices() < 1) { return; };

    SupersampleCopy(scalar_grid2, supersample_period);
    LinearInterpolate(supersample_period);
  }

  template <class GRID_CLASS, class STYPE>
  template <class STYPE2, class PTYPE>
  void SCALAR_GRID<GRID_CLASS,STYPE>::SupersampleCopy
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE2> & scalar_grid2,
   const PTYPE supersample_period)
  {
    IJK::CONSTANT<ATYPE,ATYPE> period(supersample_period);
    const DTYPE dimension = this->Dimension();
    ATYPE subgrid_axis_size[dimension];
    VTYPE axis_increment0[dimension];
    VTYPE axis_increment1[dimension];

    NTYPE numv0;
    compute_num_vertices_in_grid_facet
      (scalar_grid2, 0, numv0);


    IJK::ARRAY<VTYPE> vlist0(numv0);
    IJK::ARRAY<VTYPE> vlist1(numv0);
  
    get_vertices_in_grid_facet0(scalar_grid2, vlist0.Ptr());

    std::copy(this->AxisSize(), this->AxisSize()+this->Dimension(), 
	      subgrid_axis_size);
    subgrid_axis_size[0] = 1;

    subsample_subgrid_vertices
      (*this, 0, subgrid_axis_size, period, vlist1.Ptr());

    for (VTYPE x0 = 0; x0 < scalar_grid2.AxisSize(0); x0++) {
      VTYPE x1 = x0*supersample_period;
      for (VTYPE i = 0; i < numv0; i++) {
	VTYPE v0 = vlist0[i]+x0;
	VTYPE v1 = vlist1[i]+x1;
	STYPE s = STYPE(scalar_grid2.Scalar(v0));
	this->scalar[v1] = s;
      }
    }

  }

  template <class GRID_CLASS, class STYPE>
  template <class PTYPE>
  void SCALAR_GRID<GRID_CLASS,STYPE>::LinearInterpolate
  (const PTYPE supersample_period)
  {
    const DTYPE dimension = this->Dimension();
    ATYPE subgrid_axis_size[dimension];
    ATYPE subsample_period[dimension];
    VTYPE axis_increment[dimension];

    compute_increment(*this, axis_increment);

    for (DTYPE d = 0; d < this->Dimension(); d++) {
      for (DTYPE j = 0; j < d; j++) { subsample_period[j] = 1; };
      for (DTYPE j = d; j < this->Dimension(); j++) 
	{ subsample_period[j] = supersample_period; };
      for (DTYPE j = 0; j < this->Dimension(); j++) 
	{ subgrid_axis_size[j] = AxisSize(j); };
      subgrid_axis_size[d] = 1;

      NTYPE numv;
      compute_subsample_size
	(this->Dimension(), subgrid_axis_size, subsample_period, numv);

      IJK::ARRAY<VTYPE> vlist(numv);
      subsample_subgrid_vertices
	(*this, 0, subgrid_axis_size, subsample_period, vlist.Ptr());

      for (VTYPE x = 0; x+1 < this->AxisSize(d); x += supersample_period) {

	VTYPE inc0 = x*axis_increment[d];
	VTYPE inc1 = inc0 + supersample_period*axis_increment[d];
	for (VTYPE i = 0; i < numv; i++) {
	  VTYPE v0 = vlist[i] + inc0;
	  VTYPE v1 = vlist[i] + inc1;

	  VTYPE v2 = v0;
	  for (VTYPE j = 1; j < supersample_period; j++) {
	    v2 += axis_increment[d];
	    STYPE s0 = this->scalar[v0];
	    STYPE s1 = this->scalar[v1];
	    this->scalar[v2] = 
	      linear_interpolate(s0, 0, s1, supersample_period, j);
	  }
	}
      }
    }
  }

  // ******************************************************
  // TEMPLATE CLASS SCALAR_GRID_WRAPPER MEMBER FUNCTIONS
  // ******************************************************

  template <class GRID_CLASS, class STYPE>
  SCALAR_GRID_WRAPPER<GRID_CLASS,STYPE>::
  SCALAR_GRID_WRAPPER
  (const DTYPE dimension, const ATYPE * axis_size, STYPE * scalar):
    SCALAR_GRID_BASE<GRID_CLASS,STYPE>(dimension, axis_size)
  { this->scalar = scalar; };

  // **************************************************
  // TEMPLATE CLASS MINMAX_BASE MEMBER FUNCTIONS
  // **************************************************

  template <class GRID_CLASS, class STYPE>
  MINMAX_BASE<GRID_CLASS,STYPE>::
  MINMAX_BASE()
  // constructor
  {
    scalar_min = NULL;
    scalar_max = NULL;
  }

  template <class GRID_CLASS, class STYPE>
  MINMAX_BASE<GRID_CLASS,STYPE>::
  MINMAX_BASE
  (const DTYPE dimension, const ATYPE * axis_size) :
    GRID_CLASS(dimension, axis_size)
  // constructor
  {
    scalar_min = new STYPE[this->NumVertices()];
    scalar_max = new STYPE[this->NumVertices()];
  }

  template <class GRID_CLASS, class STYPE>
  void MINMAX_BASE<GRID_CLASS,STYPE>::FreeAll()
  {
    delete [] scalar_min;
    delete [] scalar_max;
    scalar_min = NULL;
    scalar_max = NULL;
    GRID_CLASS::FreeAll();
  }

  template <class GRID_CLASS, class STYPE>
  template <class DTYPE2, class ATYPE2>
  void MINMAX_BASE<GRID_CLASS,STYPE>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    GRID_SIZE_TYPE numv = this->NumVertices();
    if (this->scalar_min == NULL || !CompareSize(dimension, axis_size)) {
      GRID_CLASS::SetSize(dimension, axis_size);
      if (this->scalar_min == NULL || numv != this->NumVertices()) {
	if (this->scalar_min != NULL) { delete [] this->scalar_min; };
	if (this->scalar_max != NULL) { delete [] this->scalar_max; };
	this->scalar_min = new STYPE[this->NumVertices()];
	this->scalar_max = new STYPE[this->NumVertices()];
      }
    }
  }

  template <class GRID_CLASS, class STYPE> 
  template <class DTYPE2, class ATYPE2, class VTYPE2, class NTYPE2>
  void MINMAX_BASE<GRID_CLASS,STYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  template <class GRID_CLASS, class STYPE>
  void MINMAX_BASE<GRID_CLASS,STYPE>::
  Copy(const STYPE * scalar_min, const STYPE * scalar_max)
  {
    const GRID_SIZE_TYPE numv = this->NumVertices();
    IJK::PROCEDURE_ERROR error("MINMAX_BASE::Copy");

    if (numv > 0 && (scalar_min == NULL || scalar_max == NULL)) {
      error.AddMessage("Programming error. Memory for scalar_min or scalar_max has not been allocated.");
      throw error;
    }

    std::copy(scalar_min, scalar_min+numv, this->scalar_min);
    std::copy(scalar_max, scalar_max+numv, this->scalar_max);
  }

  // **************************************************
  // TEMPLATE CLASS MINMAX_GRID MEMBER FUNCTIONS
  // **************************************************

  /// Compute minimum and maximum of each region and store in primary vertices
  template <class GRID_CLASS, class STYPE>
  void MINMAX_GRID<GRID_CLASS,STYPE>::
  OverwriteWithMinMax(const ATYPE region_edge_length)
  {
    const DTYPE dimension = this->Dimension();
    ATYPE elength[dimension];
    ATYPE axis_increment[dimension];
    IJK::PROCEDURE_ERROR error("OverwriteWithMinMax");

    if (region_edge_length <= 0) {
      error("Programming error.  Region edge length must be positive.");
      throw error;
    }

    compute_increment(*this, axis_increment);

    for (DTYPE codim = 0; codim < this->Dimension(); codim++) {
      DTYPE d = this->Dimension()-codim-1;

      for (int d2 = 0; d2 < this->dimension; d2++) {
	if (d2 < d) { elength[d2] = 1; }
	else if (d2 > d) { elength[d2] = region_edge_length; }
	else { elength[d2] = this->axis_size[d2]; }
      }

      GRID_SIZE_TYPE subsample_size;
      compute_subsample_size
	(this->dimension, this->axis_size, elength, subsample_size);

      IJK::ARRAY<VTYPE> vlist(subsample_size);
      subsample_grid_vertices(this->dimension, this->axis_size, elength,
			      vlist.Ptr());

      ATYPE num_subsample_vertices_along_axis =
	compute_subsample_size(AxisSize(d), region_edge_length);

      for (ATYPE j = 0; j < num_subsample_vertices_along_axis; j++) {
	ATYPE num_edges = region_edge_length;
	if ((j+1)*region_edge_length +1 > AxisSize(d)) {
	  num_edges = AxisSize(d) - j*region_edge_length - 1;
	}

	VTYPE facet_increment = (j*region_edge_length)*axis_increment[d];
	for (VTYPE k = 0; k < subsample_size; k++) {
	  VTYPE iv0 = vlist[k]+ facet_increment;
	  VTYPE iv1 = iv0;
	  for (ATYPE j2 = 1; j2 <= num_edges; j2++) {
	    iv1 += axis_increment[d];

	    STYPE min1 = this->scalar_min[iv1];
	    STYPE max1 = this->scalar_max[iv1];

	    if (this->scalar_min[iv0] > min1)
	      { this->scalar_min[iv0] = min1; };

	    if (this->scalar_max[iv0] < max1)
	      { this->scalar_max[iv0] = max1; };
	  }
	}
      }
    }
  }

  /// Compute minimum and maximum of each region and store in primary vertices
  template <class GRID_CLASS, class STYPE>
  void MINMAX_GRID<GRID_CLASS,STYPE>::
  OverwriteWithMinMax(const ATYPE region_edge_length,
		      const ATYPE offset_edge_length)
  // Precondition: offset_edge_length < region_edge_length
  {
    const DTYPE dimension = this->Dimension();
    ATYPE elength[dimension];
    ATYPE axis_increment[dimension];
    IJK::PROCEDURE_ERROR error("OverwriteWithMinMax");

    if (region_edge_length <= 0) {
      error("Programming error.  Region edge length must be positive.");
      throw error;
    }

    if (offset_edge_length >= region_edge_length) {
      error.AddMessage("Programming error.  Offset edge length must be less than region edge length.");
      throw error;
    }

    compute_increment(*this, axis_increment);

    for (DTYPE codim = 0; codim < this->Dimension(); codim++) {
      DTYPE d = this->Dimension()-codim-1;

      for (int d2 = 0; d2 < this->dimension; d2++) {
	if (d2 < d) { elength[d2] = 1; }
	else if (d2 > d) { elength[d2] = region_edge_length; }
	else { elength[d2] = this->axis_size[d2]; }
      }

      GRID_SIZE_TYPE subsample_size;
      compute_subsample_size
	(this->dimension, this->axis_size, elength, subsample_size);

      IJK::ARRAY<VTYPE> vlist(subsample_size);
      subsample_grid_vertices(this->dimension, this->axis_size, elength,
			      vlist.Ptr());

      ATYPE num_subsample_vertices_along_axis =
	compute_subsample_size(AxisSize(d), region_edge_length);

      for (ATYPE j = 0; j < num_subsample_vertices_along_axis; j++) {
	ATYPE istart = 0;
	if (j > 0) { istart = j*region_edge_length - offset_edge_length; };
	ATYPE iend = (j+1)*region_edge_length+offset_edge_length+1;
	if (iend > AxisSize(d)) { iend = AxisSize(d); };
	ATYPE ionto = j*region_edge_length;

	for (VTYPE k = 0; k < subsample_size; k++) {
	  VTYPE iv0 = vlist[k]+ ionto*axis_increment[d];
	  VTYPE iv1 = vlist[k]+ istart*axis_increment[d];
	  for (ATYPE j2 = istart; j2 < iend; j2++) {

	    STYPE min1 = this->scalar_min[iv1];
	    STYPE max1 = this->scalar_max[iv1];

	    if (this->scalar_min[iv0] > min1)
	      { this->scalar_min[iv0] = min1; };

	    if (this->scalar_max[iv0] < max1)
	      { this->scalar_max[iv0] = max1; };

	    iv1 += axis_increment[d];
	  }
	}
      }
    }
  }

  /// Project scalar grid onto minmax grid
  /// Compute min and max of projection
  template <class GRID_CLASS, class STYPE>
  void MINMAX_GRID<GRID_CLASS,STYPE>::
  ProjectMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length, const ATYPE facet_index, 
   const ATYPE axis_increment)
  // region_edge_length = number of grid edges per region edge
  // facet_index = index of facet.  Range [0..axis_size-1]
  // axis_size = axis size of axis (dimension-1)
  // axis_increment = axis increment of axis (dimension-1)
  // Precondition: minmax_grid.Dimension()+1 == dimension
  //               minmax_grid.AxisSize(d) == axis_size[d]
  //                   for all d < minmax_grid.Dimension()
  {
    // d_last = last dimension of scalar_grid
    const DTYPE d_last = this->Dimension();

    ATYPE num_edges = region_edge_length;
    ATYPE axis_size_d_last = axis_size[d_last];
    if (facet_index + region_edge_length + 1 > axis_size_d_last) {
      num_edges = axis_size_d_last - facet_index - 1;
    }

    VTYPE facet_increment = facet_index*axis_increment;
    for (VTYPE k = 0; k < this->NumVertices(); k++) {
      VTYPE iv0 = k + facet_increment;
      STYPE s0 = scalar[iv0];
      this->scalar_min[k] = s0;
      this->scalar_max[k] = s0;
      VTYPE iv1 = iv0;
      for (ATYPE j2 = 1; j2 <= num_edges; j2++) {
	iv1 += axis_increment;

	STYPE s1 = scalar[iv1];

	if (this->scalar_min[k] > s1) { this->scalar_min[k] = s1; };
	if (this->scalar_max[k] < s1) { this->scalar_max[k] = s1; };
      }
    }
  }

  /// Project scalar grid onto minmax grid
  /// Compute min and max of projection
  template <class GRID_CLASS, class STYPE>
  void MINMAX_GRID<GRID_CLASS,STYPE>::
  ProjectMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE istart, const ATYPE iend, const ATYPE ionto,
   const ATYPE axis_increment)
  // istart = starting hyperplane
  // iend = ending hyperplane
  // ionto = project onto hyperplane ionto
  // facet_index = index of facet.  Range [0..axis_size-1]
  // axis_size = axis size of axis (dimension-1)
  // axis_increment = axis increment of axis (dimension-1)
  // Precondition: minmax_grid.Dimension()+1 == dimension
  //               minmax_grid.AxisSize(d) == axis_size[d]
  //                   for all d < minmax_grid.Dimension()
  {
    for (VTYPE k = 0; k < this->NumVertices(); k++) {
      VTYPE iv0 = k + ionto*axis_increment;
      STYPE s0 = scalar[iv0];
      this->scalar_min[k] = s0;
      this->scalar_max[k] = s0;
      VTYPE iv1 = k + istart*axis_increment;
      for (ATYPE j2 = istart; j2 < iend; j2++) {
	STYPE s1 = scalar[iv1];

	if (this->scalar_min[k] > s1) { this->scalar_min[k] = s1; };
	if (this->scalar_max[k] < s1) { this->scalar_max[k] = s1; };

	iv1 += axis_increment;
      }
    }
  }

  template <class GRID_CLASS, class STYPE>
  void MINMAX_GRID<GRID_CLASS,STYPE>::
  ProjectMinMax
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid,
   const ATYPE region_edge_length, const ATYPE facet_index, 
   const ATYPE axis_increment)
  {
    ProjectMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
		  scalar_grid.ScalarPtrConst(), region_edge_length,
		  facet_index, axis_increment);
  }



  // **************************************************
  // TEMPLATE CLASS MINMAX_REGIONS MEMBER FUNCTIONS
  // **************************************************

  template <class GRID_CLASS, class STYPE>
  MINMAX_REGIONS<GRID_CLASS,STYPE>::
  MINMAX_REGIONS()
  // constructor
  {
    region_edge_length = 0;
  }

  template <class GRID_CLASS, class STYPE>
  MINMAX_REGIONS<GRID_CLASS,STYPE>::
  MINMAX_REGIONS
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length) :
    MINMAX_BASE<GRID_CLASS,STYPE>(dimension, axis_size)
  // constructor
  {
    this->region_edge_length = region_edge_length;
  }

  template <class GRID_CLASS, class STYPE> 
  MINMAX_REGIONS<GRID_CLASS,STYPE>::
  ~MINMAX_REGIONS<GRID_CLASS,STYPE>()
  // destructor
  {
    region_edge_length = 0;
  }

  /// Compute min and max of each region
  /// Creates a grid of regions and stores min and max for each region
  template <class GRID_CLASS, class STYPE>
  void MINMAX_REGIONS<GRID_CLASS,STYPE>::ComputeMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length)
  {
    ATYPE num_regions_along_axis[dimension];
    IJK::PROCEDURE_ERROR error("MINMAX_REGIONS::ComputeMinMax");

    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    this->region_edge_length = region_edge_length;

    for (DTYPE d = 0; d < dimension; d++) {
      num_regions_along_axis[d] = 
	compute_num_regions_along_axis(axis_size[d], region_edge_length);
    }

    SetSize(dimension, num_regions_along_axis);

    compute_region_minmax
      (dimension, axis_size, scalar, region_edge_length, 
       this->scalar_min, this->scalar_max);
  }

  /// Compute min and max of each region
  /// Creates a grid of regions and stores min and max for each region
  template <class GRID_CLASS, class STYPE>
  void MINMAX_REGIONS<GRID_CLASS,STYPE>::ComputeMinMax
  (const DTYPE dimension, const ATYPE * axis_size, const STYPE * scalar,
   const ATYPE region_edge_length, const ATYPE offset_edge_length)
  {
    ATYPE num_regions_along_axis[dimension];

    this->region_edge_length = region_edge_length;

    for (DTYPE d = 0; d < dimension; d++) {
      num_regions_along_axis[d] = 
	compute_num_regions_along_axis(axis_size[d], region_edge_length);
    }

    SetSize(dimension, num_regions_along_axis);

    compute_region_minmax
      (dimension, axis_size, scalar, region_edge_length, offset_edge_length,
       this->scalar_min, this->scalar_max);
  }

  template <class GRID_CLASS, class STYPE>
  void MINMAX_REGIONS<GRID_CLASS,STYPE>::ComputeMinMax
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid, 
   const ATYPE region_edge_length)
  {
    ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
			scalar_grid.ScalarPtrConst(), region_edge_length);
  }

  template <class GRID_CLASS, class STYPE>
  void MINMAX_REGIONS<GRID_CLASS,STYPE>::ComputeMinMax
  (const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid, 
   const ATYPE region_edge_length, const ATYPE offset_edge_length)
  {
    ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
			scalar_grid.ScalarPtrConst(), region_edge_length,
			offset_edge_length);
  }

  // **************************************************
  // TEMPLATE OUTPUT FUNCTIONS
  // **************************************************

  /// Output scalar.  Template specialization for bool and char classes.
  template <class STYPE>
  inline void output_scalar(std::ostream & out, const STYPE s)
  { out << s; }

  /// Output scalar.  Template specialization for bool and char classes.
  template <>
  inline void output_scalar(std::ostream & out, const bool s)
  { out << int(s); }

  /// Output scalar.  Template specialization for bool and char classes.
  template <>
  inline void output_scalar(std::ostream & out, const char s)
  { out << int(s); }

  /// Output scalar.  Template specialization for bool and char classes.
  template <>
  inline void output_scalar(std::ostream & out, const unsigned char s)
  { out << int(s); }

  /// Output scalar grid (for debugging purposes)
  template <class GRID_CLASS, class STYPE>
  void output_scalar_grid
  (std::ostream & out, 
   const SCALAR_GRID_BASE<GRID_CLASS,STYPE> & scalar_grid)
  {
    typedef typename GRID_CLASS::DIMENSION_TYPE DTYPE;
    typedef typename GRID_CLASS::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_CLASS::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_CLASS::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const ATYPE * axis_size = scalar_grid.AxisSize();

    using namespace std;
    
    if (scalar_grid.NumVertices() < 1) { return; };

    if (scalar_grid.Dimension() <= 1) {
      output_scalar(out, scalar_grid.Scalar(0));
      for (NTYPE i = 0; i < scalar_grid.NumVertices(); ++i) {
	out << " ";
	output_scalar(out, scalar_grid.Scalar(i));
      }
      out << endl;
    }
    else {
      NTYPE num_vertices_in_ridge =
	compute_num_vertices_in_grid_ridge(dimension, axis_size, 0, 1);
      IJK::ARRAY<VTYPE> vlist(num_vertices_in_ridge);
      get_vertices_in_grid_ridge(dimension, axis_size, 0, 1, vlist.Ptr());

      for (NTYPE i = 0; i < num_vertices_in_ridge; i++) {
	VTYPE iv = vlist[i];
	for (ATYPE j1 = 0; j1 < axis_size[1]; j1++) {
	  output_scalar(out, scalar_grid.Scalar(iv));
	  ++iv;
	  for (ATYPE j0 = 1; j0 < axis_size[0]; j0++) {
	    out << " ";
	    output_scalar(out, scalar_grid.Scalar(iv));
	    ++iv;
	  }
	  out << endl;
	}
	out << endl;
      }
    }
  }

}

#endif
