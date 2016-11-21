/// \file ijkmcube_extract.cxx
/// Subroutines for extracting isosurface mesh

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

#include "ijkmcube_extract.h"
#include "ijkmcube_util.h"
#include "ijkoctree.h"
#include "ijktable.h"

#include "ijkcoord.txx"
#include "ijkinterpolate.txx"
#include "ijkpoly.txx"
#include "ijkisopoly.txx"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKTABLE;



// **************************************************
// POLY TYPES
// **************************************************

typedef ORIENTED_CUBE<int, int, VERTEX_INDEX> MC_ORIENTED_CUBE;
typedef ORIENTED_CELL<int, int, VERTEX_INDEX> MC_ORIENTED_CELL;
typedef std::vector<MC_ORIENTED_CUBE> MC_ORIENTED_CUBE_LIST;

// **************************************************
// CLASS VERTEX INCREMENT
// **************************************************

namespace {

  class VERTEX_INCREMENT {

  protected:
    int dimension;
    VERTEX_INDEX * axis_size;
    VERTEX_INDEX * axis;
    VERTEX_INDEX * cube;
    VERTEX_INDEX num_cube_vertices;
    VERTEX_INDEX * iso;
    ISO_VERTEX_INDEX num_iso_vertices;
    int num_isov_per_gridv;

  public:
    VERTEX_INCREMENT(const AXIS_SIZE_TYPE * axis_size,
		     const ISOSURFACE_TABLE & isotable,
		     const ISOTABLE_TYPE isotable_type);
    ~VERTEX_INCREMENT();

    // get functions
    int Dimension() const { return(dimension); };
    VERTEX_INDEX NumCubeVertices() const { return(num_cube_vertices); };
    ISO_VERTEX_INDEX NumIsoVertices() const { return(num_iso_vertices); };
    int NumIsoVPerGridV() const { return(num_isov_per_gridv); };
    const AXIS_SIZE_TYPE * AxisSize() const { return(axis_size); };
    VERTEX_INDEX AxisSize(const int d) const { return(axis_size[d]); };
    VERTEX_INDEX Axis(const int d) const { return(axis[d]); };
    VERTEX_INDEX Cube(const int i) const { return(cube[i]); };
    VERTEX_INDEX Iso(const int i) const { return(iso[i]); };
    const VERTEX_INDEX * AxisSizePtrConst() const { return(axis_size); };
    const VERTEX_INDEX * CubePtrConst() const { return(cube); };
  };

  VERTEX_INCREMENT::VERTEX_INCREMENT
  (const AXIS_SIZE_TYPE * axis_size, 
   const ISOSURFACE_TABLE & isotable, 
   const ISOTABLE_TYPE isotable_type)
  // constructor
  {
    PROCEDURE_ERROR error("VERTEX_INCREMENT constructor");

    this->axis_size = NULL;
    axis = NULL;
    cube = NULL;
    iso = NULL;

    assert(isotable.Polyhedron().NumVertices() == 
	   compute_num_cube_vertices(isotable.Dimension()));

    dimension = isotable.Dimension();
    num_cube_vertices = isotable.Polyhedron().NumVertices();
    num_iso_vertices = isotable.NumIsosurfaceVertices();

    this->axis_size = new VERTEX_INDEX[dimension];
    axis = new VERTEX_INDEX[dimension];
    cube = new VERTEX_INDEX[num_cube_vertices];
    iso = new VERTEX_INDEX[num_iso_vertices];

    for (int d = 0; d < dimension; d++) {
      this->axis_size[d] = axis_size[d];
    };

    compute_increment(dimension, axis_size, axis);
    compute_cube_vertex_increment(dimension, axis, cube);

    switch(isotable_type) {
    case BINARY:
      num_isov_per_gridv = get_num_iso_vertices_per_grid_vertex(dimension);
      compute_iso_vertex_increment(isotable, cube, iso);
      break;

    case NEP:
      num_isov_per_gridv = get_num_nep_iso_vertices_per_grid_vertex(dimension);
      compute_nep_iso_vertex_increment(isotable, cube, iso);
      break;

    case IVOL:
      num_isov_per_gridv = get_num_ivol_vertices_per_grid_vertex(dimension);
      compute_ivol_vertex_increment(isotable, cube, iso);
      break;

    default:
      error.AddMessage("Programming error.  Illegal isosurface table type.");
      throw error;
      break;
    }

  }

  VERTEX_INCREMENT::~VERTEX_INCREMENT()
  // destructor
  {
    dimension = 0;
    num_cube_vertices = 0;
    num_iso_vertices = 0;
    num_isov_per_gridv = 0;

    delete [] axis_size;
    delete [] axis;
    delete [] cube;
    delete [] iso;

    axis_size = NULL;
    axis = NULL;
    cube = NULL;
    iso = NULL;
  }
    
}


// **************************************************
// CLASS FACET VERTEX INCREMENT
// **************************************************

namespace {

  class FACET_VERTEX_INCREMENT:public VERTEX_INCREMENT {

  protected:
    int orth_axis;    // axis orthogonal to the facet
    int side;         // indicate side of cube
    // side = 0. Facet orthogonal to orth_axis with min orth_axis coordinates.
    // side = 1. Facet orthogonal to orth_axis with max orth_axis coordinates.

    VERTEX_INDEX * facet_iso;
    VERTEX_INDEX * facet_vertex;
    VERTEX_INDEX * cube_vertex_index;

    void ComputeFacetIso(const ISOSURFACE_TABLE & isotable,
			 const int num_cube_vertices);

  public:
    FACET_VERTEX_INCREMENT
    (const AXIS_SIZE_TYPE * axis_size, const ISOSURFACE_TABLE & isotable,
     const ISOTABLE_TYPE isotable_type, const int orth_axis, const int side);
    ~FACET_VERTEX_INCREMENT();

    // get functions
    VERTEX_INDEX FacetIso(const int i) const { return(facet_iso[i]); };
    VERTEX_INDEX FacetVertex(const int i) const { return(facet_vertex[i]); };
    VERTEX_INDEX CubeVertexIndex(const int i) const
    { return(cube_vertex_index[i]); };
    int OrthogonalAxis() const { return(orth_axis); };
    int Side() const { return(side); };

    // disable member function Iso()
    VERTEX_INDEX Iso(const int i) const;
  };

  FACET_VERTEX_INCREMENT::FACET_VERTEX_INCREMENT
    (const AXIS_SIZE_TYPE * axis_size, const ISOSURFACE_TABLE & isotable,
     const ISOTABLE_TYPE isotable_type, const int orth_axis, const int side):
    VERTEX_INCREMENT(axis_size, isotable, isotable_type)
  {
    const int dimension = isotable.Dimension();
    PROCEDURE_ERROR error("FACET_VERTEX_INCREMENT constructor");
    const ISO_VERTEX_INDEX num_iso_vertices = isotable.NumIsosurfaceVertices();

    if (isotable_type != NEP) {
      error.AddMessage("Programming error. FACET_VERTEX_INCREMENT can only be used with NEP isotable.");
      throw error;
    }

    const GRID_SIZE_TYPE num_cube_vertices = 
      compute_num_cube_vertices(dimension);
    const GRID_SIZE_TYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);

    this->orth_axis = orth_axis;
    this->side = side;

    facet_iso = new VERTEX_INDEX[num_iso_vertices];
    facet_vertex = new VERTEX_INDEX[num_facet_vertices];
    cube_vertex_index = new VERTEX_INDEX[num_facet_vertices];

    compute_facet_vertex_increment
      (isotable.Polyhedron(), orth_axis, side, axis_size,
       facet_vertex, cube_vertex_index);

    ComputeFacetIso(isotable, num_cube_vertices);
  }


  FACET_VERTEX_INCREMENT::~FACET_VERTEX_INCREMENT()
  {
    delete [] facet_iso;
    facet_iso = NULL;
    delete [] facet_vertex;
    facet_vertex = NULL;
    delete [] cube_vertex_index;
    cube_vertex_index = NULL;
  };

  void FACET_VERTEX_INCREMENT::ComputeFacetIso
  (const ISOSURFACE_TABLE & isotable, const int num_cube_vertices)
  {
    const int dimension = isotable.Dimension();
    VERTEX_INDEX cube_increment2[num_cube_vertices];

    const GRID_SIZE_TYPE num_facet_vertices =
      compute_num_cube_facet_vertices(dimension);

    for (int i = 0; i < num_cube_vertices; i++) 
      { cube_increment2[i] = 0; };

    for (int j = 0; j < num_facet_vertices; j++) {
      int iv = cube_vertex_index[j];
      cube_increment2[iv] = cube[iv];

      if (side == 0)
	{ cube_increment2[iv] -= axis[orth_axis]; };

      cube_increment2[iv] = cube_increment2[iv]*NumIsoVPerGridV()+dimension;
    }

    for (int k = 0; k < isotable.NumIsosurfaceVertices(); k++) {
      facet_iso[k] = 0;
      if (isotable.IsosurfaceVertex(k).Type() == ISOSURFACE_VERTEX::VERTEX) {
	VERTEX_INDEX iv = isotable.IsosurfaceVertex(k).Face();
	facet_iso[k] = cube_increment2[iv];
      }
    }
  }
}

// **************************************************
// CLASS NEP TABLE INDEX INCREMENT
// **************************************************

namespace {

  class NEP_TABLE_INDEX_INCREMENT {

  protected:
    int dimension;
    TABLE_INDEX * positive;
    TABLE_INDEX * equals;
    VERTEX_INDEX num_cube_vertices;
  
  public:
    NEP_TABLE_INDEX_INCREMENT(const int dimension);
    ~NEP_TABLE_INDEX_INCREMENT();

    // get functions
    int Dimension() const { return(dimension); };
    TABLE_INDEX Positive(const int i) const 
    { return(positive[i]); };
    TABLE_INDEX Equals(const int i) const 
    { return(equals[i]); };
  };

  NEP_TABLE_INDEX_INCREMENT::NEP_TABLE_INDEX_INCREMENT
  (const int dimension)
  {
    positive = NULL;
    equals = NULL;

    this->dimension = dimension;
    num_cube_vertices = compute_num_cube_vertices(dimension);
    
    positive = new TABLE_INDEX[num_cube_vertices];
    equals = new TABLE_INDEX[num_cube_vertices];

    TABLE_INDEX x = 1;
    for (VERTEX_INDEX j = 0; j < num_cube_vertices; j++) {
      equals[j] = x;
      positive[j] = 2*x;
      x = 3*x;
    }
  }

  NEP_TABLE_INDEX_INCREMENT::~NEP_TABLE_INDEX_INCREMENT()
  {
    delete [] positive;
    delete [] equals;
    positive = NULL;
    equals = NULL;
    dimension = 0;
    num_cube_vertices = 0;
  }
}


// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

namespace {
  // type definition
  typedef ARRAY<VERTEX_INDEX> VERTEX_ARRAY;
}

// local extraction routines
namespace {

  /// Average scalar values at vertices of cube
  inline SCALAR_TYPE average_scalar
  (const SCALAR_TYPE * scalar, const VERTEX_INDEX iv0, 
   const VERTEX_INDEX * increment, const int num_cube_vertices)
    /// Precondition: num_cube_vertices > 0
  {
    SCALAR_TYPE s = 0;
    for (int j = 0; j < num_cube_vertices; j++) {
      VERTEX_INDEX iv1 = iv0 + increment[j];
      s += scalar[iv1];
    };
    s = s/num_cube_vertices;
    
    return(s);
  }

  /// Average scalar values 
  inline SCALAR_TYPE average_scalar
  (const SCALAR_TYPE * scalar, const int num_vertices)
    /// Precondition: num_vertices > 0
  {
    SCALAR_TYPE s = 0;
    for (int j = 0; j < num_vertices; j++) 
      { s += scalar[j]; };
    s = s/num_vertices;
    
    return(s);
  }

  /// Compute scalar value of point deciding multilinear topology in cube.
  SCALAR_TYPE compute_scalar_of_cube_decider
  (const int dimension, const SCALAR_TYPE * scalar,
   const VERTEX_INDEX num_cube_vertices, const int * parity_sign)
  {
    SCALAR_TYPE b[dimension];
    SCALAR_TYPE x[dimension];

    SCALAR_TYPE a = 
      compute_multilinear_a<SCALAR_TYPE>(dimension, scalar, parity_sign);

    const double EPSILON = 0.00001;
    if (-EPSILON < a && a < EPSILON) {
      SCALAR_TYPE s = average_scalar(scalar, num_cube_vertices);
      return(s);
    }
    else {
      for (int d =  0; d < dimension; d++) {
	b[d] = compute_multilinear_b<SCALAR_TYPE>(dimension, scalar, parity_sign, d);
	x[d] = -b[d]/a;

	if (x[d] < 0.0 || x[d] > 1.0) {
	  SCALAR_TYPE s = average_scalar(scalar, num_cube_vertices);
	  return(s);
	}
      }

      SCALAR_TYPE s = 0.0;
      for (VERTEX_INDEX i = 0; i < num_cube_vertices; i++) {
	SCALAR_TYPE weight = 1.0;
	long mask = 1L;
	for (int d = 0; d < dimension; d++) {
	  if ((i & mask) == 0) { weight = weight*(1-x[d]); }
	  else { weight = weight*x[d]; }
	  mask = (mask << 1L);
	}
	s += weight*scalar[i];
      }

      return(s);
    }
  }


  /// Return true if decider value is greater than or equal to the isovalue.
  /// Decider value is scalar value of point deciding multilinear topology in cube.
  bool is_cube_decider_ge_isovalue
  (const int dimension, const SCALAR_TYPE * scalar,
   const VERTEX_INDEX num_cube_vertices, const SCALAR_TYPE isovalue,
   const int * parity_sign)
  {
    SCALAR_TYPE s = 
      compute_scalar_of_cube_decider
      (dimension, scalar, num_cube_vertices, parity_sign);

    if (s < isovalue) { return(false); }
    else { return(true); }
  }


  /// Compute signs of facet deciders.
  /// 0 = facet deciders not ambiguous.
  void compute_signs_of_facet_deciders
  (const MC_SCALAR_GRID_BASE & scalar_grid,
   const ISOSURFACE_TABLE & cube_isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign, int * decider_sign,
   int & num_positive, int & num_negative)
  {
  
    //cout<<"In facet decider!\n";
    const int dimension = cube_isotable.Dimension();
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();
    const int num_cube_facets = cube_isotable.Polyhedron().NumFacets();
    const int num_facet_vertices = num_cube_vertices/2;
    SCALAR_TYPE base_scalar[num_facet_vertices];

    num_positive = 0;
    num_negative = 0;

    for (int jf = 0; jf < num_cube_facets; jf++) {

      if (cube_ambig.IsFacetAmbiguous(it, jf)) {

	for (int k = 0; k < num_facet_vertices; k++) {
	  VERTEX_INDEX kv = cube_isotable.Polyhedron().FacetVertex(jf, k);
	  VERTEX_INDEX iv1 = iv0 + increment.Cube(kv);
	  base_scalar[k] = scalar_grid.Scalar(iv1);
	}

	bool is_ge_isovalue = 
	  is_cube_decider_ge_isovalue
	  (dimension-1, base_scalar, num_facet_vertices, isovalue,
	   parity_sign);

	if (is_ge_isovalue) { 
	  decider_sign[jf] = +1; 
	  num_positive++;
	}
	else { 
	  decider_sign[jf] = -1; 
	  num_negative++;
	};
      }
      else {
	decider_sign[jf] = 0;
      }
    }
  }


  /// *** SHOULD BE IN utility functions
  /// Return true if the point is in the unit cube.
  bool is_point_in_unit_cube
  (const int dimension, const COORD_TYPE * coord)
  {
    for (int d = 0; d < dimension; d++) {
      if (coord[d] < 0 || coord[d] > 1) 
	{ return(false); }
    }
    return(true);
  }

  /// Get scalar values of polyhedron vertices.
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

  /// Add isosurface simplex vertices.
  inline void add_iso_simplex_vertices_in_cube
  (const ISOSURFACE_TABLE & isotable, const TABLE_INDEX it,
   const ISO_VERTEX_INDEX isov0, const int is, 
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices)
  {
    for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
      int jv = isotable.SimplexVertex(it, is, j);
      ISO_VERTEX_INDEX isov = isov0 + increment.Iso(jv);
      iso_simplices.push_back(isov);
    };
  }

  /// Add isosurface simplices.
  inline void add_iso_simplices_in_cube
  (const ISOSURFACE_TABLE & isotable, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    ISO_VERTEX_INDEX isov0 = iv0*increment.NumIsoVPerGridV();
    for (int is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices_in_cube
	(isotable, it, isov0, is, increment, iso_simplices);
    };
  }

  /// Add isosurface simplices, reverse orientation
  inline void add_iso_simplices_in_cube_reverse_orient
  (const ISOSURFACE_TABLE & isotable, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    ISO_VERTEX_INDEX isov0 = iv0*increment.NumIsoVPerGridV();
    for (int is = 0; is < num_simplices; is++) {

      add_iso_simplex_vertices_in_cube
	(isotable, it, isov0, is, increment, iso_simplices);

      // reverse orientation by swapping last two vertices
      int ilast = iso_simplices.size()-1;
      std::swap(iso_simplices[ilast], iso_simplices[ilast-1]);
    };
  }

  /// Extract isosurface simplices in cube
  /// Note: Make this inline for faster execution
  inline void extract_from_cube
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
    // extract isosurface simplices in cube with primary vertex iv0
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment = vertex increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = 
    //     grid edge containing k'th vertex of simplex is.
    // num_simplices = number of simplices extracted
  {
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, iv0, increment.CubePtrConst(),
			num_cube_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, iv0, increment.CubePtrConst(), 
			     num_cube_vertices, it);

      add_iso_simplices_in_cube(isotable, it, iv0, increment, 
				iso_simplices, num_simplices);
    }
  }


  /// Extract isosurface simplices from ambiguous cube configuration
  /// Use average scalar values of vertices to resolve ambiguity
  void extract_from_ambig_cube_using_average
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
  {
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();

    SCALAR_TYPE s = average_scalar(scalar, iv0, increment.CubePtrConst(),
				   num_cube_vertices);

    if (s >= isovalue) {
      add_iso_simplices_in_cube(isotable, it, iv0, increment, 
				iso_simplices, num_simplices);
    }
    else {
      TABLE_INDEX it_complement = isotable.NumTableEntries() - it - 1;

      add_iso_simplices_in_cube_reverse_orient
	(isotable, it_complement, iv0, increment, 
	 iso_simplices, num_simplices);
    }
  }

  /// Extract isosurface simplices from ambiguous cube configuration
  /// Use average scalar values of vertices to resolve ambiguity
  void extract_from_ambig_cube_using_average
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices)
  {
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();

    SCALAR_TYPE s = average_scalar(scalar, iv0, increment.CubePtrConst(),
				   num_cube_vertices);

    if (s >= isovalue) {
      add_iso_simplices(isotable, it, iv0, increment.CubePtrConst(), 
				iso_simplices, endpoint, num_simplices);
    }
    else {
      TABLE_INDEX it_complement = isotable.NumTableEntries() - it - 1;

      add_reverse_orient_iso_simplices
	(isotable, it_complement, iv0, increment.CubePtrConst(), 
	 iso_simplices, endpoint, num_simplices);
    }
  }


  /// Extract isosurface simplices from ambiguous cube configuration
  /// @param ambiguity_sign = If true, use table entry it.  Otherwise, use complement of table entry it.
  void extract_from_ambig_cube
  (const ISOSURFACE_TABLE & isotable, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INDEX * increment,
   const bool ambiguity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices)
  {
    if (ambiguity_sign) {
      add_iso_simplices(isotable, it, iv0, increment,
			iso_simplices, endpoint, num_simplices);
    }
    else {
      TABLE_INDEX it_complement = isotable.NumTableEntries() - it - 1;

      add_reverse_orient_iso_simplices
	(isotable, it_complement, iv0, increment,
	 iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract isosurface simplices from ambiguous cube configuration
  void extract_from_ambig_cube
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE_AMBIG & cube_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices)
  {
    const int dimension = cube_isotable.Dimension();
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();
    SCALAR_TYPE cube_scalar[num_cube_vertices];

    get_poly_scalar(scalar, iv0, increment.CubePtrConst(), 
		    num_cube_vertices, cube_scalar);

    bool decider_sign = 
      is_cube_decider_ge_isovalue
      (dimension, cube_scalar, num_cube_vertices, isovalue, parity_sign);

    extract_from_ambig_cube
      (cube_isotable, it, iv0, increment.CubePtrConst(), decider_sign, 
       iso_simplices, endpoint, num_simplices);
  }

  /// Extract isosurface simplices from ambiguous polyhedron configuration
  void extract_from_ambig_poly
  (const ISOSURFACE_TABLE & isotable, const TABLE_INDEX it,
   const bool ambiguity_sign,
   const VERTEX_INDEX * vertex_index, const bool orientation,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices)
  {
    if (ambiguity_sign) {
      add_iso_simplices(isotable, it, 0, vertex_index, orientation,
			iso_simplices, endpoint, num_simplices);
    }
    else {
      TABLE_INDEX it_complement = isotable.NumTableEntries() - it - 1;

      bool orientation2 = !orientation;
      add_iso_simplices(isotable, it_complement, 0, vertex_index, 
			orientation2, iso_simplices, endpoint, num_simplices);
    }
  }

  /// Compute list of facets.
  void compute_facet_list
  (const MC_ORIENTED_CUBE & cube, MC_ORIENTED_CUBE_LIST & facet_list)
  {
    const int dimension = cube.Dimension();
    const int num_cube_vertices = cube.NumVertices();
    const bool cube_orientation = cube.Orientation();

    facet_list.clear();

    if (num_cube_vertices < 1) { return; };

    if (dimension > 0) {
      const int num_cube_facet_vertices = num_cube_vertices/2;

      for (int ifacet = 0; ifacet < 2*dimension; ifacet++) {

	MC_ORIENTED_CUBE facet(cube, ifacet);
	facet_list.push_back(facet);
      }
    }
  }

  bool compute_pyramid_base_center_sign
  (const SCALAR_TYPE * scalar, const VERTEX_INDEX * pyramid_index,
   const int num_pyramid_base_vertices, const SCALAR_TYPE isovalue)
  {
    SCALAR_TYPE s = 0;
    for (int j = 0; j < num_pyramid_base_vertices; j++) {
      s += scalar[pyramid_index[j]];
    };
    s = s/num_pyramid_base_vertices;

    if (s < isovalue) { return(false); }
    else { return(true); };
  }

  void extract_from_pyramids_in_cube
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const COORD_TYPE * apex_coord, const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info);

  /// Split cube into pyramids and extract isosurface simplices.
  void extract_from_pyramids_in_cube
  (const SCALAR_TYPE * scalar, const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INDEX new_grid_vertex,
   const SCALAR_TYPE new_scalar,
   const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices, SCALAR_INFO & ambiguous_info);

  inline void add_vertex_at_cube_center
  (const MC_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv0, 
   const VERTEX_INDEX * increment, const int num_cube_vertices,
   MC_MESH_VERTEX_LIST & new_mesh_vertices)
  {
    const int dimension = scalar_grid.Dimension();
    COORD_TYPE coord[dimension];

    scalar_grid.ComputeCoord(iv0, coord);

    VERTEX_INDEX k = new_mesh_vertices.NumVertices();

    new_mesh_vertices.Add();
    for (int d = 0; d < new_mesh_vertices.Dimension(); d++) 
      { new_mesh_vertices.SetCoord(k, d, coord[d]+0.5); };

    SCALAR_TYPE s = average_scalar
      (scalar_grid.ScalarPtrConst(), iv0, increment, num_cube_vertices);

    new_mesh_vertices.SetScalar(k, s);
  }

  inline void add_mesh_vertex
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const COORD_TYPE * offset_coord, const VERTEX_INDEX iv0, 
   const VERTEX_INDEX * increment, const int num_cube_vertices,
   MC_MESH_VERTEX_LIST & new_mesh_vertices)
  {
    const int dimension = scalar_grid.Dimension();
    SCALAR_TYPE cube_scalar[num_cube_vertices];
    COORD_TYPE coord0[dimension];
    COORD_TYPE coord2[dimension];

    get_poly_scalar(scalar_grid.ScalarPtrConst(), iv0, 
		    increment, num_cube_vertices,
		    cube_scalar);

    scalar_grid.ComputeCoord(iv0, coord0);

    VERTEX_INDEX k = new_mesh_vertices.NumVertices();

    new_mesh_vertices.Add();

    add_coord(dimension, coord0, offset_coord, coord2);
    new_mesh_vertices.SetCoord(k, coord2);

    SCALAR_TYPE s = 
      multilinear_interpolate_scalar
      (dimension, offset_coord, num_cube_vertices, cube_scalar);

    new_mesh_vertices.SetScalar(k, s);
  }


  /// Extract isosurface simplices from ambiguous cube configuration
  void extract_from_ambig_cube_adecider
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    const int num_cube_facets = 
      poly_isotable.cube.Polyhedron().NumFacets();
    int decider_sign[num_cube_facets];

    if (poly_isotable.ambig_cube.NumAmbiguousFacets(it) == 0) {
      SCALAR_TYPE s = average_scalar(scalar, iv0, increment.CubePtrConst(),
				     num_cube_vertices);
      bool ambiguity_sign = false;
      if (s >= isovalue) { ambiguity_sign = true; };
      extract_from_ambig_cube
	(poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	 ambiguity_sign, iso_simplices, endpoint, num_simplices);
    }
    else {
      int num_positive = 0;
      int num_negative = 0;

      /// *** CHECK ***
      compute_signs_of_facet_deciders
	(scalar_grid, poly_isotable.cube, poly_isotable.ambig_cube,
	 isovalue, it, iv0, increment, parity_sign, decider_sign,
	 num_positive, num_negative);

      if (num_positive > 0 && num_negative > 0) {
	// split into pyramids
	const int dimension = poly_isotable.cube.Dimension();
	const int num_grid_vertices = scalar_grid.NumVertices();

	VERTEX_INDEX k = new_mesh_vertices.NumVertices();
	VERTEX_INDEX new_gridv = k + num_grid_vertices;
	add_vertex_at_cube_center
	  (scalar_grid, iv0, increment.CubePtrConst(), num_cube_vertices,
	   new_mesh_vertices);

	SCALAR_TYPE s = new_mesh_vertices.Scalar(k);

	// *** SHOULD PASS decider_sign[] instead of recomputing it ***
	extract_from_pyramids_in_cube
	  (scalar, poly_isotable,
	   isovalue, iv0, new_gridv, s, increment, parity_sign,
	   iso_simplices, endpoint, num_simplices, ambiguous_info);
      }
      else if (num_positive > 0) {
	/// *** CHECK ***
	extract_from_ambig_cube
	  (poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	   true, iso_simplices, endpoint, num_simplices);
      }
      else {
	/// *** CHECK ***
	extract_from_ambig_cube
	  (poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	   false, iso_simplices, endpoint, num_simplices);
      }
    }
  }

  /// Extract isosurface simplices from cube with one ambiguous facet.
  void extract_from_cube_with_ambiguous_facet
  (const SCALAR_TYPE * cube_scalar,
   const VERTEX_INDEX * cube_vertex,
   const ISOSURFACE_TABLE & cube_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const FACET_INDEX jfacet,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = cube_isotable.Dimension();
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();
    const int num_facet_vertices = num_cube_vertices/2;
    SCALAR_TYPE facet_scalar[num_facet_vertices];

    for (int k = 0; k < num_facet_vertices; k++) {
      VERTEX_INDEX kv = cube_isotable.Polyhedron().FacetVertex(jfacet, k);
      facet_scalar[k] = cube_scalar[kv];
    }

    bool decider_sign = 
      is_cube_decider_ge_isovalue
      (dimension-1, facet_scalar, num_facet_vertices, isovalue,
       parity_sign);

    extract_from_ambig_poly
      (cube_isotable, it, decider_sign, cube_vertex, true,
       iso_simplices, endpoint, num_simplices);
  }


  /// Split cube into 2d+1 subcubes and extract from subcubes.
  void extract_from_subcubes_in_cube
  (const MC_SCALAR_GRID_BASE & scalar_grid,
   const ISOSURFACE_TABLE & cube_isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const COORD_TYPE * coordA,
   const COORD_TYPE * coordB,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = scalar_grid.Dimension();
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();
    COORD_TYPE minmax_coord[2][dimension];
    COORD_TYPE coord0[dimension];
    COORD_TYPE coord1[dimension];
    COORD_TYPE coord2[dimension];
    SCALAR_TYPE cube_scalar[num_cube_vertices];
    SCALAR_TYPE new_cube_scalar[num_cube_vertices];
    VERTEX_INDEX new_cube_vertex[num_cube_vertices];
    IJK::PROCEDURE_ERROR error("extract_from_subcubes_in_cube");

    for (int d = 0; d < dimension; d++) {
      COORD_TYPE c0 = coordA[d];
      COORD_TYPE c1 = coordB[d];
      minmax_coord[0][d] = std::min(c0, c1);
      minmax_coord[1][d] = std::max(c0, c1);
    }

    get_poly_scalar(scalar_grid.ScalarPtrConst(), iv0, 
		    increment.CubePtrConst(), num_cube_vertices,
		    cube_scalar);

    scalar_grid.ComputeCoord(iv0, coord0);

    VERTEX_INDEX k = new_mesh_vertices.NumVertices();

    for (int i = 0; i < num_cube_vertices; i++) {

      new_mesh_vertices.Add();
      long mask = 1L;
      for (int d = 0; d < dimension; d++) {
	if ((i & mask) == 0)  { coord1[d] = minmax_coord[0][d]; }
	else { coord1[d] = minmax_coord[1][d]; }

	mask = (mask << 1L);
      }

      add_coord(dimension, coord0, coord1, coord2);
      new_mesh_vertices.SetCoord(k+i, coord2);

      SCALAR_TYPE s = 
	multilinear_interpolate_scalar
	(dimension, coord1, num_cube_vertices, cube_scalar);

      new_mesh_vertices.SetScalar(k+i, s);
    }

    for (int jf = 0; jf < cube_isotable.Polyhedron().NumFacets(); jf++) {
      const int d = jf/2;
      const int jside = jf%2;
      long mask = (1L << d);

      for (int i = 0; i < num_cube_vertices; i++) {
	int iside = 0;
	if ((i & mask) != 0) { iside = 1; };
	if (iside == jside) {
	  // Cube vertex i is a grid vertex.
	  VERTEX_INDEX iv1 = iv0 + increment.Cube(i);
	  new_cube_vertex[i] = iv1;
	  new_cube_scalar[i] = scalar_grid.Scalar(iv1);
	}
	else {
	  // Cube vertex i is a new mesh vertex.
	  int i2;
	  if ((i & mask) == 0) { i2 = (i | mask); }
	  else { i2 = i - mask; }
	  int k2 = k+i2;
	  new_cube_scalar[i] = new_mesh_vertices.Scalar(k2);
	  new_cube_vertex[i] = scalar_grid.NumVertices() + k2;
	}
      }

      TABLE_INDEX new_it;
      compute_isotable_index
	(new_cube_scalar, isovalue, num_cube_vertices, new_it);

      if (cube_ambig.IsAmbiguous(new_it)) {
	// Only facet jf can be ambiguous.
	extract_from_cube_with_ambiguous_facet
	  (new_cube_scalar, new_cube_vertex, cube_isotable,
	   isovalue, new_it, jf, parity_sign,
	   iso_simplices, endpoint, new_mesh_vertices,
	   num_simplices, ambiguous_info);
      }
      else {
	add_iso_simplices(cube_isotable, new_it, 0, new_cube_vertex,
			  iso_simplices, endpoint, num_simplices);
      }
    }

    // Add isosurface patch in internal cube.
    for (int i = 0; i < num_cube_vertices; i++) {
      int k2 = k+i;
      new_cube_scalar[i] = new_mesh_vertices.Scalar(k2);
      new_cube_vertex[i] = scalar_grid.NumVertices() + k2;
    }

    TABLE_INDEX new_it;
    compute_isotable_index
      (new_cube_scalar, isovalue, num_cube_vertices, new_it);

    add_iso_simplices(cube_isotable, new_it, 0, new_cube_vertex,
		      iso_simplices, endpoint, num_simplices);
  }

  /// Extract isosurface simplices from 3D cube when trilinear interpolant
  /// Produces isosurface with linear topology.
  void extract_from_two_saddle_cube_3D
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const COORD_TYPE saddle_coord[2][3],
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {

    extract_from_subcubes_in_cube
      (scalar_grid, poly_isotable.cube, poly_isotable.ambig_cube,
       isovalue, it, 
       saddle_coord[0], saddle_coord[1], iv0, increment,
       parity_sign, iso_simplices, endpoint, new_mesh_vertices,
       num_simplices, ambiguous_info);
  }


  /// Extract isosurface simplices from 3D cube with one saddle.
  /// Produces isosurface with linear topology.
  void extract_from_one_saddle_cube_3D
  (const MC_SCALAR_GRID_BASE & scalar_grid,
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const COORD_TYPE saddle_coord[3],
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.cube.Dimension();
    const int num_grid_vertices = scalar_grid.NumVertices();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();

    VERTEX_INDEX k = new_mesh_vertices.NumVertices();
    VERTEX_INDEX new_gridv = k + num_grid_vertices;

    add_mesh_vertex
      (scalar_grid, saddle_coord, iv0, increment.CubePtrConst(), 
       num_cube_vertices, new_mesh_vertices);
    SCALAR_TYPE s = new_mesh_vertices.Scalar(k);

    extract_from_pyramids_in_cube
      (scalar, poly_isotable,
       isovalue, iv0, new_gridv, s, increment, parity_sign,
       iso_simplices, endpoint, num_simplices, ambiguous_info);

  }

  void compute_split_coord_3D
  (const MC_SCALAR_GRID_BASE & scalar_grid, const VERTEX_INDEX iv0,
   const VERTEX_INDEX * increment,
   const int * parity_sign,
   COORD_TYPE coord[3])
  {
    const int dimension = 3;
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_facets = 2*dimension;
    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    SCALAR_TYPE a;
    SCALAR_TYPE b[3];
    SCALAR_TYPE c[3];
    SCALAR_TYPE q[3];
    const double EPSILON = 0.00001;

    IJK::ARRAY<SCALAR_TYPE> cube_scalar(num_cube_vertices);
    get_poly_scalar(scalar, iv0, increment,
		    num_cube_vertices, cube_scalar.Ptr());

    a = compute_multilinear_a<SCALAR_TYPE>(dimension, cube_scalar.PtrConst(), parity_sign);

    if (-EPSILON < a && a < EPSILON) { a = 0; }

    for (int i = 0; i < 3; i++) {
      b[i] = compute_multilinear_b<SCALAR_TYPE>
	(dimension, cube_scalar.PtrConst(), parity_sign, i);

      int j1 = (i+1)%dimension;
      int j2 = (i+2)%dimension;
      c[i] = compute_multilinear_c<SCALAR_TYPE>
	(dimension, cube_scalar.PtrConst(), parity_sign, j1, j2);
    }

    for (int i = 0; i < 3; i++) {
      coord[i] = 0.5;

      COORD_TYPE denom = a*coord[i] + b[i];
      if (denom <= -EPSILON || denom >= EPSILON) {
	// denom is not near zero
	int j1 = (i+1)%dimension;
	int j2 = (i+2)%dimension;

	coord[j1] = -(b[j1]*coord[i] + c[j2]) / denom;
	coord[j2] = -(b[j2]*coord[i] + c[j1]) / denom;

	if (is_point_in_unit_cube(dimension, coord)) { return; }
      }
    }

    // Use center of cube.
    for (int d = 0; d < 3; d++) { coord[0] = 0.5; };
  }


  /// Extract isosurface simplices from 3D cube with zero saddles
  void extract_from_zero_saddle_cube_3D
  (const MC_SCALAR_GRID_BASE & scalar_grid,
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const SCALAR_TYPE a,
   const COORD_TYPE trilinear_midpoint[3],
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.cube.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    const int num_cube_facets = 
      poly_isotable.cube.Polyhedron().NumFacets();
    int decider_sign[num_cube_facets];
    int num_positive = 0;
    int num_negative = 0;
    IJK::PROCEDURE_ERROR error("extract_from_zero_saddle_cube_3D");


    compute_signs_of_facet_deciders
      (scalar_grid, poly_isotable.cube, poly_isotable.ambig_cube, isovalue, it,
       iv0, increment, parity_sign, decider_sign,
       num_positive, num_negative);

    if (num_positive > 0 && num_negative == 0) {
      extract_from_ambig_cube
	(poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	 true, iso_simplices, endpoint, num_simplices);
    }
    else  if (num_positive == 0 && num_negative > 0) {
      extract_from_ambig_cube
	(poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	 false, iso_simplices, endpoint, num_simplices);
    }
    else {
      if (a != 0 &&
	  is_point_in_unit_cube(dimension, trilinear_midpoint)) {

	extract_from_pyramids_in_cube
	  (scalar_grid, poly_isotable,
	   isovalue, iv0, increment, trilinear_midpoint, parity_sign, 
	   iso_simplices, endpoint, new_mesh_vertices, 
	   num_simplices, ambiguous_info);
      }
      else {
	const COORD_TYPE center[3] = { 0.5, 0.5, 0.5};

	extract_from_pyramids_in_cube
	  (scalar_grid, poly_isotable,
	   isovalue, iv0, increment, center, parity_sign, 
	   iso_simplices, endpoint, new_mesh_vertices, 
	   num_simplices, ambiguous_info);
      }

    }
  }


  /// Extract isosurface simplices from ambiguous cube configuration.
  /// Break cube at saddle points or midpoints of saddles
  ///   to approximate topology of trilinear interpolant.
  /// Does NOT guarantee topology of trilinear interpolant.
  void extract_from_ambig_cube_using_saddle
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const TABLE_INDEX it,
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    IJK::PROCEDURE_ERROR error("extract_from_ambig_cube_using_saddle");
    const int dimension = poly_isotable.cube.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    SCALAR_TYPE cube_scalar[num_cube_vertices];
    COORD_TYPE trilinear_midpoint[3];
    COORD_TYPE saddle_offset[3];
    COORD_TYPE saddle_coord[2][3];

    get_poly_scalar(scalar, iv0, increment.CubePtrConst(),
		    num_cube_vertices, cube_scalar);

    int num_saddles = 0;
    SCALAR_TYPE a = 0;
    compute_saddle_3D
      (cube_scalar, parity_sign, a, trilinear_midpoint, saddle_offset,
       num_saddles, 0.0001);

    int num_saddles_in_cube = 0;
    if (a != 0) {
      for (int i = 0; i < num_saddles; i++) {
	int j = num_saddles_in_cube;
	if (i == 0) {
	  add_coord(dimension, trilinear_midpoint, saddle_offset,
		    saddle_coord[j]);
	}
	else {
	  subtract_coord(dimension, trilinear_midpoint, saddle_offset,
			 saddle_coord[j]);
	}

	if (is_point_in_unit_cube(dimension, saddle_coord[j])) 
	  { num_saddles_in_cube++; }
	else
	  { set_coord(dimension, 0, saddle_coord[j]); };
      }
    }

    ambiguous_info.IncrementNumCubesWithSaddle(num_saddles_in_cube);

    switch(num_saddles_in_cube) {
      
    case 0:
      extract_from_zero_saddle_cube_3D
	(scalar_grid, poly_isotable, isovalue, it,
	 a, trilinear_midpoint,
	 iv0, increment, parity_sign, iso_simplices, endpoint,
	 new_mesh_vertices, num_simplices, ambiguous_info);

      break;

    case 1:
      extract_from_one_saddle_cube_3D
	(scalar_grid, poly_isotable, isovalue, it,
	 saddle_coord[0], iv0, increment, parity_sign, 
	 iso_simplices, endpoint, new_mesh_vertices, 
	 num_simplices, ambiguous_info);
      break;

    case 2:
      extract_from_two_saddle_cube_3D
	(scalar_grid, poly_isotable, isovalue, it,
	 saddle_coord, iv0, increment, parity_sign, 
	 iso_simplices, endpoint, new_mesh_vertices, 
	 num_simplices, ambiguous_info);
      break;

    default:
      error.AddMessage("Programming error. More than two saddle points reported.");
      throw error;
    }

  }

  /// Extract isosurface simplices from cube.
  /// Use cube center to resolve ambiguities.
  /// Note: Make this inline for faster execution.
  inline void extract_from_cube_using_cube_decider
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & cube_isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
    // extract isosurface simplices in cube with primary vertex iv0
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // cube_isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment = vertex increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = 
    //     grid edge containing k'th vertex of simplex is.
    // num_simplices = number of simplices extracted
  {
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly
	(scalar, isovalue, iv0, increment.CubePtrConst(), 
	 num_cube_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, iv0, increment.CubePtrConst(), 
			     num_cube_vertices, it);

      if (cube_ambig.IsAmbiguous(it) && 
	  cube_ambig.NumAmbiguousFacets(it) == 0) {
	extract_from_ambig_cube_using_average
	  (scalar, cube_isotable, isovalue, it,
	   iv0, increment, iso_simplices, num_simplices);
      }
      else {
	add_iso_simplices_in_cube
	  (cube_isotable, it, iv0, increment, 
	   iso_simplices, num_simplices);
      }
    }
  }

  /// Extract isosurface simplices from cube.
  /// Use cube center to resolve ambiguities.
  /// Note: Make this inline for faster execution.
  inline void extract_from_cube_using_cube_decider
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & cube_isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
    // extract isosurface simplices in cube with primary vertex iv0
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // cube_isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment = vertex increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
    // endpoint[] = endpoints of edges containing isosurface vertices.
    //   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
    //   vertex i.
    // num_simplices = number of simplices extracted
  {
    const int num_cube_vertices = cube_isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, iv0, increment.CubePtrConst(), 
			num_cube_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, iv0, increment.CubePtrConst(), 
			     num_cube_vertices, it);

      if (cube_ambig.IsAmbiguous(it) && 
	  cube_ambig.NumAmbiguousFacets(it) == 0) {

	extract_from_ambig_cube_using_average
	  (scalar, cube_isotable, isovalue, it, iv0, increment,
	   iso_simplices, endpoint, num_simplices);

	ambiguous_info.num_ambiguous_cubes++;
      }
      else {
	add_iso_simplices(cube_isotable, it, iv0, increment.CubePtrConst(), 
			  iso_simplices, endpoint, num_simplices);
      }

      if (num_simplices > 0) { ambiguous_info.num_non_empty_cubes++; };
    }
  }

  /// Extract isosurface simplices from pyramid.
  /// Use linear topology to resolve ambiguities.
  /// Note: Make this inline for faster execution.
  /// @param pyramid_scalar[i] = Scalar value of pyramid vertex \a i.
  /// @param pyramid_vertex[i] = Index of grid vertex correspoinding to pyramid vertex \a i.
  inline void extract_from_pyramid_linear_topology
  (const SCALAR_TYPE * pyramid_scalar, const VERTEX_INDEX * pyramid_vertex,
   const ISOSURFACE_TABLE & pyramid_isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & pyramid_ambig,
   const bool orientation, const SCALAR_TYPE isovalue,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint, 
   int & num_simplices, SCALAR_INFO & ambiguous_info)
    // extract isosurface simplices in pyramid with primary vertex iv0
    // returns list representing isosurface simplices
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
    // endpoint[] = endpoints of edges containing isosurface vertices.
    //   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
    //   vertex i.
    // num_simplices = number of simplices extracted
  {
    const int dimension = pyramid_isotable.Dimension();
    const int num_pyramid_vertices = 
      pyramid_isotable.Polyhedron().NumVertices();
    const int num_base_vertices = num_pyramid_vertices-1;

    num_simplices = 0;

    // check whether pyramid intersects isosurface
    if (intersects_poly(pyramid_scalar, isovalue, num_pyramid_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index
	(pyramid_scalar, isovalue, num_pyramid_vertices, it);

      if (pyramid_ambig.IsAmbiguous(it) && 
	  pyramid_ambig.NumAmbiguousFacets(it) == 1) {

	bool decider_sign =
	  is_cube_decider_ge_isovalue
	  (dimension-1, pyramid_scalar, num_base_vertices, isovalue, 
	   parity_sign);

	extract_from_ambig_poly
	  (pyramid_isotable, it, decider_sign, pyramid_vertex, 
	   orientation, iso_simplices, endpoint, num_simplices);

	ambiguous_info.num_ambiguous_pyramids++;
      }
      else {
	add_iso_simplices(pyramid_isotable, it, 0, pyramid_vertex, orientation,
			  iso_simplices, endpoint, num_simplices);
      }

      if (num_simplices > 0) { ambiguous_info.num_non_empty_pyramids++; }
    }
  }


  /// Split cube into pyramids and extract isosurface simplices.
  void extract_from_pyramids_in_cube
  (const SCALAR_TYPE * scalar, const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INDEX new_grid_vertex,
   const SCALAR_TYPE new_scalar,
   const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.pyramid.Dimension();
    const VERTEX_INDEX num_cube_vertices = increment.NumCubeVertices();
    const int num_cube_facet_vertices = num_cube_vertices/2;
    const int num_cube_facets = poly_isotable.cube.Polyhedron().NumFacets();
    const int num_pyramid_vertices = num_cube_facet_vertices+1;
    SCALAR_TYPE pyramid_scalar[num_pyramid_vertices];
    VERTEX_INDEX pyramid_vertex[num_pyramid_vertices];

    const VERTEX_INDEX apex = num_pyramid_vertices-1;
    pyramid_vertex[apex] = new_grid_vertex;
    pyramid_scalar[apex] = new_scalar;

    for (int i = 0 ; i < num_cube_facets; i++) {
      for (int j = 0; j < poly_isotable.cube.Polyhedron().NumFacetVertices(i);
	   j++) {
	VERTEX_INDEX jv = poly_isotable.cube.Polyhedron().FacetVertex(i, j);
	VERTEX_INDEX iv1 = iv0 + increment.Cube(jv);
	pyramid_vertex[j] = iv1;
	pyramid_scalar[j] = scalar[iv1];
      }

      // *** SHOULD BE PRECOMPUTED FOR EACH FACET.  
      // *** SHOULD NOT REALLY ON ORDER FACETS ARE LISTED.
      bool facet_orientation = false;
      if ((i + (i/2) + dimension)%2 == 1)
	{ facet_orientation = true; };

      extract_from_pyramid_linear_topology
	(pyramid_scalar, pyramid_vertex, poly_isotable.pyramid, 
	 poly_isotable.ambig_pyramid, facet_orientation,
	 isovalue, parity_sign, 
	 iso_simplices, endpoint, num_simplices, ambiguous_info);
    }

  }

  /// Split cube into pyramids and extract isosurface simplices.
  void extract_from_pyramids_in_cube
  (const MC_SCALAR_GRID_BASE & scalar_grid,   
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const COORD_TYPE * apex_coord, const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.pyramid.Dimension();
    const VERTEX_INDEX num_cube_vertices = increment.NumCubeVertices();
    const int num_cube_facet_vertices = num_cube_vertices/2;
    const int num_pyramid_vertices = num_cube_facet_vertices+1;
    SCALAR_TYPE pyramid_scalar[num_pyramid_vertices];
    VERTEX_INDEX pyramid_vertex[num_pyramid_vertices];
    const int num_grid_vertices = scalar_grid.NumVertices();
    const int num_cube_facets = poly_isotable.cube.Polyhedron().NumFacets();

    VERTEX_INDEX k = new_mesh_vertices.NumVertices();
    VERTEX_INDEX new_grid_vertex = k + num_grid_vertices;

    add_mesh_vertex
      (scalar_grid, apex_coord, iv0, increment.CubePtrConst(), 
       num_cube_vertices, new_mesh_vertices);
    SCALAR_TYPE s = new_mesh_vertices.Scalar(k);

    const VERTEX_INDEX apex = num_pyramid_vertices-1;
    pyramid_vertex[apex] = new_grid_vertex;
    pyramid_scalar[apex] = s;

    for (int i = 0 ; i < num_cube_facets; i++) {
      for (int j = 0; j < poly_isotable.cube.Polyhedron().NumFacetVertices(i);
	   j++) {
	VERTEX_INDEX jv = poly_isotable.cube.Polyhedron().FacetVertex(i, j);
	VERTEX_INDEX iv1 = iv0 + increment.Cube(jv);
	pyramid_vertex[j] = iv1;
	pyramid_scalar[j] = scalar_grid.Scalar(iv1);
      }

      // *** SHOULD BE PRECOMPUTED FOR EACH FACET.  
      // *** SHOULD NOT RELY ON ORDER FACETS ARE LISTED.
      bool facet_orientation = false;
      if ((i + (i/2) + dimension)%2 == 1)
	{ facet_orientation = true; };

      extract_from_pyramid_linear_topology
	(pyramid_scalar, pyramid_vertex, 
	 poly_isotable.pyramid, poly_isotable.ambig_pyramid, 
	 facet_orientation, isovalue, parity_sign, 
	 iso_simplices, endpoint, num_simplices, ambiguous_info);
    }

  }


  /// Extract isosurface simplices from cube using the asymptotic decider algorithm.
  /// Produces linear topology on 3D cube facets but not inside cube.
  /// Note: Make this inline for faster execution.
  inline void extract_from_cube_adecider
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
    // extract isosurface simplices in cube with primary vertex iv0
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment = vertex increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
    // endpoint[] = endpoints of edges containing isosurface vertices.
    //   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
    //   vertex i.
    // num_simplices = number of simplices extracted
  {
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly
	(scalar, isovalue, iv0, increment.CubePtrConst(), 
	 num_cube_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, iv0, increment.CubePtrConst(), 
			     num_cube_vertices, it);

      if (poly_isotable.ambig_cube.IsAmbiguous(it)) {
	extract_from_ambig_cube_adecider
	  (scalar_grid, poly_isotable,
	   isovalue, it, iv0, increment, parity_sign,
	   iso_simplices, endpoint, new_mesh_vertices, num_simplices,
	   ambiguous_info);

	ambiguous_info.num_ambiguous_cubes++;
      }
      else {
	add_iso_simplices
	  (poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	   iso_simplices, endpoint, num_simplices);
      }

      if (num_simplices > 0) { ambiguous_info.num_non_empty_cubes++; };
    }
  }


  /// Extract isosurface simplices from cube.
  /// Break cube at saddle points or midpoints of saddle points
  ///   to approximate topology of trilinear interpolant.
  /// Does NOT guarantee topology of trilinear interpolant.
  /// Note: Make this inline for faster execution.
  inline void extract_from_cube_using_saddle
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const int * parity_sign,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   int & num_simplices, SCALAR_INFO & ambiguous_info)
    // extract isosurface simplices in cube with primary vertex iv0
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment = vertex increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
    // endpoint[] = endpoints of edges containing isosurface vertices.
    //   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
    //   vertex i.
    // num_simplices = number of simplices extracted
  {
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly
	(scalar, isovalue, iv0, increment.CubePtrConst(), 
	 num_cube_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, iv0, increment.CubePtrConst(), 
			     num_cube_vertices, it);

      if (poly_isotable.ambig_cube.IsAmbiguous(it)) {
	extract_from_ambig_cube_using_saddle
	  (scalar_grid, poly_isotable,
	   isovalue, it, iv0, increment, parity_sign,
	   iso_simplices, endpoint, new_mesh_vertices, num_simplices,
	   ambiguous_info);

	ambiguous_info.num_ambiguous_cubes++;
      }
      else {
	add_iso_simplices
	  (poly_isotable.cube, it, iv0, increment.CubePtrConst(), 
	   iso_simplices, endpoint, num_simplices);
      }

      if (num_simplices > 0) { ambiguous_info.num_non_empty_cubes++; };
    }
  }

  /// Extract isosurface simplices in cube
  /// Note: Make this inline for faster execution
  inline void extract_from_cube_nep
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices)
    // extract isosurface simplices in cube with primary vertex iv0
    //   using isotable which differentiates between scalar values
    //   less than (negative), equals or greater than (positive) the isovalue
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment[] = vertex increments
    // table_index_increment[] = table index increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = 
    //     grid edge containing k'th vertex of simplex is.
    // num_simplices = number of simplices extracted
  {
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, iv0, increment.CubePtrConst(), 
			num_cube_vertices)) {

      // compute cube index
      TABLE_INDEX it = 0;
      for (int j = 0; j < increment.NumCubeVertices(); j++) {
	VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
	if (scalar[iv1] > isovalue) {
	  it += table_index_increment.Positive(j);
	}
	else if (scalar[iv1] == isovalue) {
	  it += table_index_increment.Equals(j);
	};
      }

      add_iso_simplices_in_cube(isotable, it, iv0, increment, 
				iso_simplices, num_simplices);
    }
  }

  /// Extract isosurface simplices in cube
  /// Note: Make this inline for faster execution
  inline void extract_from_cube_nep
  (const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   int & num_simplices, std::vector<VERTEX_INDEX> & in_facet_cube)
    // extract isosurface simplices in cube with primary vertex iv0
    //   using isotable which differentiates between scalar values
    //   less than (negative), equals or greater than (positive) the isovalue
    // returns list representing isosurface simplices
    // scalar[] = array of scalar values
    //   point (x0,x1,x2,...) has scalar value
    //     scalar[x0 + x1*grid_length[0] + 
    //                   x2*grid_length[0]*grid_length[1] + ...]
    // isotable = hypercube isosurface table for given dimension
    // isovalue = isosurface scalar value
    // iv0 = primary cube vertex (cube vertex with lowest coordinates)
    // increment[] = vertex increments
    // table_index_increment[] = table index increments
    // iso_simplices[] = vector of isosurface simplex vertices
    //   iso_simplices[dimension*is+k] = 
    //     grid edge containing k'th vertex of simplex is.
    // num_simplices = number of simplices extracted
    // in_facet_cube = list of cubes whose isosurface patches lie in a grid facet
  {
    // check whether cube intersects isosurface
    bool lt_flag = true;;
    for (int j = 0; j < increment.NumCubeVertices(); j++) {
      VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
      if (scalar[iv1] >= isovalue) {
	lt_flag = false;
	break;
      };
    };
    if (lt_flag) {
      num_simplices = 0;
      return;
    }

    bool gt_flag = true;;
    for (int j = 0; j < increment.NumCubeVertices(); j++) {
      VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
      if (scalar[iv1] <= isovalue) {
	gt_flag = false;
	break;
      };
    };
    if (gt_flag) {
      num_simplices = 0;
      return;
    }

    // compute cube index
    TABLE_INDEX it = 0;
    for (int j = 0; j < increment.NumCubeVertices(); j++) {
      VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
      if (scalar[iv1] > isovalue) {
	it += table_index_increment.Positive(j);
      }
      else if (scalar[iv1] == isovalue) {
	it += table_index_increment.Equals(j);
      };
    }

    if (isotable.IsInFacet(it))
      {
	num_simplices = 0;
	in_facet_cube.push_back(iv0);
	return;
      }


    add_iso_simplices_in_cube(isotable, it, iv0, increment, 
			      iso_simplices, num_simplices);
  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube
	  (scalar, isotable, isovalue, iv, increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };
  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_and_cube_info_binary
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & cube_list,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube
	  (scalar, isotable, isovalue, iv, increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { 
	  num_mixed_cubes++; 
	  for (int k = 0; k < num_simplices; k++) 
	    { cube_list.push_back(iv); }
	}
      }
    }
  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_isopatch_from_mesh_poly
	  (scalar, isotable, isovalue, iv, increment.CubePtrConst(),
	   iso_simplices, endpoint, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };
  }


  // newly added 
  // Mark cubes as crossed by the isosurface for which isoval is in mu1-delta1 and mu1+delta1
  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   VERTEX_INDEX & num_mixed_cubes, const SCALAR_TYPE * scalar_1)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_isopatch_from_mesh_poly
	  (scalar, isotable, isovalue, iv, increment.CubePtrConst(),
	   iso_simplices, endpoint, num_simplices,scalar_1);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };
    
  }




  /// Extract iso simplices using the asymptotic decider algorithm.
  /// Produces linear topology on 3D cube facets but not inside cube.
  /// Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_adecider
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.cube.Dimension();
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);
    const VERTEX_INDEX num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    int parity_sign[num_cube_vertices];

    int num_simplices = 0;
    ambiguous_info.Clear();

    // Compute parity_sign[].
    for (VERTEX_INDEX i = 0; i < num_cube_vertices; i++) {
      int p = compute_parity(dimension, i);
      parity_sign[i] = convert_parity_to_sign(p);
    }

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_adecider
	  (scalar_grid, poly_isotable, isovalue, iv, 
	   increment, parity_sign,
	   iso_simplices, endpoint, new_mesh_vertices,
	   num_simplices, ambiguous_info);
      };
    }

  }

  /// Extract iso simplices.
  /// Break cubes at saddle points and saddle midpoints to approximate
  ///   topology of trilinear interpolant.
  /// Does NOT guarantee linear topology.
  void extract_iso_simplices_saddle
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   MC_MESH_VERTEX_LIST & new_mesh_vertices,
   SCALAR_INFO & ambiguous_info)
  {
    const int dimension = poly_isotable.cube.Dimension();
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);
    const VERTEX_INDEX num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    int parity_sign[num_cube_vertices];

    int num_simplices = 0;
    ambiguous_info.Clear();

    // Compute parity_sign[].
    for (VERTEX_INDEX i = 0; i < num_cube_vertices; i++) {
      int p = compute_parity(dimension, i);
      parity_sign[i] = convert_parity_to_sign(p);
    }

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_using_saddle
	  (scalar_grid, poly_isotable, isovalue, iv, 
	   increment, parity_sign,
	   iso_simplices, endpoint, new_mesh_vertices,
	   num_simplices, ambiguous_info);
      };
    }

  }


  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary_cube_decider
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_using_cube_decider
	  (scalar, isotable, cube_ambig, isovalue, iv, increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    }
  }


  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_binary_cube_decider
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   SCALAR_INFO & ambiguous_info)
  {
    const int dimension = isotable.Dimension();
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);
    const VERTEX_INDEX num_cube_vertices = 
      isotable.Polyhedron().NumVertices();
    int parity_sign[num_cube_vertices];

    int num_simplices = 0;
    ambiguous_info.Clear();

    // Compute parity_sign[].
    for (VERTEX_INDEX i = 0; i < num_cube_vertices; i++) {
      int p = compute_parity(dimension, i);
      parity_sign[i] = convert_parity_to_sign(p);
    }

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_using_cube_decider
	  (scalar, isotable, cube_ambig, isovalue, iv, increment,
	   iso_simplices, endpoint, num_simplices, ambiguous_info);
      };
    }

  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_nep
  (const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   VERTEX_INDEX & num_mixed_cubes)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);

    int num_simplices = 0;
    num_mixed_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_nep
	  (scalar, isotable, isovalue, iv, increment, table_index_increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
    };
  }

  void extract_iso_patches_in_facets
  (const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   const std::vector<VERTEX_INDEX> & in_facet_cube,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   MCUBE_INFO & mcube_info)
  // extract isosurface patches which lie on grid facets
  {
    const int dimension = isotable.Dimension();
    const int numf = isotable.Polyhedron().NumFacets();
    VERTEX_INDEX facet_increment[numf];
    std::vector<VERTEX_INDEX> facet_list;
    std::vector<VERTEX_INDEX> sorted_facet_list;
    std::vector<VERTEX_INDEX> facet_list2;

    compute_facet_increment(isotable.Polyhedron(), increment.AxisSize(),
			    facet_increment);

    for (int icube = 0; icube < in_facet_cube.size(); icube++) {
      // compute cube index
      VERTEX_INDEX iv0 = in_facet_cube[icube];
      TABLE_INDEX it = 0;
      for (int j = 0; j < increment.NumCubeVertices(); j++) {
	VERTEX_INDEX iv1 = iv0 + increment.Cube(j);
	if (scalar[iv1] > isovalue) {
	  it += table_index_increment.Positive(j);
	}
	else if (scalar[iv1] == isovalue) {
	  it += table_index_increment.Equals(j);
	};
      }

      VERTEX_INDEX jf = isotable.ContainingFacet(it);
      VERTEX_INDEX kf = iv0*dimension + facet_increment[jf];

      facet_list.push_back(kf);
      sorted_facet_list.push_back(kf);
    }

    mcube_info.nep.num_in_facet_cubes = facet_list.size();

    std::sort(sorted_facet_list.begin(), sorted_facet_list.end());

    mcube_info.nep.num_dup_iso_patches = 0;
    int j = 0;
    while (j < sorted_facet_list.size()) {

      if (j+1 < sorted_facet_list.size() && 
	  sorted_facet_list[j] == sorted_facet_list[j+1]) {
	// skip sorted_facet_list[j] and sorted_facet_list[j+1]
	j = j+2;
	mcube_info.nep.num_dup_iso_patches++;
      }
      else {
	facet_list2.push_back(sorted_facet_list[j]);
	j++;
      };
    }

    int num_simplices = 0;
    for (int icube = 0; icube < in_facet_cube.size(); icube++) {

      VERTEX_INDEX jf = facet_list[icube];
      if (binary_search(facet_list2.begin(), facet_list2.end(), jf)) {
	VERTEX_INDEX iv0 = in_facet_cube[icube];
	extract_from_cube_nep
	  (scalar, isotable, isovalue, iv0, increment, table_index_increment,
	   iso_simplices, num_simplices);
      };

      if (num_simplices > 0) { mcube_info.scalar.num_non_empty_cubes++; }
    }

  }

  // Encapsulating this code in a separate function produces faster run time
  void extract_iso_simplices_nep_del_dup
  // extract nep iso simplices deleting duplicate isosurface patches
  (const SCALAR_TYPE * scalar, const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const VERTEX_INCREMENT & increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   MCUBE_INFO & mcube_info)
  {
    const AXIS_SIZE_TYPE axis_size0 = increment.AxisSize(0);
    std::vector<VERTEX_INDEX> in_facet_cube;

    int num_simplices = 0;
    mcube_info.scalar.num_non_empty_cubes = 0;

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j]; 
	   iv < facet_vlist[j] + axis_size0-1; iv++) {

	extract_from_cube_nep
	  (scalar, isotable, isovalue, iv, increment, table_index_increment,
	   iso_simplices, num_simplices, in_facet_cube);

	if (num_simplices > 0) { mcube_info.scalar.num_non_empty_cubes++; }
      };
    };

    // handle list of cubes whose isosurface patches lie in facets
    extract_iso_patches_in_facets
      (scalar, isotable, isovalue, increment, table_index_increment,
       in_facet_cube, iso_simplices, mcube_info);
  }

  void extract_iso_simplices_in_facet
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int iv0,
   const int num_cube_facet_vertices, const int orth_dir, const int side, 
   const FACET_VERTEX_INCREMENT & facet_vertex_increment,
   const NEP_TABLE_INDEX_INCREMENT & table_index_increment,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_simplices)
  {
    const int dimension = isotable.Dimension();
    const int num_cube_vertices = isotable.Polyhedron().NumVertices();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();

    num_simplices = 0;

    int num_equals = 0;
    for (int k = 0; k < num_cube_facet_vertices; k++) {
      int iv = iv0 + facet_vertex_increment.FacetVertex(k);
      if (scalar[iv] > isovalue)
	{ return; }
      else if (scalar[iv] == isovalue)
	{ num_equals++; };
    }

    // if there are not enough vertices labelled '=' to form a simplex,
    //   there is no isosurface simplex in the facet
    if (num_equals < isotable.NumVerticesPerSimplex())
      return;

    TABLE_INDEX it = 0;
    for (int k = 0; k < num_cube_facet_vertices; k++) {
      int iv = iv0 + facet_vertex_increment.FacetVertex(k);

      if (scalar[iv] == isovalue) {
	int j = facet_vertex_increment.CubeVertexIndex(k);
	it += table_index_increment.Equals(j);
      }
    }

    num_simplices = isotable.NumSimplices(it);
    ISO_VERTEX_INDEX isov0 = iv0*facet_vertex_increment.NumIsoVPerGridV();
    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
	int jv = isotable.SimplexVertex(it, is, j);
	ISO_VERTEX_INDEX isov = isov0 + facet_vertex_increment.FacetIso(jv);

	iso_simplices.push_back(isov);
      }
    }
  }


  template <class DATASTRUCT>
  void extract_iso_simplices_from_datastruct
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   MCUBE_INFO & mcube_info, DATASTRUCT & datastruct)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

    clock_t t0 = clock();

    const GRID_SIZE_TYPE num_cubes = 
      scalar_grid.ComputeNumCubes();

    VERTEX_ARRAY vlist(num_cubes);
    VERTEX_INDEX vlist_length = 0;
    get_mixed_cubes(dimension, axis_size, datastruct, isovalue,
		    vlist.Ptr(), vlist_length);

    // initialize output
    iso_simplices.clear();

    extract_iso_simplices_from_list
      (scalar_grid, isotable, isovalue, vlist.PtrConst(), vlist_length, 
       iso_simplices, mcube_info);

    clock_t t1 = clock();
    mcube_info.time.extract = clock2seconds(t1-t0);

    assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  }

  template <class DATASTRUCT>
  void extract_iso_simplices_from_datastruct_nep
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const NEP_ISOSURFACE_TABLE & isotable,
   const SCALAR_TYPE isovalue, const int nep_num_dup,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   MCUBE_INFO & mcube_info, DATASTRUCT & datastruct)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

    clock_t t0 = clock();

    const GRID_SIZE_TYPE num_cubes = scalar_grid.ComputeNumCubes();

    VERTEX_ARRAY vlist(num_cubes);
    VERTEX_INDEX vlist_length = 0;
    get_mixed_cubes(dimension, axis_size, datastruct, isovalue,
		    vlist.Ptr(), vlist_length);

    // initialize output
    iso_simplices.clear();

    extract_iso_simplices_from_list_nep
      (scalar_grid, isotable, isovalue, nep_num_dup,
       vlist.PtrConst(), vlist_length, true, iso_simplices, mcube_info);

    clock_t t1 = clock();
    mcube_info.time.extract = clock2seconds(t1-t0);

    assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  }

}


/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0(dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  extract_iso_simplices_binary
    (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, num_mixed_cubes);

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_and_cube_info
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & cube_list,
 MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// cube_list[] = list of cubes containing simplices
//   cube_list[i] = cube containing i'th simplex
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  extract_iso_simplices_and_cube_info_binary
    (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, cube_list, num_mixed_cubes);

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices.
void IJKMCUBE::extract_iso_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
  const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
  std::vector<VERTEX_INDEX> & endpoint, MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
// endpoint[] = endpoints of edges containing isosurface vertices.
//   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
//   vertex i.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  extract_iso_simplices_binary
    (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, num_mixed_cubes);

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}


// Newly added
// Consider mean grid (scalar_grid) and delta grid (scalar_grid_1) for determining cell configuration
// Do mean or most probable
/// Extract isosurface simplices.
void IJKMCUBE::extract_iso_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
  const SCALAR_TYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
  std::vector<VERTEX_INDEX> & endpoint, MCUBE_INFO & mcube_info, const MC_SCALAR_GRID_BASE & scalar_grid_1, std::vector<ISO_VERTEX_INDEX> & iso_simplices_1,
  std::vector<VERTEX_INDEX> & endpoint_1, MCUBE_INFO & mcube_info_1)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
// endpoint[] = endpoints of edges containing isosurface vertices.
//   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
//   vertex i.
{

  //cout<<"came here!\n";
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  // Newly added (for delta_grid)
 
  const int dimension_1 = scalar_grid_1.Dimension();
  const AXIS_SIZE_TYPE * axis_size_1 = scalar_grid_1.AxisSize();
  const SCALAR_TYPE * scalar_1 = scalar_grid_1.ScalarPtrConst();

  assert(scalar_grid_1.Dimension() == isotable.Dimension());

  // if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
  //   { throw error; };

  // clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0_1;
  compute_num_cubes_in_grid_facet0
    (dimension_1, axis_size_1, num_cubes_in_facet0_1);

  VERTEX_ARRAY facet_vlist_1(num_cubes_in_facet0_1);
  get_cubes_in_grid_facet0(dimension_1, axis_size_1, facet_vlist_1.Ptr());

  // initialize output
  iso_simplices_1.clear();
  endpoint_1.clear();
  VERTEX_INDEX num_mixed_cubes_1 = 0;

  VERTEX_INCREMENT increment_1(axis_size_1, isotable, BINARY);

  // Extract cell configuration corresponding to mean data since we are only passing mean grid (scalar)
   extract_iso_simplices_binary
    (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, num_mixed_cubes);

  // Pass mean (scalar) and delta (scalar_1) grids for extracting iso_simplices
  // Extract most probable cell configuration for cells of a grid using uniform density   
  /*extract_iso_simplices_binary
    (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, num_mixed_cubes, scalar_1);*/

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices from grid cubes.
void IJKMCUBE::extract_iso_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue,
 const ISOSURFACE_TOPOLOGY isosurface_topology,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & endpoint, 
 MC_MESH_VERTEX_LIST & new_mesh_vertices,
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == poly_isotable.cube.Dimension());

  if (!check_isotable_encoding
      (poly_isotable.cube, ISOSURFACE_TABLE::BINARY, error)) 
    { throw error; };

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      isosurface_topology == SADDLE_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY) {

    if (!check_isotable_encoding
	(poly_isotable.pyramid, ISOSURFACE_TABLE::BINARY, error)) 
      { throw error; };
  }

  if (isosurface_topology == SADDLE_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY) {

    if (!check_isotable_encoding
	(poly_isotable.simplex, ISOSURFACE_TABLE::BINARY, error)) 
      { throw error; };
  }

  clock_t t0 = clock();

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  SCALAR_INFO ambiguous_info(dimension);
  switch(isosurface_topology) {

  case ISOTABLE_TOPOLOGY:
    extract_iso_simplices
      (scalar_grid, poly_isotable.cube, isovalue, iso_simplices, endpoint, 
       mcube_info);
    //cout<<"I am in isotable topology!";
    break;

  case CUBE_DECIDER_TOPOLOGY:
    //cout<<"I am here in cube decider topology!"; 
    extract_iso_simplices_cube_decider
      (scalar_grid, poly_isotable.cube, poly_isotable.ambig_cube,
       isovalue, iso_simplices, endpoint, mcube_info);      
    break;

  case ASYMPTOTIC_DECIDER_TOPOLOGY:
    //cout<<"I am in asymptotic decider topology!"; 
    extract_iso_simplices_adecider
      (scalar_grid, poly_isotable, isovalue, 
       iso_simplices, endpoint, new_mesh_vertices, mcube_info);
    break;

  case SADDLE_TOPOLOGY:
    // *** SADDLE TOPOLOGY IS STILL EXPERIMENTAL ***
    // *** SADDLE TOPOLOGY MAY CREATE CRACKS IN THE ISOSURFACE ***
    // ***   IN CERTAIN EXTREME CASES ***
    /* DEACTIVATE SADDLE TOPOLOGY
    extract_iso_simplices_saddle
      (scalar_grid, poly_isotable, isovalue, 
       iso_simplices, endpoint, new_mesh_vertices, mcube_info);
    */
    error.AddMessage("Programming error. SADDLE_TOPOLOGY not implemented.");
    break;

  case LINEAR_TOPOLOGY:
    // NOT IMPLEMENTED
    error.AddMessage("Programming error. LINEAR_TOPOLOGY not implemented.");
    throw error;
    break;

  default:
    error.AddMessage("Programming error. Unknown isosurface topology.");
    throw error;
  };

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%poly_isotable.cube.NumVerticesPerSimplex() == 0);
}


// newly added
/// Extract isosurface simplices from grid cubes.
void IJKMCUBE::extract_iso_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
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
 MCUBE_INFO & mcube_info_1)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == poly_isotable.cube.Dimension());

  if (!check_isotable_encoding
      (poly_isotable.cube, ISOSURFACE_TABLE::BINARY, error)) 
    { throw error; };

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      isosurface_topology == SADDLE_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY) {

    if (!check_isotable_encoding
	(poly_isotable.pyramid, ISOSURFACE_TABLE::BINARY, error)) 
      { throw error; };
  }

  if (isosurface_topology == SADDLE_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY) {

    if (!check_isotable_encoding
	(poly_isotable.simplex, ISOSURFACE_TABLE::BINARY, error)) 
      { throw error; };
  }

  clock_t t0 = clock();

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  SCALAR_INFO ambiguous_info(dimension);

  switch(isosurface_topology) {

  case ISOTABLE_TOPOLOGY:
   

      
/*    extract_iso_simplices
      (scalar_grid, poly_isotable.cube, isovalue, iso_simplices, endpoint, 
       mcube_info);*/

 // newly added to incorporate case of isoval lying in mu1-delta1 and mu1+delta1
     extract_iso_simplices
      (scalar_grid, poly_isotable.cube, isovalue, iso_simplices, endpoint, 
       mcube_info, scalar_grid_1, iso_simplices_1, endpoint_1, 
       mcube_info_1);  


    break;

  case CUBE_DECIDER_TOPOLOGY:    
    extract_iso_simplices_cube_decider
      (scalar_grid, poly_isotable.cube, poly_isotable.ambig_cube,
       isovalue, iso_simplices, endpoint, mcube_info);
    break;

  case ASYMPTOTIC_DECIDER_TOPOLOGY:
    extract_iso_simplices_adecider
      (scalar_grid, poly_isotable, isovalue, 
       iso_simplices, endpoint, new_mesh_vertices, mcube_info);
    break;

  case SADDLE_TOPOLOGY:
    // *** SADDLE TOPOLOGY IS STILL EXPERIMENTAL ***
    // *** SADDLE TOPOLOGY MAY CREATE CRACKS IN THE ISOSURFACE ***
    // ***   IN CERTAIN EXTREME CASES ***
    /* DEACTIVATE SADDLE TOPOLOGY
    extract_iso_simplices_saddle
      (scalar_grid, poly_isotable, isovalue, 
       iso_simplices, endpoint, new_mesh_vertices, mcube_info);
    */
    error.AddMessage("Programming error. SADDLE_TOPOLOGY not implemented.");
    break;

  case LINEAR_TOPOLOGY:
    // NOT IMPLEMENTED
    error.AddMessage("Programming error. LINEAR_TOPOLOGY not implemented.");
    throw error;
    break;

  default:
    error.AddMessage("Programming error. Unknown isosurface topology.");
    throw error;
  };

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%poly_isotable.cube.NumVerticesPerSimplex() == 0);
}





/// Extract isosurface simplices using the asymptotic decider.
void IJKMCUBE::extract_iso_simplices_adecider
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & endpoint, 
 MC_MESH_VERTEX_LIST & new_mesh_vertices,
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == poly_isotable.cube.Dimension());

  if (!check_isotable_encoding
      (poly_isotable.cube, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  SCALAR_INFO ambiguous_info(dimension);

  VERTEX_INCREMENT increment(axis_size, poly_isotable.cube, BINARY);

  extract_iso_simplices_adecider
    (scalar_grid, poly_isotable, isovalue, 
     facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, new_mesh_vertices,
     ambiguous_info);
  mcube_info.scalar = ambiguous_info;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%poly_isotable.cube.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices.
/// Break cubes at saddle points to approximate linear topology.
void IJKMCUBE::extract_iso_simplices_saddle
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & endpoint, 
 MC_MESH_VERTEX_LIST & new_mesh_vertices,
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices_saddle");

  assert(scalar_grid.Dimension() == poly_isotable.cube.Dimension());
  assert(scalar_grid.Dimension() == poly_isotable.pyramid.Dimension());

  if (!check_isotable_encoding
      (poly_isotable.cube, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };
  if (!check_isotable_encoding
      (poly_isotable.pyramid, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  SCALAR_INFO ambiguous_info(dimension);

  VERTEX_INCREMENT increment(axis_size, poly_isotable.cube, BINARY);

  extract_iso_simplices_saddle
    (scalar_grid, poly_isotable, isovalue, 
     facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, new_mesh_vertices,
     ambiguous_info);
  mcube_info.scalar = ambiguous_info;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%poly_isotable.cube.NumVerticesPerSimplex() == 0);
}


/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_cube_decider
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  extract_iso_simplices_binary_cube_decider
    (scalar, isotable, cube_ambig, isovalue, facet_vlist.PtrConst(), 
     num_cubes_in_facet0, increment, iso_simplices, num_mixed_cubes);

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}


/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_cube_decider
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & isotable,
 const ISOSURFACE_TABLE_AMBIG_INFO & cube_ambig,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & endpoint, MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
// endpoint[] = endpoints of edges containing isosurface vertices.
//   (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
//   vertex i.
{ 
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BINARY, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  endpoint.clear();
  SCALAR_INFO ambiguous_info(dimension);

  VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

  extract_iso_simplices_binary_cube_decider
    (scalar, isotable, cube_ambig, isovalue, 
     facet_vlist.PtrConst(), num_cubes_in_facet0,
     increment, iso_simplices, endpoint, ambiguous_info);

  mcube_info.scalar = ambiguous_info;

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);

  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, 
const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  clock_t t0 = clock();

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  // initialize output
  iso_simplices.clear();
  mcube_info.scalar.num_non_empty_cubes = 0;

  NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
  VERTEX_INCREMENT increment(axis_size, isotable, NEP);

  switch (nep_num_dup) {

  case 0:
    extract_iso_simplices_nep_del_dup
      (scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
       increment, nep_table_index_increment, iso_simplices, mcube_info);
    break;

  case 1:
    // NOT YET IMPLEMENTED
    // DEFAULTS TO case 2

  case 2:
    {
      VERTEX_INDEX num_mixed_cubes;
      extract_iso_simplices_nep
	(scalar, isotable, isovalue, facet_vlist.PtrConst(), num_cubes_in_facet0,
	 increment, nep_table_index_increment, iso_simplices, num_mixed_cubes);


      mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

      GRID_SIZE_TYPE num_cube_facet_vertices = 
	compute_num_cube_facet_vertices(dimension);
      extract_iso_simplices_nep_boundary
	(scalar_grid, isotable, isovalue, num_cube_facet_vertices,
	 iso_simplices, num_mixed_cubes);
      mcube_info.nep.num_non_empty_boundary_facets = num_mixed_cubes;
      break;
    }

  default:
    error.AddMessage("Programming error.  Illegal value ", nep_num_dup,
		     " for nep_num_dup.");
    error.AddMessage("   Value should be 0, 1 or 2.");
    throw error;
  }

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);
}

// Extract nep isosurface vertices lying on the grid boundary.
void IJKMCUBE::extract_iso_simplices_nep_boundary
(const MC_SCALAR_GRID_BASE & scalar_grid, 
const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int num_cube_facet_vertices,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_non_empty_cubes)
{
  const int dimension = isotable.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  num_non_empty_cubes = 0;


  NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);

  VERTEX_INDEX max_num_cubes_in_facet = 0;
  for (int d = 0; d < dimension; d++) {

    GRID_SIZE_TYPE num_cubes_in_facet;
    compute_num_cubes_in_grid_facet
      (dimension, axis_size, d, num_cubes_in_facet);

    if (max_num_cubes_in_facet < num_cubes_in_facet)
      { max_num_cubes_in_facet = num_cubes_in_facet; }
  }

  VERTEX_ARRAY vlist(max_num_cubes_in_facet);

  for (int d = 0; d < dimension; d++) {
    for (int side = 0; side < 2; side++) {

      FACET_VERTEX_INCREMENT facet_vertex_increment
	(axis_size, isotable, NEP, d, side);

      GRID_SIZE_TYPE num_cubes_in_facet;
      compute_num_cubes_in_grid_facet
	(dimension, axis_size, d, num_cubes_in_facet);

      get_cubes_in_grid_facet(dimension, axis_size, d, side, vlist.Ptr());

      for (int j = 0; j < num_cubes_in_facet; j++) {
	int num_simplices = 0;
	extract_iso_simplices_in_facet
	  (scalar_grid, isotable, isovalue, vlist[j], 
	   num_cube_facet_vertices, d, side,
	   facet_vertex_increment, nep_table_index_increment,
	   iso_simplices, num_simplices);

	if (num_simplices > 0) { num_non_empty_cubes++; }
      }
    }
  }
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_from_list
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const VERTEX_INDEX * vlist,
 const VERTEX_INDEX num_cubes, std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// mcube_info = mcube information
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices_from_list");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;
  int num_simplices = 0;

  if (isotable.Encoding() == ISOSURFACE_TABLE::BINARY) {

    VERTEX_INCREMENT increment(axis_size, isotable, BINARY);

    for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

      extract_from_cube
	(scalar, isotable, isovalue, vlist[j], increment,
	 iso_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };

  }
  else if (isotable.Encoding() == ISOSURFACE_TABLE::BASE3) {

    NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
    VERTEX_INCREMENT increment(axis_size, isotable, NEP);

    for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

      extract_from_cube_nep
	(scalar, isotable, isovalue, vlist[j], increment, 
	 nep_table_index_increment, iso_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };

  }
  else {
    error.AddMessage("Illegal isotable encoding.");
    throw error;
  }

  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;
}

/// Extract isosurface simplices
void IJKMCUBE::extract_iso_simplices_from_list_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int num_dup,
 const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
 const bool extract_from_boundary,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info)
// extract isosurface mesh
// returns list representing isosurface simplices
// scalar_grid = scalar grid data
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// extract_from_boundary = extract simplices lying in grid boundary if true
// iso_simplices[] = vector of isosurface simplex vertices
//   iso_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// mcube_info = nep information
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  PROCEDURE_ERROR error("extract_iso_simplices_from_list_nep");

  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_nep_num_dup(num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  // initialize output
  iso_simplices.clear();
  VERTEX_INDEX num_mixed_cubes = 0;
  int num_simplices = 0;

  NEP_TABLE_INDEX_INCREMENT nep_table_index_increment(dimension);
  VERTEX_INCREMENT increment(axis_size, isotable, NEP);

  switch(num_dup) {

  case 0:
    {
      std::vector<VERTEX_INDEX> in_facet_cube;

      for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

	extract_from_cube_nep
	  (scalar, isotable, isovalue, vlist[j], increment, 
	   nep_table_index_increment, iso_simplices, num_simplices,
	   in_facet_cube);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
      mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;

      // handle list of cubes whose isosurface patches lie in facets
      extract_iso_patches_in_facets
	(scalar, isotable, isovalue, increment, nep_table_index_increment,
	 in_facet_cube, iso_simplices, mcube_info);

      break;
    }

  case 1:
    // NOT YET IMPLEMENTED
    // DEFAULTS TO case 2

  case 2:
    {
      for (VERTEX_INDEX j = 0; j < num_cubes; j++) {

	extract_from_cube_nep
	  (scalar, isotable, isovalue, vlist[j], increment, 
	   nep_table_index_increment, iso_simplices, num_simplices);

	if (num_simplices > 0) { num_mixed_cubes++; }
      };
      mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;
      break;
    }
  }

  if (extract_from_boundary) {
    GRID_SIZE_TYPE num_cube_facet_vertices =
      compute_num_cube_facet_vertices(dimension);

    VERTEX_INDEX num_mixed_cubes;
    extract_iso_simplices_nep_boundary
      (scalar_grid, isotable, isovalue, num_cube_facet_vertices,
       iso_simplices, num_mixed_cubes);
    mcube_info.nep.num_non_empty_boundary_facets = num_mixed_cubes;
  }

}


void IJKMCUBE::extract_iso_simplices_from_octree
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info,
 const IJKOCTREE::OCTREE & octree)
{
  extract_iso_simplices_from_datastruct
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, octree);
}

void IJKMCUBE::extract_iso_simplices_from_minmax
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 MCUBE_INFO & mcube_info,
 const IJKMCUBE::MINMAX_REGIONS & minmax)
{
  extract_iso_simplices_from_datastruct
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, minmax);
}

void IJKMCUBE::extract_iso_simplices_from_octree_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
 const IJKOCTREE::OCTREE & octree)
{
  extract_iso_simplices_from_datastruct_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, 
     mcube_info, octree);
}

void IJKMCUBE::extract_iso_simplices_from_minmax_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
 const MINMAX_REGIONS & minmax)
{
  extract_iso_simplices_from_datastruct_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, 
     mcube_info, minmax);
}


// **************************************************
// EXTRACT INTERVAL VOLUME MESH
// **************************************************


namespace {

/// Extract interval volume simplices in cube
/// Note: Make this inline for faster execution
inline void extract_ivol_simplices_in_cube
(const SCALAR_TYPE * scalar, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE lower_isovalue, const SCALAR_TYPE upper_isovalue,
 const VERTEX_INDEX iv0, const VERTEX_INCREMENT & increment,
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices, 
 int & num_simplices)
// extract isosurface simplices in cube with primary vertex iv0
// returns list representing isosurface simplices
// scalar[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// isotable = hypercube isosurface table for given dimension
// iv0 = primary cube vertex (cube vertex with lowest coordinates)
// increment[] = vertex increments
// ivol_simplices[] = vector of interval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of simplex is.
// num_simplices = number of simplices in interval volume patch
{
  TABLE_INDEX it = 0;
  long k = 1;
  for (int j = 0; j < increment.NumCubeVertices(); j++) {
    int iv1 = iv0 + increment.Cube(j);
    if (scalar[iv1] >= upper_isovalue) {
      it += 2*k;
    }
    else if (scalar[iv1] >= lower_isovalue) {
      it += k;
    };
    k = 3*k;
  };

  add_iso_simplices_in_cube(isotable, it, iv0, increment, 
			    ivol_simplices, num_simplices);
};

}

/// Extract interval volume simplices
void IJKMCUBE::extract_ivol_simplices
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
 VERTEX_INDEX & num_mixed_cubes)
// extract isosurface mesh
// returns list representing isosurface simplices
// dimension = volume dimension
// isotable = hypercube isosurface table for given dimension
// scalar_grid[] = array of scalar values
//   point (x0,x1,x2,...) has scalar value
//     scalar_grid[x0 + x1*grid_length[0] + 
//                   x2*grid_length[0]*grid_length[1] + ...]
// grid_length[i] = grid dimension i
//   # grid points = 
//      grid_length[0]*grid_length[1]*grid_length[2]*...
// isovalue = isosurface scalar value
// ivol_simplices[] = inteval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of interval volume simplex is
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  const VERTEX_INDEX num_ivolv = isotable.NumIsosurfaceVertices();
  PROCEDURE_ERROR error("extract_ivol_simplices");

  assert(isotable.Dimension() == isotable.SimplexDimension());
  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  // initialize output
  ivol_simplices.clear();
  num_mixed_cubes = 0;
  int num_simplices = 0;

  SCALAR_TYPE lower_isovalue = isovalue0;
  SCALAR_TYPE upper_isovalue = isovalue1;
  if (lower_isovalue > upper_isovalue) {
    std::swap(lower_isovalue, upper_isovalue);
  };

  // get cubes in facet 0
  GRID_SIZE_TYPE num_cubes_in_facet0;
  compute_num_cubes_in_grid_facet0
    (dimension, axis_size, num_cubes_in_facet0);

  VERTEX_ARRAY facet_vlist(num_cubes_in_facet0);
  get_cubes_in_grid_facet0(dimension, axis_size, facet_vlist.Ptr());

  VERTEX_INCREMENT increment(axis_size, isotable, IVOL);

  for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
    for (VERTEX_INDEX iv = facet_vlist[j];
	 iv < facet_vlist[j] + axis_size[0]-1; iv++) {

      extract_ivol_simplices_in_cube
	(scalar, isotable, lower_isovalue, upper_isovalue, iv,
	 increment, ivol_simplices, num_simplices);

      if (num_simplices > 0) { num_mixed_cubes++; }
    };
  }
}

/// *** NOT USED/TESTED ***
/// Extract interval volume simplices
void IJKMCUBE::extract_ivol_simplices_from_list
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
 std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
 VERTEX_INDEX & num_mixed_cubes)
// extract isosurface mesh
// returns list representing isosurface simplices
// isotable = hypercube isosurface table for given dimension
// isovalue = isosurface scalar value
// vlist[] = list of primary cube vertices 
//   (vertices with lowest coord in cubes)
// num_cubes = number of cubes
// ivol_simplices[] = inteval volume simplex vertices
//   ivol_simplices[dimension*is+k] = 
//     grid edge containing k'th vertex of interval volume simplex is
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  const VERTEX_INDEX num_ivolv = isotable.NumIsosurfaceVertices();
  PROCEDURE_ERROR error("extract_ivol_simplices_from_list");

  assert(isotable.Dimension() == isotable.SimplexDimension());
  assert(scalar_grid.Dimension() == isotable.Dimension());

  if (!check_isotable_encoding(isotable, ISOSURFACE_TABLE::BASE3, error))
    { throw error; };

  // initialize output
  ivol_simplices.clear();
  num_mixed_cubes = 0;
  int num_simplices = 0;

  SCALAR_TYPE lower_isovalue = isovalue0;
  SCALAR_TYPE upper_isovalue = isovalue1;
  if (lower_isovalue > upper_isovalue) {
    std::swap(lower_isovalue, upper_isovalue);
  };

  VERTEX_INCREMENT increment(axis_size, isotable, IVOL);

  for (VERTEX_INDEX j = 0; j < num_cubes; j++) {
    extract_ivol_simplices_in_cube
      (scalar, isotable, lower_isovalue, upper_isovalue, vlist[j],
       increment, ivol_simplices, num_simplices);

    if (num_simplices > 0) { num_mixed_cubes++; }
  }

}


// **************************************************
// CLASS MULTI-RESOLUTION POLYHEDRON
// **************************************************

/// Class multires simplex
/// Vertex indices represent the grid vertex index in a 2x2x...x2 grid
class MULTIRES_POLY:public MC_ORIENTED_CELL {

protected:
  bool is_splittable;       
  VERTEX_INDEX split_point;       // grid vertex which splits current simplex
  VERTEX_INDEX prev_split_point0;
    // grid vertex which split previous simplex creating current simplex.
    // prev_split_point0 is a vertex of the current simplex
  VERTEX_INDEX prev_split_point1; // split point before prev_split_point0
  int base_dimension;             // dimension of the cubic base;

  void Init();

public:
  MULTIRES_POLY(const int dimension, const int num_vertices, 
		const bool orientation);
  MULTIRES_POLY(const MC_ORIENTED_CELL & cell);  // copy constructor
  MULTIRES_POLY(const MULTIRES_POLY & copy);     // copy constructor

  // get functions
  VERTEX_INDEX SplitPoint() const { return(split_point); };
  VERTEX_INDEX PrevSplitPoint0() const { return(prev_split_point0); };
  VERTEX_INDEX PrevSplitPoint1() const { return(prev_split_point1); };
  int BaseDimension() const { return(base_dimension); };
  bool IsSplittable() const { return(is_splittable); };

  // set functions
  void SetSplitPoint(const VERTEX_INDEX split_point)
  { this->split_point = split_point; };
  void SetPrevSplitPoint0(const VERTEX_INDEX p)
  { this->prev_split_point0 = p; };
  void SetPrevSplitPoint1(const VERTEX_INDEX p)
  { this->prev_split_point1 = p; };
  void SetBaseDimension(const int base_dimension)
  { this->base_dimension = base_dimension; };
  void SetIsSplittable(const bool flag)
  { this->is_splittable = flag; };
  void CopySplit(const MULTIRES_POLY & poly);  // copy split data

  // compute midpoint vertex
  // Precondition: Vertex(i0) and Vertex(i1) are even
  int MidPoint(const int i0, const int i1)
  { return((Vertex(i0)+Vertex(i1))/2); }
};

typedef std::vector<MULTIRES_POLY> MULTIRES_POLY_LIST;

MULTIRES_POLY::MULTIRES_POLY
(const int dimension, const int num_vertices, const bool orientation):
  MC_ORIENTED_CELL(dimension, num_vertices, orientation)
// constructor
{
  Init();
}

void MULTIRES_POLY::Init()
{
  split_point = 0;
  prev_split_point0 = 0;
  prev_split_point1 = 0;
  base_dimension = 0;
  is_splittable = false;
}

MULTIRES_POLY::MULTIRES_POLY(const MC_ORIENTED_CELL & cell):
  MC_ORIENTED_CELL(cell)
// copy constructor
{
  Init();
}

MULTIRES_POLY::MULTIRES_POLY(const MULTIRES_POLY & poly):
  MC_ORIENTED_CELL(poly)
// copy constructor
{
  CopySplit(poly);
}

void MULTIRES_POLY::CopySplit(const MULTIRES_POLY & poly)   
// copy split data
{
  SetSplitPoint(poly.SplitPoint());
  SetPrevSplitPoint0(poly.PrevSplitPoint0());
  SetPrevSplitPoint1(poly.PrevSplitPoint1());
  SetBaseDimension(poly.BaseDimension());
  SetIsSplittable(poly.IsSplittable());
}


/// set cube1 to be the i'th subcube of cube0 
void compute_subcube(const MC_ORIENTED_CUBE & cube0,
		     const int i,
		     MC_ORIENTED_CUBE & cube1)
{
  PROCEDURE_ERROR error("compute_subcube");

  if (cube0.Dimension() != cube1.Dimension()) {
    error.AddMessage("Programming error.  Cube dimensions do not match.");
    throw error;
  }

  for (int j = 0; j < cube0.NumVertices(); j++) {
    VERTEX_INDEX iv1 = (cube0.Vertex(i) + cube0.Vertex(j))/2;
    cube1.Set(j, iv1);
  }
  cube1.SetOrientation(cube0.Orientation());
}

// **************************************************
// MULTIRES SUBDIVISION ROUTINES
// **************************************************

void compute_multires_poly
(const MC_ORIENTED_CUBE & cube, MULTIRES_POLY_LIST & plist)
// Precondition: cube vertices are the corners of a 3x3x...x3 grid
{
  const int dimension = cube.Dimension();
  const int num_cube_vertices = cube.NumVertices();
  const bool cube_orientation = cube.Orientation();

  if (num_cube_vertices < 1) { return; };
  const VERTEX_INDEX cube_center = (cube.Vertex(0) + cube.LastVertex())/2;

  if (dimension > 0) {
    const int num_cube_facet_vertices = num_cube_vertices/2;
    MULTIRES_POLY poly(dimension, num_cube_facet_vertices+1,
		       cube_orientation);

    for (int ifacet = 0; ifacet < 2*dimension; ifacet++) {

      MC_ORIENTED_CUBE facet(cube, ifacet);
      MULTIRES_POLY_LIST facet_plist;
      compute_multires_poly(facet, facet_plist);

      const VERTEX_INDEX facet_center = 
	(facet.Vertex(0) + facet.LastVertex())/2;

      poly.Join(facet, cube_center);
      poly.SetPrevSplitPoint0(cube_center);
      poly.SetPrevSplitPoint1(cube_center);
      poly.SetBaseDimension(facet.Dimension());
      if (facet.Dimension() > 0) { 
	poly.SetIsSplittable(true); 
	poly.SetSplitPoint(facet_center);
      }
      else { poly.SetIsSplittable(false); }
      plist.push_back(poly);

      for (int j = 0; j < facet_plist.size(); j++) {
	poly.Join(facet_plist[j], cube_center);
	poly.CopySplit(facet_plist[j]);
	if (facet_plist[j].PrevSplitPoint0() == 
	    facet_plist[j].PrevSplitPoint1())
	  { poly.SetPrevSplitPoint1(cube_center); };
	plist.push_back(poly);
      }

      if (facet.Dimension() > 1) {
	MC_ORIENTED_CUBE subfacet(facet.Dimension(), facet.Orientation());

	for (int j = 0; j < num_cube_facet_vertices; j++) {
	  compute_subcube(facet, j, subfacet);
	  poly.Join(subfacet, cube_center);
	  poly.SetPrevSplitPoint0(facet_center);
	  poly.SetPrevSplitPoint1(cube_center);
	  poly.SetBaseDimension(facet.Dimension());
	  poly.SetIsSplittable(false);
	  plist.push_back(poly);
	}
      }
    }
  }

}

void compute_multires_polyA
(const MC_ORIENTED_CUBE & cube, const VERTEX_INDEX iv0,
 MULTIRES_POLY_LIST & plist)
// compute multires polytopes of type A using diagonal from iv0
//   to vertex opposite iv0
// Precondition: cube vertices are the corners of a 3x3x...x3 grid
{
  const int dimension = cube.Dimension();
  const int num_cube_vertices = cube.NumVertices();
  const bool cube_orientation = cube.Orientation();
  long mask;

  if (num_cube_vertices < 1) { return; };
  const VERTEX_INDEX cube_center = (cube.Vertex(0) + cube.LastVertex())/2;

  int i0 = cube.FindIndex(iv0);
  int i1 = cube.OppositeIndex(i0);
  VERTEX_INDEX iv1 = cube.Vertex(i1);

  // Type A polyhedron have dimension 2 or greater
  if (dimension > 1) {
    const int num_cube_facet_vertices = num_cube_vertices/2;
    MULTIRES_POLY poly(dimension, num_cube_facet_vertices+1,
		       cube_orientation);

    mask = 1L;
    for (int ifacet = 0; ifacet < cube.NumFacets(); ifacet++) {

      bool bit0 = ((i0&mask) != 0);
      bool bit1 = ((ifacet%2) != 0);

      if (bit0 == bit1) {
	MC_ORIENTED_CUBE facet(cube, ifacet);
	MULTIRES_POLY_LIST facet_plist;
	compute_multires_polyA(facet, iv0, facet_plist);

	for (int j = 0; j < facet_plist.size(); j++) {
	  poly.Join(facet_plist[j], iv1);
	  poly.CopySplit(facet_plist[j]);
	  poly.SetSplitPoint(cube_center);
	  plist.push_back(poly);
	}

	int j0 = facet.FindIndex(iv0);
	int j1 = facet.OppositeIndex(j0);
	VERTEX_INDEX prev_split_point0 = facet.Vertex(j1);

	poly.Join(facet, iv1);
	poly.SetSplitPoint(cube_center);
	poly.SetPrevSplitPoint0(prev_split_point0);
	poly.SetPrevSplitPoint1(iv1);
	poly.SetBaseDimension(facet.Dimension());
	poly.SetIsSplittable(false);
	plist.push_back(poly);
      }

      if ((ifacet%2) == 1) { mask = (mask << 1L); };
    }
  }

}

void compute_multires_polyB
(const MC_ORIENTED_CUBE & cube, MULTIRES_POLY_LIST & plist)
// Precondition: cube vertices are the corners of a 3x3x...x3 grid
{
  const int dimension = cube.Dimension();
  const int num_cube_vertices = cube.NumVertices();
  const bool cube_orientation = cube.Orientation();

  if (num_cube_vertices < 1) { return; };
  const VERTEX_INDEX cube_center = (cube.Vertex(0) + cube.LastVertex())/2;

  if (dimension > 1) {
    const int num_cube_facet_vertices = num_cube_vertices/2;
    MULTIRES_POLY poly(dimension, num_cube_facet_vertices+1,
		       cube_orientation);

    for (int ifacet = 0; ifacet < 2*dimension; ifacet++) {

      MC_ORIENTED_CUBE facet(cube, ifacet);
      MULTIRES_POLY_LIST facet_plist;
      compute_multires_polyB(facet, facet_plist);

      const VERTEX_INDEX facet_center = 
	(facet.Vertex(0) + facet.LastVertex())/2;

      poly.Join(facet, cube_center);
      poly.SetPrevSplitPoint0(cube_center);
      poly.SetPrevSplitPoint1(cube_center);
      poly.SetBaseDimension(facet.Dimension());
      if (facet.Dimension() > 0) { 
	poly.SetIsSplittable(true); 
	poly.SetSplitPoint(facet_center);
      }
      else { poly.SetIsSplittable(false); }
      plist.push_back(poly);

      for (int j = 0; j < facet_plist.size(); j++) {
	poly.Join(facet_plist[j], cube_center);
	poly.CopySplit(facet_plist[j]);
	if (facet_plist[j].PrevSplitPoint0() == 
	    facet_plist[j].PrevSplitPoint1())
	  { poly.SetPrevSplitPoint1(cube_center); };
	plist.push_back(poly);
      }
    }
  }

}

void convert_vertex_to_increment
(const VERTEX_INDEX * increment, MULTIRES_POLY & poly)
{
  for (int i = 0; i < poly.NumVertices(); i++) {
    VERTEX_INDEX k = poly.Vertex(i);
    poly.Set(i, increment[k]);
  }
  VERTEX_INDEX split_point = poly.SplitPoint();
  poly.SetSplitPoint(increment[split_point]);
  VERTEX_INDEX prev_split_point0 = poly.PrevSplitPoint0();
  poly.SetPrevSplitPoint0(increment[prev_split_point0]);
  VERTEX_INDEX prev_split_point1 = poly.PrevSplitPoint1();
  poly.SetPrevSplitPoint1(increment[prev_split_point1]);
}

void convert_vertex_to_increment
(const VERTEX_INDEX * increment, MULTIRES_POLY_LIST & plist)
{
  for (int i = 0; i < plist.size(); i++) 
    { convert_vertex_to_increment(increment, plist[i]); }
}


// **************************************************
// MULTIRES SUBROUTINES
// **************************************************

namespace {

  /// Return length of max cube which fits in grid
  AXIS_SIZE_TYPE compute_max_cube_length
  (const int dimension, const AXIS_SIZE_TYPE * axis_size)
  {
    const AXIS_SIZE_TYPE * min_element_ptr = 
      std::min_element(axis_size, axis_size+dimension);

    AXIS_SIZE_TYPE min_axis_size = *min_element_ptr;

    if (min_axis_size < 1) { return(0); };

    AXIS_SIZE_TYPE k = 1;
    while(2*k + 1 < std::numeric_limits<AXIS_SIZE_TYPE>::max()) {
      if (min_axis_size < 2*k+1) { break; }
      k = 2*k;
    }

    // **** WHY k+1 (# vertices)?  (Why not k (# edges)?)
    AXIS_SIZE_TYPE length = k+1;
    return(length);
  }

  AXIS_SIZE_TYPE compute_max_cube_length(const IJKMCUBE::MC_GRID & grid)
  {
    const int dimension = grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = grid.AxisSize();
    AXIS_SIZE_TYPE length = compute_max_cube_length(dimension, axis_size);
    return(length);
  }

  /// set vertices of the multires cube
  void set_multires_cube_vertices(MC_ORIENTED_CUBE & cube)
  {
    if (cube.NumVertices() == 0) { return; };

    cube.Set(0, 0);
    int prev_num_set = 1;
    int increment = 1;

    for (int d = 0; d < cube.Dimension(); d++) {
      for (int i = 0; i < prev_num_set; i++) {
	int iv = cube.Vertex(i) + 2*increment;

	cube.Set(i+prev_num_set, iv);
      }
      prev_num_set = 2*prev_num_set;
      increment = 3*increment;
    }

    if (prev_num_set != cube.NumVertices()) {
      PROCEDURE_ERROR error("set_multires_cube_vertices");
      error.AddMessage("Programming error.  Set ", prev_num_set, 
		       " cube vertices.");
      error.AddMessage("Number of cube vertices should be ", 
		       cube.NumVertices(), ".");
      throw error;
    }
  }

  VERTEX_INDEX get_num_alternate_cubes_in_facet0
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE cube_length, const VERTEX_INDEX first_indicator)
  // get number of alternate cubes which lie in facet 0
  // cubes are listed by vertex with lowest coordinates in cube
  // first_indicator = indicates first cube.  In range [0..(2^dimension-1)].
  //    if k'th bit == 1, then k'th coord of first cube equals cube_length
  {
    if (dimension < 1) { return(0); };

    // check that first cube exists
    if ((first_indicator & 1L) == 0) {
      if (axis_size[0] < cube_length + 1) { return(0); };
    }
    else {
      if (axis_size[0] < 2*cube_length + 1) { return(0); };
    }

    VERTEX_INDEX num_cubes = 1;
    long mask = 2L;
    for (int d = 1; d < dimension; d++) {
      if ((first_indicator & mask) == 0) {
	if (axis_size[d] < 1) { return(0); };
	AXIS_SIZE_TYPE k = (axis_size[d]-1+cube_length)/(2*cube_length);
	num_cubes = num_cubes * k;
      }
      else {
	if (axis_size[d] < cube_length+1) { return(0); };
	AXIS_SIZE_TYPE k = (axis_size[d]-1)/(2*cube_length);
	num_cubes = num_cubes * k;
      }
      mask = (mask << 1L);
    }
    return(num_cubes);
  }

  void get_alternate_cubes_in_facet0
  (const int dimension, const AXIS_SIZE_TYPE * axis_size, 
   const AXIS_SIZE_TYPE cube_length, const VERTEX_INDEX first_indicator,
   VERTEX_INDEX * vlist)
  // get list of alternate cubes which lie in facet 0
  // cubes are listed by vertex with lowest coordinates in cube
  // cube_length = length of each edge of cube measured in grid edges
  // first_indicator = indicates first cube.  In range [0..(2^dimension-1)].
  //    if k'th bit == 1, then k'th coord of first cube equals cube_length
  // vlist[] = list of primary vertices of alternate cubes in facet0
  //   Precondition: vlist is preallocated to size 
  //     at least number of cubes in facet0
  {
    VERTEX_INDEX axis_increment[dimension];
    VERTEX_INDEX coord0[dimension];
    VERTEX_INDEX vfirst = 0;
    PROCEDURE_ERROR error("get_alternate_cubes_in_facet0");

    if (first_indicator < 0 || first_indicator >= (1L << dimension)) {
      error.AddMessage("Illegal value ", first_indicator, 
		       " for first_indicator.");
      error.AddMessage("Start index must be in the range [0..",
		       (1L << dimension)-1, "].");
      throw error;
    };

    compute_increment(dimension, axis_size, axis_increment);

    for (int d = 0; d < dimension; d++)
      axis_increment[d] = axis_increment[d] * cube_length;

    // compute index and coordinates of first cube
    vfirst = 0;
    long mask = 1L;
    for (int d = 0; d < dimension; d++) {
      if ((first_indicator & mask) == 0) {
	coord0[d] = 0;
      }
      else {
	coord0[d] = cube_length;
	vfirst += axis_increment[d];
      }
      mask = (mask << 1L);
    }

    // check that first cube exists
    for (int d = 0; d < dimension; d++) {
      if (coord0[d] + cube_length >= axis_size[d]) {
	// grid does not contain cube with primary vertex at given coordinates
	return; 
      }
    }

    // initialize prev_num_cubes
    vlist[0] = vfirst;
    VERTEX_INDEX prev_num_cubes = 1;

    // Process axes 1,2,..., dimension-1
    VERTEX_INDEX * vcur_ptr = vlist + prev_num_cubes;
    for (int d = 1; d < dimension; d++) {
      assert(axis_size[d] >= 1);
      VERTEX_INDEX vincrement = 2*axis_increment[d];

      for (VERTEX_INDEX i = coord0[d]+2*cube_length; 
	   i + cube_length < axis_size[d]; i += 2*cube_length) {
	for (VERTEX_INDEX * vprev_ptr = vlist; 
	     vprev_ptr != vlist+prev_num_cubes; vprev_ptr++) {
	  *(vcur_ptr) = vincrement + *(vprev_ptr);
	  vcur_ptr++;
	};
	vincrement = vincrement + 2*axis_increment[d];
      }
      VERTEX_INDEX k = (axis_size[d]-1+cube_length-coord0[d])/(2*cube_length);
      prev_num_cubes = prev_num_cubes*k;
    }

    // check
    VERTEX_INDEX num_cubes = get_num_alternate_cubes_in_facet0
      (dimension, axis_size, cube_length, first_indicator);

    if (prev_num_cubes != num_cubes) {
      PROCEDURE_ERROR error("get_alternate_cubes_in_facet0");
      error.AddMessage("Programming error.  Added ", prev_num_cubes, 
		       " primary cube vertices to list.");
      error.AddMessage
	("Number of alternating cubes in facet 0 (first_indicator = ", 
	 first_indicator, ") should be ", num_cubes, ".");
      throw error;
    }

  }

  /// Extract isosurface simplices from polytope
  /// Note: Make this inline for faster execution
  inline void extract_from_poly
  (const SCALAR_TYPE * scalar, const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, 
   const VERTEX_INDEX iv0, const MULTIRES_POLY & poly,
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint, int & num_simplices)
  // extract isosurface simplices from poly
  // returns list representing isosurface simplices
  // scalar[] = array of scalar values
  //   point (x0,x1,x2,...) has scalar value
  //     scalar[x0 + x1*grid_length[0] + 
  //                   x2*grid_length[0]*grid_length[1] + ...]
  // poly_isotable = isosurface tables of cube, pyramid and simplex
  // isovalue = isosurface scalar value
  // iv0 = base for poly vertices
  // poly = polytope increments
  // iso_simplices[] = vector of isosurface simplex vertices
  //   iso_simplices[dimension*is+k] = 
  //     index of edge in endpoint[] containing k'th vertex of simplex is.
  // endpoint[] = array of edge endpoints.
  //     (endpoint[i*k], endpoint[2*i+1]) = endpoints of i'th edge.
  // num_simplices = number of simplices extracted
  {
    // *** WORKS ONLY FOR 2D and 3D ***
    int dim_diff = poly.Dimension() - poly.BaseDimension();
    if (poly.BaseDimension() == 0) { dim_diff = poly.Dimension()-1; };

    if (poly.Dimension() > 3) {
      PROCEDURE_ERROR error("extract_iso_simplices_from_poly");
      error.AddMessage("NOT IMPLEMENTED for dimension ",
		       poly.Dimension(), ".");
      throw error;
    }

    if (dim_diff == 1) {
      // *** SHOULD BE REPLACED BY extract_isopatch_from_poly
      if (poly.Orientation()) {
	extract_isopatch_from_mesh_poly
	  (scalar, poly_isotable.pyramid, isovalue, iv0, 
	   poly.VlistPtrConst(), iso_simplices, endpoint, num_simplices);
      }
      else {
	extract_isopatch_reverse_orient_from_mesh_poly
	  (scalar, poly_isotable.pyramid, isovalue, iv0, 
	   poly.VlistPtrConst(), iso_simplices, endpoint, num_simplices);
      }
    }
    else {
      if (poly.Orientation()) {
	extract_isopatch_from_mesh_poly
	  (scalar, poly_isotable.simplex, isovalue, iv0, 
	   poly.VlistPtrConst(), iso_simplices, endpoint, num_simplices);
      }
      else {
	extract_isopatch_reverse_orient_from_mesh_poly
	  (scalar, poly_isotable.simplex, isovalue, iv0, 
	   poly.VlistPtrConst(), iso_simplices, endpoint, num_simplices);
      }
    }
  };

  /// Return true if multires poly A is in multires mesh
  inline bool in_multires_gridA
  (const MULTIRES_GRID & multires_grid, 
   const VERTEX_INDEX iv, const MULTIRES_POLY & poly)
  {
    const VERTEX_INDEX iv_prev_split0 = iv+poly.PrevSplitPoint0();
    if (multires_grid.Scalar(iv_prev_split0) != MULTIRES_A) { return(false); };
    if (multires_grid.IsMultires(iv+poly.SplitPoint())) { return(false); };
    const VERTEX_INDEX iv_prev_split1 = iv+poly.PrevSplitPoint1();
    if (multires_grid.Scalar(iv_prev_split1) == MULTIRES_A) { return(false); };
    return(true);
  }

  /// Return true if multires length1 poly A is in multires mesh
  inline bool in_multires_gridA_length1
  (const MULTIRES_GRID & multires_grid, 
   const VERTEX_INDEX iv, const MULTIRES_POLY & poly)
  {
    const VERTEX_INDEX iv_prev_split0 = iv+poly.PrevSplitPoint0();
    if (multires_grid.Scalar(iv_prev_split0) != MULTIRES_A) { return(false); };
    const VERTEX_INDEX iv_prev_split1 = iv+poly.PrevSplitPoint1();
    if (multires_grid.Scalar(iv_prev_split1) == MULTIRES_A) { return(false); };
    return(true);
  }

  /// Return true if multires poly B is in multires mesh
  inline bool in_multires_gridB
  (const MULTIRES_GRID & multires_grid, 
   const VERTEX_INDEX iv, const MULTIRES_POLY & poly)
  {
    const VERTEX_INDEX iv_prev_split0 = iv+poly.PrevSplitPoint0();
    if (!multires_grid.IsMultires(iv_prev_split0)) { return(false); };
    if (multires_grid.Scalar(iv_prev_split0) == MULTIRES_A) { return(false); };
    if (multires_grid.IsMultires(iv+poly.SplitPoint())) { return(false); };
    return(true);
  }

  void extract_iso_simplices_multires_binary
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const MULTIRES_GRID & multires_grid,
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const int cube_type,
   const AXIS_SIZE_TYPE vertex_space, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   VERTEX_INDEX & num_non_empty_cubes)
  // cube_type = cube type. In range [0..num_cube_vertices-1]
  // vertex_space = number of edges between vertices
  //   Precondition: Must be positive.
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = multires_grid.AxisSize();
    const AXIS_SIZE_TYPE axis_size0 = multires_grid.AxisSize(0);
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int cube_length = vertex_space/2;
    MC_ORIENTED_CUBE cube(dimension, true);
    IJK::PROCEDURE_ERROR error("extract_iso_simplices_multires_binary");

    if (vertex_space <= 0) {
      error.AddMessage("Programming error.  Vertex space must be positive.");
      throw error;
    }
    if (vertex_space%2 != 0 || cube_length%2 != 0) {
      error.AddMessage("Programming error.  Vertex space and cube length must be even.");
      throw error;
    }

    int num_iso_simplices = 0;
    num_non_empty_cubes = 0;
    if (num_cubes_in_facet0 == 0) { return; };

    const int num_region_corners = compute_num_cube_vertices(dimension);

    ARRAY<VERTEX_INDEX> corner_increment(num_region_corners);
    compute_cubic_region_corner_increment
      (dimension, axis_size, cube_length+1, corner_increment.Ptr());

    VERTEX_INDEX multires_increment =  
      corner_increment[cube.OppositeIndex(cube_type)];
    VERTEX_INDEX split_increment = 
      corner_increment[cube.OppositeIndex(0)]/2;

    int num_region_vertices;
    compute_num_grid_vertices_in_region
      (dimension, 2, num_region_vertices);

    ARRAY<VERTEX_INDEX> region_vertex_increment(num_region_vertices);
    compute_region_vertex_increment
      (dimension, axis_size, 2, cube_length/2, region_vertex_increment.Ptr());

    set_multires_cube_vertices(cube);


    MULTIRES_POLY_LIST plistA;
    compute_multires_polyA(cube, cube.Vertex(cube_type), plistA);
    convert_vertex_to_increment(region_vertex_increment.PtrConst(), plistA);

    MULTIRES_POLY_LIST plistB;
    compute_multires_polyB(cube, plistB);
    convert_vertex_to_increment(region_vertex_increment.PtrConst(), plistB);

    GRID_COORD_TYPE coord0 = facet_vlist[0]%axis_size0;

    int num_simplices = 0;
    

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j];
	   iv + cube_length < facet_vlist[j] + axis_size0-coord0; 
	   iv = iv + vertex_space) {

	int num_iso_simplices_in_cube = 0;

	VERTEX_INDEX iv_multires = iv + multires_increment;
	VERTEX_INDEX iv_split = iv + split_increment;

	if (multires_grid.Scalar(iv_multires) == MULTIRES_A && 
	    !multires_grid.IsMultires(iv_split)) {

	  extract_isopatch_from_mesh_poly
	    (scalar, poly_isotable.cube, isovalue, iv, 
	     corner_increment.PtrConst(), 
	     iso_simplices, endpoint, num_simplices);

	  num_iso_simplices_in_cube += num_simplices;
	}
	else {

	  for (int i = 0; i < plistA.size(); i++) {

	    if (in_multires_gridA(multires_grid, iv, plistA[i])) {
	      extract_from_poly
		(scalar, poly_isotable, isovalue, iv, plistA[i],
		 iso_simplices, endpoint, num_simplices);
	    }
	  }

	  for (int i = 0; i < plistB.size(); i++) {

	    if (in_multires_gridB(multires_grid, iv, plistB[i])) {
	      extract_from_poly
		(scalar, poly_isotable, isovalue, iv, plistB[i],
		 iso_simplices, endpoint, num_simplices);

	    }
	  }
	}

	if (num_iso_simplices_in_cube > 0) { num_non_empty_cubes++; }
      }
    }
  }

  void extract_iso_simplices_multires_binary_length1
  (const MC_SCALAR_GRID_BASE & scalar_grid, 
   const MULTIRES_GRID & multires_grid,
   const POLY_ISOTABLE & poly_isotable,
   const SCALAR_TYPE isovalue, const VERTEX_INDEX * facet_vlist,
   const VERTEX_INDEX num_cubes_in_facet0,
   const int cube_type,
   const AXIS_SIZE_TYPE vertex_space, 
   std::vector<ISO_VERTEX_INDEX> & iso_simplices,
   std::vector<VERTEX_INDEX> & endpoint,
   VERTEX_INDEX & num_non_empty_cubes)
  // cube_type = cube type. In range [0..num_cube_vertices-1]
  // vertex_space = number of edges between vertices
  //   Precondition: Must be positive.
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = multires_grid.AxisSize();
    const AXIS_SIZE_TYPE axis_size0 = multires_grid.AxisSize(0);
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int num_cube_vertices = 
      poly_isotable.cube.Polyhedron().NumVertices();
    const int cube_length = vertex_space/2;
    MC_ORIENTED_CUBE cube(dimension, true);

    assert(vertex_space > 0);

    int num_iso_simplices = 0;
    num_non_empty_cubes = 0;
    if (num_cubes_in_facet0 == 0) { return; };

    VERTEX_INCREMENT increment(axis_size, poly_isotable.cube, BINARY);

    GRID_COORD_TYPE coord0 = facet_vlist[0]%axis_size0;

    VERTEX_INDEX multires_increment = 
      increment.Cube(cube.OppositeIndex(cube_type));

    int num_simplices = 0;

    MULTIRES_POLY_LIST plistA;
    compute_multires_polyA(cube, cube.Vertex(cube_type), plistA);
    convert_vertex_to_increment(increment.CubePtrConst(), plistA);

    for (VERTEX_INDEX j = 0; j < num_cubes_in_facet0; j++) {
      for (VERTEX_INDEX iv = facet_vlist[j];
	   iv + cube_length < facet_vlist[j] + axis_size0-coord0; 
	   iv = iv + vertex_space) {

	int num_iso_simplices_in_cube = 0;

       	if (intersects_poly
	    (scalar, isovalue, iv, increment.CubePtrConst(), 
	     num_cube_vertices))
	{
	  VERTEX_INDEX iv_multires = iv + multires_increment;

	  if (multires_grid.Scalar(iv_multires) == MULTIRES_A) {

	    extract_isopatch_from_mesh_poly
	      (scalar, poly_isotable.cube, isovalue, iv, 
	       increment.CubePtrConst(),
	       iso_simplices, endpoint, num_simplices);

	    num_iso_simplices_in_cube += num_simplices;

	  }
	  else if (multires_grid.IsMultires(iv_multires)) {

	    for (int i = 0; i < plistA.size(); i++) {

	      if (in_multires_gridA_length1(multires_grid, iv, plistA[i])) {

		extract_from_poly
		  (scalar, poly_isotable, isovalue, iv, plistA[i],
		   iso_simplices, endpoint, num_simplices);

	      }
	    }

	  }
	}

	if (num_iso_simplices_in_cube > 0) { num_non_empty_cubes++; }
      }
    }
  }


}

void IJKMCUBE::multires_extract
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const MULTIRES_GRID & multires_grid,
 const POLY_ISOTABLE & poly_isotable, const SCALAR_TYPE isovalue,
 std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 std::vector<VERTEX_INDEX> & endpoint, MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  AXIS_SIZE_TYPE max_region_length =
    compute_max_cube_length(multires_grid);
  IJK::PROCEDURE_ERROR error("multires_extract");

  clock_t t0 = clock();

  if (!multires_grid.IsProcessed()) {
    error.AddMessage("Programming error.  Call multires_grid.Process() before calling multires_extract.");
    throw error;
  }
  mcube_info.scalar.num_non_empty_cubes = 0;

  AXIS_SIZE_TYPE length = 1;
  while (length <= max_region_length) {

    // allocate facet_vlist
    VERTEX_INDEX num_full_regions_in_facet0;
    compute_num_full_regions_in_grid_facet
      (dimension, axis_size, length, 0, num_full_regions_in_facet0);
    ARRAY<VERTEX_INDEX> facet_vlist(num_full_regions_in_facet0);

    const GRID_SIZE_TYPE num_cube_vertices =
      compute_num_cube_vertices(dimension);

    for (int i = 0; i < num_cube_vertices; i++) {
      VERTEX_INDEX num_non_empty_regions = 0;

      VERTEX_INDEX num_regions = 
	get_num_alternate_cubes_in_facet0(dimension, axis_size, length, i);
      get_alternate_cubes_in_facet0
	(dimension, axis_size, length, i, facet_vlist.Ptr());

      if (length > 1) {

	extract_iso_simplices_multires_binary
	  (scalar_grid, multires_grid, poly_isotable, isovalue, 
	   facet_vlist.PtrConst(), num_regions, i, 
	   2*length, iso_simplices, endpoint, num_non_empty_regions);
      }
      else {

	extract_iso_simplices_multires_binary_length1
	  (scalar_grid, multires_grid, poly_isotable, isovalue, 
	   facet_vlist.PtrConst(), num_regions, i, 
	   2, iso_simplices, endpoint, num_non_empty_regions);
      }

      mcube_info.scalar.num_non_empty_cubes += num_non_empty_regions;
    }

    length = length*2;
  }

  clock_t t1 = clock();
  mcube_info.time.extract = clock2seconds(t1-t0);
}

// **************************************************
// PROCESS MULTIRES_GRID
// **************************************************

namespace {

  bool is_cube_contained
  (const int dimension,const AXIS_SIZE_TYPE * axis_size,
   const GRID_COORD_TYPE * coord, const AXIS_SIZE_TYPE cube_length)
  {
    for (int d = 0; d < dimension; d++) {
      if (coord[d] + cube_length >= axis_size[d])
	return(false);
    }

    return(true);
  }

  bool is_cube_vertex_contained
  (const int dimension,const AXIS_SIZE_TYPE * axis_size,
   const GRID_COORD_TYPE * coord, const AXIS_SIZE_TYPE cube_length,
   const VERTEX_INDEX iv)
  // iv = index of cube  vertex
  {

    VERTEX_INDEX j = iv;
    for (int d = 0; d < dimension; d++) {
      if ((j % 2) == 1) {
	if (coord[d]+cube_length >= axis_size[d]) 
	  return(false);
      }
      else {
	if (coord[d] >= axis_size[d]) 
	  return(false);
      }
      j = (j >> 1L);
    }

    return(true);
  }

  AXIS_SIZE_TYPE compute_initial_cube_length
  (const int dimension, const AXIS_SIZE_TYPE * axis_size)
  {
    AXIS_SIZE_TYPE cube_length = 1;
    if (dimension  > 0) { 

      const AXIS_SIZE_TYPE * max_element_ptr = 
	std::max_element(axis_size, axis_size+dimension);

      AXIS_SIZE_TYPE max_axis_size = *max_element_ptr;

      while (cube_length+1 < max_axis_size) 
	{ cube_length = 2*cube_length; }
    }

    return(cube_length);
  }

  VERTEX_INDEX get_num_centers_in_facet
  (const int k, const int dimension, const AXIS_SIZE_TYPE * axis_size,
   const AXIS_SIZE_TYPE cube_length, const int ifacet)
  // get number of k-face centers in grid facet ifacet
  // if N[k,d] = # k-face centers in R^d, then
  //    N[k,d] = N[k,d-1] x ceil(axis_size[d-1])/2) +
  //             N[k-1,d-1] x floor(axis_size[d-1]/2)
  {
    assert(cube_length%2 == 0);

    VERTEX_INDEX num_centers = 0;
    VERTEX_INDEX num_lower_dim_centers = 0;
    VERTEX_INDEX num_lower_dim_facet_centers = 0;
    AXIS_SIZE_TYPE half_cube_length = cube_length/2;

    if (k > dimension) { return(0); }
    if (k < 0) { return(0); }
    if (dimension <= 0) { return(1); }
    if (dimension == ifacet+1) {
      num_centers = get_num_centers_in_facet(k, dimension-1, axis_size, 
					     cube_length, ifacet);
      return(num_centers);
    }

    num_lower_dim_centers = 
      get_num_centers_in_facet(k, dimension-1, axis_size, 
			       cube_length, ifacet);

    num_lower_dim_facet_centers = 
      get_num_centers_in_facet(k-1, dimension-1, axis_size, 
			       cube_length, ifacet);

    AXIS_SIZE_TYPE num_even = ((axis_size[dimension-1]-1)/cube_length)+1;
    AXIS_SIZE_TYPE num_odd = (axis_size[dimension-1]-1)/cube_length;

    num_centers = num_lower_dim_centers*num_even +
      num_lower_dim_facet_centers*num_odd;

    return(num_centers);
  }

  void compute_face_centers
  (const int k, const int dimension, const AXIS_SIZE_TYPE * axis_size,
   const AXIS_SIZE_TYPE cube_length, 
   const int ifacet, const AXIS_SIZE_TYPE * axis_increment, VERTEX_INDEX * vlist)
  // compute face centers of all k-faces in facet ifacet
  // cube_length = number of grid edges in cube edge
  //   Precondition: cube_length is even
  // vlist[] = list of face centers
  //   Precondition: vlist is preallocated to size 
  //     at least number of k-faces in facet ifacet
  {
    PROCEDURE_ERROR error("compute_face_centers");

    if (cube_length%2 != 0) {
      error.AddMessage("Programming error.  Illegal cube_length ",
		       cube_length, ".");
      error.AddMessage("Cube length must be even.");
      throw error;
    }

    AXIS_SIZE_TYPE half_cube_length = cube_length/2;

    if (k > dimension) { return; }

    if (k == 0 && dimension == 0) {
      vlist[0] = 0;
      return;
    }

    if (ifacet+1 == dimension) {

      compute_face_centers(k, dimension-1, axis_size, cube_length,
			   ifacet, axis_increment, vlist);
      return;
    }

    VERTEX_INDEX vlist_length = 0;
    if (0 <= k && k < dimension) {
      const VERTEX_INDEX num_lower_dim_centers = 
	get_num_centers_in_facet(k, dimension-1, axis_size, 
				 cube_length, ifacet);
      std::vector<VERTEX_INDEX> lower_dim_vlist;
      lower_dim_vlist.reserve(num_lower_dim_centers);

      compute_face_centers(k, dimension-1, axis_size, cube_length,
			   ifacet, axis_increment, &(lower_dim_vlist[0]));

      for (int i = 0; i < axis_size[dimension-1]; i += cube_length) {
	for (int j = 0; j < num_lower_dim_centers; j++) {
	  vlist[vlist_length] = lower_dim_vlist[j]+i*axis_increment[dimension-1];
	  vlist_length++;
	}
      }
    }

    if (0 < k) {
      const VERTEX_INDEX num_lower_dim_facet_centers = 
	get_num_centers_in_facet(k-1, dimension-1, axis_size, 
				 cube_length, ifacet);
      std::vector<VERTEX_INDEX> lower_dim_facet_vlist;
      lower_dim_facet_vlist.reserve(num_lower_dim_facet_centers);

      compute_face_centers(k-1, dimension-1, axis_size, cube_length,
			   ifacet, axis_increment, &(lower_dim_facet_vlist[0]));

      for (int i = half_cube_length; 
	   i+half_cube_length < axis_size[dimension-1]; i += cube_length) {
	for (int j = 0; j < num_lower_dim_facet_centers; j++) {
	  vlist[vlist_length] = 
	    lower_dim_facet_vlist[j]+i*axis_increment[dimension-1];
	  vlist_length++;
	}
      }
    }

    VERTEX_INDEX num_centers = 
      get_num_centers_in_facet(k, dimension, axis_size, cube_length, ifacet);

    if (num_centers != vlist_length) {
      error.AddMessage("Programming error.  Computed ", 
		       vlist_length, " face centers.");
      error.AddMessage("Number of ", k, "-face centers in facet ", ifacet,
		       " is ", num_centers, ".");
      throw error;
    }
  }

  void set_multires_vertex
  (const VERTEX_INDEX iv, const VERTEX_INDEX * cube_vertex_increment, 
   const VERTEX_INDEX num_cube_vertices, const MULTIRES_VERTEX_TYPE vtype,
   MULTIRES_GRID & multires_grid)
  {
    const VERTEX_INDEX ilast = num_cube_vertices-1;
    VERTEX_INDEX iv0 = iv - (cube_vertex_increment[ilast]/2);
    for (int i = 0; i < num_cube_vertices; i++) {
      VERTEX_INDEX iv2 = iv0 + cube_vertex_increment[i];
      multires_grid.Set(iv2, vtype);
    }
  }

  struct CUBIC_REGION {
    VERTEX_INDEX primary;
    AXIS_SIZE_TYPE length;
  };

  /// Partition grid into cubic regions and mark region corners 
  /// as multires vertices
  void partition_multires(MULTIRES_GRID & multires_grid)
  {
    const int dimension = multires_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = multires_grid.AxisSize();
    GRID_COORD_TYPE coord[dimension];

    std::vector<CUBIC_REGION> stack;

    VERTEX_INDEX num_cube_vertices = compute_num_cube_vertices(dimension);
    IJK::ARRAY<VERTEX_INDEX> cube_increment(num_cube_vertices);

    // initialize stack
    const AXIS_SIZE_TYPE containing_region_length =
      compute_initial_cube_length(dimension, axis_size);
    CUBIC_REGION region = {0, containing_region_length};
    stack.push_back(region);

    while (stack.size() > 0) {

      int itop = stack.size()-1;
      region = stack[itop];
      stack.pop_back();

      compute_coord(region.primary, dimension, axis_size, coord);
      compute_cubic_region_corner_increment
	(dimension, axis_size, region.length+1, cube_increment.Ptr());

      if (is_cube_contained(dimension, axis_size, coord, region.length)) {
	// mark region corners as multires vertices
	for (int i = 0; i < num_cube_vertices; i++) {
	  VERTEX_INDEX iv = region.primary + cube_increment[i];
	  multires_grid.Set(iv, MULTIRES_B);
	}
      }
      else if (region.length > 1) {

	const AXIS_SIZE_TYPE half_region_length = region.length/2;
	for (int i = 0; i < num_cube_vertices; i++) {
	  if (is_cube_vertex_contained
	      (dimension, axis_size, coord, half_region_length, i)) {
	    const VERTEX_INDEX iv = region.primary + cube_increment[i]/2;

	    CUBIC_REGION new_grid_region = {iv, half_region_length};
	    stack.push_back(new_grid_region);
	  }
	}

      }
    }

  }

  void propagate_multires(MULTIRES_GRID & multires_grid)
  {
    const int dimension = multires_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = multires_grid.AxisSize();
    AXIS_SIZE_TYPE axis_increment[dimension];

    const AXIS_SIZE_TYPE max_cube_length = 
      compute_max_cube_length(multires_grid);

    VERTEX_INDEX num_cube_vertices = compute_num_cube_vertices(dimension);
    ARRAY<VERTEX_INDEX> cube_vertex_increment(num_cube_vertices);

    compute_increment(dimension, axis_size, axis_increment);

    int num_region_vertices;
    compute_num_grid_vertices_in_region
      (dimension, 2, num_region_vertices);
    ARRAY<VERTEX_INDEX> region_vertex_increment(num_region_vertices);

    std::vector<VERTEX_INDEX> center_list;

    AXIS_SIZE_TYPE length = 2;
    while (length <= max_cube_length) {

      AXIS_SIZE_TYPE half_length = length/2;

      compute_cubic_region_corner_increment
	(dimension, axis_size, length+1, cube_vertex_increment.Ptr());
      compute_region_vertex_increment
	(dimension, axis_size, 2, length/2, region_vertex_increment.Ptr());

      // propagate from type A cube centers to face vertices
      VERTEX_INDEX num_cube_centers =
	get_num_centers_in_facet(dimension-1, dimension, axis_size, 
				 length, 0);
      
      center_list.reserve(num_cube_centers);

      compute_face_centers(dimension-1, dimension, axis_size, length,
			   0, axis_increment, &(center_list[0]));

      for (VERTEX_INDEX i = 0; i < num_cube_centers; i++) {
	for (AXIS_SIZE_TYPE j = half_length; 
	     j + half_length < axis_size[0]; j += length) {

	  VERTEX_INDEX iv = center_list[i] + j*axis_increment[0];

	  if (multires_grid.Scalar(iv) == MULTIRES_A) {
	    set_multires_vertex
	      (iv, region_vertex_increment.PtrConst(), num_region_vertices, 
	       MULTIRES_B, multires_grid);
	  }
	}
      }

      // propagate from proper face centers
      for (int k = 1; k < dimension; k++) {
	for (int ifacet = 0; ifacet < dimension; ifacet++) {
	  VERTEX_INDEX num_centers =
	    get_num_centers_in_facet(k, dimension, axis_size, length, ifacet);

	  center_list.reserve(num_centers);

	  compute_face_centers(k, dimension, axis_size, length,
			       ifacet, axis_increment, &(center_list[0]));


	  VERTEX_INDEX half_cube_increment = 
	    axis_increment[ifacet]*half_length;

	  for (VERTEX_INDEX i = 0; i < num_centers; i++) {

	    // propagate forward
	    for (AXIS_SIZE_TYPE j = 0; j + half_length < axis_size[ifacet];
		 j += length) {

	      VERTEX_INDEX iv = center_list[i] + j*axis_increment[ifacet];

	      if (multires_grid.IsMultires(iv)) {
		VERTEX_INDEX iv2 = iv+half_cube_increment;
		multires_grid.Set(iv2, MULTIRES_B);
	      }
	    }

	    // propagate backward
	    for (AXIS_SIZE_TYPE j = length; j < axis_size[ifacet];
		 j += length) {
	      VERTEX_INDEX iv = center_list[i] + j*axis_increment[ifacet];

	      if (multires_grid.IsMultires(iv)) {
		VERTEX_INDEX iv2 = iv-half_cube_increment;
		multires_grid.Set(iv2, MULTIRES_B);
	      }
	    }
	  }
	}
      }

      // propagate from cube centers
      center_list.reserve(num_cube_centers);

      compute_face_centers(dimension-1, dimension, axis_size, length,
			   0, axis_increment, &(center_list[0]));

      for (VERTEX_INDEX i = 0; i < num_cube_centers; i++) {
	for (AXIS_SIZE_TYPE j = half_length; 
	     j + half_length < axis_size[0]; j += length) {

	  VERTEX_INDEX iv = center_list[i] + j*axis_increment[0];

	  if (multires_grid.IsMultires(iv)) {
	    set_multires_vertex
	      (iv, cube_vertex_increment.PtrConst(), num_cube_vertices, 
	       MULTIRES_A, multires_grid);
	  }
	}
      }

      length = 2*length;
    }
  }

  /// Type A multires vertices are multires vertices at the center of face f
  /// such that every vertex of every subface of f is a type A multires vertex
  void compute_typeA_multires_vertices(MULTIRES_GRID & multires_grid)
  {
    const int dimension = multires_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = multires_grid.AxisSize();
    AXIS_SIZE_TYPE axis_increment[dimension];
    PROCEDURE_ERROR error("compute_typeA_multires_vertices");

    const AXIS_SIZE_TYPE max_cube_length = 
      compute_max_cube_length(multires_grid);

    VERTEX_INDEX num_cube_vertices = compute_num_cube_vertices(dimension);
    ARRAY<VERTEX_INDEX> cube_vertex_increment(num_cube_vertices);

    compute_increment(dimension, axis_size, axis_increment);

    std::vector<VERTEX_INDEX> center_list;

    // initialize all multires vertices to MULTIRES_A
    multires_grid.SetAllMultires(MULTIRES_A);

    AXIS_SIZE_TYPE length = 2;
    while (length <= max_cube_length) {

      AXIS_SIZE_TYPE half_length = length/2;

      compute_cubic_region_corner_increment
	(dimension, axis_size, length+1, cube_vertex_increment.Ptr());

      // propagate from proper face centers
      for (int k = 1; k < dimension; k++) {
	for (int ifacet = 0; ifacet < dimension; ifacet++) {
	  VERTEX_INDEX num_centers =
	    get_num_centers_in_facet(k, dimension, axis_size, length, ifacet);

	  center_list.reserve(num_centers);

	  compute_face_centers(k, dimension, axis_size, length,
			       ifacet, axis_increment, &(center_list[0]));

	  VERTEX_INDEX half_cube_increment = 
	    axis_increment[ifacet]*half_length;

	  for (VERTEX_INDEX i = 0; i < num_centers; i++) {

	    // propagate forward
	    for (AXIS_SIZE_TYPE j = 0; j + half_length < axis_size[ifacet];
		 j += length) {

	      VERTEX_INDEX iv = center_list[i] + j*axis_increment[ifacet];
	      if (multires_grid.Scalar(iv) != MULTIRES_A) {
		VERTEX_INDEX iv2 = iv+half_cube_increment;
		if (multires_grid.Scalar(iv2) == MULTIRES_A) 
		  { multires_grid.Set(iv2, MULTIRES_B); };
	      }
	    }

	    // propagate backward
	    for (AXIS_SIZE_TYPE j = length; j < axis_size[ifacet];
		 j += length) {
	      VERTEX_INDEX iv = center_list[i] + j*axis_increment[ifacet];
	      if (multires_grid.Scalar(iv) != MULTIRES_A) {
		VERTEX_INDEX iv2 = iv-half_cube_increment;
		if (multires_grid.Scalar(iv2) == MULTIRES_A) 
		  { multires_grid.Set(iv2, MULTIRES_B); };
	      }
	    }
	  }
	}
      }

      length = 2*length;
    }
  }

}

// **************************************************
// MULTIRES_GRID
// **************************************************

void MULTIRES_GRID::Init()
{
  is_processed = false;
  SetAll(NOT_MULTIRES);
  SetCorners2Multires();
}

void MULTIRES_GRID::Process()
{
  IJK::PROCEDURE_ERROR error("MULTIRES_GRID::Process()");

  if (IsProcessed()) {
    error.AddMessage("Programming error.  Process can be called only once on MULTIRES_GRID.");
    throw error;
  }

  partition_multires(*this);

  propagate_multires(*this);
  compute_typeA_multires_vertices(*this);

  is_processed = true;
}

void MULTIRES_GRID::SetAllMultires(const MULTIRES_VERTEX_TYPE vtype)
  // set all multires vertices to vtype
{
  for (VERTEX_INDEX iv = 0; iv < NumVertices(); iv++) {
    if (IsMultires(iv)) { scalar[iv] = vtype; }
  }
}

void MULTIRES_GRID::SetCorners2Multires() 
  // set corners to multires vertices
{
  SetCorners(MULTIRES_B);
}

void MULTIRES_GRID::SetSubsample2Multires(const AXIS_SIZE_TYPE scale)
  // set subsampled vertices to multires vertices
{
  SetSubsample(scale, MULTIRES_B);
}

void MULTIRES_GRID::SetRegion2Multires
(const GRID_BOX & box, const AXIS_SIZE_TYPE scale)
  // set scalar values of vertices in region to multires vertices
{
  SetRegion(box, scale, MULTIRES_B);
}

void MULTIRES_GRID::SetRegion2Multires
(const std::vector<GRID_BOX> & box_list, const AXIS_SIZE_TYPE scale)
  // set scalar values of vertices in list of regions to multires vertices
{
  SetRegion(box_list, scale, MULTIRES_B);
}

