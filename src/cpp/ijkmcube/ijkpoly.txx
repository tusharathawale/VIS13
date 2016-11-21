/// \file ijkpoly.txx
/// Polytope template classes and routines.
/// Version 0.1.0

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

#ifndef _IJKPOLY_
#define _IJKPOLY_

#include "ijk.txx"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASSES CELL AND ORIENTED_CELL
  // **************************************************

  /// Cell template.
  /// @param DTYPE Dimension type.
  /// @param NTYPE Number type.
  /// @param VTYPE Vertex type.
  template <class DTYPE, class NTYPE, class VTYPE> class CELL {

  protected:
    DTYPE dimension;        ///< Dimension of space containing cell.
    NTYPE num_vertices;     ///< Number of vertices.
    VTYPE * vlist;          ///< List of cell vertices.

    void Init(const DTYPE dimension, const NTYPE num_vertices);
    void Copy(const CELL & cell);
    void FreeAll();

  public:
    CELL(const DTYPE dimension, const NTYPE num_vertices)  ///< Constructor.
      { Init(dimension, num_vertices); };
    ~CELL() { FreeAll(); };
    CELL(const CELL & cell)            ///< Copy constructor. 
      { Copy(cell); };       
    const CELL & operator = (const CELL & right);

    // ******** set functions
    void Set(const NTYPE i, const VTYPE iv)  ///< Set i'th vertex to iv.
      { vlist[i] = iv; };
    void Join(const CELL & cell, const VTYPE iv);

    // ******** get functions
    DTYPE Dimension () const  ///< Return dimension of space containing cell.
      { return(dimension); }; 
    VTYPE Vertex(const int i) const  ///< Return i'th cell vertex.
      { return(vlist[i]); };
    NTYPE NumVertices() const        ///< Return number of cell vertices.
      { return(num_vertices); };
    VTYPE * VlistPtrConst() const    ///< Return const pointer to vertices.
      { return(vlist); }; 

    /// Return last vertex in cell.
    /// Precondition: Number of cell vertices > 0.
    VTYPE LastVertex() const { return(vlist[num_vertices-1]); };
  };

  /// Oriented cell template
  template <class DTYPE, class NTYPE, class VTYPE> 
  class ORIENTED_CELL:public CELL<DTYPE,NTYPE,VTYPE> {
 
  protected:
    bool orientation;   ///< Cell orientation.

    /// Initialize cell.
    void Init(const DTYPE dimension, const NTYPE num_vertices,
	      const bool orientation);

  public:

    /// Constructor.
    ORIENTED_CELL(const DTYPE dimension, const NTYPE num_vertices, 
		  const bool orientation):
      CELL<DTYPE,NTYPE,VTYPE>(dimension, num_vertices)
    { Init(dimension, num_vertices, orientation); };

    ORIENTED_CELL(const ORIENTED_CELL & cell);     ///< Copy constructor.
    const ORIENTED_CELL & operator =               /// Copy assignment.
      (const ORIENTED_CELL & right);  

    // ******** set functions
    void Join (const ORIENTED_CELL & cell, const VTYPE iv);
    void SetOrientation         /// Set cell orientation.
    (const bool orientation) { this->orientation = orientation; };
    void ReverseOrientation()   /// Reverse cell orientation.
      { orientation = !orientation; }

    // ******** get functions
    bool Orientation() const    /// Return cell orientation.
      { return(orientation); }
  };

  // *******************************************************
  // TEMPLATE CLASSES ORIENTED_CUBE AND ORIENTED_SIMPLEX
  // *******************************************************

  /// Oriented cube template
  template <class DTYPE, class NTYPE, class VTYPE> 
  class ORIENTED_CUBE:public ORIENTED_CELL<DTYPE,NTYPE,VTYPE> {

  public:

    /// Constructor.
    ORIENTED_CUBE(const DTYPE dimension, const bool orientation):
      ORIENTED_CELL<DTYPE,NTYPE,VTYPE>
    (dimension, (1L<<dimension), orientation) {};

    /// Constructor.
    ORIENTED_CUBE(const ORIENTED_CUBE & cube, const NTYPE facet_index);

    // ******** set functions
    void Reflect(const DTYPE d);

    // ******** get functions
    NTYPE NumFacets() const     ///< Return number of cube facets.
      { return(2*this->Dimension()); };
    NTYPE FindIndex(const VTYPE iv0) const;
    NTYPE OppositeIndex(const NTYPE i0) const;
  };

  /// Oriented simplex template
  template <class DTYPE, class NTYPE, class VTYPE> 
  class ORIENTED_SIMPLEX:public ORIENTED_CELL<DTYPE,NTYPE,VTYPE> {

  public:
    /// Constructor.
    ORIENTED_SIMPLEX (const DTYPE dimension, const bool orientation) : 
      ORIENTED_CELL<DTYPE,NTYPE,VTYPE>
    (dimension, dimension+1, orientation) {};
  };

  // **************************************************
  // TEMPLATE CELL MEMBER FUNCTIONS
  // **************************************************

  /// Initialize cell.
  template <class DTYPE, class NTYPE, class VTYPE>
  void CELL<DTYPE,NTYPE,VTYPE>::
  Init(const DTYPE dimension, const NTYPE num_vertices)
  {
    this->dimension = dimension;
    this->num_vertices = num_vertices;
    vlist = new VTYPE[num_vertices];

    // initialize to 0,..., NumVertices()-1
    for (NTYPE i = 0; i < NumVertices(); i++) 
      { vlist[i] = i; }
  }

  /// Copy cell.
  template <class DTYPE, class NTYPE, class VTYPE>
  void CELL<DTYPE,NTYPE,VTYPE>::
  Copy(const CELL<DTYPE,NTYPE,VTYPE> & cell)
  {
    Init(cell.Dimension(), cell.NumVertices());
    for (NTYPE i = 0; i < NumVertices(); i++) {
      Set(i, cell.Vertex(i));
    }
  }

  /// Copy assignment.
  template <class DTYPE, class NTYPE, class VTYPE>
  const CELL<DTYPE,NTYPE,VTYPE> & CELL<DTYPE,NTYPE,VTYPE>::
  operator = (const CELL<DTYPE,NTYPE,VTYPE> & right)
  {
    if (&right != this) {
      FreeAll();
      Copy(right);
    }
  }

  /// Set the current cell to be the join of cell c and vertex iv.
  template <class DTYPE, class NTYPE, class VTYPE>
  void CELL<DTYPE,NTYPE,VTYPE>::
  Join(const CELL<DTYPE,NTYPE,VTYPE> & c, const VTYPE iv)
  {
    if (this->Dimension() != c.Dimension()+1 ||
	this->NumVertices() != c.NumVertices()+1) {
      FreeAll();
      Init(c.Dimension()+1, c.NumVertices()+1);
    }

    for (NTYPE i = 0; i < c.NumVertices(); i++) 
      { Set(i, c.Vertex(i)); };
    Set(this->NumVertices()-1, iv);
  }

  /// Free all memory.
  template <class DTYPE, class NTYPE, class VTYPE>
  void CELL<DTYPE,NTYPE,VTYPE>::FreeAll()
  {
    dimension = 0;
    num_vertices = 0;
    delete [] vlist;
    vlist = NULL;
  }

  // **************************************************
  // TEMPLATE ORIENTED_CELL MEMBER FUNCTIONS
  // **************************************************

  /// Initialize cell orientation.
  template <class DTYPE, class NTYPE, class VTYPE>
  void ORIENTED_CELL<DTYPE,NTYPE,VTYPE>::
  Init(const DTYPE dimension, const NTYPE num_vertices, 
       const bool orientation)
  {
    this->orientation = orientation;
  }

  /// Copy constructor.
  template <class DTYPE, class NTYPE, class VTYPE>
  ORIENTED_CELL<DTYPE,NTYPE,VTYPE>::ORIENTED_CELL
  (const ORIENTED_CELL<DTYPE,NTYPE,VTYPE> & cell) :
    CELL<DTYPE,NTYPE,VTYPE>(cell)
  {
    this->orientation = cell.Orientation();
  }

  /// Copy assignment.
  template <class DTYPE, class NTYPE, class VTYPE>
  const ORIENTED_CELL<DTYPE,NTYPE,VTYPE> & 
  ORIENTED_CELL<DTYPE,NTYPE,VTYPE>::
  operator = (const ORIENTED_CELL<DTYPE,NTYPE,VTYPE> & right)
  {
    if (&right != this) {
      this->CELL<DTYPE,NTYPE,VTYPE>::operator = (right);
      this->orientation = right.Orientation();
    }
  }

  /// Set the current cell to be the join of cell c and vertex iv.
  template <class DTYPE, class NTYPE, class VTYPE>
  void ORIENTED_CELL<DTYPE,NTYPE,VTYPE>::
  Join(const ORIENTED_CELL<DTYPE,NTYPE,VTYPE> & c, const VTYPE iv)
  {
    CELL<DTYPE,NTYPE,VTYPE>::Join(c, iv);
    SetOrientation(c.Orientation());
  }

  // **************************************************
  // TEMPLATE ORIENTED_CUBE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <class DTYPE, class NTYPE, class VTYPE>
  ORIENTED_CUBE<DTYPE,NTYPE,VTYPE>::
  ORIENTED_CUBE
  (const ORIENTED_CUBE<DTYPE,NTYPE,VTYPE> & cube, const NTYPE facet_index):
    ORIENTED_CELL<DTYPE,NTYPE,VTYPE>
  (cube.Dimension()-1, cube.NumVertices()/2, cube.Orientation())
  {
    NTYPE k = 0;
    NTYPE facet_parity = facet_index%2;
    NTYPE facet_dir = facet_index/2;
    long M1 = (1L << facet_dir);
    long M2 = (M1 << 1L);

    for (NTYPE i = 0; i < cube.NumVertices(); i++) {
      if (facet_parity == ((i%M2)/M1)) {
	Set(k, cube.Vertex(i));
	k++;
      }
    }

    if ((facet_parity+facet_dir+this->Dimension())%2 == 1) 
      { this->ReverseOrientation(); };

    if (k != this->NumVertices()) {
      IJK::PROCEDURE_ERROR error("CUBE_FACE");
      error.AddMessage("Error constructing cube facet ", facet_index, ".");
      throw error;
    }
  }

  /// Reflect cube across hyperplane orthogonal to axis d.
  /// Swaps vertices of two facets orthogonal to d.
  template <class DTYPE, class NTYPE, class VTYPE>
  void ORIENTED_CUBE<DTYPE,NTYPE,VTYPE>::
  Reflect(const DTYPE d)
  {
    long M1 = (1L << d);
    long M2 = (M1 << 1L);

    for (NTYPE i = 0; i < this->NumVertices(); i++) {
      if ((i%M2) < M1) { std::swap(this->vlist[i], this->vlist[i+M1]); }
    }

    this->ReverseOrientation();
  }

  /// Return k where vertex iv0 is k'th cube vertex.
  /// Precondition: iv0 is a cube vertex.
  template <class DTYPE, class NTYPE, class VTYPE>
  NTYPE ORIENTED_CUBE<DTYPE,NTYPE,VTYPE>::
  FindIndex(const VTYPE iv0) const
  {
    NTYPE i0 = 0;
    while (i0 < this->NumVertices()) {
      if (Vertex(i0) == iv0) { return(i0); };
      i0++;
    }

    IJK::PROCEDURE_ERROR error("ORIENTED_CUBE::FindIndex");
    error.AddMessage("Vertex ", iv0, " is not one of the cube vertices.");
    throw error;
  }

  /// Return k where k'th cube vertex is opposite i'th cube vertex.
  /// Precondition: 0 <= i0 < number of cube vertices
  template <class DTYPE, class NTYPE, class VTYPE>
  NTYPE ORIENTED_CUBE<DTYPE,NTYPE,VTYPE>::
  OppositeIndex(const NTYPE i0) const
  { return(this->NumVertices()-i0-1); }

}

#endif
