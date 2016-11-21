/// \file ijkcoord.txx
/// ijk templates for coordinate arithmetic
/// Version 0.1.0

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

#ifndef _IJKCOORD_
#define _IJKCOORD_

/// Coordinate arithmetic functions.
namespace IJK {

  /// Set all vertex coordinates to \a c.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param c = Scalar constant.  Set vertex coordinates to \a c.
  /// @param coord[] = Output coordinates.
  template <class DTYPE, class STYPE, class CTYPE>
  void set_coord(const DTYPE dimension, const STYPE c, CTYPE coord)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord[d] = c; };
  }

  /// Copy \a coord0[] to \a coord1[].
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Output coordinates.
  template <class DTYPE, class CTYPE0, class CTYPE1>
  void copy_coord(const DTYPE dimension, 
		  const CTYPE0 coord0, const CTYPE1 coord1)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]; };
  }

  /// Copy \a coord0[] to \a coord1[].
  /// Faster algorithm when \a coord0[] and \a coord1[] are both type float.
  template <class DTYPE>
  void copy_coord(const DTYPE dimension, 
		  const float * coord0, const float * coord1)
  {
    std::copy(coord0, coord0+dimension, coord1);
  }

  /// Copy \a coord0[] to \a coord1[].
  /// Faster algorithm when \a coord0[] and \a coord1[] are both type double.
  template <class DTYPE>
  void copy_coord(const DTYPE dimension, 
		  const double * coord0, const double * coord1)
  {
    std::copy(coord0, coord0+dimension, coord1);
  }

  /// Add \a coord0[] to \a coord1[].
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param coord2 = Output coordinate equal to (\a coord0[] + \a coord1[]).
  template <class DTYPE, class CTYPE0, class CTYPE1, class CTYPE2>
  void add_coord(const DTYPE dimension, const CTYPE0 coord0,
		 const CTYPE1 coord1, CTYPE2 coord2)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord2[d] = coord0[d] + coord1[d]; };
  }

  /// Subtract \a coord1[] from \a coord0[].
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param coord2 = Output coordinate equal to (\a coord1[] - \a coord1[]).
  template <class DTYPE, class CTYPE0, class CTYPE1, class CTYPE2>
  inline void subtract_coord
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1,
   CTYPE2 coord2)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord2[d] = coord0[d] - coord1[d]; };
  }

  /// Compare \a coord0[] to \a coord1[].
  /// Return true if coordinates are the same.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  template <class DTYPE, class CTYPE0, class CTYPE1>
  bool compare_coord(const DTYPE dimension, 
		     const CTYPE0 coord0, const CTYPE1 coord1)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] != coord1[d]) { return(false); }; }

    return(true);
  }

  /// Compute the midpoint of two coordinates.
  template <class DTYPE, class CTYPE0, class CTYPE1, class CTYPE2>
  inline void compute_midpoint
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1,
   CTYPE2 midpoint_coord)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { midpoint_coord[d] = (coord0[d] + coord1[d])/2.0; };
  }

  /// Return true if all elements of coord0[] are less than or equal
  ///   to all elements of coord1[].
  template <class DTYPE, class CTYPE0, class CTYPE1>
  bool is_less_than_or_equal_to
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] > coord1[d]) { return(false); }; };

    return(true);
  }

}

#endif
