/// \file ijkmcube_datastruct.cxx
/// ijkmcube data structures

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

#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"
#include "ijkmcube_util.h"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKTABLE;


// **************************************************
// CLASS POLY_ISOTABLE
// **************************************************

void POLY_ISOTABLE::ComputeAmbiguityInformation()
{
  if (cube.IsTableAllocated())
    { ambig_cube.ComputeAmbiguityInformation(cube); };

  if (pyramid.IsTableAllocated())
    { ambig_pyramid.ComputeAmbiguityInformation(pyramid); };
}

// **************************************************
// MARCHING CUBES ISOSURFACE CLASS
// **************************************************

void MC_ISOSURFACE::Clear()
{
  simplex_vert.clear();
  vertex_coord.clear();
}

// **************************************************
// CLASS MC_DATA_FLAGS
// **************************************************

// Initialize MC_DATA_FLAGS
void MC_DATA_FLAGS::Init()
{
  use_minmax = false;
  use_octree = false;
  use_nep = false;
  use_list = false;
  use_multires = false;
  interval_volume_flag = false;
  cube_containing_simplex_flag = false;
  nep_num_dup = 2;
  snap_flag = false;
  snap_value = 0.5;
  isosurface_topology = ISOTABLE_TOPOLOGY;
  interpolation_type = LINEAR_INTERPOLATION;

  // Represent edge by a single integer identifier.
  edge_representation = EDGE_ID;  

  num_resolution_levels = 2;
}

// **************************************************
// CLASS MC_DATA
// **************************************************

// Initialize MC_DATA
void MC_DATA::Init()
{
  minmax = NULL;
  octree = NULL;
  snap_minmax = NULL;
  snap_octree = NULL;
  multires_grid = NULL;

  is_scalar_grid_set = false;
}

void MC_DATA::FreeAll()
{
  delete minmax;
  minmax = NULL;
  use_minmax = false;
  delete octree;
  octree = NULL;
  use_octree = false;

  is_scalar_grid_set = false;
}


// Copy scalar grid
void MC_DATA::CopyScalarGrid(const MC_SCALAR_GRID_BASE & scalar_grid2)
{
  scalar_grid.Copy(scalar_grid2);
  is_scalar_grid_set = true;
}

// Subsample scalar grid
void MC_DATA::SubsampleScalarGrid
(const MC_SCALAR_GRID_BASE & scalar_grid2, const int subsample_resolution)
{
  scalar_grid.Subsample(scalar_grid2, subsample_resolution);
  is_scalar_grid_set = true;
}

// Supersample scalar grid
void MC_DATA::SupersampleScalarGrid
(const MC_SCALAR_GRID_BASE & scalar_grid2, const int supersample_resolution)
{
  scalar_grid.Supersample(scalar_grid2, supersample_resolution);
  is_scalar_grid_set = true;
}

// Copy, subsample or supersample scalar grid.
void MC_DATA::SetScalarGrid
(const MC_SCALAR_GRID_BASE & scalar_grid2, 
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution)
{
  PROCEDURE_ERROR error("MC_DATA::SetGrid");

  if (flag_subsample && flag_supersample) {
    error.AddMessage("Scalar grid cannot both be subsampled and supersampled.");
    throw error;
  }
  
  if (flag_subsample) {
    // subsample grid
    SubsampleScalarGrid(scalar_grid2, subsample_resolution);
  }
  else if (flag_supersample) {
    // supersample grid
    SupersampleScalarGrid(scalar_grid2, supersample_resolution);
  }
  else {
    CopyScalarGrid(scalar_grid2);
  };

}

// Set octree
void MC_DATA::SetOctree()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  IJK::PROCEDURE_ERROR error("MC_DATA::SetOctree");

  if (use_octree || octree != NULL) {
    error.AddMessage("Programming error. Octree can only be set once.");
    throw error;
  }

  if (use_minmax || minmax != NULL) {
    error.AddMessage("Programming error. Minmax regions are set.");
    throw error;
  }

  if (!is_scalar_grid_set) {
    error.AddMessage("Programming error.  Scalar grid must be set before octree is created.");
    throw error;
  }

  octree = new IJKOCTREE::OCTREE(dimension, axis_size);
  octree->SetMinMax(scalar_grid.ScalarPtrConst());

  use_octree = true;
}

// Set snap octree
void MC_DATA::SetSnapOctree()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  IJK::PROCEDURE_ERROR error("MC_DATA::SetSnapOctree");

  if (use_octree || snap_octree != NULL) {
    error.AddMessage("Programming error. Octree can only be set once.");
    throw error;
  }

  if (use_minmax || snap_minmax != NULL) {
    error.AddMessage("Programming error. Minmax regions are set.");
    throw error;
  }

  if (!is_scalar_grid_set) {
    error.AddMessage("Programming error.  Scalar grid must be set before octree is created.");
    throw error;
  }

  snap_octree = new SNAP_OCTREE(dimension, axis_size);
  snap_octree->SetMinMax(scalar_grid.ScalarPtrConst());

  use_octree = true;
}

// Set minmax
void MC_DATA::SetMinmaxRegions(const int region_length)
{
  IJK::PROCEDURE_ERROR error("MC_DATA::SetMinmaxRegions");

  if (use_minmax || minmax != NULL) {
    error.AddMessage("Programming error. Minmax regions can only be set once.");
    throw error;
  }

  if (use_octree || octree != NULL) {
    error.AddMessage("Programming error. Octree is set.");
    throw error;
  }

  if (!is_scalar_grid_set) {
    error.AddMessage("Programming error.  Scalar grid must be set before minmax regions are created.");
    throw error;
  }

  minmax = new IJKMCUBE::MINMAX_REGIONS();
  minmax->ComputeMinMax(scalar_grid, region_length);

  use_minmax = true;
}

// Set snap minmax
void MC_DATA::SetSnapMinmaxRegions(const int region_length)
{
  IJK::PROCEDURE_ERROR error("MC_DATA::SetSnapMinmaxRegions");

  if (use_minmax || snap_minmax != NULL) {
    error.AddMessage("Programming error. Minmax regions can only be set once.");
    throw error;
  }

  if (use_octree || snap_octree != NULL) {
    error.AddMessage("Programming error. Octree is set.");
    throw error;
  }

  if (!is_scalar_grid_set) {
    error.AddMessage("Programming error.  Scalar grid must be set before minmax regions are created.");
    throw error;
  }

  snap_minmax = new SNAP_MINMAX();
  snap_minmax->ComputeMinMax(scalar_grid, region_length);

  use_minmax = true;
}

// Set high resolution regions
void MC_DATA::SetHighResolutionRegions
(const std::vector<GRID_BOX> & high_resolution_regions)
{
  this->high_resolution_regions.clear();

  for (VERTEX_INDEX i = 0; i < high_resolution_regions.size(); i++) {
    this->high_resolution_regions.push_back(high_resolution_regions[i]);
  }

}

// Set multires grid
void MC_DATA::SetMultiresGrid(const int num_resolution_levels)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  IJK::PROCEDURE_ERROR error("MC_DATA::SetMultires");

  if (!is_scalar_grid_set) {
    error.AddMessage("Programming error.  Scalar grid must be set before multires grid is created.");
    throw error;
  }

  this->num_resolution_levels = num_resolution_levels;

  multires_grid = new MULTIRES_GRID(dimension, axis_size);

  int high_res_scale = 1;
  int scale = high_res_scale;
  for (int ilevel = 1; ilevel < num_resolution_levels; ilevel++)
    { scale = scale * 2; }

  multires_grid->SetSubsample2Multires(scale);
  multires_grid->SetRegion2Multires(high_resolution_regions, high_res_scale);
  multires_grid->Process();

  use_multires = true;
}

/// Set mesh edge representation
void MC_DATA::SetEdgeRepresentation(const MESH_EDGE_REPRESENTATION r)
{
  edge_representation = r;
}

/// Turn NEP (negative-equals-positive) on
void MC_DATA::SetNEPOn(const int nep_num_dup)
{
  use_nep = true;
  this->nep_num_dup = nep_num_dup;
}


/// Turn NEP (negative-equals-positive) on
void MC_DATA::SetNEPOff()
{
  use_nep = false;
}

/// Set interval_volume_flag
void MC_DATA::SetIntervalVolumeFlag(const bool flag)
{
  interval_volume_flag = flag;
}

/// Set use_list
void MC_DATA::SetUseList(const bool flag)
{
  use_list = flag;
}

/// Set cube_containing_simplex_flag
void MC_DATA::SetCubeContainingSimplexFlag(const bool flag)
{
  cube_containing_simplex_flag = flag;
}

/// Turn snap on
void MC_DATA::SetSnapOff()
{
  snap_flag = false;
}

/// Turn snap on
void MC_DATA::SetSnapOn(const SNAP_TYPE snap_value, const int nep_num_dup)
{
  snap_flag = true;
  this->snap_value = snap_value;
  this->nep_num_dup = nep_num_dup;
}

/// Set desired isosurface topology
void MC_DATA::SetIsosurfaceTopology
(const ISOSURFACE_TOPOLOGY isosurface_topology)
{
  this->isosurface_topology = isosurface_topology;
}

/// Set type of interpolation
void MC_DATA::SetInterpolationType
(const INTERPOLATION_TYPE interpolation_type)
{
  this->interpolation_type = interpolation_type;
}

/// Check data structure
bool MC_DATA::Check(IJK::ERROR & error) const
{
  IJK::ERROR error2;

  if (!IsScalarGridSet()) {
    error.AddMessage("Scalar grid is not set.");
    return(false);
  }

  if (UseOctree()) {
    if ((Snap() && SnapOctree() == NULL) ||
	(!Snap() && Octree() == NULL)) {
      error.AddMessage("Octree not allocated.");
      return(false);
    }
  }

  if (Octree() != NULL) {
    if (!Octree()->Check(error2)) {
      error.AddMessage("Error in constructing octree.");
      return(false);
    }
  }

  if (SnapOctree() != NULL) {
    if (!SnapOctree()->Check(error2)) {
      error.AddMessage("Error in constructing snap octree.");
      return(false);
    }
  }

  if (UseMinmaxRegions()) {
    if ((Snap() && SnapMinmaxRegions() == NULL) ||
	(!Snap() && MinmaxRegions() == NULL)) {
      error.AddMessage("Minmax regions not allocated.");
      return(false);
    }
  }

  if (UseNEP() || Snap()) {
    if (!isotable.cube_nep.IsTableAllocated()) {
      error.AddMessage("Cube NEP (negative-equals-positive) isosurface lookup table is not allocated.");
      return(false);
    }

    if (!isotable.cube_nep.Dimension() == ScalarGrid().Dimension()) {
      error.AddMessage("Cube NEP (negative-equals-positive) isosurface lookup table dimension");
      error.AddMessage("  does not match scalar grid dimension.");
      return(false);
    }

    if (!check_isotable_encoding
	(isotable.cube_nep, ISOSURFACE_TABLE::BASE3, error)) {
      return(false);
    }
  }
  else {
    if (!isotable.cube.IsTableAllocated()) {
      error.AddMessage("Cube isosurface lookup table is not allocated.");
      return(false);
    }

    if (!isotable.cube_nep.Dimension() == ScalarGrid().Dimension()) {
      error.AddMessage("Cube isosurface lookup table dimension");
      error.AddMessage("  does not match scalar grid dimension.");
      return(false);
    }

    if (IntervalVolumeFlag()) {
      if (!check_isotable_encoding
	  (isotable.cube, ISOSURFACE_TABLE::BASE3, error)) {
	return(false);
      }
    }
    else {
      if (!check_isotable_encoding
	  (isotable.cube, ISOSURFACE_TABLE::BINARY, error)) {
	return(false);
      }
    }
  }

  if (UseMultires()) {
    if (!isotable.pyramid.IsTableAllocated()) {
      error.AddMessage("Pyramid isosurface lookup table is not allocated.");
      return(false);
    }

    if (!isotable.simplex.IsTableAllocated()) {
      error.AddMessage("Simplex isosurface lookup table is not allocated.");
      return(false);
    }
  }


  if (IsosurfaceTopology() == CUBE_DECIDER_TOPOLOGY ||
      IsosurfaceTopology() == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      IsosurfaceTopology() == SADDLE_TOPOLOGY ||
      IsosurfaceTopology() == LINEAR_TOPOLOGY) {

    if (!check_isotable_ambig_info
	(isotable.cube, isotable.ambig_cube, error)) {
      return(false);
    }
  }

  if (IsosurfaceTopology() == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      IsosurfaceTopology() == SADDLE_TOPOLOGY ||
      IsosurfaceTopology() == LINEAR_TOPOLOGY) {

    if (!isotable.pyramid.IsTableAllocated()) {
      error.AddMessage("Pyramid isosurface lookup table is not allocated.");
      return(false);
    }

    if (!check_isotable_encoding
	(isotable.pyramid, ISOSURFACE_TABLE::BINARY, error)) {
      return(false);
    }

    if (!check_isotable_ambig_info
	(isotable.pyramid, isotable.ambig_pyramid, error)) {
      return(false);
    }

  }

  if (IsosurfaceTopology() == SADDLE_TOPOLOGY ||
      IsosurfaceTopology() == LINEAR_TOPOLOGY) {

    if (!isotable.simplex.IsTableAllocated()) {
      error.AddMessage("Simplex isosurface lookup table is not allocated.");
      return(false);
    }

    if (!check_isotable_encoding
	(isotable.simplex, ISOSURFACE_TABLE::BINARY, error)) {
      return(false);
    }
  }

  return(true);
}


// **************************************************
// MCUBE TIME
// **************************************************

IJKMCUBE::MCUBE_TIME::MCUBE_TIME()
// constructor
{
  Clear();
}

void IJKMCUBE::MCUBE_TIME::Clear()
{
  preprocessing = 0.0;
  process_multires = 0.0;
  snap = 0.0;
  extract = 0.0;
  merge = 0.0;
  position = 0.0;
  total = 0.0;
}

void IJKMCUBE::MCUBE_TIME::Add(const MCUBE_TIME & mcube_time)
{
  preprocessing += mcube_time.preprocessing;
  process_multires += mcube_time.process_multires;
  snap += mcube_time.snap;
  extract += mcube_time.extract;
  merge += mcube_time.merge;
  position += mcube_time.position;
  total += mcube_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

IJKMCUBE::GRID_INFO::GRID_INFO()
{
  Clear();
}

void IJKMCUBE::GRID_INFO::Clear()
{
  num_cubes = 0;

}

void IJKMCUBE::SCALAR_INFO::Init(const int dimension)
{
  num_cubes_with_saddle = NULL;
  this->dimension = 0;
  SetDimension(dimension);

  Clear();
}

void IJKMCUBE::SCALAR_INFO::FreeAll()
{
  delete [] num_cubes_with_saddle;
  num_cubes_with_saddle = NULL;
  dimension = 0;
}

void IJKMCUBE::SCALAR_INFO::Clear()
{
  num_non_empty_cubes = 0;
  num_bipolar_edges = 0;
  num_ambiguous_cubes = 0;
  num_non_empty_pyramids = 0;
  num_ambiguous_pyramids = 0;

  for (int d = 0; d < dimension; d++)
    { num_cubes_with_saddle[d] = 0; };
}

void IJKMCUBE::SCALAR_INFO::SetDimension(const int dimension)
{
  FreeAll();

  this->dimension = dimension;
  num_cubes_with_saddle = new VERTEX_INDEX[dimension];

  Clear();
}

void IJKMCUBE::SCALAR_INFO::Copy(const SCALAR_INFO & info)
{
  Init(info.Dimension());
  num_non_empty_cubes = info.num_non_empty_cubes;
  num_bipolar_edges = info.num_bipolar_edges;
  num_ambiguous_cubes = info.num_ambiguous_cubes;
  num_non_empty_pyramids = info.num_non_empty_pyramids;
  num_ambiguous_pyramids = info.num_ambiguous_pyramids;

  for (int d = 0; d < Dimension(); d++) 
    { SetNumCubesWithSaddle(d, info.NumCubesWithSaddle(d));  }
}

/// Copy assignment.
const SCALAR_INFO &  IJKMCUBE::SCALAR_INFO::operator =
(const SCALAR_INFO & right)
{
  if (&right != this) {
    FreeAll();
    Copy(right);
  }
}

IJKMCUBE::SCALAR_INFO::~SCALAR_INFO()
{
  delete [] num_cubes_with_saddle;
  num_cubes_with_saddle = NULL;
  dimension = 0;
  Clear();
}

IJKMCUBE::MCUBE_INFO::MCUBE_INFO()
{
  Clear();
}

IJKMCUBE::MCUBE_INFO::MCUBE_INFO(const int dimension):scalar(dimension)
{
  Clear();
}

void IJKMCUBE::MCUBE_INFO::Clear()
{
  grid.Clear();
  scalar.Clear();
  time.Clear();
  nep.Clear();
}

void IJKMCUBE::NEP_INFO::Clear()
{
  num_in_facet_cubes = 0;
  num_dup_iso_patches = 0;
  num_non_empty_boundary_facets = 0;
}

// **************************************************
// MERGE DATA
// **************************************************

void IJKMCUBE::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  this->num_obj_per_grid_vertex = 
    dimension*num_obj_per_edge + num_obj_per_vertex;
  compute_num_grid_vertices(dimension, axis_size, num_vertices);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj = 
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX,MERGE_INDEX>::Init(num_obj);
}

bool IJKMCUBE::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxNumInt() < 
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}

bool IJKMCUBE::MERGE_DATA::Check
(const ISOTABLE_TYPE isotable_type, ERROR & error) const
{
  if (isotable_type == NEP) {
    if (NumObjPerVertex() < 1) {
      error.AddMessage("Programming error:  MERGE_DATA has too few object per vertex.");
      error.AddMessage
	("  NEP isosurface lookup table may position isosurface vertices ");
      error.AddMessage("  at grid vertices.");
      error.AddMessage
	("  Number of objects per vertex in MERGE_DATA should be at least 1.");

      return(false);
    }
  }

  return(true);
}


// **************************************************
// NEP_ISOSURFACE_TABLE
// **************************************************

NEP_ISOSURFACE_TABLE::NEP_ISOSURFACE_TABLE(const int d):ISOSURFACE_TABLE(d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d, d-1);
}

NEP_ISOSURFACE_TABLE::NEP_ISOSURFACE_TABLE
(const int dimension, const int simplex_dimension):
  ISOSURFACE_TABLE(dimension, simplex_dimension)
// constructor
// dimension = dimension of space containing isosurface.  Should be 2, 3 or 4.
// simplex_dimension = dimension of isosurface simplices
//        simplex_dimension equals dimension-1 for isosurfaces.
//        simplex_dimension equals dimension for interval volumes
{
  Init(dimension, simplex_dimension);
}


void NEP_ISOSURFACE_TABLE::Init
(const int dimension, const int simplex_dimension)
// initialization routine.  Called by constructors.
{
  ISOSURFACE_TABLE::Init(dimension, simplex_dimension);

  is_in_facet = NULL;
  containing_facet = NULL;
}

NEP_ISOSURFACE_TABLE::~NEP_ISOSURFACE_TABLE()
{
  FreeAll();
}

void NEP_ISOSURFACE_TABLE::SetIsInFacet(const bool flag)
// set is_in_facet[k] = flag, for all table entries k
{
  for (IJKTABLE::TABLE_INDEX it = 0; it < NumTableEntries(); it++)
    { is_in_facet[it] = flag; };
}

void NEP_ISOSURFACE_TABLE::SetNumTableEntries(const int num_table_entries)
{
  const char * procname = "NEP_ISOSURFACE_TABLE::SetNumTableEntries";

  ISOSURFACE_TABLE::SetNumTableEntries(num_table_entries);

  is_in_facet = new bool[num_table_entries];
  containing_facet = new IJKTABLE::FACET_INDEX[num_table_entries];

  if (is_in_facet == NULL || containing_facet == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for NEP isosurface table.");
}

void NEP_ISOSURFACE_TABLE::FreeAll()
{
  ISOSURFACE_TABLE::FreeAll();

  delete [] is_in_facet;
  is_in_facet = NULL;
  delete [] containing_facet;
  containing_facet = NULL;
}

// **************************************************
// CLASS SNAP_MINMAX
// **************************************************

void IJKMCUBE::SNAP_MINMAX::ComputeMinMax
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const SCALAR_TYPE * scalar, const AXIS_SIZE_TYPE region_edge_length)
{
  // Offset each region by one cube to include all snapped cubes
  const AXIS_SIZE_TYPE offset_length = 1;
  MINMAX_REGIONS::ComputeMinMax(dimension, axis_size, scalar,
				region_edge_length, offset_length);
}

void IJKMCUBE::SNAP_MINMAX::ComputeMinMax
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 AXIS_SIZE_TYPE region_edge_length)
{
  ComputeMinMax(scalar_grid.Dimension(), scalar_grid.AxisSize(),
		scalar_grid.ScalarPtrConst(), region_edge_length);
}

// **************************************************
// CLASS SNAP_OCTREE
// **************************************************

/// Offset each leaf node by one cube to include all snapped cubes
void IJKMCUBE::SNAP_OCTREE::SetMinMax(const SCALAR_TYPE * scalar)
{
  if (NumLevels() == 0) { return; };

  SNAP_MINMAX minmax;

  minmax.ComputeMinMax(Dimension(), AxisSize(), scalar, 2);

  OCTREE::SetMinMax(minmax);
}

