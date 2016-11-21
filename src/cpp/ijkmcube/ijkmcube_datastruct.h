/// \file ijkmcube_datastruct.h
/// Data structure definitions for ijkmcube.

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

#ifndef _IJKMCUBE_DATASTRUCT_
#define _IJKMCUBE_DATASTRUCT_

#include <string>

#include "ijk.txx"
#include "ijkgrid.txx"
#include "ijkmerge.txx"


#include "ijkmcube_types.h"
#include "ijktable.h"
#include "ijkoctree.h"

namespace IJKMCUBE {

  using IJKTABLE::ISOSURFACE_TABLE;
  using IJKTABLE::ISOSURFACE_TABLE_AMBIG_INFO;
  using IJKTABLE::ISOSURFACE_TABLE_POLYHEDRON;
  using IJKTABLE::TABLE_INDEX;

  class MERGE_DATA;
  class SNAP_MINMAX;
  class SNAP_OCTREE;
  class MULTIRES_GRID;

// **************************************************
// GRID DATA STRUCTURES
// **************************************************

  typedef IJK::GRID<int, AXIS_SIZE_TYPE, VERTEX_INDEX, VERTEX_INDEX> 
    MC_GRID;                  ///< Marching Cubes grid.

  typedef IJK::SCALAR_GRID_BASE<MC_GRID, SCALAR_TYPE> 
    MC_SCALAR_GRID_BASE;      ///< Marching Cubes base scalar grid.
  typedef IJK::SCALAR_GRID_WRAPPER<MC_GRID, SCALAR_TYPE>
    MC_SCALAR_GRID_WRAPPER;   ///< Marching Cubes scalar grid wrapper.
  typedef IJK::SCALAR_GRID<MC_GRID, SCALAR_TYPE> 
    MC_SCALAR_GRID;           ///< Marching Cubes scalar grid.
  typedef IJK::MINMAX_REGIONS<MC_GRID, SCALAR_TYPE> 
    MINMAX_REGIONS;           ///< Marching Cubes structures for computing minimum/maximum in regions.

// **************************************************
// NEP_ISOSURFACE_TABLE
// **************************************************

  /// NEP isosurface lookup table.
  /// Includes information about isosurface patches contained on polyhedron facets.
  class NEP_ISOSURFACE_TABLE:public ISOSURFACE_TABLE {

  protected:
    bool * is_in_facet;   // is_in_facet[k] = true if iso patch is in facet
    IJKTABLE::FACET_INDEX * containing_facet;
                          // containing_facet[k] = facet containing iso patch

    // initialization routine
    void Init(const int dimension, const int simplex_dimension);

  public:
    // constructors
    NEP_ISOSURFACE_TABLE() { Init(3,2); };
    NEP_ISOSURFACE_TABLE(const int d);
    NEP_ISOSURFACE_TABLE(const int dimension, const int simplex_dimension);

    ~NEP_ISOSURFACE_TABLE();                // destructor

    // get functions
    bool IsInFacet(const TABLE_INDEX it) const
    { return(is_in_facet[it]); };
    IJKTABLE::FACET_INDEX ContainingFacet(const TABLE_INDEX it) const
      { return(containing_facet[it]); };

    // set functions

    // set is_in_facet[k] = flag, for all table entries k
    void SetIsInFacet(const bool flag); 

    void SetIsInFacet(const TABLE_INDEX it, const bool flag)
    { is_in_facet[it] = flag; };
    void SetContainingFacet(const TABLE_INDEX it, 
			    const IJKTABLE::FACET_INDEX jf)
    { is_in_facet[it] = true; containing_facet[it] = jf; };

    virtual void SetNumTableEntries(const int num_table_entries);

    // free memory
    virtual void FreeAll();                     // free all memory
  };

// **************************************************
// CLASS POLY_ISOTABLE
// **************************************************

  /// Isosurface lookup tables for polyhedra.
  class POLY_ISOTABLE {

  public:
    /// NOTE: Not all isosurface tables or ambiguity information
    ///   is set or used.
    ISOSURFACE_TABLE cube;                     ///< Cube isotable.
    ISOSURFACE_TABLE pyramid;                  ///< Pyramid isotable.
    ISOSURFACE_TABLE simplex;                  ///< Simplex isotable.
    ISOSURFACE_TABLE_AMBIG_INFO ambig_cube;    ///< Cube ambiguity info.
    ISOSURFACE_TABLE_AMBIG_INFO ambig_pyramid; ///< Pyramid ambiguity info.

    NEP_ISOSURFACE_TABLE cube_nep;             ///< Cube NEP isotable.

  public:
    POLY_ISOTABLE() {};
    void ComputeAmbiguityInformation();
  };


// **************************************************
// MARCHING CUBES ISOSURFACE CLASS
// **************************************************

  /// Marching cubes isosurface.
  /// Representation of isosurface (or interval volume)
  ///   returned by Marching Cubes and its variants.
  class MC_ISOSURFACE {

  public:
    /// List of simplex vertices.
    VERTEX_INDEX_ARRAY simplex_vert;

    /// List of vertex coordinates.
    COORD_ARRAY vertex_coord;

    /// cube_containing_simplex[i] = cube containing i'th simplex.
    /// Not always set.
    VERTEX_INDEX_ARRAY cube_containing_simplex;

  public:
    MC_ISOSURFACE() {};

    void Clear();
  };

// **************************************************
// MARCHING CUBES INPUT DATA AND DATA STRUCTURES
// **************************************************

  /// Marching cubes flags.
  class MC_DATA_FLAGS {

  protected:
    void Init();

  public:
    bool interval_volume_flag;
    MESH_EDGE_REPRESENTATION edge_representation;
    bool use_minmax;
    bool use_octree;
    bool use_nep;
    bool use_list;            ///< Extract non-empty cubes to list.
    bool use_multires;
    bool snap_flag;
    bool cube_containing_simplex_flag;

    // control parameters
    int nep_num_dup;
    SNAP_TYPE snap_value;
    ISOSURFACE_TOPOLOGY isosurface_topology;
    INTERPOLATION_TYPE interpolation_type;

    // multiresolution parameters
    int num_resolution_levels;
    std::vector<GRID_BOX> high_resolution_regions;

  public:
    MC_DATA_FLAGS() { Init(); };
    ~MC_DATA_FLAGS() { Init(); };
  };

  /// Input data to Marching Cubes and related algorithms
  class MC_DATA:protected MC_DATA_FLAGS {

  protected:
    MC_SCALAR_GRID scalar_grid;
    MINMAX_REGIONS * minmax;
    IJKOCTREE::OCTREE * octree;
    SNAP_MINMAX * snap_minmax;
    SNAP_OCTREE * snap_octree;
    MULTIRES_GRID * multires_grid;

    // flags
    bool is_scalar_grid_set;

    void Init();
    void FreeAll();

  public:
    POLY_ISOTABLE isotable;
    MERGE_EDGES_PARAMETERS merge_edges_parameters;

  public:
    MC_DATA() { Init(); };
    ~MC_DATA() { FreeAll(); };

    // Set functions
    void CopyScalarGrid             /// Copy scalar_grid to MC_DATA
      (const MC_SCALAR_GRID_BASE & scalar_grid2);
    void SubsampleScalarGrid        /// Subsample scalar_grid.
      (const MC_SCALAR_GRID_BASE & scalar_grid2, 
       const int subsample_resolution);
    void SupersampleScalarGrid      /// Supersample scalar_grid.
      (const MC_SCALAR_GRID_BASE & scalar_grid2, 
       const int supersample_resolution);
    void SetOctree();               ///< Set and use octree.
    void SetMinmaxRegions           /// Set and use minmax.
      (const int region_length); 
    void SetSnapOctree();           ///< Set and use snap octree.
    void SetSnapMinmaxRegions       /// Set and use snap minmax regions.
      (const int region_length); 
    void SetMultiresGrid            /// Set and use multires grid.
      (const int num_resolution_levels);
    void SetEdgeRepresentation      /// Set mesh edge representation
      (const MESH_EDGE_REPRESENTATION r);
    void SetNEPOn                   /// Turn NEP (neg-equals-pos) on
      (const int nep_num_dup);
    void SetNEPOff();               ///< Turn NEP (neg-equals-pos) off
    void SetIntervalVolumeFlag      /// Set interval_volume_flag.
      (const bool flag);
    void SetSnapOn                  /// Turn snap on
      (const SNAP_TYPE snap_value, const int nep_num_dup);
    void SetSnapOff();              ///< Turn snap on
    void SetIsosurfaceTopology      /// Set desired isosurface topology.
      (const ISOSURFACE_TOPOLOGY isosurface_topology);
    void SetInterpolationType       /// Set type of interpolation.
      (const INTERPOLATION_TYPE interpolation_type);
    void SetUseList(const bool flag);  ///< Set use_list.
    void SetHighResolutionRegions   /// Set high resolution regions.
      (const std::vector<GRID_BOX> & high_resolution_regions);
    void SetCubeContainingSimplexFlag  /// Set cube_containing_simplex flag.
      (const bool flag);

    /// Copy, subsample or supersample scalar grid.
    /// Precondition: flag_subsample and flag_supersample are not both true.
    void SetScalarGrid
      (const MC_SCALAR_GRID_BASE & scalar_grid2, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution);

    // Get functions
    bool IsScalarGridSet() const     /// Return true if scalar grid is set.
      { return(is_scalar_grid_set); };
    bool EdgeRepresentation() const  /// Return type of edge representation.
      { return(edge_representation); };
    bool UseMinmaxRegions() const            /// Return minmax flag.
      { return(use_minmax); };
    bool UseOctree() const                   /// Return octree flag.
      { return(use_octree); };
    bool UseMultires() const                 /// Return multires flag.
      { return(use_multires); };
    const MC_SCALAR_GRID_BASE & ScalarGrid() const /// Return scalar_grid
      { return(scalar_grid); };
    const IJKOCTREE::OCTREE * Octree() const /// Return pointer to octree.
      { return(octree); };
    const MINMAX_REGIONS * MinmaxRegions() const    
      { return(minmax); };          ///< Return pointer to minmax.
    const SNAP_OCTREE * SnapOctree() const 
      { return(snap_octree); };     ///< Return pointer to snap_octree.
    const SNAP_MINMAX * SnapMinmaxRegions() const    
      { return(snap_minmax); };     ///< Return pointer to minmax.
    const MULTIRES_GRID * MultiresGrid() const
      { return(multires_grid); };   ///< Return pointer to multires_grid.
    bool UseNEP() const             /// Return NEP flag.
      { return(use_nep); };
    bool UseList() const            /// Return use_list flag.
      { return(use_list); };
    int NEPNumDup() const           /// Return nep_num_dup
      { return(nep_num_dup); };
    bool IntervalVolumeFlag() const /// Return interval_volume_flag
      { return(interval_volume_flag); };
    bool Snap() const               /// Return true if snap is on.
      { return(snap_flag); };
    SNAP_TYPE SnapValue() const     /// Return snap_value
      { return(snap_value); };

    /// Return isosurface topology.
    ISOSURFACE_TOPOLOGY IsosurfaceTopology() const
      { return(isosurface_topology); };     

    /// Return interpolation type.
    INTERPOLATION_TYPE InterpolationType() const
      { return(interpolation_type); };      

    /// Return flag indicating if algorithm should return the cube 
    /// containing each isosurface simplex.
    bool CubeContainingSimplexFlag() const
      { return(cube_containing_simplex_flag); };

    /// Check data structure
    bool Check(IJK::ERROR & error) const;
  };


// **************************************************
// CLASS MULTIRES_GRID
// **************************************************

  typedef IJK::SCALAR_GRID<MC_GRID, MULTIRES_VERTEX_TYPE>
    MULTIRES_GRID_BASE;   ///< Marching Cubes base multiresolution grid.

  /// Multi-resolution grid.
  class MULTIRES_GRID:public MULTIRES_GRID_BASE {

  protected:
    bool is_processed;
    void Init();

  public:
    MULTIRES_GRID(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      IJK::SCALAR_GRID<MC_GRID, MULTIRES_VERTEX_TYPE>
      (dimension, axis_size) 
      { Init(); };

    // get functions
    inline bool IsMultires(const VERTEX_INDEX iv) const
      { return(Scalar(iv) != NOT_MULTIRES); };

    // set functions
    void Process();
    void SetAllMultires(const MULTIRES_VERTEX_TYPE vtype);
      // set all multires vertices to vtype
    void SetCorners2Multires(); 
       // set corners to multires vertices
    void SetSubsample2Multires(const AXIS_SIZE_TYPE scale); 
      // set subsampled vertices to multires vertices
    void SetRegion2Multires(const GRID_BOX & box, const AXIS_SIZE_TYPE scale);
      // set scalar values of vertices in region to multires vertices
    void SetRegion2Multires(const std::vector<GRID_BOX> & box_list, 
			    const AXIS_SIZE_TYPE scale);
      // set scalar values of vertices in list of regions to multires vertices

    // get functions
    bool IsProcessed() const { return(is_processed); };
  };


// **************************************************
// MCUBE TIME
// **************************************************

  /// Marching cubes time.
  /// Uses system clock() function to determine time.  
  /// System clock() function should return cpu time, 
  ///   but may return wall time.
  class MCUBE_TIME {

  public:
    // all times are in seconds
    float preprocessing;  
      // time to create data structure for faster isosurface extraction
    float process_multires;  // time to process multires grid
    float snap;       // time to snap grid scalar values 
    float extract;    // time to extract isosurface mesh
    float merge;      // time to merge identical vertices
    float position;   // time to position isosurface vertices
    float total;      // extract_time+merge_time+position_time

    MCUBE_TIME();
    void Clear();
    void Add(const MCUBE_TIME & mcube_time);
  };

// **************************************************
// GRID INFO
// **************************************************

  /// Regular grid information.
  class GRID_INFO {

  public:
    GRID_INFO();                    ///< Constructor.

    VERTEX_INDEX num_cubes;         ///< Number of grid cubes.

    void Clear();                   ///< Clear all data.
  };

// **************************************************
// SCALAR INFO
// **************************************************

  /// Scalar grid information.
  class SCALAR_INFO {

  protected:
    int dimension;
    VERTEX_INDEX * num_cubes_with_saddle; ///< Number of ambiguous cubes with i saddles.

    void Copy(const SCALAR_INFO & info);  ///< Copy scalar info.
    void Init(const int dimension);       ///< Initialize scalar info.
    void FreeAll();                       ///< Free all memory.

  public:
    SCALAR_INFO() { Init(3); };
    SCALAR_INFO(const int dimension) { Init(dimension); };
    ~SCALAR_INFO();
    SCALAR_INFO(const SCALAR_INFO & info) ///< Copy constructor. 
      { Copy(info); };       
    const SCALAR_INFO & operator =        ///< Copy assignment.
      (const SCALAR_INFO & right);

    VERTEX_INDEX num_non_empty_cubes;     ///< Number of cubes containing an isosurface simplex.

    VERTEX_INDEX num_bipolar_edges; ///< Number of bipolar edges.
    //   (number of edges with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)

    VERTEX_INDEX num_ambiguous_cubes; ///< Number of ambiguous cubes.
    VERTEX_INDEX num_non_empty_pyramids; ///< Number of non empty pyramids.
    VERTEX_INDEX num_ambiguous_pyramids; ///< Number of ambiguous pyramids.


    // Set functions.
    void SetDimension(const int dimension); ///< Set dimension.
    void SetNumCubesWithSaddle(const int i, const int n)
      { num_cubes_with_saddle[i] = n; };
    void IncrementNumCubesWithSaddle(const int i)
      { num_cubes_with_saddle[i]++; };


    // Get functions.
    int Dimension() const { return(dimension); };
    int MaxNumSaddles() const { return(dimension-1); };

    /// Return number of cubes with i saddles.
    int NumCubesWithSaddle(const int i) const
      { return(num_cubes_with_saddle[i]); };

    void Clear();     // clear all data
  };

// **************************************************
// MCUBE INFO
// **************************************************

  /// NEP information.
  /// Information on extracting isosurfaces using NEP (negative-equals-positive)
  ///   isosurface lookup table.
  class NEP_INFO {
  public:
    
    VERTEX_INDEX num_in_facet_cubes;
    // number of cubes whose isosurface patch lies in a cube facet

    VERTEX_INDEX num_dup_iso_patches;
    // number of duplicate isosurface patches
    // each pair of duplicate patches is counted only once

    VERTEX_INDEX num_non_empty_boundary_facets;
    // number of facets containing an isosurface patch

    NEP_INFO() { Clear(); };
    void Clear();     // clear all data
  };


  /// Marching cubes information.
  /// Statistical and timing information from the Marching Cubes algorithm.
  class MCUBE_INFO {

  public:
    GRID_INFO grid;
    SCALAR_INFO scalar;
    MCUBE_TIME time;
    NEP_INFO nep;

    MCUBE_INFO();
    MCUBE_INFO(const int dimension);

    void Clear();     // clear all data
  };

// **************************************************
// MERGE DATA
// **************************************************

/// Internal data structure for merge_identical_vertices 
  class MERGE_DATA: public IJK::INTEGER_LIST<MERGE_INDEX, MERGE_INDEX> {

  protected:
    MERGE_INDEX num_edges;             ///< Number of edges.
    MERGE_INDEX num_vertices;          ///< Number of vertices.
    MERGE_INDEX num_obj_per_vertex;    ///< Number of objects per vertex.
    MERGE_INDEX num_obj_per_edge;      ///< Number of objects per edge.
    MERGE_INDEX num_obj_per_grid_vertex; ///< Number of objects per grid vertex.
    MERGE_INDEX vertex_id0;            ///< First vertex identifier.

    /// Initialize.
    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size,
	      const MERGE_INDEX num_obj_per_vertex,
	      const MERGE_INDEX num_obj_per_edge);

  public:
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size)
      { Init(dimension, axis_size, 0, 1); };
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size,
	       const MERGE_INDEX num_obj_per_vertex, 
	       const MERGE_INDEX num_obj_per_edge)
      { Init(dimension, axis_size, num_obj_per_vertex, num_obj_per_edge); };

    // get functions
    MERGE_INDEX NumEdges() const        /// Number of edges.
      { return(num_edges); };
    MERGE_INDEX NumVertices() const     /// Number of vertices.
      { return(num_vertices); };
    MERGE_INDEX NumObjPerVertex() const { return(num_obj_per_vertex); };
    MERGE_INDEX NumObjPerEdge() const { return(num_obj_per_edge); };
    MERGE_INDEX NumObjPerGridVertex() const
      { return(num_obj_per_grid_vertex); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv) const { return(vertex_id0 + iv); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv, const MERGE_INDEX j) const 
      { return(vertex_id0 + iv + j * num_vertices); };
    MERGE_INDEX EdgeIdentifier         /// Edge identifier.
      (const MERGE_INDEX ie) const { return(ie); };
    MERGE_INDEX EdgeIdentifier         /// edge identifier
      (const MERGE_INDEX ie, const MERGE_INDEX j) const 
    { return(ie + j * num_edges); };

    /// Get first endpoint of edge containing isosurface vertex isov.
    inline VERTEX_INDEX GetFirstEndpoint(const MERGE_INDEX isov) const
      { return(isov/NumObjPerGridVertex()); };

    /// Get direction of edge containing isosurface vertex isov.
    inline MERGE_INDEX GetEdgeDir(const MERGE_INDEX isov) const
      { return(isov%NumObjPerGridVertex()); };

    bool Check(IJK::ERROR & error) const;     ///< Check allocated memory.
    bool Check     ///  Check data structure for given isotable type.
      (const ISOTABLE_TYPE type, IJK::ERROR & error) const;
  };

  /// Merge data structure for interval volume vertices.
  /// Vertices are identified by a single integer.
  class IVOL_MERGE_DATA: public MERGE_DATA {
  public:
    IVOL_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 2) {};
  };

  /// Merge data structure for isosurface vertices.
  /// Vertices are identified by a single integer.
  class ISO_MERGE_DATA: public MERGE_DATA {
  public:
    ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 0, 1) {};
  };

  /// Merge data structure for NEP isosurface vertices.
  /// Vertices are identified by a single integer.
  class NEP_ISO_MERGE_DATA: public MERGE_DATA {
  public:
    NEP_ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 1) {};
  };

// **************************************************
// MESH VERTEX LIST
// **************************************************

  /// List of mesh vertices.
  template <class DTYPE, class CTYPE, class STYPE>
  class MESH_VERTEX_LIST {

  protected:
    DTYPE dimension;        ///< Dimension.

    /// Coordinates of mesh vertices.
    /// coord[i*dimension+j] = j'th coordinate of i'th mesh vertex.
    std::vector<CTYPE> coord;
    std::vector<STYPE> scalar;   ///< scalar[i] = scalar value of i'th mesh vertex.

  public:
    MESH_VERTEX_LIST(const DTYPE dim) { this->dimension = dim; };

    // Get functions.
    DTYPE Dimension() const { return(dimension); } ///< Return dimension.

    /// Return j'th coordinate of i'th mesh vertex.
    CTYPE Coord(const int i, const int j) const
    { return(coord[i*Dimension()+j]); };

    /// Copy coordinates of \a i'th mesh vertex into array \a coord[].
    /// @pre coord[] is preallocated to length at least \a dimension.
    template <class CTYPE2>
      void CopyCoord(const int i, CTYPE2 * coord2) const
      {
	const CTYPE * coord_begin = &(coord[i*Dimension()]);
	std::copy(coord_begin, coord_begin+Dimension(), coord2);
      };

    /// Return scalar value of i'th mesh vertex.
    STYPE Scalar(const int i) const { return(scalar[i]); };

    int NumVertices() const     ///< Return number of vertices in list.
    { return(scalar.size()); };

    // Set functions.

    void Add()                 ///< Add a vertex to the list.
    {
      const int k = coord.size();
      coord.resize(k+Dimension());
      scalar.push_back(0);
    }

    void SetScalar(const int i, const STYPE s)  ///< Set scalar of vertex i.
    { scalar[i] = s; };

    /// Set j'th coordinate of vertex i.
    void SetCoord(const int i, const int j, const CTYPE c)
    { coord[i*Dimension()+j] = c; };

    /// Set coordinates of vertex i
    template <class CTYPE2>
      void SetCoord(const int i, CTYPE2 * coord2)
      {
	for (DTYPE d = 0; d < Dimension(); d++)
	  { SetCoord(i, d, coord2[d]); };
      };

  };

  typedef MESH_VERTEX_LIST<int, COORD_TYPE, SCALAR_TYPE> 
    MC_MESH_VERTEX_LIST;    ///< Marching Cubes mesh vertex list.

// **************************************************
// CLASS SNAP_MINMAX
// **************************************************

  /// MINMAX_REGIONS for snapMC
  class SNAP_MINMAX:public MINMAX_REGIONS {

  public:
    SNAP_MINMAX(){};   ///< Constructor.

    /// Compute min and max of each region.
    /// Offset each region by one cube to include all snapped cubes.
    void ComputeMinMax
    (const int dimension, const AXIS_SIZE_TYPE * axis_size,
     const SCALAR_TYPE * scalar, const AXIS_SIZE_TYPE region_edge_length);

    /// Compute min and max of each region.
    /// Offset each region by one cube to include all snapped cubes.
    void ComputeMinMax
    (const MC_SCALAR_GRID_BASE & scalar_grid,
     AXIS_SIZE_TYPE region_edge_length);
  };

// **************************************************
// CLASS SNAP_OCTREE
// **************************************************

  /// OCTREE for snapMC
  class SNAP_OCTREE:public IJKOCTREE::OCTREE {

  public:
    /// Constructor
    SNAP_OCTREE(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      IJKOCTREE::OCTREE(dimension, axis_size) {};

    /// Offset each leaf node by one cube to include all snapped cubes
    void SetMinMax(const SCALAR_TYPE * scalar);
  };

};

#endif
