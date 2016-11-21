/// \file ijktable.h
/// Class containing a table of isosurface patches in a given polyhedron.
/// All 2^numv +/- patterns are stored in the table 
///   where numv = # polyhedron vertices.
/// Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2007, 2006, 2003, 2001 Rephael Wenger

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

#ifndef _IJKTABLE_
#define _IJKTABLE_

#include <iostream>

#include "ijk.txx"

/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKTABLE {

  typedef unsigned char 
    ISOSURFACE_VERTEX_INDEX;  ///< Index of isosurface vertex.
  typedef unsigned char EDGE_INDEX;    ///< Index of edge.
  typedef unsigned char FACET_INDEX;   ///< Index of facet.
  typedef int TABLE_INDEX;    ///< Index of entry in isosurface lookup table.
  typedef int FACET;          ///< Bits representing vertices in facet.
  typedef int FACET_SET;      ///< Bits representing set of facets.

  const int NO_VERTEX = -1;

//**************************************************
// ISOSURFACE TABLE POLYHEDRON
//**************************************************

/// Isosurface table polyhedron.
class ISOSURFACE_TABLE_POLYHEDRON {

 protected:
  int dimension;         ///< Polyhedron dimension.
  int num_vertices;      ///< Number of polyhedron vertices.
  int num_edges;         ///< Number of polyhedron edges.
  int num_facets;        ///< Number of polyhedron facets.
  int * vertex_coord;    ///< Polyhedron vertex coordinates.
  int * edge_endpoint;   ///< Polyhedron edge endpoints.
  int * num_facet_vertices; ///< Number of vertices of each facet.
  int ** facet_vertex_list; ///< List of vertices in each facet.
  FACET * facet;         ///< Polyhedron facets.
  void Init();           ///< Initialize.
  void FreeFacets();     ///< Free all facet arrays.

 public:
  ISOSURFACE_TABLE_POLYHEDRON(const int d);  ///< Constructor
  ~ISOSURFACE_TABLE_POLYHEDRON();            ///< Destructor
  ISOSURFACE_TABLE_POLYHEDRON
    (const ISOSURFACE_TABLE_POLYHEDRON & init);  ///< Copy constructor.
  const ISOSURFACE_TABLE_POLYHEDRON & operator = 
    (const ISOSURFACE_TABLE_POLYHEDRON &);  ///< Assignment.

  /// @name Get Functions
  //@{
  int Dimension() const   ///< Polyhedron dimension.
    { return(dimension); };
  int NumVertices() const ///< Number of polyhedron vertices.
    { return(num_vertices); };
  int NumEdges() const    ///< Number of polyhedron edges.
    { return(num_edges); };
  int NumFacets() const   ///< Number of polyhedron facets.
    { return(num_facets); };
  int NumFacetVertices    ///< Number of facet vertices of facet \a jf.
    (const FACET_INDEX jf) const
    { return(num_facet_vertices[jf]); };
  int VertexCoord         ///< \a ic'th vertex coordinate of vertex \a iv.
    (const int iv, const int ic) const
    { return(vertex_coord[iv*dimension + ic]); };
  int EdgeEndpoint        ///< \a j'th endpoint of edge \a ie. \a j = 0 or 1.
    (const EDGE_INDEX ie, const int j) const
    { return(edge_endpoint[int(ie)*2 + j]); };
  int MidpointCoord       ///< \a ic'th coordinate of midpoint of edge \a ie.
    (const EDGE_INDEX ie, const int ic) const;
  FACET Facet             ///< Bits representing vertices in facet \a jf.
    (const FACET_INDEX jf) const
    { return(facet[jf]); };
  bool IsVertexInFacet    ///< Return true if vertex \a iv is in facet \a jf.
    (const FACET_INDEX jf, const int iv) const
    { return(facet[jf] & ((1L) << iv)); };
  int FacetVertex         ///< Return \a k'th vertex in facet \a jf.
    (const FACET_INDEX jf, const int k) const
    { return(facet_vertex_list[jf][k]); };
  //@}

  /// @name Set Functions
  //@{
  void SetDimension(const int d);      ///< Set polyhedron dimension.
  void SetNumVertices(const int numv); ///< Set number of polyhedron vertices.
  void SetNumEdges(const int nume);    ///< Set number of polyhedron edges.
  void SetNumFacets(const int numf);   ///< Set number of polyhedron facets.

  /// Set number of polyhedron vertices, edges and facets.
  void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };

  /// Set \a ic'th coordinate of vertex \a iv.
  /// @pre SetNumVertices or SetSize must be called before SetVertexCoord.
  void SetVertexCoord 
    (const int iv, const int ic, const int coord);

  /// Set endpoints of edge \a ie.
  /// @pre SetNumEdges or SetSize must be called before SetEdge.
  void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1);

  /// Set number of vertices in facet \a jf.
  /// @pre SetNumFacets or SetSize must be called before SetNumFacetVertices.
  void SetNumFacetVertices 
    (const FACET_INDEX jf, const int numv);

  /// Set \a k'th facet vertex of facet \a jf to vertex \a iv.
  /// @pre SetNumFacetVertices(jf) must be called before SetFacetVertex.
  void SetFacetVertex(const FACET_INDEX jf, const int k, const int iv);
  //@}

  /// @name Memory Management Functions
  //@{
  void FreeAll();             ///< Free all memory.
  //@}

  /// @name Check Functions
  //@{
  bool CheckDimension() const;
  bool Check(IJK::ERROR & error_msg) const;
  //@}

  /// @name Generate Polyhedron
  //@{

  /// Generate a square, cube or hypercube.
  void GenCube(const int cube_dimension); 

  /// Generate a triangle, tetrahedron or simplex.      
  void GenSimplex(const int simplex_dimension);

  /// Generate a pyramid over a square, cube or hypercube base.
  void GenPyramid(const int pyramid_dimension);

  //@}

};

typedef ISOSURFACE_TABLE_POLYHEDRON * ISOSURFACE_TABLE_POLYHEDRON_PTR;


//**************************************************
// ISOSURFACE VERTEX
//**************************************************

/// Isosurface vertex class.
class ISOSURFACE_VERTEX {

 public:
  typedef enum {VERTEX, EDGE, FACET, POINT} ISOSURFACE_VERTEX_TYPE;
  typedef float COORD_TYPE;

 protected:
  ISOSURFACE_VERTEX_TYPE vtype;
  int face;
  int num_coord;
  COORD_TYPE * coord;
  std::string label;
  bool is_label_set;

 public:
  ISOSURFACE_VERTEX();        // constructor
  ~ISOSURFACE_VERTEX();       // destructor

  // Get Functions

  /// Return isosurface vertex type.
  ISOSURFACE_VERTEX_TYPE Type() const { return(vtype); };

  /// Return index of face (vertex, edge, facet) containing isosurface.
  /// Valid only if vertex type is VERTEX, EDGE or FACET.
  int Face() const { return(face); };

  /// Return coordinate of isosurface vertex.
  /// Valid only if vertex type is POINT.
  COORD_TYPE Coord(const int d) const { return(coord[d]); };

  /// Return number of isosurface vertex coordinates.
  /// Valid only if vertex type is POINT.
  int NumCoord() const { return(num_coord); };

  /// Return label of isosurface vertex.
  /// Used for extending vertex types.
  std::string Label() const { return(label); };

  /// Return true if label is set.
  bool IsLabelSet() const { return(is_label_set); };

  // Set Functions
  void SetType(const ISOSURFACE_VERTEX_TYPE t) { vtype = t; };
  void SetFace(const int index) { face = index; };
  void SetNumCoord(const int numc);
  void SetCoord(const int ic, const COORD_TYPE c) 
  { coord[ic] = c; };
  void SetLabel(const std::string & s) 
  { label = s; is_label_set = true; };
};


// **************************************************
// ISOSURFACE TABLE
// **************************************************

/// Isosurface lookup table.
/// Stores isosurface patches for each configuration 
///   of +/- labels at polyhedron vertices.
class ISOSURFACE_TABLE {

 protected:

  /// Entry in the isosurface lookup table.
  class ISOSURFACE_TABLE_ENTRY {

  public:
    int num_simplices;
    ISOSURFACE_VERTEX_INDEX * simplex_vertex_list;
    ISOSURFACE_TABLE_ENTRY();                  // constructor
    ~ISOSURFACE_TABLE_ENTRY();                 // destructor

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };


 public:
  /// Configuration encodings.
  /// Standard Marching Cubes lookup table with "-/+" vertex labels 
  ///   uses binary encoding for "-/+".
  /// Interval Volume lookup table uses base 3 encoding for "-/*/+"
  ///   where '*' represents a scalar value between the two isovalues.
  /// NEP (negative-equals-positive) isosurface lookup tables
  ///   with "-/=/+" labels uses base 3 encoding for "-/=/+".
  typedef enum {BINARY, BASE3, NONSTANDARD} ENCODING;
  // Programmer's note: Update standard_encoding_name[] in ijktable.cxx
  //   whenever ENCODING is changed.

  /// Index of entry in isosurface lookup table.
  /// Define within ISOSURFACE_TABLE for use in templates.
  typedef IJKTABLE::TABLE_INDEX TABLE_INDEX;    

 protected:
  ENCODING encoding;           // type of encoding
  std::string encoding_name;   // string storing encoding name

  ISOSURFACE_TABLE_POLYHEDRON polyhedron;  ///< Mesh polyhedron.
  int simplex_dimension;                   ///< Simplex dimension.
  ISOSURFACE_VERTEX * isosurface_vertex;   ///< Array of isosurface vertex descriptors.
  int num_isosurface_vertices;    ///< Number of vertices in array isosurface_vertex[].
  ISOSURFACE_TABLE_ENTRY * entry; ///< Array of isosurface table entries.
  long num_table_entries;         ///< Number of entries in table.

  /// Maximum number of vertices allowed for table polyhedron.
  int max_num_vertices; 

  bool is_table_allocated;  ///< True, if array num_table_entries[] is allocated.

  /// Check if isosurface vertices are allocated.
  /// Throw error if not enough isosurface vertices are allocated.
  void CheckIsoVerticesAlloc
    (const char * procname, const int vstart, const int numv);

  /// Initialization routine.
  void Init(const int dimension, const int simplex_dimension);

 public:
  /// @name Constructors and destructors.
  //@{
  ISOSURFACE_TABLE();
  ISOSURFACE_TABLE(const int d);
  ISOSURFACE_TABLE(const int dimension, const int simplex_dimension);

  ~ISOSURFACE_TABLE();                ///< Destructor

  //@}

  /// @name Get Functions
  //@{

  /// Return table encoding.
  ENCODING Encoding() const { return(encoding); };

  /// Return string for table encoding.
  std::string EncodingName() const { return(encoding_name); };

  /// Return polyhedron dimension.
  int Dimension() const { return(polyhedron.Dimension()); };

  /// Return isosurface simplex dimension.
  int SimplexDimension() const { return(simplex_dimension); };

  /// Return number of vertices in each isosurface simplex.
  int NumVerticesPerSimplex() const { return(SimplexDimension()+1); };

  /// Return number of isosurface vertices in polyhedron.
  int NumIsosurfaceVertices() const { return(num_isosurface_vertices); };

  /// Return number of lookup table entries.
  int NumTableEntries() const { return(num_table_entries); };

  /// Access isosurface table polyhedron.
  const ISOSURFACE_TABLE_POLYHEDRON & Polyhedron() const
    { return(polyhedron); };

  /// Access i'th isosurface vertex.
  const ISOSURFACE_VERTEX & IsosurfaceVertex(const int i) const
    { return(isosurface_vertex[i]); }; 

  /// Return number of simplices in isosurface patch for table entry \a it.
  int NumSimplices(const TABLE_INDEX it) const
    { return(entry[it].num_simplices); }; 

  /// Return \a k'th vertex of isosurface simplex \a is, table entry \a it.
  /// @param it = Index of table entry.
  /// @param is = Simplex \a is of table entry \a it.
  /// @param k = Return \a k'th vertex of simplex \a is.
  ISOSURFACE_VERTEX_INDEX SimplexVertex
    (const TABLE_INDEX it, const int is, const int k) const
    { return(entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k]); };

  /// Return maximum number of polyhedron vertices permitted in any table.
  /// Note: Even tables for polyhedra of this size are probably impossible 
  ///   to compute/store.
  int MaxNumVertices() const { return(max_num_vertices); };

  /// Return true if table memory is allocated.
  bool IsTableAllocated() const
    { return(is_table_allocated); };
  //@}

  /// Return standard string for the encoding.
  static std::string StandardEncodingName(const ENCODING encoding);

  /// @name Set Polyhedron Functions
  //@{
  void SetDimension(const int d) { polyhedron.SetDimension(d); };
  void SetNumPolyVertices(const int numv) 
    { polyhedron.SetNumVertices(numv); };
  void SetNumPolyEdges(const int nume) { polyhedron.SetNumEdges(nume); };
  void SetNumPolyFacets(const int numf) { polyhedron.SetNumFacets(numf); };
  void SetPolySize(const int numv, const int nume, const int numf)
    { SetNumPolyVertices(numv); SetNumPolyEdges(nume); 
      SetNumPolyFacets(numf); };
  void SetPolyVertexCoord(const int iv, const int ic, const int coord)
    { polyhedron.SetVertexCoord(iv, ic, coord); };
  // Note: SetNumPolyVertices or SetPolySize must be called before 
  //   SetPolyVertexCoord
  void SetPolyEdge(const int ie, const int iv0, const int iv1)
    { polyhedron.SetEdge(ie, iv0, iv1); };
  // Note: SetNumPolyEdges or SetPolySize must be called before SetPolyEdge
  void SetPolyNumFacetVertices(const int jf, const int numv)
    { polyhedron.SetNumFacetVertices(jf, numv); }
  void SetPolyFacetVertex(const int jf, const int k, const int iv)
    { polyhedron.SetFacetVertex(jf, k, iv); };
  // Note: SetPolyNumFacetVertices must be called before SetPolyFacetVertex
  void Set(const ISOSURFACE_TABLE_POLYHEDRON & polyhedron)
    { this->polyhedron = polyhedron; };
  //@}

  /// @name Set Isosurface Vertices Functions
  //@{
  void SetNumIsosurfaceVertices(const int num_vertices);
  void SetIsoVertexType(const int i, 
			const ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE t) 
  { isosurface_vertex[i].SetType(t); };
  void SetIsoVertexFace(const int i, const int index) 
  { isosurface_vertex[i].SetFace(index); };
  void SetIsoVertexNumCoord(const int i, const int numc)
  { isosurface_vertex[i].SetNumCoord(numc); };
  void SetIsoVertexCoord(const int i, 
			 const int ic, const ISOSURFACE_VERTEX::COORD_TYPE c) 
  { isosurface_vertex[i].SetCoord(ic, c); };
  void SetIsoVertexLabel(const int i, const std::string & s) 
  { isosurface_vertex[i].SetLabel(s); };

  // store polyhedron vertices, edges or faces as isosurface vertices
  void StorePolyVerticesAsIsoVertices(const int vstart);
  void StorePolyEdgesAsIsoVertices(const int vstart);
  void StorePolyFacetsAsIsoVertices(const int vstart);
  //@}

  /// @name Set Isosurface Table Functions
  //@{
  void SetSimplexDimension(const int d) { this->simplex_dimension = d; };
  void SetEncoding(const ENCODING encoding);
  void SetBinaryEncoding() { SetEncoding(BINARY); };
  void SetBase3Encoding() { SetEncoding(BASE3); };
  void SetNonstandardEncoding(const std::string & name);
  virtual void SetNumTableEntries(const int num_table_entries);
  void SetNumSimplices(const TABLE_INDEX it, const int nums);
  void SetSimplexVertex(const TABLE_INDEX it, const int is, 
			const int iv, const ISOSURFACE_VERTEX_INDEX isov);
  //@}

  /// @name Generate Polyhedron Functions
  //@{
  void GenCube(const int cube_dimension) 
    { polyhedron.GenCube(cube_dimension); };
  // Note: Cubes of dimension > 4 will have too many vertices
  void GenSimplex(const int simplex_dimension) 
    { polyhedron.GenSimplex(simplex_dimension); };
  void GenPyramid(const int pyramid_dimension) 
    { polyhedron.GenPyramid(pyramid_dimension); };
  //@}

  /// @name Check Functions
  //@{
  bool CheckDimension(const int d) const;
  bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
  bool CheckTable(IJK::ERROR & error_msg) const;
  bool Check(IJK::ERROR & error_msg) const;
  //@}

  /// @name Memory Management Functions
  //@{
  virtual void FreeAll();                     /// Free all memory.
  //@}

};

typedef ISOSURFACE_TABLE * ISOSURFACE_TABLE_PTR;

// **************************************************
// ISOSURFACE EDGE TABLE
// **************************************************

/// Isosurface edge table.
/// Store list of edges containing isosurface vertices.
class ISOSURFACE_EDGE_TABLE:public ISOSURFACE_TABLE {

 protected:

  /// Entry in isosurface edge table.
  class ISOSURFACE_EDGE_TABLE_ENTRY {

  public:
    int num_edges;
    EDGE_INDEX * edge_endpoint_list;

    ISOSURFACE_EDGE_TABLE_ENTRY();                  // constructor
    ~ISOSURFACE_EDGE_TABLE_ENTRY();                 // destructor

    bool Check(IJK::ERROR & error_msg) const;
    void FreeAll();                            // free all memory
  };

 protected:
  ISOSURFACE_EDGE_TABLE_ENTRY * edge_entry;

  // initialization routine
  void Init(const int d);

 public:
  ISOSURFACE_EDGE_TABLE(const int d);      // constructor
  ~ISOSURFACE_EDGE_TABLE();                // destructor

  // get functions
  int NumEdges(const TABLE_INDEX it) const 
  { return(edge_entry[it].num_edges); };
  EDGE_INDEX EdgeEndpoint
    (const TABLE_INDEX it, const int ie, const int iend) const
    // it = table entry index. ie = edge index. iend = endpoint index (0 or 1)
  { return(edge_entry[it].edge_endpoint_list[2*ie+iend]); };

  // set isosurface table functions
  virtual void SetNumTableEntries(const int num_table_entries);

  // generate edge lists
  void GenEdgeLists();

  // check functions
  bool CheckTable(IJK::ERROR & error_msg) const;
  bool Check(IJK::ERROR & error_msg) const;

  // free memory
  virtual void FreeAll();                     // free all memory
};

typedef ISOSURFACE_EDGE_TABLE * ISOSURFACE_EDGE_TABLE_PTR;

//**************************************************
// AMBIGUITY INFORMATION FOR ISOSURFACE TABLE
//**************************************************

 /// Isosurface table with ambiguity information.
 class ISOSURFACE_TABLE_AMBIG_INFO {

 protected:
   long num_table_entries;         ///< Number of table entries
   bool * is_ambiguous;            ///< True for ambiguous configurations.
   FACET_INDEX * 
     num_ambiguous_facets;         ///< Number of ambiguous facts.
   FACET_SET * ambiguous_facet;    ///< k'th bit is 1 if facet k is ambiguous

   void Init();                    ///< Initialization routine.
   void Alloc(const long num_table_entries);  ///< Allocate memory.
   void FreeAll();                 ///< Free all memory.

 public:
    
   // constructors
   ISOSURFACE_TABLE_AMBIG_INFO() { Init(); };
   ~ISOSURFACE_TABLE_AMBIG_INFO();                // destructor

   // get functions
   bool IsAmbiguous(const TABLE_INDEX it) const
     { return(is_ambiguous[it]); };
   FACET_INDEX NumAmbiguousFacets(const TABLE_INDEX it) const
     { return(num_ambiguous_facets[it]); };
   // get functions
   bool IsFacetAmbiguous(const TABLE_INDEX it, const FACET_INDEX jf) const
     { return(ambiguous_facet[it] & ((1L) << jf)); };
   long NumTableEntries() const { return(num_table_entries); };

   // compute functions
   void ComputeAmbiguityInformation
     (const ISOSURFACE_TABLE & isotable);

 };

//**************************************************
// ISOSURFACE TABLE WITH AMBIGUITY INFORMATION
//**************************************************

 /// Isosurface table with ambiguity information.
 class ISOSURFACE_TABLE_AMBIG:public ISOSURFACE_TABLE {

 protected:
   bool * is_ambiguous;            ///< True for ambiguous configurations.
   FACET_INDEX * 
     num_ambiguous_facets;         ///< Number of ambiguous facts.
   FACET_SET * ambiguous_facet;    ///< k'th bit is 1 if facet k is ambiguous

   void Init();                    ///< Initialization routine.
   void FreeAll();                 ///< Free all memory.

   /// Compute ambiguity information for all table entries.
   bool ComputeAmbiguous(const int * vertex_sign) const;

   /// Compute number of ambiguous facets.
   FACET_INDEX ComputeNumAmbiguousFacets(const int * vertex_sign) const;

   /// Compute set of ambiguous facets.
   void ComputeAmbiguousFacets
     (const int * vertex_sign, 
      FACET_SET & facet_set, FACET_INDEX & num_ambiguous_facets) const;

 public:
    
   // constructors
   ISOSURFACE_TABLE_AMBIG():ISOSURFACE_TABLE() { Init(); };
   ISOSURFACE_TABLE_AMBIG(const int d):ISOSURFACE_TABLE(d)
     { Init(); };
   ISOSURFACE_TABLE_AMBIG(const int dimension, const int simplex_dimension):
     ISOSURFACE_TABLE(dimension, simplex_dimension)
     { Init(); };

   ~ISOSURFACE_TABLE_AMBIG();                // destructor

   // get functions
   bool IsAmbiguous(const TABLE_INDEX it) const
     { return(is_ambiguous[it]); };
   FACET_INDEX NumAmbiguousFacets(const TABLE_INDEX it) const
     { return(num_ambiguous_facets[it]); };
   // get functions
   bool IsFacetAmbiguous(const TABLE_INDEX it, const FACET_INDEX jf) const
     { return(ambiguous_facet[it] & ((1L) << jf)); };

   // set functions
   virtual void SetNumTableEntries(const int num_table_entries);

   // compute functions
   void ComputeAmbiguityInformation() const;

 };

//**************************************************
// UTILITY FUNCTIONS
//**************************************************

// calculate number of entries required in ISOSURFACE_TABLE
unsigned long calculate_num_entries(const int num_vert, const int num_colors);

// convert integer to base "base"
void convert2base(const unsigned long ival, const unsigned int base, 
		  int * digit, const unsigned int max_num_digits);


//**************************************************
// ROUTINES FOR GENERATING POLYHEDRA
//**************************************************

void generate_prism(const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
		    ISOSURFACE_TABLE_POLYHEDRON & prism);

//**************************************************
// AMBIGUITY ROUTINES
//**************************************************

/// Return true if isosurface topology is ambiguous
/// @param vertex_sign[i] = Sign of isosurface vertex i. (0 or 1);
bool is_poly_ambiguous
(const ISOSURFACE_TABLE_POLYHEDRON & poly, const int * vertex_sign);

/// Return true if facet jf is ambiguous.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
bool is_facet_ambiguous
(const ISOSURFACE_TABLE_POLYHEDRON & poly,
 const FACET_INDEX jf, const int * vertex_sign);

/// Return number of vertices connected by edges to iv with same sign as iv.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
int compute_num_connected
(const ISOSURFACE_TABLE_POLYHEDRON & poly,
 const int iv, const int * vertex_sign);

/// Return number of vertices connected by edges in facet jf to vertex iv with same sign as iv.
/// @param jf = Facet index.
/// @param iv = Vertex index.  Precondition: Vertex iv is in facet jf.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
int compute_num_connected_in_facet
(const ISOSURFACE_TABLE_POLYHEDRON & poly,
 const FACET_INDEX jf, const int iv, const int * vertex_sign);

/// Compute ambiguous facets.
/// @param vertex_sign[i] = Sign of isosurface vertex i.
void compute_ambiguous_facets
  (const ISOSURFACE_TABLE_POLYHEDRON & poly, const int * vertex_sign,
   FACET_SET & facet_set, FACET_INDEX & num_ambiguous_facets);

};

#endif
