// ijk output templates
// output OpenInventor .iv file
// output Geomview file
// output Fig file

#ifndef _IJKIO_
#define _IJKIO_

#include <iostream>
#include <vector>
#include <string>

#include "ijk.txx"

namespace IJK {

  //******************************************
  // Output OpenInventor .iv file
  //******************************************

  template <class CTYPE, class VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // output OpenInventor .iv file
    // out = output stream
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    IJK::ERROR error("ijkoutIV");

    if (dim != 3)
      throw error("Illegal dimension.  OpenInventor files are only for dimension 3.");

    out << "#Inventor V2.1 ascii" << std::endl;
    out << std::endl;

    out << "Separator {" << std::endl;

    // set vertex ordering to clockwise to turn on two-sided lighting
    out << "  ShapeHints {" << std::endl;
    out << "    vertexOrdering CLOCKWISE" << std::endl;
    out << "  }" << std::endl;
    out << std::endl;

    out << "  IndexedFaceSet {" << std::endl;

    out << "    vertexProperty VertexProperty {" << std::endl;
    out << "      vertex [" << std::endl;

    // output vertex coordinates
    out << std::endl << "# vertex coordinates" << std::endl;
    for (int i = 0; i < numv; i++) {
      for (int d = 0; d < dim; d++) {
	out << coord[dim*i+d];
	if (d < dim-1) { out << " "; }
	else {	
	  if (i < numv-1) { out << "," << std::endl; };
	};
      };
    };

    out << " ]" << std::endl;
    out << "    }" << std::endl;

    out << "    coordIndex [" << std::endl;
    // output triangle vertices
    out << std::endl << "# triangle vertices" << std::endl;
    for (int it = 0; it < numt; it++) {
      for (int d = 0; d < dim; d++) {
	out << tri[dim*it+d] << ",";
      };
      out << "-1";
      if (it < numt-1) {
	out << "," << std::endl;
      };
    };
    out << " ]" << std::endl;

    out << "  }" << std::endl;

    out << "}" << std::endl;
  }

  template <class CTYPE, class VTYPE> void ijkoutIV
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * tri, const int numt)
    // output OpenInventor .iv format to standard outpu
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    ijkoutIV(std::cout, dim, coord, numv, tri, numt); 
  }

  template <class CTYPE, class VTYPE> void ijkoutIV
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(out, dim, &(coord[0]), coord.size()/dim,
	     &(simplex_vert[0]), simplex_vert.size()/dim);
  }

  template <class CTYPE, class VTYPE> void ijkoutIV
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutIV(dim, &(coord[0]), coord.size()/dim,
	     &(simplex_vert[0]), simplex_vert.size()/dim);
  }	

  //******************************************
  // Geomview OFF file
  //******************************************

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF files
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    IJK::ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
	out << coord[iv*dim + d];
	if (d < dim-1) { out << " "; }
	else { out << std::endl; };
      }
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
	out << simplex_vert[is*numv_per_simplex + iv];
	if (iv < numv_per_simplex-1) { out << " "; }
	else { out << std::endl; };
      };

    };

  }

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex,
   const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  // output Geomview OFF format to standard output
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, numv, 
	      simplex_vert, nums);
  }

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  // set numv_per_simplex to dim
  {
    ijkoutOFF(out, dim, dim, coord, numv, simplex_vert, nums);
  }

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const CTYPE * coord, const int numv,
   const VTYPE * simplex_vert, const int nums)
  // output Geomview OFF format to standard output
  // set numv_per_simplex to dim
  {
    ijkoutOFF(dim, dim, coord, numv, simplex_vert, nums);
  }


  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  // output Geomview OFF files
  // vector input format
  {
    ijkoutOFF(out, dim, numv_per_simplex, &(coord[0]), coord.size()/dim,
	      &(simplex_vert[0]), simplex_vert.size()/numv_per_simplex);
  }

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (std::ostream & out, const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  // output Geomview OFF files
  // vector input format
  // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(out, dim, numv_per_simplex, coord, simplex_vert);
  }


  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const int numv_per_simplex, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  // output Geomview OFF format to standard output
  // vector input format
  {
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  template <class CTYPE, class VTYPE> void ijkoutOFF
  (const int dim, const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert)
  // output Geomview OFF format to standard output
  // vector input format
  // set numv_per_simplex to dim
  {
    int numv_per_simplex = dim;
    ijkoutOFF(std::cout, dim, numv_per_simplex, coord, simplex_vert);
  }

  template <class T, class colorT> void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
    // output Geomview OFF files
    // out = output stream
    // dim = dimension
    // numv_per_simplex = num vertices per simplex
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
    // front_color = array of front colors, 4 entries (RGBA) per vertex
    // back_color = array of backface colors, 4 entries (RGBA) per vertex
    //              May be NULL.
  {
    IJK::ERROR error("ijkoutColorVertOFF");

	//std::cout<<"ijkcolorvert!!\n";	

    assert(front_color != NULL);

    out << "C";
    if (back_color != NULL && false) out << "C";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
	out << coord[iv*dim + d];
	if (d+1 < dim) out << " ";
      }
      out << "  ";
  //    std::cout<<"front color:";	
      for (int ic = 0; ic < 4; ic++) {
	out << front_color[4*iv+ic];
//	std::cout<<front_color[4*iv+ic]<<" ";
	if (ic < 3) out << " ";
      }
  //     std::cout<<"\n";	
  //     std::cout<<"back color:";	
      out << "  ";
      if (back_color != NULL && false) {
	for (int ic = 0; ic < 4; ic++) {
	  out << back_color[4*iv+ic];
	  std::cout<<back_color[4*iv+ic]<<" ";
	  if (ic < 3) out << " ";
	}
//	 std::cout<<"\n";	
      };
      out << std::endl;
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
	out << simplex_vert[is*numv_per_simplex + iv];
	if (iv < numv_per_simplex-1) { out << " "; }
	else { out << std::endl; };
      };
    };
  }

  template <class T, class colorT> void ijkoutColorVertOFF
  (const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
  {	
    ijkoutColorVertOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  template <class CTYPE, class VTYPE, class colorT> void ijkoutColorVertOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const colorT * front_color, const colorT * back_color)
  {
	 
    ijkoutColorVertOFF
      (out, dim, numv_per_simplex, &(coord[0]), coord.size()/dim,
       &(simplex_vert[0]), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  template <class T, class colorT> void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
    // output Geomview OFF files
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
    // front_color = array of front colors, 4 entries (RGBA) per face
    // back_color = array of backface colors, 4 entries (RGBA) per face
    //              May be NULL.
  {
    IJK::ERROR error("ijkoutColorFacesOFF");

    assert(front_color != NULL);

    out << "D";
    if (back_color != NULL) out << "D";
    if (dim == 3) { out << "OFF" << std::endl; }
    else if (dim == 4) { out << "4OFF" << std::endl;}
    else {
      out << "nOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
	out << coord[iv*dim + d];
	if (d+1 < dim) out << " ";
      }
      out << std::endl;
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
	out << simplex_vert[is*numv_per_simplex + iv];
	if (iv < numv_per_simplex-1) { out << " "; }
      };
      out << "  ";
      for (int ic = 0; ic < 4; ic++) {
	out << front_color[4*is+ic];
	if (ic < 3) out << " ";
      }
      out << "  ";
      if (back_color != NULL) {
	for (int ic = 0; ic < 4; ic++) {
	  out << back_color[4*is+ic];
	  if (ic < 3) out << " ";
	}
      };
      out << std::endl;
    };

  }

  template <class T, class colorT> void ijkoutColorFacesOFF
  (const int dim, const int numv_per_simplex, 
   const T * coord, const int numv,
   const int * simplex_vert, const int nums,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorFacesOFF
      (std::cout, dim, numv_per_simplex, coord, numv, simplex_vert, nums,
       front_color, back_color);
  }

  template <class CTYPE, class VTYPE, class colorT> void ijkoutColorFacesOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,   
   const std::vector<CTYPE> & coord,
   const std::vector<VTYPE> & simplex_vert,
   const colorT * front_color, const colorT * back_color)
  {
    ijkoutColorFacesOFF
      (out, dim, numv_per_simplex, &(coord[0]), coord.size()/dim,
       &(simplex_vert[0]), simplex_vert.size()/numv_per_simplex,
       front_color, back_color);
  }

  template <class CTYPE, class NTYPE, class VTYPE> void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const CTYPE * coord, const NTYPE * normal, const int numv,
   const VTYPE * simplex_vert, const int nums)
    // output Geomview OFF files with vertex normal information
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    IJK::ERROR error("ijkoutOFF");

    if (dim == 3) { out << "NOFF" << std::endl; }
    else if (dim == 4) { out << "N4OFF" << std::endl;}
    else {
      out << "NnOFF" << std::endl;
      out << dim << std::endl;
    };

    out << numv << " " << nums << " " << 0 << std::endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
	out << coord[iv*dim + d] << " ";
      }
      for (int d = 0; d < dim; d++) {
	out << normal[iv*dim+d];
	if (d < dim-1) { out << " "; }
	else { out << std::endl; };
      }
    };
    out << std::endl;

    for (int is = 0; is < nums; is++) {
      out << numv_per_simplex << " ";
      for (int iv = 0; iv < numv_per_simplex; iv++) {
	out << simplex_vert[is*numv_per_simplex + iv];
	if (iv < numv_per_simplex-1) { out << " "; }
	else { out << std::endl; };
      };

    };

  }

  template <class CTYPE, class NTYPE, class VTYPE> void ijkoutNormalsOFF
  (std::ostream & out, const int dim, const int numv_per_simplex,
   const std::vector<CTYPE> & coord,
   const std::vector<NTYPE> & normal,
   const std::vector<VTYPE> & simplex_vert)
  {
    ijkoutNormalsOFF
      (out, dim, numv_per_simplex,
       &(coord[0]), &(normal[0]), coord.size()/dim,
       &(simplex_vert[0]),
       simplex_vert.size()/numv_per_simplex);
  }

  template <class T> void ijkinOFF
  (std::istream & in, int & dim, int & mesh_dim,
	T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF files
    // in = input stream
    // dim = dimension
    // mesh_dim = mesh dimension
    //   Precondition: All simplices have the same dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    std::string header_line;
    int nume;

    IJK::ERROR error("ijkinOFF");

    if (!in.good()) {
      error.AddMessage("Error reading from input stream in.");
      throw error;
    }

    coord = NULL;
    simplex_vert = NULL;
    mesh_dim = 0;            // default mesh dimension

    // read input
    header_line = "";
    while (header_line == "" && !in.eof())
      in >> header_line;

    if (header_line == "OFF") {
      dim = 3;
    }
    else if (header_line == "4OFF") {
      dim = 4;
    }
    else if (header_line == "nOFF") {
      in >> dim;
    }
    else {
      std::string errmsg =  
	"Illegal Geomview .off file header: " + header_line;
      throw error(errmsg);
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nums;
    in >> nume;

    if (!in.good()) {
      error.AddMessage("Error reading number of vertices and polyhedra from input stream in.");
      throw error;
    }

    coord = new T[numv*dim];
    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
	{ in >> coord[iv*dim + d]; };

      if (!in.good()) {
	error.AddMessage("Error reading coordinates of vertex ", iv,
			 " from input stream in.");
	throw error;
      }
    }
      
    if (!in.good()) {
      error.AddMessage("Error reading vertex coordinates from input stream in.");
      throw error;
    }

    if (nums > 0) {
      // use first simplex to set mesh dimension
      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert > 0) { 
	mesh_dim = num_simplex_vert-1; 
      }
      else { mesh_dim = 0; };

      simplex_vert = new int[nums*num_simplex_vert];

      // read in first simplex
      for (int d = 0; d < num_simplex_vert; d++)
	in >> simplex_vert[d];

      // read in remaining simplices
      int nvert;
      for (int is = 1; is < nums; is++) {

	in >> nvert;
	if (nvert != num_simplex_vert) {
	  delete [] coord;
	  delete [] simplex_vert;
	  coord = NULL;
	  simplex_vert = NULL;
	  throw error("Input simplices do not all have the same dimension.");
	}
	for (int d = 0; d < num_simplex_vert; d++)
	  in >> simplex_vert[is*num_simplex_vert + d];
      };
    };
  }

  template <class T> void ijkinOFF
  (int & dim, int & mesh_dim,
   T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF format to standard input
    // in = input stream
    // dim = dimension
    // mesh_dim = mesh dimension
    //   Precondition: All simplices have the same dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkinOFF(std::cin, dim, mesh_dim, coord, numv, simplex_vert, nums);
  }


  template <class T> void ijkinOFF
  (std::istream & in, int & dim, T * & coord, int & numv,
   int * & simplex_vert, int & nums)
    // input Geomview OFF files
    // in = input stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    std::string header_line;
    int nume;

    IJK::ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;

    if (!in.good()) {
      throw error("Error: Corrupted input stream. Unable to read input.");
    }

    // read input
    header_line = "";
    while (header_line == "" && !in.eof())
      in >> header_line;

    if (header_line == "OFF") {
      dim = 3;
    }
    else if (header_line == "4OFF") {
      dim = 4;
    }
    else if (header_line == "nOFF") {
      in >> dim;
    }
    else {
      std::string errmsg =  
	"Illegal Geomview .off file header: " + header_line;
      throw error(errmsg);
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nums;
    in >> nume;

    coord = new T[numv*dim];
    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
	in >> coord[iv*dim + d];
    };

    simplex_vert = new int[nums*dim];
    for (int is = 0; is < nums; is++) {

      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert != dim) {
	delete [] coord;
	delete [] simplex_vert;
	coord = NULL;
	simplex_vert = NULL;
	throw error("Wrong number of vertices in geomview polytope list.");
      }
      for (int d = 0; d < dim; d++)
	in >> simplex_vert[is*dim + d];
    };
  }

  template <class T> void ijkinOFF
  (int & dim, T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF format to standard input
    // in = input stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkinOFF(std::cin, dim, coord, numv, simplex_vert, nums);
  }

  //******************************************
  // Fig file
  //******************************************

  template <class CTYPE, class VTYPE, class SCALE_TYPE> void ijkoutFIG
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord,
   const bool flag_polyline = false)
    // output .fig file
    // out = output stream
    // dim = dimension.  Must be 2.
    // coord[2*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // seg[2*j+k] = k'th vertex index of segment j (k < dim)
    // nums = number of line segments
    // flag_polyline: if true, join segments to form FIG polylines
    //                if false, output segments individually
  {
    IJK::ERROR error("ijkoutOFF");
    const int cap_style = 1;             // round cap

    if (dim != 2)
      throw error("Illegal dimension.  Fig files are only for dimension 2.");

    // FIG header
    out << "#FIG 3.2" << std::endl;

    out << "Landscape" << std::endl;
    out << "Center" << std::endl;
    out << "Metric" << std::endl;
    out << "Letter" << std::endl;
    out << "100.00" << std::endl;
    out << "Single" << std::endl;
    out << "-2" << std::endl;
    out << "1200 2" << std::endl;

    if (!flag_polyline) {
      // output each line segment separately line segments

      for (int i = 0; i < nums; i++) {
	out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
	    << cap_style << " -1 0 0 2" << std::endl;
	out << "    ";
	for (int j = 0; j < 2; j++) {
	  int iv = seg[2*i+j];
	  for (int k = 0; k < 2; k++) {
	    int x = int(scale_coord*coord[2*iv+k]);
	    out << "  " << x;
	  }
	}
	out << std::endl;
      }
    }
    else {
      ijkoutFIGpolyline(out, dim, coord, numv, seg, nums, scale_coord);
    }
  }

  // local namespace
  namespace {
    inline int index_other_endpoint(int k)
    {
      int kseg = k/2;
      int j = (k+1)%2;
      return(2*kseg+j);
    }
  };

  template <class CTYPE, class VTYPE, class SCALE_TYPE> void ijkoutFIGpolyline
  (std::ostream & out, const int dim, const CTYPE * coord, const int numv,
   const VTYPE * seg, const int nums, const SCALE_TYPE scale_coord)
    // output FIG polygonal lines
  {
    using std::vector;
    const int cap_style = 1;             // round cap
    IJK::ARRAY<bool> is_edge_processed(nums, false);

    // degree of vertex iv
    IJK::ARRAY<int> degree(numv, 0);

    // next[2*i] = next element (in next[]) around vertex 0 of segment i
    // next[2*i+1] = next element (in next[]) around vertex 1 of segment i
    IJK::ARRAY<int> next(2*nums);

    // first[iv] = first element in next[] incident on vertex iv
    // last[iv] = last element in next[] incident on vertex iv
    IJK::ARRAY<int> first(numv);
    IJK::ARRAY<int> last(numv);

    // set degree[], first[], last[], next[]
    for (int i = 0; i < nums; i++)
      for (int j = 0; j < 2; j++) {
	VTYPE iv = seg[2*i+j];
	if (degree[iv] == 0) {
	  first[iv] = 2*i+j;
	  last[iv] = 2*i+j;
	  next[2*i+j] = 2*i+j;
	}
	else {
	  next[2*i+j] = next[last[iv]];
	  next[last[iv]] = 2*i+j;
	  last[iv] = 2*i+j;
	}
	degree[iv]++;
      }

    vector<VTYPE> vlist;
    for (int i = 0; i < nums; i++) {
      if (is_edge_processed[i]) { continue; };

      vlist.clear();
      VTYPE iv0 = seg[2*i];
      VTYPE iv1 = seg[2*i+1];
      vlist.push_back(iv0);
      vlist.push_back(iv1);
      is_edge_processed[i] = true;

      iv0 = iv1;
      int k = next[2*i+1];
      k = index_other_endpoint(k);
      while (!is_edge_processed[k/2]) {
	iv1 = seg[k];
	vlist.push_back(iv1);
	is_edge_processed[k/2] = true;
	iv0 = iv1;
	k = next[k];
	k = index_other_endpoint(k);
      }

      out << "2 1 0 1 0 7 50 -1 -1 0.000 0 "
	  << cap_style << " -1 0 0 " << vlist.size() << std::endl;

      out << "    ";
      for (int j = 0; j < vlist.size(); j++) {
	int iv = vlist[j];
	for (int k = 0; k < 2; k++) {
	  int x = int(scale_coord*coord[2*iv+k]);
	  out << "  " << x;
	}
      }
      out << std::endl;
    }

  }

}

#endif
