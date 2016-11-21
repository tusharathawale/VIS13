/// \file ijkmcubeIO.cxx
/// IO routines for ijkmcube

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
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ijkmcubeIO.h"
#include "ijkmcube_util.h"
#include "ijkIO.txx"
#include "ijkxitIO.h"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKTABLE;

using namespace std;

// **************************************************
// PARSE COMMAND LINE
// **************************************************

// local namespace
namespace {

typedef enum
  {REGION_PARAM, OCTREE_PARAM, LIST_PARAM, NEP_PARAM, NEP_DUP_PARAM,
   IVOL_PARAM, SNAP_PARAM, HIGHRES_PARAM, 
   TOPOLOGY_PARAM, INTERPOLATE_PARAM, NUMRES_LEVELS_PARAM,
   SUBSAMPLE_PARAM, SUPERSAMPLE_PARAM, EDGE1_PARAM, EDGE2_PARAM,
   LOG_INTERVAL_PARAM,
   COLOR_ALTERNATING_PARAM, 
   HELP_PARAM, OFF_PARAM, IV_PARAM, 
   DIR_PARAM, OUTPUT_FILENAME_PARAM, STDOUT_PARAM, 
   NOWRITE_PARAM, SILENT_PARAM, TIME_PARAM, UNKNOWN_PARAM} PARAMETER;
const char * parameter_string[] = 
  {"-region", "-octree", "-list", 
   "-nep", "-nep_dup", "-ivol", "-snap", 
   "-highres", "-topology", "-interpolate", "-numres_levels", 
   "-subsample", "-supersample",
   "-edge1", "-edge2", "-log_interval",
   "-color_alternating",
   "-help", "-off", "-iv", "-dir", 
   "-o", "-stdout",
   "-nowrite", "-s", "-time", "-unknown"};

PARAMETER get_parameter_token(char * s)
// convert string s into parameter token
{
  for (int i = 0; i < int(UNKNOWN_PARAM); i++)
    if (strcmp(parameter_string[i], s) == 0)
      return(PARAMETER(i));
  return(UNKNOWN_PARAM);
}

void get_box(const char * s, IJK::BOX<GRID_COORD_TYPE> & box)
{
  istringstream coord_string;
  vector<GRID_COORD_TYPE> coord;

  string s2 = s;
  // remove trailing blanks from s2
  size_t pos = 0;
  for (size_t i = 0; i < s2.length(); i++) {
    if (!isspace(s2[i])) { pos = i+1; }
  }
  if (pos < s2.length()) { s2.erase(pos); };

  coord_string.str(s2);
  while (coord_string.good()) {
    GRID_COORD_TYPE c;
    coord_string >> c;
    coord.push_back(c);
  }

  if (coord_string.fail() && !coord_string.eof()) {
    cerr << "Error reading -highres coordinates: "
	 << "\"" << s << "\"" << endl;
    cerr << "  Non-numeric character in coordinate string." << endl;
    exit(600);
  }

  if (coord.size() == 1) {
    cerr << "Error in -highres coordinates: "
	 << "\"" << s << "\"" << endl;
    cerr << "  Coordinate list must be contained in quotation marks." << endl;
    cerr << "  Number of coordinates must be at least two." << endl;
    exit(605);
  }

  if (coord.size()%2 == 1) {
    cerr << "Error in -highres coordinates: "
	 << "\"" << s << "\"" << endl;
    cerr << "  Number of coordinates must be even." << endl;
    exit(610);
  }

  int box_dim = coord.size()/2;
  box.SetDimension(box_dim);
  for (int d = 0; d < box_dim; d++) {
    GRID_COORD_TYPE minc = coord[d];
    GRID_COORD_TYPE maxc = coord[d+box_dim];
    if (minc > maxc) {
      cerr << "Error in -highres coordinates: "
	   << "\"" << s << "\"" << endl;
      cerr << "  Minimum coordinate " << minc
	   << " (position " << d << ")"
	   << " is greater than maximum coordinate " << maxc
	   << " (position " << d+box_dim << ")." << endl;
      cerr << "  List all minimum coordinates before maximum coordinates."
	   << endl;
      exit(620);
    }

    box.SetMinCoord(d, coord[d]);
    box.SetMaxCoord(d, coord[d+box_dim]);
  }

}

ISOSURFACE_TOPOLOGY get_topology(char * s)
// convert string s into parameter token
{
  ISOSURFACE_TOPOLOGY topology = ISOTABLE_TOPOLOGY;

  string str = s;

  if (str == "isotable") 
    { topology = ISOTABLE_TOPOLOGY; }
  else if (str == "adecider") 
    { topology = ASYMPTOTIC_DECIDER_TOPOLOGY; }
  else if (str == "cube_decider") 
    { topology = CUBE_DECIDER_TOPOLOGY; }
  else if (str == "linear") 
    { topology = LINEAR_TOPOLOGY; }
  else {
    cerr << "Error in input parameter -topology.  Illegal isosurface topology: " 
	 << s << "." << endl;
    exit(1020);
  }

  return(topology);
}

INTERPOLATION_TYPE get_interpolation_type(char * s)
// convert string s into parameter token
{
  INTERPOLATION_TYPE type = LINEAR_INTERPOLATION;

  if (strcmp(s, "linear") == 0) 
    { type = LINEAR_INTERPOLATION; }
  else if (strcmp(s, "multilinear") == 0) 
    { type = MULTILINEAR_INTERPOLATION; }
  else if (strcmp(s, "alpha_uncertainty") == 0) 
    { type = ALPHA_UNCERTAINTY; }
  else {
    cerr << "Error in input parameter -interpolate.  Illegal interpolation type: " 
	 << s << "." << endl;
    exit(1030);
  }

  return(type);
}

bool check_highres
(const int dimension, const IO_INFO & io_info, IJK::ERROR & error)
  // check that highres options have correct number of parameters
  // dimension = grid dimension
{
  using std::string;

  for (int i = 0; i < io_info.high_resolution_regions.size(); i++) {
    if (io_info.high_resolution_regions[i].Dimension() != dimension) {
      string msg = string("Error in option: ") + io_info.high_resolution_option[i];
      error.AddMessage(msg);
      error.AddMessage
	("Option -highres should be followed by a string containing ",
	 2*dimension, " coordinates.");
      return(false);
    }
  }
  return(true);
}

};

void IJKMCUBE::parse_command_line(int argc, char **argv, IO_INFO & io_info)
// parse command line
// control parameters, followed by one or more isovalues, 
// followed by input file name
{
  if (argc == 1) { usage_error(); };

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {
    PARAMETER param = get_parameter_token(argv[iarg]);
    if (param == UNKNOWN_PARAM) break;

    switch(param) {
    case REGION_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.region_length);
      io_info.use_minmax = true;
      break;

    case OCTREE_PARAM:
      io_info.use_octree = true;
      break;

    case LIST_PARAM:
      io_info.use_list = true;
      break;

    case NEP_PARAM:
      // use isotable which differentiates between scalar values less
      //   than the isovalue (negative), equals, and greater than
      //   the isovalue (positive)
      io_info.use_nep = true;  
      break;

    case NEP_DUP_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.nep_num_dup);
      break;	

    case IVOL_PARAM:
      io_info.interval_volume_flag = true;
      break;

    case SNAP_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &io_info.snap_value);
      io_info.snap_flag = true;
      break;

    case TOPOLOGY_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.isosurface_topology = get_topology(argv[iarg]);
      break;

    case INTERPOLATE_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.interpolation_type = get_interpolation_type(argv[iarg]);
      break;

    case SUBSAMPLE_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.subsample_resolution);
      io_info.flag_subsample = true;
      break;

    case SUPERSAMPLE_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.supersample_resolution);
      io_info.flag_supersample = true;
      break;

    case HIGHRES_PARAM:
      {
	IJK::BOX<GRID_COORD_TYPE> box;
	string s = argv[iarg];
	iarg++;
	if (iarg >= argc) usage_error();
	get_box(argv[iarg], box);
	io_info.high_resolution_regions.push_back(box);
	s = s + " " + argv[iarg];
	io_info.high_resolution_option.push_back(s);
	io_info.use_multires = true;
	break;
      }

    case NUMRES_LEVELS_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.num_resolution_levels);
      break;

    case EDGE1_PARAM:
      io_info.edge_representation = EDGE_ID;
      break;

    case EDGE2_PARAM:
      io_info.edge_representation = EDGE_ENDPOINT_PAIR;
      break;

    case LOG_INTERVAL_PARAM:
      int log2_interval;
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &log2_interval);
      io_info.merge_edges_parameters.SetLog2Interval(log2_interval);
      break;

    case COLOR_ALTERNATING_PARAM:
      io_info.flag_color_alternating = true;
      break;

    case DIR_PARAM: 
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.isotable_directory = argv[iarg];
      break;

    case OFF_PARAM:
      io_info.output_format = OFF;
      break;

    case IV_PARAM:
      io_info.output_format = IV;
      break;

    case OUTPUT_FILENAME_PARAM:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.output_filename = argv[iarg];
      break;

    case STDOUT_PARAM:
      io_info.use_stdout = true;
      break;

    case NOWRITE_PARAM:
      io_info.nowrite_flag = true;
      break;

    case SILENT_PARAM:
      io_info.flag_silent = true;
      break;

    case TIME_PARAM:
      io_info.report_time_flag = true;
      break;

    case HELP_PARAM:
      help();
      break;
    };

    iarg++;
  };

  // remaining parameters should be list of isovalues followed
  // by input file name

  // check for more parameter tokens
  for (int j = iarg; j < argc; j++) {
    if (get_parameter_token(argv[j]) != UNKNOWN_PARAM) {
      // argv[iarg] is not an isovalue
      cerr << "Error. Illegal parameter: " << argv[iarg] << endl;
      usage_error();
    }
  }

  if (iarg+2 > argc) {
    cerr << "Error.  Missing input isovalue or input file name." << endl;
    usage_error();
  };

  // store isovalues
  for (int j = iarg; j+1 < argc; j++) {
    io_info.isovalue_string.push_back(argv[j]);
    SCALAR_TYPE value;

    istringstream input_string(argv[j]);
    input_string >> value;

    if (input_string.fail()) {
      cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue." 
	   << endl;
      usage_error();
    };

    io_info.isovalue.push_back(value);
  }

  io_info.input_filename = argv[argc-1];

  if (io_info.use_octree && io_info.use_minmax) {
    cerr << "Error.  Can't use both -region and -octree parameters.";
    usage_error();
  }

  if (io_info.interval_volume_flag && io_info.isovalue.size() < 2) {
    cerr << "Error.  Need at least two isovalues for interval volume generation." << endl;
    exit(225);
  }

  if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
    cerr << "Error.  Subsample resolution must be an integer greater than 1."
	 << endl;
    exit(230);
  };

  if (io_info.output_filename != NULL && io_info.use_stdout) {
    cerr << "Error.  Can't use both -o and -stdout parameters."
	 << endl;
    exit(230);
  };

  if (io_info.snap_flag && 
      (io_info.snap_value < 0.0 || io_info.snap_value > 0.5)) {
    cerr << "Error.  Illegal snap value " << io_info.snap_value << "."
	 << endl;
    cerr << "        Snap value must be in range [0.0, 0.5]." << endl;
    exit(550);
  };

  if (io_info.flag_subsample && io_info.use_multires) {
    cerr << "Error.  Can't use both -subsample and -highres parameters."
	 << endl;
    exit(555);
  }

  if (io_info.flag_subsample && io_info.flag_supersample) {
    cerr << "Error.  Can't use both -subsample and -supersample parameters."
	 << endl;
    exit(555);
  }
}

// Check input information/flags.
bool IJKMCUBE::check_input
(const IO_INFO & io_info, 
 const MC_SCALAR_GRID_BASE & scalar_grid,
 IJK::ERROR & error)
{
  if (io_info.isotable_directory == "") {
    error.AddMessage("Error.  Unknown isotable directory.");
    error.AddMessage
      ("  Use -dir {isotable_directory} argument or set environment variable IJK_ISOTABLE_DIR.");
    return(false);
  }

  if (io_info.interval_volume_flag) {
    // Construct interval volume
    if (io_info.isovalue.size() > 2 && io_info.use_stdout) {
      error.AddMessage
	("Error.  Cannot use stdout for more than one interval volume.");
      return(false);
    }

    if (io_info.isovalue.size() > 2 && io_info.output_filename != NULL) {
	error.AddMessage
	  ("Error.  Cannot specify output file for more than one interval volume.");
	return(false);
    }
  }
  else {
    // Construct isosurface
    if (io_info.isovalue.size() > 1 && io_info.use_stdout) {
      error.AddMessage
	("Error.  Cannot use stdout for more than one isovalue.");
      return(false);
    }

    if (io_info.isovalue.size() > 1 && io_info.output_filename != NULL) {
      error.AddMessage
	("Error.  Cannot specify output file for more than one isovalue.");
      return(false);
    }
  }

  if (!check_highres(scalar_grid.Dimension(), io_info, error)) 
    { return(false); }

  return(true);
}

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

void IJKMCUBE::read_nrrd_file
(const char * input_filename, MC_SCALAR_GRID & scalar_grid, NRRD_INFO & nrrd_info,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  int dimension = 0;
  Nrrd *nin;

  // get scalar field data from nrrd file
  nin = nrrdNew();
  if (nrrdLoad(nin, input_filename, NULL)) {
    char *err = biffGetDone(NRRD);
    cerr << "Error reading: " << input_filename << endl;
    cerr << "  Error: " << err << endl;
    exit(35);
  };
  dimension = nin->dim;

  if (dimension < 1) {
    cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
    exit(20);
  };

  IJK::ARRAY<AXIS_SIZE_TYPE> axis_size(dimension);
  size_t size[NRRD_DIM_MAX];
  double grid_spacing[NRRD_DIM_MAX];

  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size); 
  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSpacing, grid_spacing); 

  nrrd_info.grid_spacing.clear();
  for (int d = 0; d < dimension; d++) { 
    axis_size[d] = size[d]; 
    nrrd_info.grid_spacing.push_back(grid_spacing[d]);
  }

  scalar_grid.SetSize(dimension, axis_size.PtrConst());
  nrrd2scalar(nin, scalar_grid.ScalarPtr());

  nrrdNuke(nin);

  nrrd_info.dimension = dimension;

  io_time.read_nrrd_time = wall_time.getElapsed();
}

// **************************************************
// READ ISOSURFACE LOOKUP TABLE(S)
// **************************************************

void IJKMCUBE::read_cube_isotable
(const int dimension, const ISOTABLE_TYPE isotable_type,
 const string & isotable_directory,
 ISOSURFACE_TABLE & isotable, IO_TIME & io_time)
{
  read_isosurface_table
    (dimension, "cube", isotable_type, isotable_directory, 
     isotable, io_time);
}

void IJKMCUBE::read_poly_isotable
(const std::string & isotable_directory, MC_DATA & mc_data, 
 IO_TIME & io_time)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  PROCEDURE_ERROR error("read_poly_isotable");

  const ISOTABLE_TYPE isotable_type = get_isotable_type(mc_data);

  if (mc_data.UseNEP() || mc_data.Snap()) {
    read_cube_isotable(dimension, isotable_type, isotable_directory, 
		       mc_data.isotable.cube_nep, io_time);
    set_in_facet_iso_patches(mc_data.isotable.cube_nep);
  }
  else {
    read_cube_isotable(dimension, isotable_type, isotable_directory,
		       mc_data.isotable.cube, io_time);
  }

  ISOSURFACE_TOPOLOGY isosurface_topology = mc_data.IsosurfaceTopology();
  bool flag_multires = mc_data.UseMultires();

  if (!check_dimension(mc_data.isotable.cube, mc_data.ScalarGrid(), error))
    { throw(error); };

  if (isosurface_topology != ISOTABLE_TOPOLOGY) {
    IJK::ERROR warning;

    if (!check_cube_isotable_fits_topology
	(mc_data.isotable.cube, isosurface_topology, warning)) {
      cerr << endl;
      cerr << "*** Warning: ";
      warning.Print(cerr);
      cerr << "  Continuing..." << endl;
      cerr << endl;
    }
  }

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY ||
      flag_multires) {

    read_isosurface_table
      (dimension, "pyramid", isotable_type, isotable_directory,
       mc_data.isotable.pyramid, io_time);

    if (!check_dimension(mc_data.isotable.pyramid, 
			 mc_data.ScalarGrid(), error))
      { throw(error); };
  }

  if (flag_multires) {

    read_isosurface_table
      (dimension, "simplex", isotable_type, isotable_directory,
       mc_data.isotable.simplex, io_time);

    if (!check_dimension(mc_data.isotable.simplex, 
			 mc_data.ScalarGrid(), error))
      { throw(error); };
  }

  if (isosurface_topology != ISOTABLE_TOPOLOGY) 
    { mc_data.isotable.ComputeAmbiguityInformation(); }
}

void IJKMCUBE::read_isosurface_table_from_file
(const string & isotable_filename, const string & isotable_directory,
 ISOSURFACE_TABLE & isotable, const int nrrd_dimension)
{
  string isotable_pathname = 
    isotable_directory + "/" + isotable_filename;

  ifstream isotable_file(isotable_pathname.c_str(), ios::in);
  if (!isotable_file) {
    isotable_file.clear();
    isotable_file.open(isotable_filename.c_str(), ios::in);
  };

  if (!isotable_file) {
    cerr << "Unable to obtain isosurface table file "
	 << isotable_pathname << "." << endl;
    exit(30);
  };

  try {
    IJKXIO::read_xit(isotable_file, isotable);
  }
  catch(...) {
    cerr << "Error reading file: " << isotable_filename << "." << endl;
    throw;
  };

  isotable_file.close();

  ERROR error_msg;
  if (!isotable.Check(error_msg)) {
    cerr << "Warning: Data structure inconsistency in isosurface table "
	 << isotable_pathname << "." << endl;
    for (int i = 0; i < error_msg.NumMessages(); i++) 
      cerr << error_msg.Message(i) << endl;
    cerr << "  Attempting to continue..." << endl << endl;
  }

}

void IJKMCUBE::read_isosurface_table
(const int dimension, const char * poly_name,
 const ISOTABLE_TYPE isotable_type,
 const string & isotable_directory,
 ISOSURFACE_TABLE & isotable, IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  wall_time.getElapsed();

  string filename = 
    get_isotable_filename(isotable_type, dimension, poly_name);

  read_isosurface_table_from_file
    (filename, isotable_directory, isotable, dimension);

  io_time.read_table_time = wall_time.getElapsed();
}


ISOTABLE_TYPE IJKMCUBE::get_isotable_type(const MC_DATA & mc_data)
{
  ISOTABLE_TYPE isotable_type;
  if (mc_data.UseNEP() || mc_data.Snap()) { isotable_type = NEP;  }
  else if (mc_data.IntervalVolumeFlag()) { isotable_type = IVOL; }
  else { isotable_type = BINARY; };

  return(isotable_type);
}

// Get default isotable directory and get directory from environment.
void IJKMCUBE::get_isotable_directory(std::string & isotable_directory)
{
#ifdef IJK_ISOTABLE_DIR
  isotable_directory = std::string(IJK_ISOTABLE_DIR);
#endif

  const char * envir_isotable_dir = getenv("IJK_ISOTABLE_DIR");
  if (envir_isotable_dir != NULL) {
    isotable_directory = std::string(envir_isotable_dir);
  };
}

// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  output_isosurface(output_info, mc_data, 
		    mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
		    mcube_info, io_time);
}

// newly added
void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time, float* f_color)
{
  output_isosurface(output_info, mc_data, 
		    mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
		    mcube_info, io_time,f_color);
}



void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
	
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }

  if (!output_info.nowrite_flag) 
    { write_mesh(output_info, vertex_coord, slist, io_time); }
}

// newly added
void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time, float* f_color)
{
	
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }

  if (!output_info.nowrite_flag) 
    { //write_mesh(output_info, vertex_coord, slist, io_time); 

	write_mesh(output_info, vertex_coord, slist, io_time, f_color); 
   }
}


void IJKMCUBE::output_nep_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  output_nep_isosurface
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     mcube_info, io_time);
}

void IJKMCUBE::output_nep_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_nep_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }

  if (!output_info.nowrite_flag) 
    { write_mesh(output_info, vertex_coord, slist, io_time); }
}

void IJKMCUBE::output_snapmc_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const SNAP_INFO & snap_info, IO_TIME & io_time)
{
  output_snapmc_isosurface
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     snap_info, io_time);
}

void IJKMCUBE::output_snapmc_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const SNAP_INFO & snap_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_snap_info(output_info, mc_data, vertex_coord, slist, snap_info);
  }

  if (!output_info.nowrite_flag) 
    { write_mesh(output_info, vertex_coord, slist, io_time); }
}

void IJKMCUBE::output_isosurface_color
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  output_isosurface_color
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     front_color, back_color, mcube_info, io_time);
}

void IJKMCUBE::output_isosurface_color
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  if (!output_info.use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }
  
  if (!output_info.nowrite_flag) {
    write_mesh_color
      (output_info, vertex_coord, slist, front_color, back_color, io_time);
  }
}

void IJKMCUBE::output_isosurface_color_alternating
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info, IO_TIME & io_time)
{
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const VERTEX_INDEX nums = 
    mc_isosurface.simplex_vert.size()/numv_per_simplex;

  IJK::ARRAY<COLOR_TYPE> front_color(4*nums);
  IJK::ARRAY<COLOR_TYPE> back_color(4*nums);
  set_color_alternating
    (mc_data.ScalarGrid(), mc_isosurface.cube_containing_simplex,
     front_color.Ptr());
  set_color_alternating
    (mc_data.ScalarGrid(), mc_isosurface.cube_containing_simplex,
     back_color.Ptr());

  output_isosurface_color
    (output_info, mc_data, mc_isosurface.vertex_coord,
     mc_isosurface.simplex_vert, 
     front_color.PtrConst(), back_color.PtrConst(),
     mcube_info, io_time);
}

// **************************************************
// WRITE_MESH
// **************************************************

void IJKMCUBE::write_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_simplex;
  const bool use_stdout = output_info.use_stdout;

  ofstream output_file;
  ERROR error_mcube("write_mesh");

  string ofilename = output_info.output_filename;
  switch (output_info.output_format) {

  case OFF:
    if (!use_stdout) {
      output_file.open(ofilename.c_str(), ios::out);  	
 
	int numv = vertex_coord.size()/dimension;
	cout<<"Number of vertices:"<<numv<<"\n";
	 COLOR_TYPE front_color[4*numv];
	 COLOR_TYPE back_color[4*numv];

	const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
  	const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
  	const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

	COLOR_TYPE* color_ptr_f = front_color;
	COLOR_TYPE* color_ptr_b = back_color;

	for(int l=0; l< 4*numv; l++)
	{
		if(l%4 == 0)
		{
			color_ptr_f[l] = 0.5;
			color_ptr_b[l] = 0.5;
		}	
		else if(l%4 == 3)
		{
			color_ptr_f[l] = 1;
			color_ptr_b[l] = 1;
		}	
		else
		{
			color_ptr_f[l] = 0.5;
			color_ptr_b[l] = 0.5;
		}

	}
	
	 ijkoutColorVertOFF(output_file, dimension, numv_per_simplex,
			  vertex_coord, slist, front_color, back_color);

	      //ijkoutOFF(output_file, dimension, numv_per_simplex,
		//	vertex_coord, slist);
      output_file.close();
    }
    else {
      ijkoutOFF(dimension, numv_per_simplex, vertex_coord, slist);
    };
    break;

  case IV:
    if (dimension == 3) {
      if (!use_stdout) {
	output_file.open(ofilename.c_str(), ios::out);
	ijkoutIV(output_file, dimension, vertex_coord, slist);
	output_file.close();
      }
      else {
	ijkoutOFF(dimension, vertex_coord, slist);
      }
    }
    else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
    break;

  default:
    throw error_mcube("Illegal output format.");
    break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}


// Newly added!
void IJKMCUBE::write_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist, float* f_color)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_simplex;
  const bool use_stdout = output_info.use_stdout;

  ofstream output_file;
  ERROR error_mcube("write_mesh");

  string ofilename = output_info.output_filename;

  switch (output_info.output_format) {

  case OFF:
    if (!use_stdout) {
      output_file.open(ofilename.c_str(), ios::out);  
	

	int numv = vertex_coord.size()/dimension;
	cout<<"Number of vertices:"<<numv<<"\n";
        // COLOR_TYPE back_color[4*numv];

	 float *back_color = new float[4*numv]; 
     	 for(int l=0; l< 4*numv; l++)
	 {
		back_color[l] = f_color[l];	
	 }

   /*       for (int ct = 0; ct < numv; ct++) {
	cout<<"\n";
	cout<< f_color[4*ct + 0]<<" ";	
 	cout<< f_color[4*ct + 1]<<" " ;
	cout<< f_color[4*ct + 2]<<" " ;
	cout<< f_color[4*ct + 3]<<" " ;		
    }*/


	/* COLOR_TYPE front_color[4*numv];
	

	const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
  	const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
  	const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

	COLOR_TYPE* color_ptr_f = front_color;
	COLOR_TYPE* color_ptr_b = back_color;

	for(int l=0; l< 4*numv; l++)
	{
		if(l%4 == 0)
		{
			color_ptr_f[l] = 0.5;
			color_ptr_b[l] = 0.5;
		}	
		else if(l%4 == 3)
		{
			color_ptr_f[l] = 1;
			color_ptr_b[l] = 1;
		}	
		else
		{
			color_ptr_f[l] = 0;
			color_ptr_b[l] = 0;
		}

	}*/
	
	 ijkoutColorVertOFF(output_file, dimension, numv_per_simplex,
			  vertex_coord, slist, f_color, back_color);

	      //ijkoutOFF(output_file, dimension, numv_per_simplex,
		//	vertex_coord, slist);
      output_file.close();
    }
    else {
      ijkoutOFF(dimension, numv_per_simplex, vertex_coord, slist);
    };
    break;

  case IV:
    if (dimension == 3) {
      if (!use_stdout) {
	output_file.open(ofilename.c_str(), ios::out);
	ijkoutIV(output_file, dimension, vertex_coord, slist);
	output_file.close();
      }
      else {
	ijkoutOFF(dimension, vertex_coord, slist);
      }
    }
    else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
    break;

  default:
    throw error_mcube("Illegal output format.");
    break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}






void IJKMCUBE::write_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_mesh(output_info, vertex_coord, slist);

  io_time.write_time += wall_time.getElapsed();
}


// newly added
void IJKMCUBE::write_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 IO_TIME & io_time, float * f_color)
{
  ELAPSED_TIME wall_time;

  write_mesh(output_info, vertex_coord, slist, f_color);

  io_time.write_time += wall_time.getElapsed();
}



void IJKMCUBE::write_mesh_color
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_simplex;
  const bool use_stdout = output_info.use_stdout;
  
  

  ofstream output_file;
  ERROR error_mcube("write_mesh_color");

  string ofilename = output_info.output_filename;

  switch (output_info.output_format) {

  case OFF:
    if (!use_stdout) {
      output_file.open(ofilename.c_str(), ios::out);
      ijkoutColorFacesOFF(output_file, dimension, numv_per_simplex,
			  vertex_coord, slist, front_color, back_color);
      output_file.close();
    }
    else {
      ijkoutColorFacesOFF(std::cout, dimension, numv_per_simplex,
			  vertex_coord, slist, front_color, back_color);
    };
    break;

  case IV:
    if (dimension == 3) {
      if (!use_stdout) {
	output_file.open(ofilename.c_str(), ios::out);
	ijkoutIV(output_file, dimension, vertex_coord, slist);
	output_file.close();
      }
      else {
	ijkoutOFF(dimension, vertex_coord, slist);
      }
    }
    else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
    break;

  default:
    throw error_mcube("Illegal output format.");
    break;
  }

  if (!use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}

void IJKMCUBE::write_mesh_color
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_mesh_color(output_info, vertex_coord, slist, 
		   front_color, back_color);

  io_time.write_time += wall_time.getElapsed();
}

// **************************************************
// RESCALE ROUTINES
// **************************************************

namespace {

  void grow_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
    for (int i = 0; i < vertex_coord.size(); i++) {
      vertex_coord[i] = scale * vertex_coord[i];
    };
  }

  void shrink_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
    for (int i = 0; i < vertex_coord.size(); i++) {
      vertex_coord[i] = vertex_coord[i]/scale;
    };
  }

  bool unit_spacing(const std::vector<COORD_TYPE> & spacing)
    // return true if spacing not defined or spacing along all axes equals 1.0
  {
    for (int d = 0; d < spacing.size(); d++) {
      if (!AIR_EXISTS(spacing[d])) { return(true); }
      else if (spacing[d] != 1.0) { return(false); };
    }

    return(true);
  }

  void rescale_coord(const std::vector<COORD_TYPE> & grid_spacing,
		     std::vector<COORD_TYPE> & vertex_coord)
  {
    const int dimension = grid_spacing.size();

    if (unit_spacing(grid_spacing)) { return; }

    const VERTEX_INDEX numv = vertex_coord.size()/dimension;
    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dimension; d++) {
	vertex_coord[iv*dimension+d] *= grid_spacing[d];
      }
    };
  }

};

/// Rescale subsampled/supersampled vertex coordinates.
/// Also rescale to reflect grid spacing.
void IJKMCUBE::rescale_vertex_coord
(const OUTPUT_INFO & output_info, vector<COORD_TYPE> & vertex_coord)
{
  const int grow_factor = output_info.grow_factor;
  const int shrink_factor = output_info.shrink_factor;
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (grow_factor <= 0) {
    error.AddMessage("Illegal grow factor ", grow_factor, ".");
    error.AddMessage("  Grow factor must be a positive integer");
  }

  if (shrink_factor <= 0) {
    error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
    error.AddMessage("  Shrink factor must be a positive integer");
  }

  if (output_info.dimension != output_info.grid_spacing.size()) {
    error.AddMessage("Size of grid spacing array does not equal volume dimension.");
    error.AddMessage("  Grid spacing array has ", 
		     output_info.grid_spacing.size(), " elements.");
    error.AddMessage("  Volume dimension = ", output_info.dimension, ".");
  }

  if (output_info.grow_factor != 1) 
    { grow_coord(output_info.grow_factor, vertex_coord); };

  if (output_info.shrink_factor != 1) 
    { shrink_coord(output_info.shrink_factor, vertex_coord); };

  rescale_coord(output_info.grid_spacing, vertex_coord);
}

/// Rescale subsampled/supersampled vertex coordinates.
/// Also rescale to reflect grid spacing.
void IJKMCUBE::rescale_vertex_coord
(const int grow_factor, const int shrink_factor,
 const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord)
{
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (grow_factor <= 0) {
    error.AddMessage("Illegal grow factor ", grow_factor, ".");
    error.AddMessage("  Grow factor must be a positive integer");
  }

  if (shrink_factor <= 0) {
    error.AddMessage("Illegal shrink factor ", shrink_factor, ".");
    error.AddMessage("  Shrink factor must be a positive integer");
  }

  if (vertex_coord.size() == 0) { return; };

  if (grid_spacing.size() < 1) {
    error.AddMessage("Illegal size ", grid_spacing.size(), 
		     " of array grid spacing.");
    error.AddMessage("Size must equal vertex dimension.");
    throw error;
  }

  if (grow_factor != 1) 
    { grow_coord(grow_factor, vertex_coord); };

  if (shrink_factor != 1) 
    { shrink_coord(shrink_factor, vertex_coord); };

  rescale_coord(grid_spacing, vertex_coord);
}

// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

void IJKMCUBE::report_num_cubes
(const MC_GRID & full_scalar_grid, const IO_INFO & io_info, 
 const MC_DATA & mc_data)
{
  const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
  const int num_cubes_in_mc_data = 
    mc_data.ScalarGrid().ComputeNumCubes();

  if (!io_info.use_stdout && !io_info.flag_silent) {

    if (io_info.flag_subsample) {
      // subsampled grid
      cout << num_grid_cubes << " grid cubes.  "
	   << num_cubes_in_mc_data << " subsampled grid cubes." << endl;
    }
    else if (io_info.flag_supersample) {
      // supersample grid
      cout << num_grid_cubes << " grid cubes.  "
	   << num_cubes_in_mc_data << " supersampled grid cubes." << endl;
    }
    else {
      // use full_scalar_grid
      cout << num_grid_cubes << " grid cubes." << endl;
    }
  }

}

void IJKMCUBE::report_iso_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const MCUBE_INFO & mcube_info)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.num_vertices_per_simplex;
  const int isosurface_topology = mc_data.IsosurfaceTopology();

  const char * indent4 = "    ";
  string grid_element_name = "cubes";
  if (dimension == 2) { grid_element_name = "squares"; };

  VERTEX_INDEX numv = (vertex_coord.size())/dimension;
  VERTEX_INDEX nums = (slist.size())/numv_per_simplex;
  VERTEX_INDEX num_grid_cubes = mcube_info.grid.num_cubes;
  VERTEX_INDEX num_non_empty_cubes = mcube_info.scalar.num_non_empty_cubes;

  float percent = 0.0;
  if (num_grid_cubes > 0)
    { percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
  int ipercent = int(100*percent);
  cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
       << numv << " isosurface vertices.  "
       << nums << " isosurface simplices." << endl;
  cout << indent4 << num_non_empty_cubes
       << " (" << ipercent << "%) non-empty " << grid_element_name << "." << endl;

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY) {
    cout << indent4 << mcube_info.scalar.num_ambiguous_cubes 
	 << " ambiguous " << grid_element_name << "." << endl;
    if (dimension > 2) {
      cout << indent4 << mcube_info.scalar.num_non_empty_pyramids << " non-empty pyramids." << endl;
      cout << indent4 << mcube_info.scalar.num_ambiguous_pyramids << " ambiguous_pyramids." << endl;
    }
  }

/* saddle topology not implemented
  if (mcube_info.scalar.Dimension() == 3 && 
      isosurface_topology == SADDLE_TOPOLOGY) {
    cout << indent4 << mcube_info.scalar.NumCubesWithSaddle(1)
	 << " cubes with 1 saddle." << endl;
    cout << indent4 << mcube_info.scalar.NumCubesWithSaddle(2)
	 << " cubes with 2 saddles." << endl;
  }
*/

}

void IJKMCUBE::report_nep_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const MCUBE_INFO & mcube_info)
{
  report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);

  if (!output_info.use_stdout && !output_info.flag_silent) {
    cout << "    " << mcube_info.nep.num_non_empty_boundary_facets
	 << " non-empty boundary cube facets." << endl;

    if (mc_data.NEPNumDup() != 2) {
      cout << "    " << mcube_info.nep.num_in_facet_cubes
	   << " cubes with isosurface patches contained in a facet." << endl;
      cout << "    " << mcube_info.nep.num_dup_iso_patches
	   << " duplicate isosurface patches eliminated." << endl;
    }
  }
}

void IJKMCUBE::report_snap_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const SNAP_INFO & snap_info)
{
  report_nep_info(output_info, mc_data, vertex_coord, slist, snap_info);

  if (!output_info.use_stdout && !output_info.flag_silent) {
    cout << "    " << snap_info.num_snapped_iso_vertices
	 << " snapped isosurface vertices." << endl;
  }
}

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

void IJKMCUBE::report_mcube_time
(const IO_INFO & io_info, const MCUBE_TIME & mcube_time, 
 const char * mesh_type_string)
{
  cout << "CPU time to run Marching Cubes: " 
       << mcube_time.total << " seconds." << endl;
  if (io_info.use_octree) {
    cout << "    Time to create octree: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.use_minmax) {
    cout << "    Time to create min/max regions: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.snap_flag) {
    cout << "    Time to snap grid scalar values: "
	 << mcube_time.snap << " seconds." << endl;
  }
  if (io_info.use_list) {
    cout << "    Time to create list of non-empty cubes: "
	 << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.use_multires) {
    cout << "    Time to preprocess multi-resolution mesh: "
	 << mcube_time.process_multires << " seconds." << endl;
  }
  cout << "    Time to extract " << mesh_type_string << " triangles: "
       << mcube_time.extract << " seconds." << endl;
  cout << "    Time to merge identical "
       << mesh_type_string << " vertices: " 
       << mcube_time.merge << " seconds." << endl;
  cout << "    Time to position "
       << mesh_type_string << " vertices: "
       << mcube_time.position << " seconds." << endl;
}


void IJKMCUBE::report_time
(const IO_INFO & io_info, const IO_TIME & io_time, 
 const MCUBE_TIME & mcube_time, const double total_elapsed_time)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * INTERVAL_VOLUME_STRING = "interval volume";
  const char * mesh_type_string = NULL;
  
  if (!io_info.interval_volume_flag) {
    mesh_type_string = ISOSURFACE_STRING;
  }
  else {
    mesh_type_string = INTERVAL_VOLUME_STRING;
  };

  cout << "Time to read file " << io_info.input_filename << ": "
       << io_time.read_nrrd_time << " seconds." << endl;

  cout << "Time to read " << mesh_type_string << " lookup tables: "
       << io_time.read_table_time << " seconds." << endl;

  report_mcube_time(io_info, mcube_time, mesh_type_string);
  if (!io_info.nowrite_flag) {
    cout << "Time to write "
	 << mesh_type_string << ": " 
	 << io_time.write_time << " seconds." << endl;
  };
  cout << "Total elapsed time: " << total_elapsed_time
       << " seconds." << endl;
}

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

// local namespace
namespace {

void usage_msg(std::ostream & out)
{
  out << "Usage: ijkmcube [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
}

void options_msg()
{
  cerr << "OPTIONS:" << endl;
  cerr << "[-octree] [-region L] [-list]" << endl;
  cerr << "[-nep] [-snap D] [-ivol]" << endl;
  cerr << "[-highres \"region coordinates\"] [-numres_levels N]" << endl;
  cerr << "[-topology {isotable|cube_decider|adecider}]" << endl;
  cerr << "[-interpolate {linear|multilinear}]" << endl;
  cerr << "[-subsample S] [-supersample S]" << endl;
  cerr << "[-color_alternating]" << endl;
  cerr << "[-off|-iv] [-dir {isotable_directory}] [-o {output_filename}] [-stdout]" 
       << endl; 
  cerr << "[-help] [-s] [-nowrite] [-time]" << endl;
}

};

void IJKMCUBE::usage_error()
{
  usage_msg(cerr);
  options_msg();
  exit(10);
}

void IJKMCUBE::help()
{
  usage_msg(cout);
  cout << endl;
  cout << "ijkmcube - Marching cubes isosurface generation algorithm." << endl;
  cout << endl;
  cout << "OPTIONS:" << endl;

  cout << "  -octree: Create octree for faster isosurface extraction."
       << endl;
  cout << "  -region L: Preprocess into regions of dimension LxLx... for faster"
       << endl;
  cout << "              isosurface extraction.  Each region has a minimum and maximum" << endl;
  cout << "              scalar value."
       << endl;
  cout << "  -nep:  Use isotable which differentiates scalar values less than (negative)," << endl;
  cout << "         equal to, and greater than (positive) the isovalue." << endl;
  cout << "  -snap D: Snap isosurface vertices within distance D of grid vertices." << endl;
  cout << "           Value D must be in range [0.0, 0.5]." << endl;
  cout << "  -ivol: Generate interval volume." << endl;
  cout << "  -subsample S: Subsample grid at every S vertices." << endl;
  cout << "                S must be an integer greater than 1." << endl;
  cout << "  -highres \"region coordinates\":" << endl;
  cout << "          High resolution region in multiresolution isosurface."
       << endl;
  cout << "          Coordinates must be surrounded by quotation marks."
       << endl;
  cout << "          List coordinates of lowest point in region followed by coordinates of highest point." << endl;
  cout << "  -topology T:  Construct isosurface with topology type T." 
       << endl;
  cout << "     Topology types:" << endl;
  cout << "       isotable: Topology based on isosurface lookup table."
       << endl;
  cout << "       cube_decider: Topology based on resolving cube ambiguities." << endl;
  cout << "       adecider: Asymptotic decider topology." << endl;
  cout << "  -interpolate {linear|multilinear}: Use linear or multilinear interpolation." << endl;
  cout << "          Multilinear interpolation only works with asymptotic decider or linear topology." << endl;
  cout << "  -list:  Preprocess by creating list of mixed cubes." << endl;
  cout << "  -color_alternating: Color alternating cubes.  Works only with -off." << endl;
  cout << "  -off: Output in geomview OFF format. (Default.)" << endl;
  cout << "  -iv: Output in OpenInventor .iv format." << endl;
  cout << "  -dir {isotable_directory}: Directory containing appropriate isosurface table." << endl;
  cout << "  -o {output_filename}: Write isosurface to file {output_filename}." << endl;
  cout << "  -stdout: Write isosurface to standard output." << endl;
  cout << "  -nowrite: Don't write isosurface." << endl;
  cout << "  -time: Output running time." << endl;
  cout << "  -s: Silent mode." << endl;
  cout << "  -help: Print this help message." << endl;
  exit(20);
}


// **************************************************
// CLASS IO_INFO
// **************************************************

/// IO information
void IJKMCUBE::IO_INFO::Init()
{
  isovalue.clear();
  isovalue_string.clear();
  input_filename = NULL;
  output_filename = NULL;
  isotable_directory = "";
  output_format = OFF;
  report_time_flag = false;
  use_stdout = false;
  nowrite_flag = false;
  flag_silent = false;
  flag_subsample = false;
  subsample_resolution = 2;
  flag_supersample = false;
  supersample_resolution = 2;
  flag_color_alternating = false;  // color simplices in alternating cubes
  region_length = 1;
};

// **************************************************
// class OUTPUT_INFO
// **************************************************

void IJKMCUBE::OUTPUT_INFO::Init()
{
  output_filename = "";
  dimension = 3;
  num_vertices_per_simplex = 3;
  isovalue[0] = 0;
  isovalue[1] = 0;
  nowrite_flag = false;
  use_stdout = false;
  flag_silent = false;
  output_format = OFF;
  grow_factor = 1;
  shrink_factor = 1;
  grid_spacing.resize(3,1);
}

namespace {

  void split_string(const string & s, const char c,
		    string & prefix, string & suffix)
    // split string at last occurrence of character c into prefix and suffix
  {
    string::size_type i = s.rfind(c);
    if (i == string::npos) {
      prefix = s;
      suffix = "";
    }
    else {
      if (i > 0) { prefix = s.substr(0,i); }
      else { prefix = ""; };

      if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
      else { suffix = ""; };
    }
  }

  string construct_output_filename
  (const IO_INFO & io_info, const int i)
  {
    // create output filename
    string fname = string(io_info.input_filename);

    // remove path from file name
    string prefix, suffix;
    split_string(fname, '/', prefix, suffix);
    if (suffix != "") { fname = suffix; }

    string ofilename;

    // construct output filename
    split_string(fname, '.', prefix, suffix);
    if (suffix == "nrrd" || suffix == "nhdr") { ofilename = prefix; }
    else { ofilename = string(io_info.input_filename); }

    if (!io_info.interval_volume_flag) {
      ofilename += string(".") + string("isov=") + io_info.isovalue_string[i];
    }
    else {
      ofilename += string(".") + string("ivol=") + io_info.isovalue_string[i]
	+ "-" + io_info.isovalue_string[i+1];
    }

    switch (io_info.output_format) {
    case OFF: 
      ofilename += ".off";
      break;

    case IV:
      ofilename += ".iv";
      break;
    }

    return(ofilename);
  }

}

// **************************************************
// SET ROUTINES
// **************************************************

void IJKMCUBE::set_mc_data
(const IO_INFO & io_info, MC_DATA & mc_data, MCUBE_TIME & mcube_time)
{
  const bool use_stdout = io_info.use_stdout;
  const bool flag_silent = io_info.flag_silent;
  PROCEDURE_ERROR error("set_mc_data");

  if (!mc_data.IsScalarGridSet()) {
    error.AddMessage("Programming error. Scalar field must be set before set_mc_data is called.");
    throw error;
  }
  
  // Set data structures in mc data  
  clock_t t0 = clock();

  if (io_info.use_octree) {
    if (io_info.snap_flag)
      { mc_data.SetSnapOctree(); }
    else
      { mc_data.SetOctree(); }
  }
  else if (io_info.use_minmax) {
    if (io_info.snap_flag) 
      { mc_data.SetSnapMinmaxRegions(io_info.region_length); }
    else
      { mc_data.SetMinmaxRegions(io_info.region_length); };
  }

  // Set flags in mc data  
  if (io_info.use_nep) {
    mc_data.SetNEPOn(io_info.nep_num_dup);
  };

  if (io_info.interval_volume_flag)
    { mc_data.SetIntervalVolumeFlag(true); };

  if (io_info.snap_flag) 
    { mc_data.SetSnapOn(io_info.snap_value, io_info.nep_num_dup); }

  if (io_info.use_multires) {
    mc_data.SetHighResolutionRegions(io_info.high_resolution_regions);
    mc_data.SetMultiresGrid(io_info.num_resolution_levels);
  }

  if (io_info.use_octree || io_info.use_minmax || io_info.use_multires) {
    clock_t t1 = clock();
    mcube_time.preprocessing += float(t1-t0)/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;
  }

  if (io_info.use_list) { mc_data.SetUseList(true); }

  if (io_info.flag_color_alternating) 
    { mc_data.SetCubeContainingSimplexFlag(true); };

  mc_data.SetEdgeRepresentation(io_info.edge_representation);
  mc_data.SetIsosurfaceTopology(io_info.isosurface_topology);
  mc_data.SetInterpolationType(io_info.interpolation_type);

  mc_data.merge_edges_parameters = io_info.merge_edges_parameters;
}

void IJKMCUBE::set_io_info
(const NRRD_INFO & nrrd_info, IO_INFO & io_info)
{
  io_info.grid_spacing.clear();
  for (int d = 0; d < nrrd_info.dimension; d++) {
    io_info.grid_spacing.push_back(nrrd_info.grid_spacing[d]);
  }
}

void IJKMCUBE::set_output_info
(const ISOSURFACE_TABLE & isotable, const IO_INFO & io_info, 
 const int i, OUTPUT_INFO & output_info)
{
	
  output_info.dimension = isotable.Dimension();
  output_info.num_vertices_per_simplex = isotable.NumVerticesPerSimplex();
  output_info.nowrite_flag = io_info.nowrite_flag;
  output_info.use_stdout = io_info.use_stdout;
  output_info.flag_silent = io_info.flag_silent;

  output_info.grow_factor = 1;
  if (io_info.flag_subsample) 
    { output_info.grow_factor = io_info.subsample_resolution; }

  output_info.shrink_factor = 1;
  if (io_info.flag_supersample) 
    { output_info.shrink_factor = io_info.supersample_resolution; }

  output_info.grid_spacing.clear();
  for (int j = 0; j < io_info.grid_spacing.size(); j++)
    { output_info.grid_spacing[j] = io_info.grid_spacing[j]; }

  output_info.output_format = io_info.output_format;
  output_info.isovalue[0] = io_info.isovalue[i];
  if (i+1 < io_info.isovalue.size()) 
    { output_info.isovalue[1] = io_info.isovalue[i+1]; };

  if (io_info.output_filename != NULL) {
    output_info.output_filename = string(io_info.output_filename);
  }
  else {
    output_info.output_filename = 
      construct_output_filename(io_info, i);
  }
}

void IJKMCUBE::set_color_alternating
(const MC_GRID & grid, const vector<VERTEX_INDEX> & cube_list, 
 COLOR_TYPE * color)
{
  const int dimension = grid.Dimension();
  GRID_COORD_TYPE coord[dimension];

  const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
  const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
  const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

  VERTEX_INDEX icube = 0;
  int parity = 0;
  COLOR_TYPE * color_ptr = color;
  for (int i = 0; i < cube_list.size(); i++) {
    int new_cube = cube_list[i];
    if (icube != new_cube) {
      icube = new_cube;
      grid.ComputeCoord(icube, coord);
      int sum = 0;
      for (int d = 0; d < dimension; d++) 
	{ sum += coord[d]; }
      parity = sum%2;
    }

    if (parity == 0) 
      { std::copy(red, red+3, color_ptr); }
    else
      { std::copy(blue, blue+3, color_ptr); }

    // set opacity
    color_ptr[3] = 1.0;
    color_ptr += 4;
  }
  
}

// **************************************************
// NRRD INFORMATION
// **************************************************

NRRD_INFO::NRRD_INFO()
{
  Clear();
}

NRRD_INFO::~NRRD_INFO()
{
  Clear();
}

void NRRD_INFO::Clear()
{
  dimension = 0;
}
