/// \file ijkmcubeIO.h
/// IO classes and routines for ijkmcube.

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

#ifndef _IJKMCUBEIO_
#define _IJKMCUBEIO_

#include <string>

#include "ijk.txx"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"
#include "ijksnapmc.h"
#include "ijkNrrd.h"

namespace IJKMCUBE {

// **************************************************
// TYPE DEFINITIONS
// **************************************************

  typedef float COLOR_TYPE;           /// Color type.

  using IJKSNAPMC::SNAP_INFO;

// **************************************************
// NRRD INFORMATION
// **************************************************

  /// NRRD information.
  class NRRD_INFO {

  public:
    int dimension;                        ///< Volume dimension.
    COORD_ARRAY grid_spacing;             ///< Grid spacing.

  public:
    NRRD_INFO();
    ~NRRD_INFO();

    void Clear();
  };
		       
// **************************************************
// OUTPUT_FORMAT
// **************************************************

  typedef enum { OFF, IV } OUTPUT_FORMAT;   ///< Output format.

// **************************************************
// IO INFORMATION
// **************************************************

  /// IO information
  class IO_INFO:public MC_DATA_FLAGS {

  protected:
    void Init();

  public:
    SCALAR_ARRAY isovalue;        ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    COORD_ARRAY grid_spacing;
    char * input_filename;
    char * output_filename;
    std::string isotable_directory;
    //std::string isotable_directory_1;
    OUTPUT_FORMAT output_format;
    bool report_time_flag;
    bool use_stdout;
    bool nowrite_flag;
    bool flag_silent;
    bool flag_subsample;
    int subsample_resolution;
    bool flag_supersample;
    int supersample_resolution;
    bool flag_color_alternating;  ///< Color simplices in alternating cubes
    int region_length;

    /// Merge edges parameters.
    MERGE_EDGES_PARAMETERS merge_edges_parameters;

    /// List of high resolution arguments,
    ///   e.g., "-highres {coord list}".
    std::vector<std::string> high_resolution_option;

  public:
    IO_INFO() { Init(); };
    ~IO_INFO() { Init(); };
  };

// **************************************************
// OUTPUT INFORMATION
// **************************************************

  /// Output information.
  class OUTPUT_INFO {

  protected:
    void Init();

  public:
    std::string output_filename;
    int dimension;
    int num_vertices_per_simplex;
    SCALAR_TYPE isovalue[2];
    bool nowrite_flag;
    bool use_stdout;
    bool flag_silent;
    OUTPUT_FORMAT output_format;
    COORD_ARRAY grid_spacing;
    int grow_factor;
    int shrink_factor;

    OUTPUT_INFO() { Init(); };
    ~OUTPUT_INFO() { Init(); };
  };

// **************************************************
// TIMING FUNCTIONS/CLASSES
// **************************************************

  /// Elapsed CPU time.
  class ELAPSED_CPU_TIME {

  protected:
    clock_t t;

  public:
    ELAPSED_CPU_TIME() { t = clock(); };

    clock_t getElapsed() {
      clock_t old_t = t;
      t = clock();
      return(t - old_t);
    };
  };

  /// Elapsed wall time.
  class ELAPSED_TIME {

  protected:
    time_t t;

  public:
    ELAPSED_TIME() { time(&t);  };

    double getElapsed() {
      time_t old_t = t;
      time(&t);
      return(difftime(t,old_t));
    };
  };

  /// IO time.
  struct IO_TIME {
    double read_table_time; ///< Wall time to read isosurface lookup table.
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };

// **************************************************
// PARSE COMMAND LINE
// **************************************************

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, IO_INFO & io_info);

  /// Check input information in io_info
  bool check_input
    (const IO_INFO & io_info, 
     const MC_SCALAR_GRID_BASE & scalar_grid,
     IJK::ERROR & error);

// **************************************************
// READ NEARLY RAW RASTER DATA (nrrd) FILE
// **************************************************

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, MC_SCALAR_GRID & scalar_grid, 
     NRRD_INFO & nrrd_info, IO_TIME & io_time);

// **************************************************
// READ ISOSURFACE LOOKUP TABLE(S)
// **************************************************

  /// Read cube isosurface lookup table.
  void read_cube_isotable
    (const int dimension, const ISOTABLE_TYPE isotable_type,
     const std::string & isotable_directory,
     ISOSURFACE_TABLE & cube_isotable, IO_TIME & io_time);

  /// Read cube, pyramid and simplex isosurface lookup tables.
  /// Only reads pyramid and simplex tables if needed.
  void read_poly_isotable
    (const std::string & isotable_directory, MC_DATA & mc_data, 
     IO_TIME & io_time);

  /// Determine isosurface table type from io_info.
  ISOTABLE_TYPE get_isotable_type(const MC_DATA & mc_data);

  /// Get default isotable directory and get directory from environment.
  void get_isotable_directory(std::string & isotable_directory);

  /// Read isosurface lookup table from file.
  void read_isosurface_table_from_file
    (const std::string & isotable_filename, 
     const std::string & isotable_directory,
     ISOSURFACE_TABLE & isotable, const int nrrd_dimension);

  /// Read isosurface lookup table.
  /// Construct lookup table file name from dimension, poly_name and isotable_type.
  void read_isosurface_table
    (const int dimension, const char * poly_name,
     const ISOTABLE_TYPE isotable_type,
     const std::string & isotable_directory,
     ISOSURFACE_TABLE & isotable, IO_TIME & io_time);

// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

  void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

  // newly added
    void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time, float* f_color);

  void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);
 
  // newly added
  void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time, float* f_color);

  void output_nep_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

  void output_nep_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

  void output_snapmc_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const SNAP_INFO & snap_info, IO_TIME & io_time);

  void output_snapmc_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const SNAP_INFO & snap_info, IO_TIME & io_time);

  void output_isosurface_color
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

  void output_isosurface_color
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

  void output_isosurface_color_alternating
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info, IO_TIME & io_time);

// **************************************************
// RESCALE ROUTINES
// **************************************************

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
  void rescale_vertex_coord
    (const OUTPUT_INFO & output_info, std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale vertex coordinates by grow and shrink factor and by grid_spacing.
  /// Precondition: grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
    (const int grow_factor, const int shrink_factor,
     const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);

// **************************************************
// WRITE_MESH
// **************************************************

  void write_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist);

   // newly added
   void write_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, float* f_color);


  void write_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time); 

  // Newly added
  void write_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time, float* f_color); 
 

  void write_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  void write_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     IO_TIME & io_time);

// **************************************************
// SET ROUTINES
// **************************************************

  /// Set mc_data based on io_info.
  /// Precondition: Scalar field in mc_data must be set before
  ///   this routines is called.
  void set_mc_data
    (const IO_INFO & io_info, MC_DATA & mc_data, MCUBE_TIME & mcube_time);

  /// Copy nrrd_info into io_info.
  void set_io_info
    (const NRRD_INFO & nrrd_info, IO_INFO & io_info);

  /// Set output_info based on isotable, io_info and isovalue index i.
  void set_output_info
    (const ISOSURFACE_TABLE & isotable, const IO_INFO & io_info, 
     const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
    (const MC_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
     COLOR_TYPE * color);

// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

  void report_num_cubes
    (const MC_GRID & full_grid, const IO_INFO & io_info, 
     const MC_DATA & mc_data);

  void report_iso_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const MCUBE_INFO & mcube_info);

  void report_nep_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const MCUBE_INFO & mcube_info);

  void report_snap_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const SNAP_INFO & snap_info);

// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

  void report_mcube_time
    (const IO_INFO & io_info, const MCUBE_TIME & mcube_time, 
     const char * mesh_type_string);

  void report_time
    (const IO_INFO & io_info, const IO_TIME & io_time, 
     const MCUBE_TIME & mcube_time, const double total_elapsed_time);

// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

  void usage_error();
  void help();
};

#endif
