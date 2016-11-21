/// \file ijkmcube.cxx
/// generate isosurface from scalar field
/// Version 0.2.9

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009, 2008, 2007, 2006, 2003, 2001 Rephael Wenger

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


#include <iostream>

#include "ijkmcubeIO.h"
#include "string.h"


using namespace IJK;
using namespace IJKMCUBE;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_isosurface
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time, const IO_INFO & io_info_1, const MC_DATA & mc_data_1,
 MCUBE_TIME & mcube_time_1, IO_TIME & io_time_1);
void construct_interval_volume
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);
  string isotable_filename;

  MCUBE_TIME mcube_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;


  //change	
  MCUBE_TIME mcube_time_1;
  IO_TIME io_time_1 = {0.0, 0.0, 0.0};
  IO_INFO io_info_1;
  IJK::ERROR error_1;	
 

  try {

    std::set_new_handler(memory_exhaustion);

    std::string isotable_directory;
    get_isotable_directory(isotable_directory);
    io_info.isotable_directory = isotable_directory;

    parse_command_line(argc, argv, io_info);

    // Read mu grid (abc.nrrd) and delta grid (2_abc.nrrd)
    MC_SCALAR_GRID full_scalar_grid;
    NRRD_INFO nrrd_info;
    read_nrrd_file
      (io_info.input_filename, full_scalar_grid,  nrrd_info, io_time);

    if (!check_input(io_info, full_scalar_grid, error)) 
      { throw(error); };

    // copy nrrd_info into io_info
    set_io_info(nrrd_info, io_info);

    // set mc datastructures and flags
    MC_DATA mc_data;

    // Note: mc_data.SetScalarGrid must be called before set_mc_data.
    mc_data.SetScalarGrid
      (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
       io_info.flag_supersample, io_info.supersample_resolution);
    set_mc_data(io_info, mc_data, mcube_time);
    report_num_cubes(full_scalar_grid, io_info, mc_data);

    // Note: All flags in mc_data should be set before calling 
    //       read isosurface lookup tables.
    // read isosurface lookup tables
    read_poly_isotable(isotable_directory, mc_data, io_time);


    char* c = io_info.input_filename;

    // change	
    char buffer [100];

    char* st = new char[2];
    strcpy(st,"2_");  
    sprintf (buffer, strcat(st,c));  
    cout<<"buffer:"<<buffer<<"\n";
   

    MC_SCALAR_GRID full_scalar_grid_1;
    NRRD_INFO nrrd_info_1;
    // set mc datastructures and flags
    MC_DATA mc_data_1;

    INTERPOLATION_TYPE interp = mc_data.InterpolationType();	
   
    // alpha uncertainty interpolation
    if (interp == 2) 
    {
    	read_nrrd_file
      	(buffer, full_scalar_grid_1,  nrrd_info_1, io_time_1);    

    //if (!check_input(io_info_1, full_scalar_grid_1, error_1)) 
    //  { throw(error_1); };

    // copy nrrd_info into io_info
    set_io_info(nrrd_info_1, io_info_1);  

    // Note: mc_data.SetScalarGrid must be called before set_mc_data.
    mc_data_1.SetScalarGrid
      (full_scalar_grid_1, io_info_1.flag_subsample, io_info_1.subsample_resolution,
       io_info_1.flag_supersample, io_info_1.supersample_resolution);
    set_mc_data(io_info_1, mc_data_1, mcube_time_1);   
    report_num_cubes(full_scalar_grid_1, io_info_1, mc_data_1);	
    } 


    // construct isosurface or interval volume
    if (mc_data.IntervalVolumeFlag()) {	
      construct_interval_volume(io_info, mc_data, mcube_time, io_time);
    }
    else {	
      // change
      construct_isosurface(io_info, mc_data, mcube_time, io_time, io_info_1, mc_data_1, mcube_time_1, io_time_1);         		
    }

    
    if (io_info.report_time_flag) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, mcube_time, total_elapsed_time);
    };
  }   
 
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

void construct_isosurface
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time, const IO_INFO & io_info_1, const MC_DATA & mc_data_1,
 MCUBE_TIME & mcube_time_1, IO_TIME & io_time_1 )
{         
   // Array to store colors proportional to variance
   float* f_color;
   float** q  =  &f_color;
   q = NULL;    	
   float* numv;	 

  // mean grid
  const int dimension = mc_data.ScalarGrid().Dimension();
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

  // delta grid
  const int dimension_1 = mc_data_1.ScalarGrid().Dimension();
  const int numv_per_simplex_1 = mc_data_1.isotable.cube.NumVerticesPerSimplex();
  const int num_cubes_1 = mc_data_1.ScalarGrid().ComputeNumCubes();   
  
  io_time.write_time = 0;
  //change
  io_time_1.write_time = 0; 

  for (int i = 0; i < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue = io_info.isovalue[i];

    MC_ISOSURFACE mc_isosurface;
    MCUBE_INFO mcube_info(dimension);

    // mean grid
    mcube_info.grid.num_cubes = num_cubes;
    SNAP_INFO snap_info(dimension);
    snap_info.grid.num_cubes = num_cubes;

    // delta grid
    MC_ISOSURFACE mc_isosurface_1;
    MCUBE_INFO mcube_info_1(dimension_1);
    mcube_info_1.grid.num_cubes = num_cubes_1;   

    if (mc_data.Snap()) {

      if (io_info.use_list) {
	std::vector<VERTEX_INDEX> cube_list;

	float preprocessing_time;
	IJKSNAPMC::get_nonempty_snap_cubes
	  (mc_data.ScalarGrid(), mc_data.isotable.cube_nep, 
	   isovalue, cube_list, preprocessing_time);

	mcube_time.preprocessing += preprocessing_time;
	mcube_time.total += preprocessing_time;
      
	snapMC(mc_data, isovalue, cube_list, mc_isosurface, snap_info);
      }
      else {
	snapMC(mc_data, isovalue, mc_isosurface, snap_info);
      }

      mcube_time.Add(snap_info.time);
    }
    else {	
      
      // marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info);

      // Pass mean (mc_data) and delta (mc_data_1) grids to marching cubes function and get colors in f_color array       
      marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info, mc_data_1, mc_isosurface_1, mcube_info_1, f_color, numv);
      mcube_time.Add(mcube_info.time);      	
    }

    OUTPUT_INFO output_info;
    set_output_info(mc_data.isotable.cube, io_info, i, output_info);

    VERTEX_INDEX nums = 
      mc_isosurface.simplex_vert.size()/numv_per_simplex;

    int grow_factor = 1;
    int shrink_factor = 1;
    if (io_info.flag_subsample) 
      { grow_factor = io_info.subsample_resolution; }
    if (io_info.flag_supersample) 
      { shrink_factor = io_info.supersample_resolution; }

    rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
			 mc_isosurface.vertex_coord);

    if (mc_data.Snap()) {
      output_snapmc_isosurface
	(output_info, mc_data, mc_isosurface, snap_info, io_time);
    }
    else if (mc_data.UseNEP()) {
      output_nep_isosurface
	(output_info, mc_data, mc_isosurface, mcube_info, io_time);
    }
    else if (io_info.flag_color_alternating &&
	     mc_isosurface.cube_containing_simplex.size() == nums) {
      output_isosurface_color_alternating
	(output_info, mc_data, mc_isosurface, mcube_info, io_time);
    }
    else {
	INTERPOLATION_TYPE interpolation_type = mc_data.InterpolationType();	
      // linear interpolation
      if (interpolation_type == 0) 
      {		      		         
      		output_isosurface(output_info, mc_data, mc_isosurface, mcube_info, io_time);	
      		 cout<<"In linear interpolation!\n"; 
      }
      // alpha uncertainty
      else if (interpolation_type == 2)  	
		// Output isosurface (.off format) with colors as obtained in f_color 
       		output_isosurface(output_info, mc_data, mc_isosurface, mcube_info, io_time, f_color);
    }
  }
}

void construct_interval_volume
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;
  for (int i = 0; i+1 < io_info.isovalue.size(); i++) {

    MC_ISOSURFACE mc_ivol;
    MCUBE_INFO mcube_info(dimension);
    mcube_info.grid.num_cubes = num_cubes;
    SNAP_INFO snap_info(dimension);
    snap_info.grid.num_cubes = num_cubes;

    MCVol(mc_data, io_info.isovalue[i], io_info.isovalue[i+1], 
	  mc_ivol, mcube_info);

    mcube_time.Add(mcube_info.time);

    OUTPUT_INFO output_info;
    set_output_info(mc_data.isotable.cube, io_info, i, output_info);

    VERTEX_INDEX nums = 
      mc_ivol.simplex_vert.size()/numv_per_simplex;

    int grow_factor = 1;
    int shrink_factor = 1;
    if (io_info.flag_subsample) 
      { grow_factor = io_info.subsample_resolution; }
    if (io_info.flag_supersample) 
      { shrink_factor = io_info.supersample_resolution; }

    rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
			 mc_ivol.vertex_coord);

    output_isosurface
      (output_info, mc_data, mc_ivol, mcube_info, io_time);
  }

}

void memory_exhaustion()
{
   cerr << "Error: Out of memory.  Terminating program." << endl;
   exit(10);
}

