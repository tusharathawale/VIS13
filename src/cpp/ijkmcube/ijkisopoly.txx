/// \file ijkisopoly.txx
/// ijk templates for extracting an isosurface patch from a polyhedron
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

#ifndef _IJKISOPOLY_
#define _IJKISOPOLY_

#include "ijk.txx"

namespace IJK {

  // **************************************************
  // SUBROUTINES FOR EXTRACTING ISOSURFACE PATCH
  // **************************************************

  /// Return true if isosurface intersects polyhedron
  template <class STYPE, class STYPE2, class VTYPE, class ITYPE, 
	    class NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const VTYPE iv0, const ITYPE * increment, const NTYPE num_poly_vertices)
  // return true if isosurface intersects polyhedron
  // Precondition: num_poly_vertices > 0
  {
    VTYPE iv1 = iv0 + increment[0];
    if (scalar[iv1] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
	iv1 = iv0 + increment[j];
	if (scalar[iv1] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[iv1] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
	iv1 = iv0 + increment[j];
	if (scalar[iv1] < isovalue) { return(true); };
      }
    }

    return(false);
  }

  
  // Newly added
  // Return true if isosurface intersects polyhedron i.e. when isovalue is in the range [mu-delta,mu+delta]
  template <class STYPE, class STYPE2, class VTYPE, class ITYPE, 
	    class NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const VTYPE iv0, const ITYPE * increment, const NTYPE num_poly_vertices, const STYPE * scalar_1)
  // return true if isosurface intersects polyhedron
  // Precondition: num_poly_vertices > 0
  {
   
    // First find the probability of corner value is greater than isovalue
    // and less than isovalue for each of the 8 corners.
     

    double low[8];
    double high[8];
    double lowest = 500;
    double highest = -10;
    double prob_x_greater_than_c = 1;
    double prob_x_less_than_c = 1;
    double afCubeValue[8];
    double vari[8];
    double positive[8];
    double negative[8];


    VTYPE iv1;
    
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        // mean
        afCubeValue[j] = scalar[iv1];
        // variance
        vari[j] = scalar_1[iv1];     
     }      

    for(int i=0; i<8; i++)
    {
        low[i] = afCubeValue[i] - vari[i];
        high[i] = afCubeValue[i] + vari[i];
        if(low[i] < lowest)
              lowest = low[i];
        if(high[i] > highest)
              highest = high[i];
    }

    if((isovalue > lowest) && (isovalue < highest))
    {

   			for(int i=0; i<8; i++)
        {
              if(isovalue > high[i])
                    prob_x_greater_than_c = 0;

              else if(isovalue < low[i])
                    prob_x_greater_than_c = 1;
                    
              else
                    prob_x_greater_than_c =   (double)((high[i] - isovalue)/(high[i]-low[i]));    
  
			 			  positive[i] =  prob_x_greater_than_c;

        }

        for(int i=0; i<8; i++)
        {
             if(isovalue > high[i])
                   prob_x_less_than_c = 1;

             else if(isovalue < low[i])
                   prob_x_less_than_c = 0;
                   
             else
                   prob_x_less_than_c = (double)((isovalue - low[i])/(high[i]-low[i]));

             negative[i] =  prob_x_less_than_c;
         }
               
    }
	
	int configuration[256][8];
  double crossprob[256];
  int conf_ind = 0;
  
		for (int a0 = 0; a0 <= 1; a0++)
		  for (int a1 = 0; a1 <= 1; a1++)
			 	for (int a2 = 0; a2 <= 1; a2++)
			    for (int a3 = 0; a3 <= 1; a3++)
	          for (int a4 = 0; a4 <= 1; a4++)                                                    
	            for (int a5 = 0 ; a5 <= 1; a5++)
        		    for (int a6 = 0 ; a6 <= 1; a6++)
 				 	        for (int a7 = 0 ; a7 <= 1; a7++)
				          {
			            	// save configuration
										// 0 negative , 1 positive
										configuration[conf_ind][0] = a0;						
										configuration[conf_ind][1] = a1;						
										configuration[conf_ind][2] = a2;						
										configuration[conf_ind][3] = a3;						
										configuration[conf_ind][4] = a4;						
										configuration[conf_ind][5] = a5;						
										configuration[conf_ind][6] = a6;						
										configuration[conf_ind][7] = a7;
					
										// save corrensponding probability	
										double temp_prob = 1;
										if (configuration[conf_ind][0] == 0)
											temp_prob = temp_prob*negative[0];

										if (configuration[conf_ind][0] == 1)
                      temp_prob = temp_prob*positive[0];	
                                              
                    if (configuration[conf_ind][1] == 0)
                      temp_prob = temp_prob*negative[1];

                    if (configuration[conf_ind][1] == 1)
                      temp_prob = temp_prob*positive[1];

										if (configuration[conf_ind][2] == 0)
                      temp_prob = temp_prob*negative[2];
                                                
                    if (configuration[conf_ind][2] == 1)
                      temp_prob = temp_prob*positive[2];

										if (configuration[conf_ind][3] == 0)
                      temp_prob = temp_prob*negative[3];
                                                
                    if (configuration[conf_ind][3] == 1)
                      temp_prob = temp_prob*positive[3];

										if (configuration[conf_ind][4] == 0)
                      temp_prob = temp_prob*negative[4];
                                                
                  	if (configuration[conf_ind][4] == 1)
                      temp_prob = temp_prob*positive[4];
                      
										if (configuration[conf_ind][5] == 0)
                      temp_prob = temp_prob*negative[5];
                                                
                    if (configuration[conf_ind][5] == 1)
                      temp_prob = temp_prob*positive[5];
                      
										if (configuration[conf_ind][6] == 0)
                      temp_prob = temp_prob*negative[6];
                                                
                    if (configuration[conf_ind][6] == 1)
                      temp_prob = temp_prob*positive[6];

										if (configuration[conf_ind][7] == 0)
                      temp_prob = temp_prob*negative[7];
                                                
                    if (configuration[conf_ind][7] == 1)
                      temp_prob = temp_prob*positive[7];

										crossprob[conf_ind] = temp_prob;

										conf_ind++;

          			}

 // if crossing probability is greater than 0 for any configuration return true
 
 int intersectflag = 0;
 for(int i=0; i<256 ; i++)
 {
	if(crossprob[i] > 0.001)
	{
		intersectflag = 1; 
		break;
	}	
 }

 if(intersectflag)
	return true;
 else
	return false;	
    
  }


  /// Return true if isosurface intersects polyhedron
  template <class STYPE, class STYPE2, class NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const NTYPE num_poly_vertices)
  // return true if isosurface intersects polyhedron
  // Precondition: num_poly_vertices > 0
  {
    if (scalar[0] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
	if (scalar[j] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[0] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
	if (scalar[j] < isovalue) { return(true); };
      }
    }

    return(false);
  }

  /// Compute isosurface table index of polyhedron with primary vertex iv0
  template <class STYPE, class STYPE2, class VTYPE, class INC_TYPE, 
	    class NTYPE, class ITYPE>
  inline void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const VTYPE iv0, const INC_TYPE * increment,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
	VTYPE iv1 = iv0 + increment[j];
	if (scalar[iv1] >= isovalue) {
	  it = (it | (1L << j));
	};
    };
  }
  


   // newly added
   // Find cell configuration corresponding to maximum cell-crossing probability
   /// Compute isosurface table index of polyhedron with primary vertex iv0
  template <class STYPE, class STYPE2, class VTYPE, class INC_TYPE, 
	    class NTYPE, class ITYPE>
  inline void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const VTYPE iv0, const INC_TYPE * increment,
   const NTYPE num_poly_vertices, ITYPE & it, const STYPE * scalar_1)
  {
  	double low[8];
    double high[8];
    double lowest = 1000000;
    double highest = -1000000;
    double prob_x_greater_than_c = 1;
    double prob_x_less_than_c = 1;
    double afCubeValue[8];
    double vari[8];
    double positive[8];
    double negative[8];


    VTYPE iv1;
    
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        // mean
        afCubeValue[j] = scalar[iv1];
        // variance
        vari[j] = scalar_1[iv1];     
     }      

    for(int i=0; i<8; i++)
    {
        low[i] = afCubeValue[i] - vari[i];
        high[i] = afCubeValue[i] + vari[i];
        if(low[i] < lowest)
              lowest = low[i];
        if(high[i] > highest)
              highest = high[i];
    }

  

   			for(int i=0; i<8; i++)
        {
              if(isovalue > high[i])
                    prob_x_greater_than_c = 0;

              else if(isovalue < low[i])
                    prob_x_greater_than_c = 1;
                    
              else
                    prob_x_greater_than_c =   (double)((high[i] - isovalue)/(high[i]-low[i]));    
  
			 			  positive[i] =  prob_x_greater_than_c;

        }

        for(int i=0; i<8; i++)
        {
             if(isovalue > high[i])
                   prob_x_less_than_c = 1;

             else if(isovalue < low[i])
                   prob_x_less_than_c = 0;
                   
             else
                   prob_x_less_than_c = (double)((isovalue - low[i])/(high[i]-low[i]));

             negative[i] =  prob_x_less_than_c;
         }
               
    
	
	int configuration[256][8];
  double crossprob[256];
  int conf_ind = 0;
  // configuration with maximum probability
  double max_prob = 0;
  // Find and store cell configuration with maximum probability
  double max_prob_conf[8];
 
		for (int a0 = 0; a0 <= 1; a0++)
		  for (int a1 = 0; a1 <= 1; a1++)
			 	for (int a2 = 0; a2 <= 1; a2++)
			    for (int a3 = 0; a3 <= 1; a3++)
	          for (int a4 = 0; a4 <= 1; a4++)                                                    
	            for (int a5 = 0 ; a5 <= 1; a5++)
        		    for (int a6 = 0 ; a6 <= 1; a6++)
 				 	        for (int a7 = 0 ; a7 <= 1; a7++)
				          {
			            	// save configuration
										// 0 negative , 1 positive
										configuration[conf_ind][0] = a0;						
										configuration[conf_ind][1] = a1;						
										configuration[conf_ind][2] = a2;						
										configuration[conf_ind][3] = a3;						
										configuration[conf_ind][4] = a4;						
										configuration[conf_ind][5] = a5;						
										configuration[conf_ind][6] = a6;						
										configuration[conf_ind][7] = a7;
					
										// save corrensponding probability	
										double temp_prob = 1;
										// Major bug
										//max_prob = 0;
										
										if (configuration[conf_ind][0] == 0)
											temp_prob = temp_prob*negative[0];

										if (configuration[conf_ind][0] == 1)
                      temp_prob = temp_prob*positive[0];	
                                              
                    if (configuration[conf_ind][1] == 0)
                      temp_prob = temp_prob*negative[1];

                    if (configuration[conf_ind][1] == 1)
                      temp_prob = temp_prob*positive[1];

										if (configuration[conf_ind][2] == 0)
                      temp_prob = temp_prob*negative[2];
                                                
                    if (configuration[conf_ind][2] == 1)
                      temp_prob = temp_prob*positive[2];

										if (configuration[conf_ind][3] == 0)
                      temp_prob = temp_prob*negative[3];
                                                
                    if (configuration[conf_ind][3] == 1)
                      temp_prob = temp_prob*positive[3];

										if (configuration[conf_ind][4] == 0)
                      temp_prob = temp_prob*negative[4];
                                                
                  	if (configuration[conf_ind][4] == 1)
                      temp_prob = temp_prob*positive[4];
                      
										if (configuration[conf_ind][5] == 0)
                      temp_prob = temp_prob*negative[5];
                                                
                    if (configuration[conf_ind][5] == 1)
                      temp_prob = temp_prob*positive[5];
                      
										if (configuration[conf_ind][6] == 0)
                      temp_prob = temp_prob*negative[6];
                                                
                    if (configuration[conf_ind][6] == 1)
                      temp_prob = temp_prob*positive[6];

										if (configuration[conf_ind][7] == 0)
                      temp_prob = temp_prob*negative[7];
                                                
                    if (configuration[conf_ind][7] == 1)
                      temp_prob = temp_prob*positive[7];

										crossprob[conf_ind] = temp_prob;
										
										if(temp_prob > max_prob)
										{
												max_prob = temp_prob;
												max_prob_conf[0] = a0;						
											  max_prob_conf[1] = a1;						
												max_prob_conf[2] = a2;						
												max_prob_conf[3] = a3;						
												max_prob_conf[4] = a4;						
												max_prob_conf[5] = a5;						
												max_prob_conf[6] = a6;						
												max_prob_conf[7] = a7;											
										}
										
										conf_ind++;

          			}

  						 it = 0;  

    					 for (int j = 0; j < num_poly_vertices; j++) 
     					 {
									//VTYPE iv1 = iv0 + increment[j];
									if (max_prob_conf[j] == 1) 
									{
	  								it = (it | (1L << j));	   
									}	
     					 }
   
  }


  
  /// Compute isosurface table index of polyhedron
  template <class STYPE, class STYPE2, class NTYPE, class ITYPE>
  void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
	if (scalar[j] >= isovalue) {
	  it = (it | (1L << j));
	};
    };
  }

  /// Add isosurface simplex vertices.
  template <class ISOTABLE_TYPE, class INDEX_TYPE, 
	    class VTYPE, class ITYPE, class STYPE, 
	    class VTYPE2, class VTYPE3>
  inline void add_iso_simplex_vertices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv_primary, const STYPE is, 
   const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint)
  {
    for (VTYPE j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
      VTYPE jv = isotable.SimplexVertex(it, is, j);
      VTYPE je = isotable.IsosurfaceVertex(jv).Face();
      VTYPE poly_v0 = isotable.Polyhedron().EdgeEndpoint(je, 0);
      VTYPE poly_v1 = isotable.Polyhedron().EdgeEndpoint(je, 1);

      VTYPE iv0 = iv_primary + increment[poly_v0];
      VTYPE iv1 = iv_primary + increment[poly_v1];
      if (iv0 > iv1) { std::swap(iv0, iv1); };

      iso_simplices.push_back(endpoint.size()/2);
      endpoint.push_back(iv0);
      endpoint.push_back(iv1);
    };

  }

  /// Add isosurface simplices.
  template <class ISOTABLE_TYPE, class INDEX_TYPE, 
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
			       increment, iso_simplices, endpoint);
    };
  }

  /// Add isosurface simplices and reverse orientation.
  template <class ISOTABLE_TYPE, class INDEX_TYPE, 
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void add_reverse_orient_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
			       increment, iso_simplices, endpoint);

      // reverse orientation by swapping last two vertices
      NTYPE ilast = iso_simplices.size()-1;
      std::swap(iso_simplices[ilast], iso_simplices[ilast-1]);
    };
  }

  /// Add isosurface simplices.
  template <class ISOTABLE_TYPE, class INDEX_TYPE, 
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   const bool orientation,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      add_iso_simplices(isotable, it, iv0, increment, iso_simplices,
			endpoint, num_simplices);
    }
    else {
      add_reverse_orient_iso_simplices
	(isotable, it, iv0, increment, iso_simplices,
	 endpoint, num_simplices);
    }
  }

  // **************************************************
  // EXTRACT ISOSURFACE PATCH FROM A POLYHEDRON
  // **************************************************

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param mesh_scalar = Array of scalar values of all mesh vertices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iv0 = Index of first vertex in polyhedron.
  /// @param increment = k'th polyhedron vertex has index iv0+increment[k].
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Numver of simplices added to isosurface.
  template <class STYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void extract_isopatch_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
			num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
			     num_poly_vertices, it);

      add_iso_simplices(isotable, it, iv0, increment, 
			iso_simplices, endpoint, num_simplices);
    }
  }

  // Newly added
  // Mark cubes to be crossed by the isosurface if the isovalue lies in the range [mu-delta,mu+delta]
  // mean (mesh_scalar) and delta (mesh_scalar_1) are taken as input
  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param mesh_scalar = Array of scalar values of all mesh vertices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iv0 = Index of first vertex in polyhedron.
  /// @param increment = k'th polyhedron vertex has index iv0+increment[k].
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Numver of simplices added to isosurface.
  template <class STYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void extract_isopatch_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices, const STYPE * mesh_scalar_1)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;    
    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface i.e. isovalue is in the range [mu-delta,mu+delta]
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
			num_poly_vertices, mesh_scalar_1)) {
			//cout<<"In intersectpoly!!\n";

      // Find the cell configuaration corresponding to maximum cell-crossing probability      		
      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
			     num_poly_vertices, it,mesh_scalar_1);

      add_iso_simplices(isotable, it, iv0, increment, 
			iso_simplices, endpoint, num_simplices);
    }
  }

 


  /// Extract reverse orientation isosurface simplices from polyhedron.
  /// Note: Make this inline for faster execution.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param mesh_scalar = Array of scalar values of all mesh vertices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iv0 = Index of first vertex in polyhedron.
  /// @param increment = k'th polyhedron vertex has index iv0+increment[k].
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <class STYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE, class ITYPE, class VTYPE2, class VTYPE3,
	    class NTYPE>
  inline void extract_isopatch_reverse_orient_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
			num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
			     num_poly_vertices, it);

      add_reverse_orient_iso_simplices
	(isotable, it, iv0, increment, 
	 iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <class STYPE, class VTYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE2, class VTYPE3, class NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_iso_simplices(isotable, it, 0, poly_vertex,
			iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract reverse orientation isosurface simplices from polyhedron.
  /// Note: Make this inline for faster execution.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isovalue = Isovalue.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <class STYPE, class VTYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE2, class VTYPE3, class NTYPE>
  inline void extract_isopatch_reverse_orient_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polyhedron().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_reverse_orient_iso_simplices
	(isotable, it, 0, poly_vertex, iso_simplices, endpoint, num_simplices);
    }
  }

  /// Extract isosurface simplices from polyhedron.
  /// Returns list of isosurface simplex vertices and list of endpoints of grid edges containing simplex vertices.
  /// @param scalar = Array of polyhedron scalar values.
  /// @param poly_vertex = Array of polyhedron vertex indices.
  /// @param isotable = Isosurface lookup table.
  /// @param isovalue = Isovalue.
  /// @param orientation = Orientation of isosurface simplices.
  /// @param iso_simplices = List of vertices of isosurface simplices.
  /// @param endpoint = List of endpoints of edges containing isosurface vertices.
  ///      (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
  /// @param num_simplices = Number of simplices added to isosurface.
  template <class STYPE, class VTYPE, class ISOTABLE_TYPE, class STYPE2,
	    class VTYPE2, class VTYPE3, class NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   const bool orientation, std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      extract_isopatch_from_poly
	(scalar, poly_vertex, isotable, isovalue, 
	 iso_simplices, endpoint, num_simplices);
    }
    else {
      extract_isopatch_reverse_orient_from_poly
	(scalar, poly_vertex, isotable, isovalue, 
	 iso_simplices, endpoint, num_simplices);
    }
  }

}

#endif
