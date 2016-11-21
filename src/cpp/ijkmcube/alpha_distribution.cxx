//---------------------------------------------------------------------------------------
/*
	Title : Uncertainty Quantification in Linear Interpolation for Isosurface Extraction	
       Authors: Tushar Athawale, Alireza Entezari
        Date  : Jun 27, 2013.
*/
//---------------------------------------------------------------------------------------

#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "alpha_distribution.h"

# define infinity 10000
# define minus_infinity -10000
#define _USE_MATH_DEFINES

using namespace std;

// p is pdf over -infinity to infinity range. Function returns part of pdf just over range 0 and 1.
piecewise alpha_distribution :: getPdfOver0To1(piecewise p)
{

	piecewise p01;

	int  k=0;       

	// obtain number of pieces in the piecewise function 
	int piece_count = p.numPieces;

	// a1<p1<b1<p2<a2  .. so if 2 pieces, then 3 limits
	int limit_count = piece_count + 1;       

	// Determine limits that fall in [0,1]

	for (int i=0; i<(limit_count-1); i++) 
	{
		//if both limits are below 0 or above 1 don't do anything. skip to next limit
     
		// if any limits of piece lie in 0-1 range   
		if (((p.limits[i] >= 0) && (p.limits[i] <= 1)) || ((p.limits[i+1] >= 0) && (p.limits[i+1] <= 1)))
		{
			p01.limits[k] = p.limits[i];
			p01.limits[k+1] = p.limits[i+1];
			p01.pc[k] = p.pc[i];
			k++;
		}

		// if a piece limits completely overlap [0,1] such that 1 limit is less than 0 and other is greater than 1  
		if ((p.limits[i] < 0) && (p.limits[i+1] > 1))
		{
			p01.limits[k] = 0;
			p01.limits[k+1] = 1;
			p01.pc[k] = p.pc[i];
			k++; 
		}	  

	} 

	// if limit is less than 0 then set it to 0. That can happen for only 1st entry of the lm_in_0_1[10][2] array
	// if limit is greater than 1 then set it to 1. That can happen for only last entry of the lm_in_0_1[10][2] array

	if (p01.limits[0] < 0)
		p01.limits[0] = 0;
	if (p01.limits[k] > 1)
		p01.limits[k] = 1;

	p01.numPieces = k;

	return p01;    

}

double alpha_distribution :: getExpectedValOverSinglePiece(double k1, double k2, double k3, double a, double b)
{
	double expectation = 0;
	double epsilon = 1e-10;

	double c1 = (double)(k2/k3)*log((double)((b)/(a))); 
	double c2 = (double)(k1/k3)*log((double)((b-1)/(a-1))); 
	
	if(fabs(k1)<epsilon && fabs(k2)<epsilon)
		expectation = 0;

	else if(fabs(k1)<epsilon && fabs(k2)>epsilon)
		expectation = c1;

	else if(fabs(k1)>epsilon && fabs(k2)<epsilon)
		expectation = c2 - (double)(k1/k3)*((double)(1 / (b - 1))) + (double)(k1/k3)*((double)(1 / (a - 1)));
	
	else 
		expectation = c1 + c2 - (double)(k1/k3)*((double)(1 / (b - 1))) + (double)(k1/k3)*((double)(1 / (a - 1)));

	return expectation;
}

double alpha_distribution :: getEdgeCrossingProbability(double k1, double k2, double k3, double a, double b)
{
	double crossingProb = 0;
	double epsilon = 1e-10; 
	if((fabs(k1)<epsilon) && (fabs(k2)<epsilon))
		crossingProb = 0;
	else if((fabs(k1)<epsilon) && (fabs(k2)>epsilon))
	{
		crossingProb = (double)(k2/k3)*((double)(1/a)-(double)(1/b));
	}
	else if((fabs(k1)>epsilon) && (fabs(k2)<epsilon))
	{
		crossingProb = (double)(k1/k3)*((double)(1/(1-b))-(double)(1/(1-a)));	
	}
	else
	{ 			
		crossingProb = (double)(1/k3)*(k2*((double)(1/a) - (double)(1/b)) + k1*((double)(1/(a-1))- (double)(1/(b-1))));
	}
	return crossingProb;
}

double alpha_distribution :: getSecondMoment(double k1, double k2, double k3, double a, double b)
{
	double secondMoment = 0;
	double epsilon = 1e-10; 
	
	if(fabs(k1)<epsilon && fabs(k2)<epsilon)
		secondMoment = 0;

	else if(fabs(k1)<epsilon && fabs(k2)>epsilon)
		secondMoment = (double)(k2/k3)*(b-a);

	else if(fabs(k1)>epsilon && fabs(k2)<epsilon)
		secondMoment = (double)(k1/k3)*(b-a) + (double)(k1/k3)*((double)(1/(a-1))-(double)(1/(b-1))) + (double)((2*k1)/k3)*log((double)((b-1)/(a-1)));
	else 
		secondMoment = ((double)(k1/k3) + (double)(k2/k3))*(b-a) + (double)(k1/k3)*((double)(1/(a-1))-(double)(1/(b-1))) + (double)((2*k1)/k3)*log((double)((b-1)/(a-1)));

	return secondMoment;

}

void alpha_distribution :: Compute0To1(piecewise p, double* expt, double* cross_prob, double* var)
{

	piecewise pdf01 = getPdfOver0To1(p);

	// compute expectation (pdf which is unnormalized over range 0 to 1)
	double expectation_over_0_1 = 0;	
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		expectation_over_0_1 += getExpectedValOverSinglePiece(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	// compute crossing probability 
	double crossProb_0_1 = 0;	
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		crossProb_0_1 += getEdgeCrossingProbability(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	// compute variance 
	double sec_0_1 = 0;
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		sec_0_1 += getSecondMoment(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	*expt = expectation_over_0_1;
	*cross_prob = crossProb_0_1;
	*var = (double)(((sec_0_1*crossProb_0_1) - (double)(expectation_over_0_1*expectation_over_0_1))/(crossProb_0_1*crossProb_0_1));
}

// Create piece by setting piece coefficients
// Set k1, k2, k3, k4 for piece type 0, Set k1 if piece is of type 1 or 2
piecewise alpha_distribution :: setPiece(piecewise P, int pieceIndex, double k1, double k2)
{                
	P.pc[pieceIndex].k1 = k1; 
	P.pc[pieceIndex].k2 = k2;
	P.pc[pieceIndex].k3 = denominator;
	return P; 
}

// 1) nonOverlapping Intervals
piecewise alpha_distribution :: nonOverlapping(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P; 
	int pieceID;

	// mu1-delta1 <= c <= mu1+delta1
	if ((c >= (mu1-delta1)) && (c <= (mu1+delta1))) 
	{
		// Vertex Order : SPQR		
		
		// specify number of pieces
		P.numPieces = 3;
		pieceID = 0;

		// set pieces

		// Set piece 0 coefficients                    
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;           	

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
		pieceID++;          	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, -f4, f2);    
		pieceID++;                 	

		// set limits
		P.limits[0] = slope_OS;  
		P.limits[1] = slope_OP;
		P.limits[2] = slope_OQ; 
		P.limits[3] = slope_OR;

	}    

	// mu1+delta1 < c <= m2-delta2           
	else if ((c > (mu1+delta1)) && (c <= (mu2-delta2))) 
	{
		
		// Vertex Order : PQSR		

		if(slope_OS > slope_OQ)
		{
			// specify number of pieces
			P.numPieces = 3;
			pieceID = 0; 


			// set pieces

			// Set piece 0 coefficients                    
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;           	

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));   
			pieceID++;          	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, -f4, f2);    
			pieceID++;                 	

			// set limits
			P.limits[0] = slope_OP;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OR;

		}

		else if (slope_OQ >= slope_OS)
		{			
			// Vertex Order : PSQR		

			// specify number of pieces
			P.numPieces = 3;
			pieceID = 0; 


			// set pieces

			// Set piece 0 coefficients                    
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;           	

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
			pieceID++;          	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, -f4, f2);    
			pieceID++;               

			// set limits
			P.limits[0] = slope_OP;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
		}                
	}

	// mu2-delta2 < c <= mu2+delta2
	else if ((c > (mu2-delta2)) && (c <= (mu2+delta2)))
	{				
				// Vertex Order : PQRS		

				// specify number of pieces
				P.numPieces = 3;
				pieceID = 0; 

				// set pieces

				// Set piece 0 coefficients                    
				P = setPiece(P, pieceID, f3, -f1);        
				pieceID++;           	

				// Set piece 1 coefficients
				P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));   
				pieceID++;          	

				// Set piece 2 coefficients
				P = setPiece(P, pieceID, f4, -f1);    
				pieceID++;       

				// set limits
				P.limits[0] = slope_OP;  
				P.limits[1] = slope_OQ;
				P.limits[2] = slope_OR; 
				P.limits[3] = slope_OS;
  
	}

	return P; 
}


// 2) Overlapping Intervals
piecewise alpha_distribution :: overlapping(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu1-delta1 <= c <= mu2-delta2
	if ((c >= mu1-delta1) && (c <= mu2-delta2))
	{

		// Vertex Order : PQRS		

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 

		// set pieces

		// Set piece 0 coefficients                        
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, -f4, f2);        
		pieceID++;    					      	

		// Set piece 3 coefficients                        
		P = setPiece(P, pieceID, 0, 0);   
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OP;
		P.limits[2] = slope_OQ; 
		P.limits[3] = slope_OR;
		P.limits[4] = slope_OS;
		P.limits[5] = infinity;
	}     

	// mu2-delta2 < c <= mu1+delta1 
	else if ((c > mu2-delta2) && (c <= mu1+delta1))
	{
		if (slope_OS < slope_OQ)
		{
			// Vertex Order : PSQR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);   
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f3, f1);        
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;       

			// set limits

			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}

		else if (slope_OQ <= slope_OS) 
		{

			// Vertex Order : PQSR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);   
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f4, f2);        
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}		
	}

	// case mu1+delta1 < c <= mu2+delta2
	else if ((c > mu1+delta1) && (c <= mu2+delta2))
	{
			// Vertex Order : SPQR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, -f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, 0, 0);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));                  
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, -f1);        
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;     
	}	

	return P;	
}

// 3a) Contained Intervals
piecewise alpha_distribution :: containedA(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu1-delta1 <= c <= mu2-delta2 
	if ((c >= (mu1-delta1)) && (c <= (mu2-delta2))) 
	{

		// Vertex Order : QRSP		

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 


		// set pieces

		// Set piece 0 coefficients       
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)), 0);                  
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, -f4, f2);  
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, 0, 0);        
		pieceID++;    					      	

		// Set piece 3 coefficients          
		P = setPiece(P, pieceID, -f4, f1);                         
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)), 0);   
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OQ;
		P.limits[2] = slope_OR; 
		P.limits[3] = slope_OS;
		P.limits[4] = slope_OP;
		P.limits[5] = infinity;
	}	   

	// mu2-delta2 < c <= mu2+delta2 
	else if ((c > (mu2-delta2)) && (c <= (mu2+delta2)))
	{
		if ((slope_OS <= slope_OQ) && (slope_OR <= slope_OP))
		{

			// Vertex Order : SQRP	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OP;
			P.limits[5] = infinity;
		}

		else if ((slope_OS >= slope_OQ) && (slope_OR <= slope_OP))
		{
	
			// Vertex Order : QSRP
	
			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OP;
			P.limits[5] = infinity;
		}

		else if ((slope_OS <= slope_OQ) && (slope_OR >= slope_OP))
		{

			// Vertex Order : SQPR	
	
			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}

		else if ((slope_OS >= slope_OQ) && (slope_OR >= slope_OP)) 
		{
		
			// Vertex Order : QSPR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     	


			// set limits

			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}	
	}	

	// mu2+delta2 < c <= mu1+delta1 
	else if ((c > (mu2+delta2)) && (c <= (mu1+delta1)))
	{

		// Vertex Order : SPQR

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 

		// set pieces

		// Set piece 0 coefficients       
		P = setPiece(P, pieceID, (double)(4*delta2*(c-mu2)), 0);                  
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, -f3, f1);  
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, 0, 0);        
		pieceID++;    					      	

		// Set piece 3 coefficients          
		P = setPiece(P, pieceID, -f3, f2);                         
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(c-mu2)), 0);   
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OS;
		P.limits[2] = slope_OP; 
		P.limits[3] = slope_OQ;
		P.limits[4] = slope_OR;
		P.limits[5] = infinity;

	}   
	return P;	
}

// 3b) Contained Intervals
piecewise alpha_distribution :: containedB(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu2-delta2 <= c <= mu1-delta1
	if ((c >= (mu2-delta2)) && (c <= (mu1-delta1)))
	{		
			// Vertex Order : PQRS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients       
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(mu1-c)));                  
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, -f2);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0, 0);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, f4, -f2);                         
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(mu1-c)));
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;	
	}      

	// mu1-delta1 < c <= mu1+delta1
	else if ((c > (mu1-delta1)) && (c <= (mu1+delta1)))
	{
		if ((slope_OP >= slope_OR) && (slope_OQ >= slope_OS))
		{

			// Vertex Order : RPSQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
		}

		else if ((slope_OP <= slope_OR) && (slope_OQ >= slope_OS))
		{

			// Vertex Order : PRSQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     
	
			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OR; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
		}

		else if ((slope_OP >= slope_OR) && (slope_OQ <= slope_OS)) 
		{

			// Vertex Order : RPQS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;
		}

		else if ((slope_OP <= slope_OR) && (slope_OQ <= slope_OS)) 
		{

			// Vertex Order : PRQS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     
	
			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OR; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;     						                              	
		}
	}       

	// mu1+delta1 < c <= mu2+delta2
	else if ((c > (mu1+delta1)) && (c <= mu2+delta2))
	{	
			// Vertex Order : RSPQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients       
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));                  
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, -f1);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0, 0);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, f3, -f1);                         
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
	}	 	   
	return P;	
}

// Get Distribution of Z = (c-X1)/(X2-X1) assuming endpoints have independent uniform noise
piecewise alpha_distribution::alpha_pdf(double mu1, double delta1, double mu2, double delta2, double c)
{
	piecewise P;         

	// Assume mu2 is always greater than mu1
	if (mu1 > mu2)
	{
		mean_temp = mu1;
		mu1 = mu2;
		mu2 = mean_temp;

		delta_temp = delta1;
		delta1 = delta2;
		delta2 = delta_temp;
	}   

	// Precomputation of slopes and variables which are part of final PDF    
	if ((mu2 - mu1 + delta2 - delta1) != 0)
		slope_OP = (double)((c - mu1 - delta1) / (mu2 - mu1 + delta2 - delta1));

	if ((mu2 - mu1 + delta2 + delta1) != 0) 
		slope_OQ = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta2 + delta1)); 

	if ((mu2 - mu1 - delta2 - delta1) != 0)
		slope_OS = (double)((c - mu1 - delta1) / (mu2 - mu1 - delta2 - delta1));		

	if ((mu2 - mu1 + delta1 - delta2) != 0)
		slope_OR = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta1 - delta2));     	


	// Pre-define constants
	
	f1 = (mu1 + delta1 -c)*(mu1 + delta1 -c);
	f2 = (mu1 - delta1 -c)*(mu1 - delta1 -c);
	f3 = (mu2 + delta2 -c)*(mu2 + delta2 -c);
	f4 = (mu2 - delta2 -c)*(mu2 - delta2 -c);
	denominator = 8*delta1*delta2;


	// 1) Non-overlapping Intervals
	if ((mu2 - delta2) >= (mu1 + delta1))
	{            
		P = nonOverlapping(mu1, delta1, mu2, delta2, c); 
	} 

	// 2) Overlapping Intervals
	else if (((mu2-delta2) <= (mu1+delta1)) && ((mu2-delta2) > (mu1-delta1)) && ((mu2+delta2) > (mu1+delta1)))	
	{
		P = overlapping(mu1, delta1, mu2, delta2, c);          
	}	

	// 3a) Contained Intervals
	else if (((mu2-delta2) < (mu1+delta1)) && ((mu2-delta2) >= (mu1-delta1)) && ((mu2+delta2) <= (mu1+delta1)) && ((mu2+delta2) > (mu1-delta1))) 	     	
	{
		if((mu2+delta2) == (mu1+delta1))
		{
			if(c < mu1+delta1)
				slope_OP = infinity;					
		}		 
		P = containedA(mu1, delta1, mu2, delta2, c);          
	}	

	// 3b) Contained Intervals
	else if ((mu2-delta2) <= (mu1-delta1))
	{
		if((mu2-delta2) == (mu1-delta1))
		{
			slope_OR = -infinity;					
		}		

		P = containedB(mu1, delta1, mu2, delta2, c);            
	}		 

	return P; 	
} 

/*
int main()
{
	//mu1: 2.85714 delta1: 57.1429 mu2: 11.4286 delta2: 114.286
	
	// ratio distribution assuming data is uniformly distributed
	alpha_distribution* d = new alpha_distribution();

	// ratio distribution assuming data is sampled from kde
	alpha_distribution* kde_d = new alpha_distribution();*/
	
/*	// Allocate memory for mu and delta arrays
	double* mu1 = new double[2];
	double* mu2 = new double[2];
	double* delta1 = new double[2];
	double* delta2 = new double[2];

	// kernel 1
	mu1[0] = 5;
	delta1[0] = 3;
	mu1[1] = 9;
	delta1[1] = 1;

	// kernel 2
	mu2[0] = 23;
	delta2[0] = 3;
	mu2[1] = 27;
	delta2[1] = 1;*/

	/*piecewise p1, p2;

	// Assuming data is uniformly distributed
	p1 = d->alpha_pdf(8,7, 9, 4, 12.3);
	p2 = d->getPdfOver0To1(p1);

	double expected, crossProb, var;
	d->Compute0To1(p2, &expected, &crossProb, &var);
	cout<<"Expected Value is:"<<expected<<"\n";
	cout<<"Crossing Probability is:"<<crossProb<<"\n";
	cout<<"Variance is:"<<var<<"\n";

	//p1 = d->alpha_pdf(mu1[0], delta1[0], mu2[0], delta2[0], 10);*/
	/*cout<<p1.numPieces<<"\n";
	for(int i=0;i<=p1.numPieces;i++)
		cout<<p1.limits[i]<<" ";
	cout<<"\n";

	cout<<p2.numPieces<<"\n";
	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.limits[i]<<" ";
	cout<<"\n";*/

	/*// Assuming data is kde sampled
	p2 = kde_d->kde_alpha_pdf(mu1, delta1, mu2, delta2, 10);
	cout<<p2.numPieces<<"\n";
	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.limits[i]<<" ";
	cout<<"\n";

	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.pc[i].type<<" ";
	cout<<"\n ";*/

//	return 0;
//}
