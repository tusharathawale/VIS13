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

# define infinity 10000
# define minus_infinity -10000

typedef struct inverse_quadratic
{  
	// Each piece of form (k1*z^2 + k2*(1-z)^2)/(k3*z^2*(1-z)^2), 
	double k1;
	double k2;
	double k3;
}piece;

typedef struct fn
{         
	// number of pieces in a piecewise function
	int numPieces;
	// maximum number of pieces supported(10)
	piece pc[10];
	// store limits  
	double limits[12]; 

}piecewise;

using namespace std;

class alpha_distribution{

	private :
		double slope_OP, slope_OQ, slope_OR, slope_OS, f1, f2, f3, f4, denominator, mean_temp, delta_temp;

	public :
		// Provide piece type(from 0,1,2), actual piece coefficients calculated using alpha_pdf, and limits a(lower) and b(higher) of the integration
		double getSecondMoment(double k1, double k2, double k3, double a, double b);

		// Provide piece type(from 0,1,2), actual piece coefficients calculated using alpha_pdf, and limits a(lower) and b(higher) of the integration
		double getEdgeCrossingProbability(double k1, double k2, double k3, double a, double b);

		// Provide piece type(from 0,1,2), actual piece coefficients calculated using alpha_pdf, and limits a(lower) and b(higher) of the integration
		double getExpectedValOverSinglePiece(double k1, double k2, double k3, double a, double b);

		// p is pdf over -infinity to infinity range. Function returns part of pdf just over range 0 and 1.
		piecewise getPdfOver0To1(piecewise p);

		// Computes expected value, crossing probability and variance over [0,1]
		void Compute0To1(piecewise p, double* expt, double* cross_prob, double* var);

		// use k1 only if piece is of type 1 or 2
		piecewise setPiece(piecewise P, int pieceIndex, double k1, double k2);

		// Piecewise function returned assuming data is sampled from a kernel density estimation
		piecewise kde_alpha_pdf(double* mu1, double* delta1, double* mu2, double* delta2, double c);

		// 1) Non-overlapping Intervals
		piecewise nonOverlapping(double mu1, double delta1, double mu2, double delta2, double c);

		// 2) Overlapping Intervals
		piecewise overlapping(double mu1, double delta1, double mu2, double delta2, double c);

		// 3a) Contained Intervals
		piecewise containedA(double mu1, double delta1, double mu2, double delta2, double c); 

		// 3b) Contained Intervals
		piecewise containedB(double mu1, double delta1, double mu2, double delta2, double c);

		// Piecewise function returned assuming data is sampled from uniform distribution
		piecewise alpha_pdf(double mu1, double delta1, double mu2, double delta2, double c);

};

