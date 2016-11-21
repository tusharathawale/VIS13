#include "uniform_kernel_polynomial.h"
#include <cfloat>
#include "math.h"


// Issues faced : If z value is 0 or 1, the 3 quantities, i.e., pdf, expected and second moment evaluate to NAN because of division by 0.
// Therefore, when z is 0 or 1, z is fluctuated by small amount 1.0e-5. 

// Expressions for expected and second moment contain term log(fabs(-0.1e1 + z)) which gives NAN,  since -0.1e1 + z is a negative quantity (as z is in the range [0,1]) and log of negative numeber in C is
// NAN. however it can be derived log(-z) = log(z) + i*pi. These imaginary parts cancel out. Thus, every term log(fabs(-0.1e1 + z)) is replaced with log(fabs(-0.1e1 + z)).

// there are total 12 basic pieces so piece number will vary from  1 to 12 (These basic pieces vary depending upon the parallelogram)

// very sensitive. can make computations unstable so if some case fails try for 1e-2. For climate ensembel tweak is 1e-1.
double u_tweak = 1.0e-3;
double u_epsilon = 1.0e-4;


// Check later for DBL_EPSILON working/not

double uniform_kernel_polynomial :: PQRS_integrate_piece_value(double z, int piecenum)
{
// Expressions for pdf(z)
	
 if(fabs(z) < u_epsilon)	
	{
		z = u_tweak;
	}
	if(fabs(1-z) < u_epsilon)	
	{
		z = 1 - u_tweak;
	}


	switch (piecenum)
	{
		
		case 1:		
				return -pow(c - mu1 - delta1 - mu2 * z + mu1 * z - delta2 * z + delta1 * z, 0.2e1) / z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;

				
		case 2:
		
				return  pow(-mu2 * z - delta2 * z + c - mu1 + mu1 * z + delta1 - delta1 * z, 0.2e1) / z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;

		case 3:
			
				return -pow(c - mu1 + delta1 - mu2 * z + mu1 * z - delta1 * z + delta2 * z, 0.2e1) / z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;


		case 4:
			
				return pow(-mu2 * z + delta2 * z + c - mu1 + mu1 * z - delta1 + delta1 * z, 0.2e1) / z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;


		case 5:
			
				return -(-delta2 * z + c - mu1 + mu1 * z - mu2 * z) / delta2 / z / 0.2e1;

				
		case 6:
		
				return (-mu2 * z + mu1 * z + delta1 * z + c - mu1 - delta1) / delta1 / (-0.1e1 + z) / 0.2e1;
				
		case 7:
			
				return -pow(-mu2 - delta2 + c, 0.2e1) * z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;

		case 8:
			
				return -pow(-mu2 + delta2 + c, 0.2e1) * z / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;
			
		case 9:
		
				return -(c - mu1 + delta1) * (c * z + mu1 * z - delta1 * z - 0.2e1 * mu2 * z - 0.2e1 * delta2 * z + c - mu1 + delta1) / delta1 / delta2 / z / 0.8e1;

		case 10:
		
				return -(c - mu1 - delta1) * (c * z + mu1 * z + delta1 * z - 0.2e1 * mu2 * z + 0.2e1 * delta2 * z + c - mu1 - delta1) / delta1 / delta2 / z / 0.8e1;
				 
		case 11:
			
				return -(mu2 * mu2 * z + delta2 * delta2 * z - 0.2e1 * mu2 * z * delta2 + 0.4e1 * mu1 * delta2 * z - 0.4e1 * delta1 * delta2 * z - 0.2e1 * c * delta2 * z - 0.2e1 * c * mu2 * z + c * c * z + 0.4e1 * delta2 * delta1 - 0.4e1 * mu1 * delta2 + 0.4e1 * c * delta2) / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;

		case 12:
		
				return -(mu2 * mu2 * z + delta2 * delta2 * z + 0.2e1 * mu2 * z * delta2 - 0.4e1 * mu1 * delta2 * z - 0.4e1 * delta1 * delta2 * z + 0.2e1 * c * delta2 * z - 0.2e1 * c * mu2 * z + c * c * z + 0.4e1 * delta2 * delta1 + 0.4e1 * mu1 * delta2 - 0.4e1 * c * delta2) / delta1 / delta2 / (-0.1e1 + z) / 0.8e1;

		
		default:
			cout<<"This basic piece doesn't exist"; 
	}

}

// there are total 12 basic pieces  so piece number will vary from  1 to 12
double uniform_kernel_polynomial :: PQRS_expected_piece_value(double z, int piecenum)
{

  if(fabs(z) < u_epsilon)	
	{
		z = u_tweak;
	}
	if(fabs(1-z) < u_epsilon)	
	{
		z = 1 - u_tweak;
	}

	switch (piecenum)
	{
		// Expressions for z.pdf(z) and not pdf(z)
		case 1:
			return 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 + 0.1e1 / delta2 * log(z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 - 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 - delta1 / delta2 * log(z) / 0.8e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;

		case 2:
			return -0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 + 0.1e1 / delta2 * log(z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 - 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 + delta1 / delta2 * log(z) / 0.8e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;
			
		case 3:
			return -0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 - 0.1e1 / delta2 * log(z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 + 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 - delta1 / delta2 * log(z) / 0.8e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;
			
		case 4:
			return 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 - 0.1e1 / delta2 * log(z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 + 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 + delta1 / delta2 * log(z) / 0.8e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;
			
		case 5:
			return 0.1e1 / delta2 * log(z) * c / 0.2e1 - 0.1e1 / delta2 * log(z) * mu1 / 0.2e1;
			
		case 6:
			return -0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.2e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.2e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1;
			
		case 7:
			return 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;
			
		case 8:
			return -0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;
			
		case 9:
			return 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 + 0.1e1 / delta2 * log(z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 - 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 + delta1 / delta2 * log(z) / 0.8e1;
			
		case 10:
			return 0.1e1 / delta1 / delta2 * log(z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(z) * c * mu1 / 0.4e1 - 0.1e1 / delta2 * log(z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 * log(z) * mu1 * mu1 / 0.8e1 + 0.1e1 / delta2 * log(z) * mu1 / 0.4e1 + delta1 / delta2 * log(z) / 0.8e1;
			
		case 11:
			return -0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;

		case 12:
			return 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.8e1;

	
		default:
			cout<<"This basic piece doesn't exist"; 
	}
}

// there are total 12 basic pieces  so piece number will vary from  1 to 12
double uniform_kernel_polynomial :: PQRS_second_moment_piece_value(double z, int piecenum)
{

	if(fabs(z) < u_epsilon)	
	{
		z = u_tweak;
	}
	if(fabs(1-z) < u_epsilon)	
	{
		z = 1 - u_tweak;
	}

	
	switch (piecenum)
	{
		// Expressions for z^2.pdf(z) and not pdf(z)
		case 1:
			return -delta1 / delta2 * z / 0.8e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 * mu1 * mu1 * z / 0.8e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 + 0.1e1 / delta1 * z * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 - 0.1e1 / delta2 * mu1 * z / 0.4e1 + 0.1e1 / delta1 / delta2 * c * mu1 * z / 0.4e1 - 0.1e1 / delta1 * c * z / 0.4e1 + 0.1e1 / delta2 * c * z / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1;
			
		case 2:
			return delta1 / delta2 * z / 0.8e1 - 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * mu1 * mu1 * z / 0.8e1 - 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * z * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 - 0.1e1 / delta2 * mu1 * z / 0.4e1 - 0.1e1 / delta1 / delta2 * c * mu1 * z / 0.4e1 + 0.1e1 / delta1 * c * z / 0.4e1 + 0.1e1 / delta2 * c * z / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1;
			
		case 3:
			return -delta1 / delta2 * z / 0.8e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 * mu1 * mu1 * z / 0.8e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * z * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta2 * mu1 * z / 0.4e1 + 0.1e1 / delta1 / delta2 * c * mu1 * z / 0.4e1 + 0.1e1 / delta1 * c * z / 0.4e1 - 0.1e1 / delta2 * c * z / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1;

		case 4:
			return delta1 / delta2 * z / 0.8e1 - 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * mu1 * mu1 * z / 0.8e1 - 0.1e1 / delta1 * delta2 * z / 0.8e1 + 0.1e1 / delta1 * z * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta2 * mu1 * z / 0.4e1 - 0.1e1 / delta1 / delta2 * c * mu1 * z / 0.4e1 - 0.1e1 / delta1 * c * z / 0.4e1 - 0.1e1 / delta2 * c * z / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1;
			
		case 5:
			return (c - mu1) / delta2 * z / 0.2e1;
			
		case 6:
			return 0.1e1 / delta1 * z * mu2 / 0.2e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.2e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 - 0.1e1 / delta1 * c * z / 0.2e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.2e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c;
			
		case 7:
			return -0.1e1 / delta1 * c * z / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 + 0.1e1 / delta1 * z * mu2 / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * c * c * z / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1;
			
		case 8:
			return 0.1e1 / delta1 * c * z / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 - 0.1e1 / delta1 * z * mu2 / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * c * c * z / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1;

		case 9:
			return 0.1e1 / delta1 / delta2 * pow(c - mu1 + delta1, 0.2e1) * z / 0.8e1;

		case 10:
			return 0.1e1 / delta1 / delta2 * pow(c - mu1 - delta1, 0.2e1) * z / 0.8e1;
	
		case 11:
			return 0.1e1 / delta1 * c * z / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 - 0.1e1 / delta1 * z * mu2 / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * c * c * z / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1;
			
		case 12:
			return -0.1e1 / delta1 * c * z / 0.4e1 + 0.1e1 / delta1 / (-0.1e1 + z) * c / 0.4e1 - 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * c / 0.2e1 + 0.1e1 / delta1 / delta2 * z * mu2 * mu2 / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * mu2 * mu2 / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * mu2 * mu2 / 0.4e1 + 0.1e1 / delta1 * delta2 * z / 0.8e1 - 0.1e1 / delta1 * delta2 / (-0.1e1 + z) / 0.8e1 + 0.1e1 / delta1 * delta2 * log(fabs(-0.1e1 + z)) / 0.4e1 + 0.1e1 / delta1 * z * mu2 / 0.4e1 - 0.1e1 / delta1 / (-0.1e1 + z) * mu2 / 0.4e1 + 0.1e1 / delta1 * log(fabs(-0.1e1 + z)) * mu2 / 0.2e1 - 0.1e1 / delta1 / delta2 * c * mu2 * z / 0.4e1 + 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * mu2 / 0.4e1 - 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * mu2 / 0.2e1 + 0.1e1 / delta1 / delta2 * c * c * z / 0.8e1 - 0.1e1 / delta1 / delta2 / (-0.1e1 + z) * c * c / 0.8e1 + 0.1e1 / delta1 / delta2 * log(fabs(-0.1e1 + z)) * c * c / 0.4e1;

		default:
			cout<<"This basic piece doesn't exist"; 
	}
}

