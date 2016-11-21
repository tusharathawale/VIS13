#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"

using namespace std;

class uniform_kernel_polynomial{

	private: double mu1, delta1, mu2, delta2, c;

	public:
		uniform_kernel_polynomial(double m1,double d1,double m2,double d2,double isoval)
		{
			mu1 = m1;
			delta1 = d1;
			mu2 = m2;
			delta2 = d2;
			c = isoval;
		}

		double PQRS_integrate_piece_value(double z, int piecenum);
		double PQRS_expected_piece_value(double z, int piecenum);
		double PQRS_second_moment_piece_value(double z, int piecenum);

		double approx_PQRS_integrate_piece_value(double z, int piecenum);
		double approx_PQRS_expected_piece_value(double z, int piecenum);
		double approx_PQRS_second_moment_piece_value(double z, int piecenum);
};

