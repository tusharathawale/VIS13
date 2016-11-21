using namespace std;

class z_density_uniform{

	private :
		double numUniformInKde1, numUniformInKde2;
		
	public :

   double alpha_density_uniform(double m1, double d1, double  m2, double  d2, double isovalue, double*  exp, double*  var, double*  crossprob);

	double kde_z_pdf_expected(float* mu1, double h1, float* mu2, double h2, double c);
	
	double kde_z_pdf_expected_NL(float* mu1, float* mu1_wts, double h1, float* mu2, float* mu2_wts, double h2, double c);

	double kde_z_pdf_variance(float* mu1, double h1, float* mu2, double h2, double c, double expected_crossing);
	
	double kde_z_pdf_variance_NL(float* mu1, float* mu1_wts, double h1, float* mu2, float* mu2_wts, double h2, double c, double expected_crossing);
	
	// set number of uniform distributions in kde_1 
	void setNumUniformInKde1(int a);

	// get number of uniform distributions in kde_1 
	int getNumUniformInKde1();

	// set number of uniform distributions in kde_2 
	void setNumUniformInKde2(int a);

	// get number of uniform distributions in kde_2
	int getNumUniformInKde2();
		
};

