#include <math.h>
#include <limits>
#include <vector>
#include <algorithm>

#define PI 3.14159265358979323846
#define PIX1 (PI * 1.0)
#define PIX2 (PI * 2.0)
#define PIX4 (PI * 4.0)

#define n_bands 5
#define n_coeff n_bands * n_bands
#define n_samples_squared 100
#define n_samples n_samples_squared * n_samples_squared

using std::vector;
typedef unsigned long long LargeInt;
typedef double(*SH_polar_fn)(double theta, double phi);

struct Vector3d {
	double x, y, z;
};

struct SHSample {
	Vector3d sph;
	Vector3d vec;
	double coeff[n_coeff] = { 0 };
};

double Random() {
	int a = rand();
	int b = std::numeric_limits<int>::max();
	int c = RAND_MAX;
	double r = (double)a / (double)c;
	return r;
}

LargeInt Factorial(LargeInt n) {
	static LargeInt factorialTable[16] = { 0 };

	if (n <= 1) {
		factorialTable[n] = 1;
		return factorialTable[n];
	} else if (factorialTable[n] != 0) {
		return factorialTable[n];
	} else {
		LargeInt f = n * Factorial(n - 1);
		factorialTable[n] = f;
		return f;
	}
}

double P(int l, int m, double x) {
	// evaluate an Associated Legendre Polynomial P(l,m,x) at x
	double pmm = 1.0;
	if (m>0) {
		double somx2 = sqrt((1.0 - x)*(1.0 + x));
		double fact = 1.0;
		for (int i = 1; i <= m; i++) {
			pmm *= (-fact) * somx2;
			fact += 2.0;
		}
	}
	if (l == m) return pmm;
	double pmmp1 = x * (2.0*m + 1.0) * pmm;
	if (l == m + 1) return pmmp1;
	double pll = 0.0;
	for (int ll = m + 2; ll <= l; ++ll) {
		pll = ((2.0*ll - 1.0)*x*pmmp1 - (ll + m - 1.0)*pmm) / (ll - m);
		pmm = pmmp1;
		pmmp1 = pll;
	}
	return pll;
}

double K(int l, int m) {
	// renormalisation constant for SH function
	double temp = ((2.0*l + 1.0)*Factorial(l - m)) / (4.0*PI*Factorial(l + m));
	return sqrt(temp);
}

double SH(int l, int m, double theta, double phi) {
	// return a point sample of a Spherical Harmonic basis function
	// l is the band, range [0..N]
	// m in the range [-l..l]
	// theta in the range [0..Pi]
	// phi in the range [0..2*Pi]
	const double sqrt2 = sqrt(2.0);
	if (m == 0) return K(l, 0)*P(l, m, cos(theta));
	else if (m>0) return sqrt2*K(l, m)*cos(m*phi)*P(l, m, cos(theta));
	else return sqrt2*K(l, -m)*sin(-m*phi)*P(l, -m, cos(theta));
}

void SH_setup_spherical_samples(SHSample* samples, int sqrt_n_samples) {
	// fill an N*N*2 array with uniformly distributed
	// samples across the sphere using jittered stratification
	int i = 0; // array index
	double oneoverN = 1.0 / sqrt_n_samples;
	for (int a = 0; a<sqrt_n_samples; a++) {
		for (int b = 0; b<sqrt_n_samples; b++) {
			// generate unbiased distribution of spherical coords
			double x = (a + Random()) * oneoverN; // do not reuse results
			double y = (b + Random()) * oneoverN; // each sample must be random
			double theta = 2.0 * acos(sqrt(1.0 - x));
			double phi = 2.0 * PI * y;
			samples[i].sph = Vector3d{ theta, phi, 1.0 };
			// convert spherical coords to unit vector
			Vector3d vec = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) };
			samples[i].vec = vec;
			// precompute all SH coefficients for this sample
			for (int l = 0; l<n_bands; ++l) {
				for (int m = -l; m <= l; ++m) {
					int index = l*(l + 1) + m;
					double c = SH(l, m, theta, phi);
					samples[i].coeff[index] = c;
				}
			}
			++i;
		}
	}
}

double Sample_two_lights(double theta, double phi) {
	double a = std::max(0.0,  5.0 * cos(theta) - 4.0);
	double b = std::max(0.0, -4.0 * sin(theta - PI) * cos(phi - 2.5) - 3.0);
	return a + b;
}

void SH_project_polar_function(SH_polar_fn fn, const vector<SHSample>& samples, vector<double>& result) {
	const double weight = 4.0*PI;
	// for each sample
	for (int i = 0; i<n_samples; ++i) {
		double theta = samples[i].sph.x;
		double phi = samples[i].sph.y;
		for (int n = 0; n<n_coeff; ++n) {
			result[n] += fn(theta, phi) * samples[i].coeff[n];
		}
	}
	// divide the result by weight and number of samples
	double factor = weight / double(n_samples);
	for (int i = 0; i<n_coeff; ++i) {
		result[i] = result[i] * factor;
	}
}

// From: https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
double SH_project(double theta, double phi, vector<double>& c) {

	double x = sin(theta) * cos(phi);
	double y = sin(theta) * sin(phi);
	double z = cos(theta);

	double Y00 = (1.0 / 2.0)  * sqrt(  1.0 / PIX1);

	double Y10 =                sqrt(  3.0 / PIX4) * y;
	double Y11 =                sqrt(  3.0 / PIX4) * z;
	double Y12 =                sqrt(  3.0 / PIX4) * x;

	double Y20 = (1.0 / 2.0)  * sqrt( 15.0 / PIX1) * (x*y);
	double Y21 = (1.0 / 2.0)  * sqrt( 15.0 / PIX1) * (y*z);
	double Y22 = (1.0 / 4.0)  * sqrt(  5.0 / PIX1) * (-x*x - y*y + 2.0*z*z);
	double Y23 = (1.0 / 2.0)  * sqrt( 15.0 / PIX1) * (x*z);
	double Y24 = (1.0 / 4.0)  * sqrt( 15.0 / PIX1) * (x*x - y*y);

	double Y30 = (1.0 / 4.0)  * sqrt( 35.0 / PIX2) * (y*(3.0*x*x - y*y));
	double Y31 = (1.0 / 2.0)  * sqrt(105.0 / PIX1) * (x*y*z);
	double Y32 = (1.0 / 4.0)  * sqrt( 21.0 / PIX2) * (y*(4.0*z*z - x*x - y*y));
	double Y33 = (1.0 / 4.0)  * sqrt(  7.0 / PIX1) * (z*(2.0*z*z - 3.0*x*x - 3.0*y*y));
	double Y34 = (1.0 / 4.0)  * sqrt( 21.0 / PIX2) * (x*(4.0*z*z - x*x - y*y));
	double Y35 = (1.0 / 4.0)  * sqrt(105.0 / PIX1) * (z*(x*x - y*y));
	double Y36 = (1.0 / 4.0)  * sqrt( 35.0 / PIX2) * (x*(x*x - 3.0*y*y));

	double Y40 = (3.0 / 4.0)  * sqrt( 35.0 / PIX1) * (x*y*(x*x - y*y));
	double Y41 = (3.0 / 4.0)  * sqrt( 35.0 / PIX2) * (y*z*(3.0*x*x - y*y));
	double Y42 = (3.0 / 4.0)  * sqrt(  5.0 / PIX1) * (x*y*(7.0*z*z - 1.0));
	double Y43 = (3.0 / 4.0)  * sqrt(  5.0 / PIX2) * (y*z*(7.0*z*z - 3.0));
	double Y44 = (3.0 / 16.0) * sqrt(  1.0 / PIX1) * (35.0*z*z*z*z - 30.0*z*z + 3.0);
	double Y45 = (3.0 / 4.0)  * sqrt(  5.0 / PIX2) * (x*z*(7.0*z*z - 3.0));
	double Y46 = (3.0 / 8.0)  * sqrt(  5.0 / PIX1) * (x*x - y*y)*(7.0*z*z - 1.0);
	double Y47 = (3.0 / 4.0)  * sqrt( 35.0 / PIX2) * (x*z*(x*x - 3.0*y*y));
	double Y48 = (3.0 / 16.0) * sqrt( 35.0 / PIX1) * (x*x*(x*x - 3.0*y*y) - y*y*(3.0*x*x - y*y));

	double SH =
		Y00*c[0]  +
		Y10*c[1]  + Y11*c[2]  + Y12*c[3]  +
		Y20*c[4]  + Y21*c[5]  + Y22*c[6]  + Y23*c[7]  + Y24*c[8]  +
		Y30*c[9]  + Y31*c[10] + Y32*c[11] + Y33*c[12] + Y34*c[13] + Y35*c[14] + Y36*c[15] +
		Y40*c[16] + Y41*c[17] + Y42*c[18] + Y43*c[19] + Y44*c[20] + Y45*c[21] + Y46*c[22] + Y47*c[23] + Y48*c[24];

	return SH;
}

int main() {
	vector<SHSample> samples(n_samples);
	SH_setup_spherical_samples(samples.data(), n_samples_squared);

	vector<double> coefficients(n_coeff);
	SH_project_polar_function(Sample_two_lights, samples, coefficients);

	double theta = 0.0;
	double phi = 0.0;
	double valueFromFunction = Sample_two_lights(theta, phi);
	double valueFromSHProjection = SH_project(theta, phi, coefficients);

	// Print GL coeficients for gl shader
	for (int i = 0; i < (int)coefficients.size(); ++i) {
		printf("float c%i = %f;\n", i, coefficients[i]);
	}
}