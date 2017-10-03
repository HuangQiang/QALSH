#ifndef __RANDOM_H
#define __RANDOM_H


// -----------------------------------------------------------------------------
//  Functions used for generating random variables (r.v.).
// -----------------------------------------------------------------------------
float uniform(						// r.v. from Uniform(min, max)
	float min,							// min value
	float max);							// max value

// -----------------------------------------------------------------------------
float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mu,							// mean (location)
	float sigma);						// stanard deviation (scale > 0)

// -----------------------------------------------------------------------------
float cauchy(						// r.v. from Cauchy(gamma, delta)
	float gamma,						// scale factor (gamma > 0)
	float delta);						// location

// -----------------------------------------------------------------------------
float levy(							// r.v. from Levy(gamma, delta)
	float gamma,						// scale factor (gamma > 0)
	float delta);						// location

// -----------------------------------------------------------------------------
float p_stable(						// r.v. from p-satble distr.
	float p,							// p value, where p in (0,2]
	float zeta,							// symmetric factor (zeta in [-1, 1])
	float gamma,						// scale factor (gamma > 0)
	float delta);						// location


// -----------------------------------------------------------------------------
//  Functions used for calculating probability distribution function (pdf) and 
//  cumulative distribution function (cdf).
// -----------------------------------------------------------------------------
float gaussian_pdf(					// pdf of N(0, 1)
	float x);							// variable

// -----------------------------------------------------------------------------
float gaussian_cdf(					// cdf of N(0, 1) in range (-inf, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float new_gaussian_cdf(				// cdf of N(0, 1) in range [-x, x]
	float x,							// integral border (x > 0)
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float levy_pdf(						// pdf of Levy(1, 0)
	float x);							// variable

// -----------------------------------------------------------------------------
float levy_cdf(						// cdf of Levy(0, 1) in range (0, x]
	float x,							// integral border (x > 0)
	float step = 0.001f);				// step increment


// -----------------------------------------------------------------------------
//  query-oblivious and query-aware collision probability under gaussian
//  distribution, cauchy distribution and levy distribution
// -----------------------------------------------------------------------------
float orig_gaussian_prob(			// calc original gaussian probability
	float x);							// x = w / r

// -----------------------------------------------------------------------------
float new_gaussian_prob(			// calc new gaussian probability
	float x);							// x = w / (2 * r)

// -----------------------------------------------------------------------------
float orig_cauchy_prob(				// calc original cauchy probability
	float x);							// x = w / r

// -----------------------------------------------------------------------------
float new_cauchy_prob(				// calc new cauchy probability
	float x);							// x = w / (2 * r)

// -----------------------------------------------------------------------------
float orig_levy_prob(				// calc original levy probability
	float x);							// x = w / r

// -----------------------------------------------------------------------------
float new_levy_prob(				// calc new levy probability
	float x);							// x = w / (2 * r)

// -----------------------------------------------------------------------------
void orig_stable_prob(				// calc orig stable probability
	float  p,							// the p value, where p in (0, 2]
	float  zeta,						// symmetric factor (zeta in [-1, 1])
	float  ratio,						// approximation ratio
	float  radius,						// radius
	float  w,							// bucket width
	int    num,							// number of repetition
	float& p1,							// p1 = p(w / r), returned
	float& p2);							// p2 = p(w / (c * r)), returned

// -----------------------------------------------------------------------------
void new_stable_prob(				// calc new stable probability
	float  p,							// the p value, where p in (0, 2]
	float  zeta,						// symmetric factor (zeta in [-1, 1])
	float  ratio,						// approximation ratio
	float  radius,						// radius
	float  w,							// bucket width
	int    num,							// number of repetition
	float& p1,							// p1 = p(w / (2 *r)), returned
	float& p2);							// p2 = p(w / (2 * c * r)), returned


// -----------------------------------------------------------------------------
//  probability vs. w for a fixed ratio c
// -----------------------------------------------------------------------------
void prob_of_gaussian();			// curve of p1, p2 vs. w under gaussian

// -----------------------------------------------------------------------------
void prob_of_cauchy();				// curve of p1, p2 vs. w under cauchy

// -----------------------------------------------------------------------------
void prob_of_levy();				// curve of p1, p2 vs. w under levy


// -----------------------------------------------------------------------------
//  the difference (p1 - p2) vs. w for a fixed ratio
// -----------------------------------------------------------------------------
void diff_prob_of_gaussian();		// curve of p1 - p2 vs. w under gaussian

// -----------------------------------------------------------------------------
void diff_prob_of_cauchy();			// curve of p1 - p2 vs. w under cauchy

// -----------------------------------------------------------------------------
void diff_prob_of_levy();			// curve of p1 - p2 vs. w under levy


// -----------------------------------------------------------------------------
//  rho = log(1/p1) / log(1/p2) vs. w for a fixed ratio c
// -----------------------------------------------------------------------------
void rho_of_gaussian();				// curve of rho vs. w under gaussian

// -----------------------------------------------------------------------------
void rho_of_cauchy();				// curve of rho vs. w under cauchy

// -----------------------------------------------------------------------------
void rho_of_levy();					// curve of rho vs. w under levy

#endif
