#ifndef _gummocap_PSO_
#define _gummocap_PSO_

#include <iostream>
#include <gu/gurandom.h>

using namespace std;

namespace gummocap
{
namespace algorithms
{
  class PSO
     /**\brief This class represents the Particle Swarm Optimization algorithm.
     */
     {
       public:
	 struct Params {
	    Params(unsigned int noag, unsigned int dim, double maximumv, double irl, double irr, double iw, double c_1, double c_2, double cutoff, double maxiter, double maxeval, double numruns) {
		number_of_agents=noag;
		dimension=dim;
		stop_cutoff=cutoff;
		stop_maxeval=maxeval;
		stop_maxiter=maxiter;
		no_runs=numruns;
		maxv=maximumv;
		irang_l=irl;
		irang_r=irr;
		initial_weight=iw;
		c1=c_1;
		c2=c_2;
	    }
	    Params() {
		dimension=1;
		number_of_agents=1;
		no_runs=1;
		stop_cutoff=0.0;
		stop_maxiter=1.;
		stop_maxeval=1.;
		maxv=0.0;
		irang_l=0.0;
		irang_r=0.0;
		initial_weight=1.0;
		c1=2.0;
		c2=2.0;
	    }
	    unsigned int number_of_agents;	// Number of particles.
	    unsigned int dimension;		// Number of variables in the search.
	    double maxv;		// Maximum velocity: particles velocity must be in [-maxv,maxv]
	    double irang_l, irang_r;	// Initialization range: [particle initial value - irang_l, particle initial value + irang_r]
	    double stop_cutoff;		// Termination criterion: best fitness <= cutoff.
	    double stop_maxiter;	// Termination criterion: maximum number of iterations.
	    double stop_maxeval;	// Termination criterion: maximum number of evaluations.
	    double initial_weight;	// Starting value of the inertia weight
	    double c1;			// Parameters which influence the social and cognition components of the swarm behavior. 
	    double c2;
	    unsigned int no_runs;	// Number of independent runs of the algorithm.
	};
	
	/**
	*/
	PSO();
	/** Releases the memory
	*/
	virtual ~PSO();
	/**Sets the required params
	*/
	void setParams(Params  p);
	/**Initializes the algorithm
	*/
	void initialize();
	/**Initializes the algorithm taking as initial point x
	*/
	void initialize(double* x);
	/**Applies the inertia weight updating rule
	*/
	void updateInertiaWeight();
	/**Evaluates the particles
	*/
	virtual void evaluate()=0;
	/**Updates the current best particle
	*/
	void updateBest();
	/**Updates the velocity and position vectors
	*/
	void updateVelocityAndPosition();
	/**Checks the termination criteria:
	*	a) Number of iterations >= Max iterations.
	*	b) Number of evaluations >= Max evaluations.
	*	c) Best fitness		 <= Cutoff.
	* @return: True if any of termination criteria is fulfilled. False otherwise.
	*/
	bool checkTerminationCriteria(); 
	/**Runs an iteration of the algorithm
	*/
	void nextIteration();
	/**Makes a run of the algorithm taking as initial point x
	*/
	void run(double* x);
	/**Makes a run of the algorithm
	*/
	void run();
	/**Reset the best ever found solution
	*/
	void resetBestEver();
	/**Executes the algorithm
	*/
	void execute();
	/**Executes the algorithm
	*/
	void execute(double* x);
	/**Returns the best solution ever found
	*/
	double* getBest();
	/**Returns the best fitness ever found
	*/
	double getBestFitness();
	
	const static double INFINITE=1000000.0; // Infinite value
       protected:
	 Params _params;
	 bool _areParams;

	 double** _xx;			// Position vector
	 double** _vx;			// Velocity vector			
	 double** _tx;			// TX = XX + VX
	 double** _pbestx;		// Best position vector
	 double*  _pbest;		// Best fitness vector
	 double*  _pf;			// Current fitness vector
	 
	 double* _besteverx;		// Best solution ever found
	 double  _besteverf;		// Best fitness ever found
	 
	 double _weight_up;		// Current inertia weight
	 double _C;			// Desired number of inertia weight changes (No of evaluations/No of agents)
	 double _incr_c;		// Parameter of the inertia weight updating rule
	 double _c;			// Sampling variable of the inertia weight updating rule
	 
	 unsigned int _current_it;	// Current iteration
	 unsigned int _current_nevals;	// Current number of evaluations
	 unsigned int _gbest;		// Best particle
	 
  };
}
}

#endif
