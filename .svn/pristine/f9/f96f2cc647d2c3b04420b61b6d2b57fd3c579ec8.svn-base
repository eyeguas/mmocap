
#ifndef __DEHUMANPARTTRACKER__
#define __DEHUMANPARTTRACKER__

#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <omp.h>
#include "dehumantracker.h"
#include "poseevaluators.h"
#include "humanparttracker.h"
#include "dehpt_fitnessfunction.h"

namespace gummocap{

  using namespace gumocap;
  using namespace gumocap::models;
  using namespace gumocap::skeleton;
  using namespace gummocap::poseevaluators;

  
  namespace algorithms{
 
    class DEHumanPartTracker : public DEHumanTracker, public HumanPartTracker
    /**\brief This class represents a human tracker using the DE algorithm.
    * First the body is tracked and then the different limbs.
    */
    {
      public:
	struct Params {
	    Params(unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints, unsigned int population_size, unsigned int nLPC, double minbound, double maxbound, double c_r, double F, double maxeval) {
		stop_maxeval=maxeval;
		min_bound=minbound;
		max_bound=maxbound;
		pop_size=population_size;
		cr=c_r;
		f=F;
		nLPCsols=nLPC;
		nbody_joints=nbodyjoints;
		nlarm_joints=nlarmjoints;
		nrarm_joints=nrarmjoints;
		nlleg_joints=nllegjoints;
		nrleg_joints=nrlegjoints;
	    }
	    Params() {
		stop_maxeval=1.;
		min_bound=0.0;
		max_bound=0.0;
		pop_size=1;
		cr=0.9;
		f=0.5;
		nLPCsols=1;
		nbody_joints=3;
		nlarm_joints=3;
		nrarm_joints=3;
		nlleg_joints=3;
		nrleg_joints=3;
	    }
	    double stop_maxeval;	// Termination criterion: maximum number of evaluations.
	    double min_bound, max_bound;	// Minimum and maximum bound in the search: [x-min_bound,x+max_bound]
	    unsigned int pop_size;		// Population size
	    double cr;				// Crossover constant
	    double f;				// Controls the amplification of the differential variation (Mutation)
	    unsigned int nLPCsols;		// Number of Leading, Placing and Correcting solutions
	    unsigned int nbody_joints;
	    unsigned int nlarm_joints;
	    unsigned int nrarm_joints;
	    unsigned int nlleg_joints;
	    unsigned int nrleg_joints;
      
	};
	
	/**
	*/
	DEHumanPartTracker ();
	/**Releases the memory
	*/
	~DEHumanPartTracker (); 
	/**Sets the required params and initializes
	*/
	void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask, unsigned int nBodyParts);
	/**Performs an step of the algorithm in a full body.
	*/
	void step(guscene::view::ViewSet *VSet)throw(GUException);

	
      private:
	
	Params _params;
	vector<vector< FitnessFunction*> > _pffs;
    
	
    };
    
  }
  
}

#endif