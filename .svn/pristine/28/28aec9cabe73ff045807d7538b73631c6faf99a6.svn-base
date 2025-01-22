#ifndef _gummocap_DE_
#define _gummocap_DE_


#include <gsl/gsl_rng.h>
#include "deht_fitnessfunction.h"
#include "de_d40_mm.h"
#include "humantracker.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <values.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

namespace gummocap{
  
  namespace algorithms{
 
    #define abss(a)     (a<0 ? (-a) : a)

    class DEHumanTracker: public HumanTracker
    /**\brief This class represents the Diferential Evolution algorithm (DE-D40+Mm) applied to Human Tracking.
     */
    {
      public:
	struct Params {
	    Params(unsigned int population_size, unsigned int nLPC, double minbound, double maxbound, double c_r, double F, double maxeval) {
		stop_maxeval=maxeval;
		min_bound=minbound;
		max_bound=maxbound;
		pop_size=population_size;
		cr=c_r;
		f=F;
		nLPCsols=nLPC;
	    }
	    Params() {
		stop_maxeval=1.;
		min_bound=0.0;
		max_bound=0.0;
		pop_size=1;
		cr=0.9;
		f=0.5;
		nLPCsols=1;
	    }
	    double stop_maxeval;	// Termination criterion: maximum number of evaluations.
	    double min_bound, max_bound;	// Minimum and maximum bound in the search: [x-min_bound,x+max_bound]
	    unsigned int pop_size;		// Population size
	    double cr;				// Crossover constant
	    double f;				// Controls the amplification of the differential variation (Mutation)
	    unsigned int nLPCsols;		// Number of Leading, Placing and Correcting solutions
	};
	
	/**
	*/
	DEHumanTracker ();
	/**Releases the memory
	*/
	~DEHumanTracker ();
	/**Sets the required params and initializes
	*/
	void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask);
	void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, vector<bool> *DOFXYZmask);
	/**Performs an step of the algorithm: a DE search.
	*/
	void step(guscene::view::ViewSet *VSet)throw(gu::Exception);
	/**Gets the current best solution
	*/
	SkeletonPose getCurrentEstimation();
	/**
	*/
	void enableDebug(bool enable){}
      
      protected:
	  const static double PI;
	  const static double E;
	  
	  bool _areParams;
	  Params _params;
	  
	  BodyModel<HumanSkeleton> _bodyModel,_auxBodyModel;
	  vector<BodyModel<HumanSkeleton> > _auxBodyModelV;
	  
	  gumocap::skeleton::SkeletonPoseWoT _reducedPose;
	  
	  guscene::view::ViewSet *_TheViewSet;
	  
	  PoseEvaluator* _fitness_function;
	  vector<FitnessFunction*> _ffs;
	  vector<PoseEvaluator*> _evaluators;
	  
	  vector<bool> *_reducedDOFXYZmask;
	  

	  HTFitnessFunction* _ff;
	  double long _fvalue;
	  
	  double* _xbest;
	  unsigned int _N;
	  unsigned int _num_steps;
	
    };
    
  }
  
}

#endif