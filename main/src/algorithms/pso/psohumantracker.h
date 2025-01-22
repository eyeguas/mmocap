#ifndef _PSO_HUMANTRACKER_
#define _PSO_HUMANTRACKER_

#include "humantracker.h"
#include "pso.h"
#include "psohtprobleminstance.h"
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include "poseevaluators.h"

namespace gummocap
{
namespace algorithms {
  
  using namespace gumocap;
  using namespace gumocap::models;
  using namespace gumocap::skeleton;
  using namespace gummocap::poseevaluators;
 
  class PSOHumanTracker: public HumanTracker
   /**\brief This class represents a human tracker using the PSO algorithm.
   */
  { 
    public:
	struct Params {
	    Params(unsigned int noag, double maximumv, double ir, double iw, double maxeval) {
		number_of_agents=noag;
		stop_maxeval=maxeval;
		maxv=maximumv;
		irang=ir;
		initial_weight=iw;
	    }
	    Params() {
		number_of_agents=1;
		stop_maxeval=1.;
		maxv=0.0;
		irang=0.0;
		initial_weight=1.0;
	    }
	    unsigned int number_of_agents;	// Number of particles.
	    double maxv;		// Maximum velocity: particles velocity must be in [-maxv,maxv]
	    double irang;	// Initialization range: [particle initial value + irang, particle initial value + irang]
	    double stop_maxeval;	// Termination criterion: maximum number of evaluations.
	    double initial_weight;	// Starting value of the inertia weight
	};
	/**
	*/
	PSOHumanTracker ();
	/**Releases the memory
	*/
	~PSOHumanTracker ();
	/**Sets the required params and initializes
	*/
	void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask);
	void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, vector<bool> *DOFXYZmask);
	/**Performs an step of the algorithm: a PSO search.
	*/
	void step(guscene::view::ViewSet *VSet)throw(gu::Exception);
	/**Gets the current best solution
	*/
	SkeletonPose getCurrentEstimation();
	/**Enables debug mode
	*/
	void enableDebug(bool enable){}
	
    protected:
	bool _areParams;
	Params _params;
	BodyModel<HumanSkeleton> _bodyModel;
	gumocap::skeleton::SkeletonPoseWoT _reducedPose;
	guscene::view::ViewSet *_TheViewSet;

	unsigned int _num_steps;
	unsigned int _N;
	PoseEvaluator* _fitness_function;
	double* _xbest;
	double _xbestf;
    
	vector<bool> *_reducedDOFXYZmask;
	
	PSOHTProblemInstance* _pso_instance;
	static const double _c_1=2.0;   // PSO social component
	static const double _c_2=2.0;	// PSO cognition component
  };
  
}
}

#endif