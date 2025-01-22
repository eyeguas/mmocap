#ifndef _PSOHUMANPARTTRACKER_
#define _PSOHUMANPARTTRACKER_

#include "psohumantracker.h"
#include "hpsohtprobleminstance.h"
#include "humanparttracker.h"

namespace gummocap
{
namespace algorithms {
  
  class PSOHumanPartTracker: public PSOHumanTracker, public HumanPartTracker
   /**\brief This class represents a human tracker using the PSO algorithm.
    * First the body is tracked and then the different limbs.
    */
  { 
    public:
      struct Params {
      Params(unsigned int noag, double maximumv, double ir, double iw, double maxeval, unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints)
      {
       nbody_joints=nbodyjoints;
       nlarm_joints=nlarmjoints;
       nrarm_joints=nrarmjoints;
       nlleg_joints=nllegjoints;
       nrleg_joints=nrlegjoints;
       number_of_agents=noag;
       stop_maxeval=maxeval;
       maxv=maximumv;
       irang=ir;
       initial_weight=iw;
      }
      Params(){
       nbody_joints=3;
       nlarm_joints=3;
       nrarm_joints=3;
       nlleg_joints=3;
       nrleg_joints=3;
       number_of_agents=1;
       stop_maxeval=1.;
       maxv=0.0;
       irang=0.0;
       initial_weight=1.0;
      }
      unsigned int number_of_agents;	// Number of particles.
      double maxv;			// Maximum velocity: particles velocity must be in [-maxv,maxv]
      double irang;			// Initialization range: [particle initial value + irang, particle initial value + irang]
      double stop_maxeval;		// Termination criterion: maximum number of evaluations.
      double initial_weight;		// Starting value of the inertia weight
      unsigned int nbody_joints;
      unsigned int nlarm_joints;
      unsigned int nrarm_joints;
      unsigned int nlleg_joints;
      unsigned int nrleg_joints;
    };
    /**
    */
    PSOHumanPartTracker ();
    /**Releases the memory
    */
    ~PSOHumanPartTracker (); 
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask, int nBodyParts);
    /**Performs an step of the algorithm: a PSO search.
    */
    void step(guscene::view::ViewSet *VSet)throw(GUException);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();

    protected:
      
      Params _params;
            
      vector<bool>* _DOFmask;
      
      HPSOHTProblemInstance** _pso_instances;
      
      

  };
  
  
}

}

#endif