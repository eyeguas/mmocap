#ifndef _gummocap_MACMACHAINSHUMANPARTTRACKER_
#define _gummocap_MACMACHAINSHUMANPARTTRACKER_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <omp.h>
#include "macmachainshumantracker.h"
#include "humanparttrackerprobleminstance.h"
#include "poseevaluators.h"
#include "humanparttracker.h"
#include <vector>

namespace gummocap
{
using namespace std;  
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;
namespace algorithms{
class MACMAChainsHumanPartTracker : public MACMAChainsHumanTracker, public HumanPartTracker
/**\brief This class represents a human tracker using the MA-CMA-Chains algorithm.
* First the body is tracked and then the different limbs.
*/
{

  public:
    struct Params {
      Params(double maxiter, double max_ev, double min_dom, double max_dom, double ls_int, double lg_rat, unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints)
      {
       nbody_joints=nbodyjoints;
       nlarm_joints=nlarmjoints;
       nrarm_joints=nrarmjoints;
       nlleg_joints=nllegjoints;
       nrleg_joints=nrlegjoints;
       stop_maxiter=maxiter;
       min_domain=min_dom;
       max_domain=max_dom;
       ls_intensity=ls_int;
       lg_ratio=lg_rat;
       max_eval=max_ev;
      }
      Params(){
       nbody_joints=3;
       nlarm_joints=3;
       nrarm_joints=3;
       nlleg_joints=3;
       nrleg_joints=3;
       stop_maxiter=1;
       min_domain=0;
       max_domain=0;
       ls_intensity=500;
       lg_ratio=0.5;
       max_eval=10000;
      }
      unsigned int nbody_joints;
      unsigned int nlarm_joints;
      unsigned int nrarm_joints;
      unsigned int nlleg_joints;
      unsigned int nrleg_joints;
      double stop_maxiter;	// Maximum number of iterations for each step
      double min_domain;	// Search domain: [xbest-min_domain,xbest+max_domain]
      double max_domain;
      double max_eval;	// Maximum number of evaluations per iteration = max_eval * dimensionality of the problem
      double ls_intensity;	// Local search intensity in number of evaluations
      double lg_ratio;	// Local search / Global search ratio
    
    };

    /**
    */
    MACMAChainsHumanPartTracker ();
    /**Releases the memory
    */
    ~MACMAChainsHumanPartTracker (); 
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask,unsigned int nBodyParts);
    /**Performs an step of the algorithm in a full body.
    */
    void step(guscene::view::ViewSet *VSet)throw(gu::Exception);
    //void stepSeq(guscene::view::ViewSet *VSet)throw(gu::Exception);
    //void stepParallel(guscene::view::ViewSet *VSet)throw(gu::Exception);
    //void stepParallelLimbs(guscene::view::ViewSet *VSet)throw(gu::Exception);
    //void stepLimbs(unsigned int first_limb, unsigned int second_limb, tChromosomeReal xlimbs_in, tChromosomeReal xlimbs_out, unsigned int ini_first_limb, unsigned int end_first_limb, unsigned int ini_second_limb, unsigned int end_second_limb, double niterations, double numevals)throw(gu::Exception);
    /**Performs an step of the algorithm in a body part: a CMAES search.
    */
    void stepPart(tChromosomeReal* x, tChromosomeReal* xbest, double* best_fitness, unsigned int inibodypart, unsigned int endbodypart, double niterations, double maxeval)throw(gu::Exception);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();
    /**
     */
    void enableDebug(bool enable){}
  protected:

/** converts a vector in a reduced skeleton pose
    */

    void fromPartVector2PoseWoT(const tChromosomeReal &x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPoseWoT &pose){
      assert(_areParams);
      pose=_reducedPose;
      unsigned int index=3, index_x=0;
      
      if(inibodypart==_initrans){
	pose.Translation[0]=x[0];
	pose.Translation[1]=x[1];
	pose.Translation[2]=x[2];
	index_x=3;
      }
      else{
	pose.Translation[0]=(*_xbest)[0];
	pose.Translation[1]=(*_xbest)[1];
	pose.Translation[2]=(*_xbest)[2];
      }
      
      // pose[0]=root
      for(unsigned int i=0;i<pose.size();i++){
	  if((index<inibodypart)||(index>endbodypart)){
	    pose[i].x=(*_xbest)[index++];
	    pose[i].y=(*_xbest)[index++];
	    pose[i].z=(*_xbest)[index++];
	  }
	  else{
	    pose[i].x=x[index_x++];
	    pose[i].y=x[index_x++];
	    pose[i].z=x[index_x++];
	    index+=3;
	  }
      }
    }

    /** converts a vector in a full skeleton pose
    */
    
    void fromPartVector2Pose(const tChromosomeReal &x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPose &pose)
    {
      SkeletonPoseWoT redPose;
      fromPartVector2PoseWoT(x,inibodypart, endbodypart, redPose);
      SkeletonPoseWoT::fromWoT2Complete(redPose,pose);
      
    }
 
  private:

    Params _params;

};
}

} 
#endif

