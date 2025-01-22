#ifndef _gummocap_CMAESHUMANPARTTRACKER_
#define _gummocap_CMAESHUMANPARTTRACKER_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <omp.h>
#include "cmaeshumantracker.h"
#include "poseevaluators.h"
#include "humanparttracker.h"

namespace gummocap
{
  
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;
namespace algorithms{
class CMAESHumanPartTracker : public CMAESHumanTracker, public HumanPartTracker
/**\brief This class represents a human tracker using the CMAES algorithm.
* First the body is tracked and then the different limbs.
*/
{

  public:
    struct Params {
      Params(int lambd, double stddev_ini, double maxiter, double maxeval, unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints, double stddev_body,double stddev_arms,double stddev_legs,double stddev_trans, double stddev_rot)
      {std_ini=stddev_ini; stop_maxiter=maxiter;
       stop_maxeval=maxeval;
       nbody_joints=nbodyjoints;
       nlarm_joints=nlarmjoints;
       nrarm_joints=nrarmjoints;
       nlleg_joints=nllegjoints;
       nrleg_joints=nrlegjoints;
       std_body=stddev_body;
       std_arms=stddev_arms;
       std_legs=stddev_legs;
       std_trans=stddev_trans;
       std_rot=stddev_rot;
       lambda=lambd;
      }
      Params(){std_ini=0.0; stop_maxiter=1;
       stop_maxeval=1;
       nbody_joints=3;
       nlarm_joints=3;
       nrarm_joints=3;
       nlleg_joints=3;
       nrleg_joints=3;
       std_body=0.0;
       std_arms=0.0;
       std_legs=0.0;
       std_trans=0.0;
       std_rot=0.0;
       lambda=0;
      }
      double std_ini;		// Standard deviation for each point (parameter: rotation or translation) in the search
      double stop_maxiter;	// Maximum number of iterations for each step
      double stop_maxeval;	// Maximum number of evaluations for each step
      unsigned int nbody_joints;
      unsigned int nlarm_joints;
      unsigned int nrarm_joints;
      unsigned int nlleg_joints;
      unsigned int nrleg_joints;
      double std_body;
      double std_arms;
      double std_legs;
      double std_trans;
      double std_rot;
      int lambda;
    };

    /**
    */
    CMAESHumanPartTracker ();
    /**Releases the memory
    */
    ~CMAESHumanPartTracker (); 
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask, unsigned int nBodyParts);
    /**Performs an step of the algorithm in a full body.
    */
    void step(guscene::view::ViewSet *VSet)throw(gu::Exception);
    //void step2(guscene::view::ViewSet *VSet)throw(gu::Exception);
    /**Performs an step of the algorithm in a body part: a CMAES search.
    */
    /**Performs an step of the algorithm in a body part: a CMAES search.
    */
    void stepPart(double* x, double* xbest, double* best_fitness, unsigned int inibodypart, unsigned int endbodypart, double niterations, double maxevals)throw(gu::Exception);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();
    /**
     */
    void enableDebug(bool enable){}
  protected:

/** converts a vector in a reduced skeleton pose
    */

    void fromPartVector2PoseWoT(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPoseWoT &pose){
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
	pose.Translation[0]=_xbest[0];
	pose.Translation[1]=_xbest[1];
	pose.Translation[2]=_xbest[2];
      }
      
      // pose[0]=root
      for(unsigned int i=0;i<pose.size();i++){
	  if((index<inibodypart)||(index>endbodypart)){
	    pose[i].x=_xbest[index++];
	    pose[i].y=_xbest[index++];
	    pose[i].z=_xbest[index++];
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
    
    void fromPartVector2Pose(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPose &pose)
    {
      SkeletonPoseWoT redPose;
      fromPartVector2PoseWoT(x,inibodypart, endbodypart, redPose);
      SkeletonPoseWoT::fromWoT2Complete(redPose,pose);
      
    }
 
  private:
    
    Params _params;
    
};
};

}; 
#endif

