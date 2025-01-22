#ifndef _HPSO_HT_PROBLEMINSTANCE_
#define _HPSO_HT_PROBLEMINSTANCE_

#include "psohtprobleminstance.h"

namespace gummocap
{
namespace algorithms
{
  class HPSOHTProblemInstance: public PSOHTProblemInstance
  /**\brief This class represents the Hierarchical Particle Swarm Optimization algorithm applied to an instance of the human tracker problem.
     */
  {
    public:
      /**
	*/
      HPSOHTProblemInstance();
      /** Releases the memory
	*/
      ~HPSOHTProblemInstance();
      /**Sets the required params
	*/
      void setParams(BodyModel<HumanSkeleton> &sk,unsigned int globalN,unsigned int inibodypart, unsigned int endbodypart, guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, PSO::Params p, vector<bool> *DOFmask);
      /**Evaluates the particles
	*/
      void evaluate();
        /**Executes the algorithm
	*/
      void execute(double* x, guscene::view::ViewSet *VSet,unsigned int FrameId);
       /**Gets the current best solution
	*/
	SkeletonPose getCurrentEstimation();
      
    private:
      /** converts a vector in a reduced skeleton pose
    */

    void fromPartVector2PoseWoT(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPoseWoT &pose){
      assert(_areParams);
      pose=_reducedPose;
      unsigned int index=3, index_x=0;
      
      if(inibodypart==0){
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
      for(unsigned int i=1;i<pose.size();i++){
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

      
      unsigned int _inibodypart;
      unsigned int _endbodypart;
      unsigned int _N;
      double* _xbest;
  };
  
}

}

#endif