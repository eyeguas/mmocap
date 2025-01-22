#ifndef _PSOPROBLEMINSTANCE_
#define _PSOPROBLEMINSTANCE_

#include "pso.h"
#include "poseevaluators.h"
#include <omp.h>

namespace gummocap
{
namespace algorithms
{
  
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;
  
  class PSOHTProblemInstance: public PSO
  /**\brief This class represents the Particle Swarm Optimization algorithm applied to an instance of the human tracker problem.
     */
  {
    public:
      /**
	*/
      PSOHTProblemInstance();
      /** Releases the memory
	*/
      virtual ~PSOHTProblemInstance();
      /**Sets the required params
	*/
      void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, PSO::Params p, vector<bool> *DOFmask);
      void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, PSO::Params p, vector<bool> *DOFmask, vector<bool> *DOFXYZmask);
      /**Evaluates the particles
	*/
      void evaluate();
      /**Executes the algorithm
	*/
      void execute(double* x, guscene::view::ViewSet *VSet,unsigned int FrameId);
      /**Gets the current best solution
	*/
	SkeletonPose getCurrentEstimation();
     
    protected:
      /** converts a vector in a reduced skeleton pose
        */
/*
    void fromVector2PoseWoT(double const *x, SkeletonPoseWoT &pose) {
        assert(_areParams);
        pose=_reducedPose;

        pose.Translation[0]=x[0];
        pose.Translation[1]=x[1];
        pose.Translation[2]=x[2];
        int index=3;
        //int index=0;
        // pose[0]=root
        for (unsigned int i=0;i<pose.size();i++) {
            pose[i].x=x[index++];
            pose[i].y=x[index++];
            pose[i].z=x[index++];
        }
    }
*/

void fromVector2PoseWoT(double const *x, SkeletonPoseWoT &pose) {
        //assert(_areParams);
        pose=_reducedPose;

        pose.Translation[0]=x[0];
        pose.Translation[1]=x[1];
        pose.Translation[2]=x[2];
        int index=3;
        //int index=0;
        // pose[0]=root
	if(_reducedDOFXYZmask!=NULL){
        for (unsigned int i=0,j=0;i<pose.size();i++) {
	  if((*_reducedDOFXYZmask)[j++])
            pose[i].x=x[index++];
	 
	  if((*_reducedDOFXYZmask)[j++])
            pose[i].y=x[index++];
	  
	  if((*_reducedDOFXYZmask)[j++])
            pose[i].z=x[index++];
	 
        }
	}
        else{
	  for (unsigned int i=0;i<pose.size();i++) {
	    pose[i].x=x[index++];
            pose[i].y=x[index++];
            pose[i].z=x[index++];
	  }
	}
    }

     /** converts a vector in a full skeleton pose
    */

    void fromVector2Pose(double const *x, SkeletonPose &pose)
    {
        SkeletonPoseWoT redPose;
        fromVector2PoseWoT(x,redPose);
        SkeletonPoseWoT::fromWoT2Complete(redPose,pose);

    }

      
      PoseEvaluator* _fitness_function;
      vector<PoseEvaluator*> _evaluators;
      gumocap::skeleton::SkeletonPoseWoT _reducedPose;
      BodyModel<HumanSkeleton> _bodyModel,_auxBodyModel;
      vector<BodyModel<HumanSkeleton> > _auxBodyModelV;
      guscene::view::ViewSet *_TheViewSet;
      vector<bool> *_reducedDOFXYZmask;
  };
  
}

}

#endif
