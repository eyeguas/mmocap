#ifndef __HUMANPARTTRACKERPROBLEM__
#define __HUMANPARTTRACKERPROBLEM__

#include <gurealea/realea/problem.h>
#include <gurealea/realea/ea.h>
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include "poseevaluators.h"

using namespace std;
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;
using namespace realea;

namespace gummocap{
  
  namespace humantracker{
/**
 * @class HumanTrackerProblem
 *
 * @brief This class represent a human tracker problem to be resolved
 *
 * To apply the EA to new problems, a subclass of ProblemFactory must be defined, and it must
 * return an adequated initialized problem (instance of this class)
 *
 * @see Problem, ProblemFactory
 */
class HumanPartTrackerProblem: public GeneralProblem{
      public:
	 HumanPartTrackerProblem();
 /*
 *
 *
 */
      HumanPartTrackerProblem(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
 /*
 *
 *
 */
 
 ~HumanPartTrackerProblem();
 /*
 *
 *
 */
 void setParams(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
 /*
 *
 *
 */
 void setEval(PoseEvaluator* fitness_function);
 /*
 *
 *
 */   
 tFitness eval(const tChromosomeReal &sol);
      
      protected:
	PoseEvaluator* _fitness_function;
	BodyModel<HumanSkeleton> _bodyModel,_auxBodyModel;
	bool _areParams;
	gumocap::skeleton::SkeletonPoseWoT _reducedPose;
	guscene::view::ViewSet *_TheViewSet;
	unsigned int _inipart;
	unsigned int _endpart;
	tChromosomeReal* _xglobal;
	
    void fromPartVector2PoseWoT(const tChromosomeReal &x, SkeletonPoseWoT &pose){
      assert(_areParams);
      pose=_reducedPose;
      unsigned int index=3, index_x=0;
      
      if(_inipart==0){
	pose.Translation[0]=x[0];
	pose.Translation[1]=x[1];
	pose.Translation[2]=x[2];
	index_x=3;
      }
      else{
	pose.Translation[0]=(*_xglobal)[0];
	pose.Translation[1]=(*_xglobal)[1];
	pose.Translation[2]=(*_xglobal)[2];
      }
      
      // pose[0]=root
      for(unsigned int i=1;i<pose.size();i++){
	  if((index<_inipart)||(index>_endpart)){
	    pose[i].x=(*_xglobal)[index++];
	    pose[i].y=(*_xglobal)[index++];
	    pose[i].z=(*_xglobal)[index++];
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
    
    void fromPartVector2Pose(const tChromosomeReal &x, SkeletonPose &pose)
    {
      SkeletonPoseWoT redPose;
      fromPartVector2PoseWoT(x,redPose);
      SkeletonPoseWoT::fromWoT2Complete(redPose,pose);
      
    }
     };

 

  }

}

#endif