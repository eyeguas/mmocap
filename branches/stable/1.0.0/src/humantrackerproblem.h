#ifndef __HUMANTRACKERPROBLEM__
#define __HUMANTRACKERPROBLEM__

#include <gurealea/realea/problem.h>
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
class HumanTrackerProblem: public GeneralProblem{
      public:
	 HumanTrackerProblem();
 /*
 *
 *
 */
      HumanTrackerProblem(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
 /*
 *
 *
 */
 
 ~HumanTrackerProblem();
 /*
 *
 *
 */
 void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
 /*
 *
 *
 */
 void setViewSet(guscene::view::ViewSet *VSet);
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
	
	 /** converts a vector in a reduced skeleton pose
        */

    void fromVector2PoseWoT(const tChromosomeReal &x, SkeletonPoseWoT &pose) {
        assert(_areParams);
        pose=_reducedPose;

        pose.Translation[0]=x[0];
        pose.Translation[1]=x[1];
        pose.Translation[2]=x[2];
        int index=3;
        //int index=0;
        // pose[0]=root
        for (unsigned int i=1;i<pose.size();i++) {
            pose[i].x=x[index++];
            pose[i].y=x[index++];
            pose[i].z=x[index++];
        }
    }

    /** converts a vector in a full skeleton pose
    */

    void fromVector2Pose(const tChromosomeReal &x, SkeletonPose &pose)
    {
        SkeletonPoseWoT redPose;
        fromVector2PoseWoT(x,redPose);
        SkeletonPoseWoT::fromWoT2Complete(redPose,pose);

    }
    };
    
    typedef auto_ptr<HumanTrackerProblem> HTProblemPtr;
    typedef HTProblemPtr& HTProblemParamPtr;
 

  }

}

#endif