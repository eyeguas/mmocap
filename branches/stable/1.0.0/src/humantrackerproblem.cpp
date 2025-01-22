#include "humantrackerproblem.h"

namespace gummocap{
  
  namespace humantracker{
    
    HumanTrackerProblem::HumanTrackerProblem(){
	  _areParams=false;
	  _TheViewSet=NULL;
	}
 /*
 *
 *
 */
      HumanTrackerProblem::HumanTrackerProblem(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func){
	_areParams=false;
	_TheViewSet=NULL;
	setParams(sk,VSet,fitness_func);
	}
 /*
 *
 *
 */
 
 HumanTrackerProblem::~HumanTrackerProblem(){
 }
 /*
 *
 *
 */
 void HumanTrackerProblem::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func){
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose);
    setEval(fitness_func);
    _areParams=true;
 }
 /*
 *
 *
 */
  void HumanTrackerProblem::setViewSet(guscene::view::ViewSet *VSet){
    _TheViewSet=VSet;
  }
  /*
 *
 *
 */
 void HumanTrackerProblem::setEval(PoseEvaluator* fitness_function){
      _fitness_function=fitness_function;
 }
 /*
 *
 *
 */   
 tFitness HumanTrackerProblem::eval(const tChromosomeReal &sol){
	    SkeletonPose pose;
            fromVector2Pose(sol,pose );
            _bodyModel.setPose(_auxBodyModel,pose );
	    return _fitness_function->evaluate(&_auxBodyModel,_TheViewSet);
 }
      
    
  }
  
  
}