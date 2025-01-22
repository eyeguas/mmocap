#include "humanparttrackerproblem.h"

namespace gummocap{
  
  namespace humantracker{
    
    HumanPartTrackerProblem::HumanPartTrackerProblem(){
	  _areParams=false;
	  _TheViewSet=NULL;
	}
 /*
 *
 *
 */
      HumanPartTrackerProblem::HumanPartTrackerProblem(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func){
	_areParams=false;
	_TheViewSet=NULL;
	setParams(xglobal,inipart,endpart,sk,VSet,fitness_func);
	}
 /*
 *
 *
 */
 
 HumanPartTrackerProblem::~HumanPartTrackerProblem(){
 }
 /*
 *
 *
 */
 void HumanPartTrackerProblem::setParams(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func){
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose);
    setEval(fitness_func);
    _xglobal=xglobal;
    _inipart=inipart;
    _endpart=endpart;
    _areParams=true;
 }
 /*
 *
 *
 */
 void HumanPartTrackerProblem::setEval(PoseEvaluator* fitness_function){
      _fitness_function=fitness_function;
 }
 /*
 *
 *
 */   
 tFitness HumanPartTrackerProblem::eval(const tChromosomeReal &sol){
	    SkeletonPose pose;
            fromPartVector2Pose(sol,pose );
            _bodyModel.setPose(_auxBodyModel,pose );
	    return _fitness_function->evaluate(&_auxBodyModel,_TheViewSet);
 }
      
    
  }
  
  
}