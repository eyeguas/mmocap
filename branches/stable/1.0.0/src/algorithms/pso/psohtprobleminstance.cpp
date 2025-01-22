#include "psohtprobleminstance.h"
#include <omp.h>
namespace gummocap
{
namespace algorithms
{
 
  PSOHTProblemInstance::PSOHTProblemInstance(){
    _areParams=false;
    _TheViewSet=NULL;
  }
  /**
  *
  */
  PSOHTProblemInstance::~PSOHTProblemInstance(){
    if(_areParams){
      for (unsigned int i=0;i<_evaluators.size();i++)
            delete _evaluators[i];
    }
  }
  /**
  *
  */
  void PSOHTProblemInstance::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, PSO::Params p, vector<bool> *DOFmask){
    
    PSO::setParams(p);
    _bodyModel=sk;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
    _TheViewSet=VSet;
    _fitness_function=fitness_func;
    for (int i=0;i<omp_get_max_threads();i++) { 
        _evaluators.push_back( _fitness_function->makeCopy());
        _auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
    }
    _areParams=true;
  }
 
   /**
 *
 *
 */ 
  void PSOHTProblemInstance::execute(double* x, guscene::view::ViewSet *VSet,unsigned int FrameId){
      _TheViewSet=VSet;
      //indicate to evaluators that frame analysis starts
      for (unsigned int i=0;i<_evaluators.size();i++)
        _evaluators[i]->startFrame(VSet,FrameId);
      resetBestEver();
      PSO::execute(x);
      //indicate to evaluators that frame analysis is finished
      SkeletonPose bestEstimation=getCurrentEstimation();
      _bodyModel.setPose(_auxBodyModel,bestEstimation);
      for(unsigned int i=0;i<_evaluators.size();i++)
	_evaluators[i]->frameFinished(&_auxBodyModel);
  }
  /**
  *
  */
  void PSOHTProblemInstance::evaluate(){
    #pragma omp parallel for
    for(unsigned int i=0;i<_params.number_of_agents;i++){
            int threadId=omp_get_thread_num();
	    SkeletonPose pose;
            fromVector2Pose(_xx[i],pose );
            _bodyModel.setPose(_auxBodyModelV[threadId],pose );
            _pf[i]=_evaluators[threadId]->evaluate(&_auxBodyModelV[threadId],_TheViewSet);
    }
      
  }
   /**
*
*/
SkeletonPose PSOHTProblemInstance::getCurrentEstimation()
{
    SkeletonPose pose;

    fromVector2Pose(_besteverx,pose);

    return pose;

}
}

}