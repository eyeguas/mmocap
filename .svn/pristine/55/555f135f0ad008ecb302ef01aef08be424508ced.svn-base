#include "hpsohtprobleminstance.h"

namespace gummocap
{
namespace algorithms
{
  
  HPSOHTProblemInstance::HPSOHTProblemInstance(){
    
  }
  /**
  *
  */
  HPSOHTProblemInstance::~HPSOHTProblemInstance(){
    if(_areParams){
      delete[] _xbest;
    }
  }
  /**
  *
  */
  void HPSOHTProblemInstance::setParams(BodyModel<HumanSkeleton> &sk, unsigned int globalN, unsigned int inibodypart, unsigned int endbodypart, guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, PSO::Params p, vector<bool> *DOFmask){
    _N=globalN;
    _xbest=new double[globalN];
    _inibodypart=inibodypart;
    _endbodypart=endbodypart;
    PSOHTProblemInstance::setParams(sk,VSet,fitness_func,p,DOFmask);
  }
  /**
 *
 *
 */ 
  void HPSOHTProblemInstance::execute(double* x, guscene::view::ViewSet *VSet,unsigned int FrameId){
      _TheViewSet=VSet;
      for(unsigned int i=0;i<_N;i++)
	_xbest[i]=x[i];
      //indicate to evaluators that frame analysis starts
      for (unsigned int i=0;i<_evaluators.size();i++)
        _evaluators[i]->startFrame(VSet,FrameId);
      resetBestEver();
      PSO::execute(x+_inibodypart);
      //indicate to evaluators that frame analysis is finished
      SkeletonPose bestEstimation=getCurrentEstimation();
      _bodyModel.setPose(_auxBodyModel,bestEstimation);
      for(unsigned int i=0;i<_evaluators.size();i++)
	_evaluators[i]->frameFinished(&_auxBodyModel);
  }
  /**
  *
  */
  void HPSOHTProblemInstance::evaluate(){
    #pragma omp parallel for
    for(unsigned int i=0;i<_params.number_of_agents;i++){
            int threadId=omp_get_thread_num();
	    SkeletonPose pose;
            fromPartVector2Pose(_xx[i],_inibodypart,_endbodypart,pose );
            _bodyModel.setPose(_auxBodyModelV[threadId],pose );
            _pf[i]=_evaluators[threadId]->evaluate(&_auxBodyModelV[threadId],_TheViewSet);
    }
      
  }
  /**
*
*/
SkeletonPose HPSOHTProblemInstance::getCurrentEstimation()
{
    SkeletonPose pose;

    fromPartVector2Pose(_besteverx,_inibodypart,_endbodypart,pose);

    return pose;

}
 
}

}
