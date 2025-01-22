#include "psohumanparttracker.h"

namespace gummocap
{
namespace algorithms {
  
  
/**
*/
PSOHumanPartTracker::PSOHumanPartTracker ()
{
    _areParams=false;
}
/**
*/
PSOHumanPartTracker::~PSOHumanPartTracker ()
{
  if(_areParams){
    for(unsigned int i=0;i<_BODY_PARTS;i++){
      delete _pso_instances[i];
    }
    delete[] _pso_instances;
    
  }
}
/**Sets the required params and initializes
*/
void PSOHumanPartTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, int nBodyParts)
{
    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose,DOFmask);
    _N=(_reducedPose.size())*3+3;
    assert((unsigned int)_N==3*(1+_params.nbody_joints+_params.nlarm_joints+_params.nrarm_joints+_params.nlleg_joints+_params.nrleg_joints)+3);

    HumanPartTracker::setParams(nBodyParts,_params.nbody_joints,_params.nlarm_joints,_params.nrarm_joints,_params.nlleg_joints,_params.nrleg_joints);
          
    _xbest=new double[_N];
    for(unsigned int i=0;i<_N;i++){
      _xbest[i]=0.0;
    }
    
    _fitness_function=fitness_func;
   
    _pso_instances=new HPSOHTProblemInstance*[_BODY_PARTS];
    for(unsigned int i=0;i<_BODY_PARTS;i++){
      _pso_instances[i]=new HPSOHTProblemInstance();
      unsigned int ini=_ini_end_body_parts[i*2];
      unsigned int end=_ini_end_body_parts[i*2+1];
      unsigned int n=end-ini+1;
      _pso_instances[i]->setParams(_bodyModel,_N,ini,end,_TheViewSet,_fitness_function,PSO::Params(_params.number_of_agents,n,_params.maxv,_params.irang,_params.irang,_params.initial_weight,_c_1,_c_2,0.0,PSO::INFINITE,_params.stop_maxeval/_BODY_PARTS,1),DOFmask);
    }
    
    _num_steps=0;
    _areParams=true;
}
/**Sets the required params and initializes
*/
void PSOHumanPartTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, vector<bool> *DOFXYZmask,int nBodyParts)
{
    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose,DOFmask);
    //_N=gumocap::skeleton::SkeletonPoseWoT::nDOF(*sk.getSkptr(), DOFmask, DOFXYZmask);
    _reducedDOFXYZmask=new   vector<bool>();
    gumocap::skeleton::SkeletonPoseWoT::toReducedMaskXYZ(*sk.getSkptr(), DOFmask, DOFXYZmask, _reducedDOFXYZmask);
    
    _N=(_reducedPose.size())*3+3;
    //assert((unsigned int)_N==3*(1+_params.nbody_joints+_params.nlarm_joints+_params.nrarm_joints+_params.nlleg_joints+_params.nrleg_joints)+3);

    HumanPartTracker::setParams(nBodyParts,_params.nbody_joints,_params.nlarm_joints,_params.nrarm_joints,_params.nlleg_joints,_params.nrleg_joints);
          
    _xbest=new double[_N];
    for(unsigned int i=0;i<_N;i++){
      _xbest[i]=0.0;
    }
    
    _fitness_function=fitness_func;
   
    _pso_instances=new HPSOHTProblemInstance*[_BODY_PARTS];
    for(unsigned int i=0;i<_BODY_PARTS;i++){
      _pso_instances[i]=new HPSOHTProblemInstance();
      unsigned int ini=_ini_end_body_parts[i*2];
      unsigned int end=_ini_end_body_parts[i*2+1];
      unsigned int n=end-ini+1;
      _pso_instances[i]->setParams(_bodyModel,_N,ini,end,_TheViewSet,_fitness_function,PSO::Params(_params.number_of_agents,n,_params.maxv,_params.irang,_params.irang,_params.initial_weight,_c_1,_c_2,0.0,PSO::INFINITE,_params.stop_maxeval/_BODY_PARTS,1),DOFmask, _reducedDOFXYZmask);
    }
    
    _num_steps=0;
    _areParams=true;
}
/**
  *
  *
  */
  void PSOHumanPartTracker::step(guscene::view::ViewSet *VSet)throw(gu::Exception){
    
    _TheViewSet=VSet;
    
    for(unsigned int h=0;h<_BODY_PARTS;h++){
      _pso_instances[h]->execute(_xbest,_TheViewSet,_num_steps);
      _xbestf=_pso_instances[h]->getBestFitness();
      double* auxb=_pso_instances[h]->getBest();
      unsigned int ini=_ini_end_body_parts[h*2];
      unsigned int end=_ini_end_body_parts[h*2+1];
      for(unsigned int i=ini,j=0;i<=end;i++,j++)
	_xbest[i]=auxb[j];
    }
    
    _num_steps++;
    
  }
  
  /**
*
*/
SkeletonPose PSOHumanPartTracker::getCurrentEstimation()
{
  if(_BODY_PARTS==6)
    return _pso_instances[RLEG]->getCurrentEstimation();
  else
    return _pso_instances[RL_LEG]->getCurrentEstimation();

}
  
  
}

}
