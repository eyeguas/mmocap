#include "psohumantracker.h"

namespace gummocap
{
namespace algorithms
{
  /**
  *
  *
  */
  PSOHumanTracker::PSOHumanTracker()
  {
  }
  /**
  *
  *
  */
  PSOHumanTracker::~PSOHumanTracker()
  {
    if(_areParams)
      delete[] _xbest;
  }
  /**
  *
  *
  */
  void PSOHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask)
  {
    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
    _N=(_reducedPose.size())*3+3;
    _xbest=new double[_N];
    for(unsigned int i=0;i<_N;i++)
      _xbest[i]=0.0;
    _xbestf=PSO::INFINITE;
    
    _fitness_function=fitness_func;
    
    _pso_instance=new PSOHTProblemInstance();
    _pso_instance->setParams(_bodyModel,_TheViewSet,_fitness_function, PSO::Params(_params.number_of_agents, _N, _params.maxv, _params.irang, _params.irang, _params.initial_weight, _c_1, _c_2, 0.0, PSO::INFINITE, _params.stop_maxeval, 1), DOFmask);
    
    _num_steps=0;
    _areParams=true;

  }
  /**
  *
  *
  */
  void PSOHumanTracker::step(guscene::view::ViewSet *VSet)throw(GUException){
    
    _TheViewSet=VSet;
    
    _pso_instance->execute(_xbest,_TheViewSet,_num_steps);
    _xbestf=_pso_instance->getBestFitness();
    double* auxb=_pso_instance->getBest();
    for(unsigned int i=0;i<_N;i++)
      _xbest[i]=auxb[i];
    _num_steps++;
    
  }
  /**
*
*/
SkeletonPose PSOHumanTracker::getCurrentEstimation()
{

    return _pso_instance->getCurrentEstimation();

}
}

}
