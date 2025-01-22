
#include "dehumanparttracker.h"

namespace gummocap{
  
  namespace algorithms{
    
/**
*/
DEHumanPartTracker::DEHumanPartTracker ()
{
}
/**
*/
DEHumanPartTracker::~DEHumanPartTracker ()
{
  if(_areParams){
    }
}
/**Sets the required params and initializes
*/
void DEHumanPartTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, unsigned int nBodyParts)
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
    _prev_best_fitness=INFINITE;
    _fitness_function=fitness_func;
    //Create independent evaluators , one for each limb
    assert(_evaluators.size()==0);
    
        //Create independent evaluators , one for each thread
    assert(_evaluators.size()==0);
    for (int i=0;i<omp_get_max_threads();i++) { 
	_auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
	_evaluators.push_back( _fitness_function->makeCopy());
        _ffs.push_back(new HTFitnessFunction(_N,_evaluators[i],_auxBodyModelV[i],_TheViewSet,DOFmask));
    }
    for (unsigned int h=0;h<_BODY_PARTS;h++){
	_pffs.push_back(vector<FitnessFunction* >());
	for(int i=0;i<omp_get_max_threads();i++){
	  _pffs[h].push_back(new HPTFitnessFunction(_xbest,_N,_ini_end_body_parts[h*2],_ini_end_body_parts[h*2+1],_evaluators[i],_auxBodyModelV[i],_TheViewSet,DOFmask));
	}
    }
    
    _ff = new HTFitnessFunction(_N,_fitness_function,_bodyModel,_TheViewSet,DOFmask);  

    _num_steps=0;
    _areParams=true;
}
/*
*/
void DEHumanPartTracker::step(guscene::view::ViewSet *VSet)throw(GUException)
{
  _TheViewSet=VSet;
  //indicate to evaluators that frame analysis starts
  for (int i=0;i<omp_get_max_threads();i++) { 
	_evaluators[i]->startFrame(VSet,_num_steps);
    }   
  for(unsigned int h=0;h<_BODY_PARTS;h++){
    for (int i=0;i<omp_get_max_threads();i++) { 
	(dynamic_cast<HTFitnessFunction* >(_pffs[h][i]))->setViewSet(VSet);
      }
  }
  
   _prev_best_fitness=INFINITE;
    for(unsigned int h=0;h<_BODY_PARTS;h++){
	int n=_ini_end_body_parts[h*2+1]-_ini_end_body_parts[h*2]+1;
        DE_D40_MM alg(_pffs[h],n,_xbest+_ini_end_body_parts[h*2],_params.pop_size,_params.nLPCsols,_params.cr,_params.f,_params.min_bound,_params.max_bound,time(NULL));
        _fvalue = alg.run(_xbest+_ini_end_body_parts[h*2],_xbest+_ini_end_body_parts[h*2],_params.stop_maxeval/_BODY_PARTS);
    }
	      
      //indicate to evaluators that frame analysis is finished
      SkeletonPose bestEstimation=getCurrentEstimation();
      _bodyModel.setPose(_auxBodyModel,bestEstimation);
      for(int i=0;i<omp_get_max_threads();i++)
	_evaluators[i]->frameFinished(&_auxBodyModel);    
      
      _num_steps++;
}	

    
    
  }
  
}
