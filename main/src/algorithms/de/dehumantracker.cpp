
#include "dehumantracker.h"

namespace gummocap{
  
  namespace algorithms{
    
    const double DEHumanTracker::PI = acos(-1.0);
    const double DEHumanTracker::E = exp(1.0);
    
/**
*/
DEHumanTracker::DEHumanTracker ()
{
    _areParams=false;
}
/**
*/
DEHumanTracker::~DEHumanTracker ()
{
    if (_areParams) {
        delete[] _xbest;
        for (unsigned int i=0;i<_evaluators.size();i++){
            delete _evaluators[i];
	    delete _ffs[i];
	}	
	delete _ff;
    }
}
/**Sets the required params and initializes
*/
void DEHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask)
{

    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
    _N=(_reducedPose.size())*3+3;
    _xbest=new double[_N];
    for (unsigned int i=0;i<_N;i++) {
        _xbest[i]=0.0;
    }
 
    _fitness_function=fitness_func;
    //Create independent evaluators , one for each thread
    assert(_evaluators.size()==0);
    for (int i=0;i<omp_get_max_threads();i++) { 
	_auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
	_evaluators.push_back( _fitness_function->makeCopy());
        _ffs.push_back(new HTFitnessFunction(_N,_evaluators[i],_auxBodyModelV[i],_TheViewSet,DOFmask));
    }
    
    _ff = new HTFitnessFunction(_N,_fitness_function,_bodyModel,_TheViewSet,DOFmask);  
        
     
    _num_steps=0;
    _areParams=true;
}
/**Sets the required params and initializes
*/
void DEHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, vector<bool> *DOFXYZmask)
{

    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
    //_N=gumocap::skeleton::SkeletonPoseWoT::nDOF(*sk.getSkptr(), DOFmask, DOFXYZmask);
    _reducedDOFXYZmask=new   vector<bool>();
    gumocap::skeleton::SkeletonPoseWoT::toReducedMaskXYZ(*sk.getSkptr(), DOFmask, DOFXYZmask, _reducedDOFXYZmask);
    //cout << "N" << _N << endl;
    //cout << "REDUCED POSE" << _reducedPose.size() << "REDUCED DOFXYZMASK" << _reducedDOFXYZmask->size();
    _N=(_reducedPose.size())*3+3;
    
    _xbest=new double[_N];
    for (unsigned int i=0;i<_N;i++) {
        _xbest[i]=0.0;
    }
 
    _fitness_function=fitness_func;
    //Create independent evaluators , one for each thread
    assert(_evaluators.size()==0);
    for (int i=0;i<omp_get_max_threads();i++) { 
	_auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
	_evaluators.push_back( _fitness_function->makeCopy());
        //_ffs.push_back(new HTFitnessFunction(_N,_evaluators[i],_auxBodyModelV[i],_TheViewSet,DOFmask));
	_ffs.push_back(new HTFitnessFunction(_N,_evaluators[i],_auxBodyModelV[i],_TheViewSet,DOFmask,_reducedDOFXYZmask));
    }
    
    //_ff = new HTFitnessFunction(_N,_fitness_function,_bodyModel,_TheViewSet,DOFmask);
    _ff = new HTFitnessFunction(_N,_fitness_function,_bodyModel,_TheViewSet,DOFmask,_reducedDOFXYZmask);  
        
     
    _num_steps=0;
    _areParams=true;
}
/**
*
*
*/
void DEHumanTracker::step(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
  
      _TheViewSet=VSet;

      for (int i=0;i<omp_get_max_threads();i++) { 
	(dynamic_cast<HTFitnessFunction* >(_ffs[i]))->setViewSet(VSet);
	(dynamic_cast<HTFitnessFunction* >(_ffs[i]))->startFrame(VSet,_num_steps);
      }
      
      DE_D40_MM alg(_ffs,_N,_xbest,_params.pop_size,_params.nLPCsols,_params.cr,_params.f,_params.min_bound,_params.max_bound,time(NULL));
      
      _fvalue = alg.run(_xbest,_xbest,_params.stop_maxeval);

      //indicate to evaluators that frame analysis is finished
      SkeletonPose bestEstimation=getCurrentEstimation();
      _bodyModel.setPose(_auxBodyModel,bestEstimation);
      for(int i=0;i<omp_get_max_threads();i++)
	(dynamic_cast<HTFitnessFunction* >(_ffs[i]))->frameFinished(&_auxBodyModel);    
      
      _num_steps++;

  
}
/**
*
*/
SkeletonPose DEHumanTracker::getCurrentEstimation()
{
    SkeletonPose pose;

    _ff->fromVector2Pose(_xbest,pose);

    return pose;

}

    
  }
  
}
