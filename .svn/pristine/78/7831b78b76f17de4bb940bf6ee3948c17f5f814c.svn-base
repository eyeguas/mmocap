
#include "dehpt_fitnessfunction.h"

namespace gummocap{
  
  namespace algorithms{
 
     
    /**
    *
    *
    */
    HPTFitnessFunction::HPTFitnessFunction(double* xbest, unsigned int N, unsigned int inibodypart, unsigned int endbodypart, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton>& sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask) :
      HTFitnessFunction(N,fitness_function,sk,TheViewSet,DOFmask)
    {
	_inibodypart=inibodypart;
	_endbodypart=endbodypart;
	_xbest=xbest;
    }
     /**
    *
    *
    */
    HPTFitnessFunction::HPTFitnessFunction(double* xbest, unsigned int N, unsigned int inibodypart, unsigned int endbodypart, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton>& sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask,vector<bool> *DOFXYZmask) :
      HTFitnessFunction(N,fitness_function,sk,TheViewSet,DOFmask,DOFXYZmask)
    {
	_inibodypart=inibodypart;
	_endbodypart=endbodypart;
	_xbest=xbest;
    }
    /**
    *
    *
    */
    HPTFitnessFunction::~HPTFitnessFunction()
    {
    }
    /**
    *
    *
    */
    long double HPTFitnessFunction::fitness(const double *vars)
    {
      numFEs++;
      SkeletonPose pose;
      fromPartVector2Pose(vars,_inibodypart,_endbodypart,pose );
      _bodyModel.setPose(_auxBodyModel,pose );
      return _fitness_function->evaluate(&_auxBodyModel,_TheViewSet);       
    }
    /**
    *
    *
    */
    long double HPTFitnessFunction::fitness(const double *vars, int dimension)
    {
	return fitness(vars);
    }
    /**
    *
    *
    */
    int HPTFitnessFunction::getDim()
    {
	return _endbodypart-_inibodypart+1;
    }
    
    
  }
  
}