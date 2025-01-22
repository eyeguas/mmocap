
#include "deht_fitnessfunction.h"

namespace gummocap{
  
  namespace algorithms{
    
    /**
    *
    *
    */
    HTFitnessFunction::HTFitnessFunction(unsigned int N, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton>& sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask)
    {
      _fitness_function=fitness_function;
      _bodyModel=sk;
      _TheViewSet=TheViewSet;
      SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
      gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
      _reducedDOFXYZmask=NULL;
      _auxBodyModel=BodyModel<HumanSkeleton>(_bodyModel);
      _N=N;
    }
    /**
    *
    *
    */
    HTFitnessFunction::HTFitnessFunction(unsigned int N, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton>& sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask,vector<bool> *DOFXYZmask)
    {
      
      _fitness_function=fitness_function;
      _bodyModel=sk;
      _TheViewSet=TheViewSet;
      SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
      gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
      _reducedDOFXYZmask=DOFXYZmask;
      _auxBodyModel=BodyModel<HumanSkeleton>(_bodyModel);
      _N=N;
    }
    /**
    *
    *
    */
    HTFitnessFunction::~HTFitnessFunction()
    {
    }
    /**
    *
    *
    */
    long double HTFitnessFunction::fitness(const double *vars)
    {
      numFEs++;
      SkeletonPose pose;
      fromVector2Pose(vars,pose );
      _bodyModel.setPose(_auxBodyModel,pose );
      return _fitness_function->evaluate(&_auxBodyModel,_TheViewSet);       
    }
    /**
    *
    *
    */
    long double HTFitnessFunction::fitness(const double *vars, int dimension)
    {
	return fitness(vars);
    }
    /**
    *
    *
    */
    int HTFitnessFunction::getDim()
    {
	return _N;
    }
        /**
     * In this class, the functions are to be minimised.
     */
    long double HTFitnessFunction::compare(long double f1, long double f2){

        return (f2 - f1);
    }
    /**
    *
    *
    */
    ostringstream* HTFitnessFunction::getName(){
        return NULL;
    }
    /**
    *
    *
    */
    void HTFitnessFunction::setViewSet(ViewSet* VS){
       _TheViewSet=VS;
    }
    /**
    *
    *
    */
    void HTFitnessFunction::startFrame(ViewSet* VS,unsigned int idFrame)
    {
	_fitness_function->startFrame(VS,idFrame);
    }
    /**
    *
    *
    */
    void HTFitnessFunction::frameFinished(BodyModel<HumanSkeleton>* bodymodel)
    {
	_fitness_function->frameFinished(bodymodel);
    }
	


  }
  
  
}