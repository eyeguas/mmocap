#include <gu/gurandom.h>
#include <gumocap/gupainter.h>
#include <gu/gutimemark.h>
#include "cmaeshumanparttracker.h"

namespace gummocap
{
namespace algorithms
{
/**
*/
CMAESHumanPartTracker::CMAESHumanPartTracker ()
{
}
/**
*/
CMAESHumanPartTracker::~CMAESHumanPartTracker ()
{
  if(_areParams){
    }
}

/**Sets the required params and initializes
*/
void CMAESHumanPartTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask, unsigned int nBodyParts)
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
    _stdini=new double[_N];
    for(unsigned int i=_initrans;i<=_endtrans;i++){
      _xbest[i]=0.0;
      _stdini[i]=p.std_trans;
    }
    for(unsigned int i=_inirot;i<=_endrot;i++){
      _xbest[i]=0.0;
      _stdini[i]=p.std_rot;
    }
    for(unsigned int i=_inibody;i<=_endbody;i++){
      _xbest[i]=0.0;
      _stdini[i]=p.std_body;
    }
    for(unsigned int i=_inilarm;i<=_endrarm;i++){
      _xbest[i]=0.0;
      _stdini[i]=p.std_arms;
    }
    for(unsigned int i=_inilleg;i<=_endrleg;i++){
      _xbest[i]=0.0;
      _stdini[i]=p.std_legs;
    }
    _prev_best_fitness=INFINITE;
    _fitness_function=fitness_func;
    //Create independent evaluators , one for each limb
    assert(evaluators.size()==0);
    
    for (int i=0;i<omp_get_max_threads();i++) { 
        evaluators.push_back( _fitness_function->makeCopy());
        _auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
    }
    
    _num_steps=0;
    _areParams=true;
}

/*
*/
void CMAESHumanPartTracker::step(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
  _TheViewSet=VSet;
  //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<evaluators.size();i++)
        evaluators[i]->startFrame(VSet,_num_steps);

   _prev_best_fitness=INFINITE;
    for(unsigned int h=0;h<_BODY_PARTS;h++)
	stepPart(_xbest,_xbest,&_prev_best_fitness,_ini_end_body_parts[h*2],_ini_end_body_parts[h*2+1],_params.stop_maxiter/_BODY_PARTS,_params.stop_maxeval/_BODY_PARTS);
	      
    //indicate to evaluators that frame analysis is finished
    SkeletonPose bestEstimation=getCurrentEstimation();
    _bodyModel.setPose(_auxBodyModel,bestEstimation);
    for(unsigned int i=0;i<evaluators.size();i++)
      evaluators[i]->frameFinished(&_auxBodyModel);
    
    _num_steps++;
}	
/*
*/
/*
*/
/*
void CMAESHumanPartTracker::step2(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
  _TheViewSet=VSet;
  unsigned int inilegs;
  unsigned int order12[4]={8,9,10,11};
  unsigned int order6[2]={5,4};
  unsigned int* order;
  
  if(_BODY_PARTS==6){
    inilegs=4;
    order=order6;
  }
  else{
    inilegs=8;
    order=order12;
  }
  
  double* auxbest=new double[_N];
  double auxbestf;
  
  //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<evaluators.size();i++)
        evaluators[i]->startFrame(VSet,_num_steps);

   _prev_best_fitness=INFINITE;
    for(unsigned int h=0;h<inilegs;h++)
	stepPart(_xbest,_xbest,&_prev_best_fitness,_ini_end_body_parts[h*2],_ini_end_body_parts[h*2+1],_params.stop_maxiter/_BODY_PARTS,_params.stop_maxeval/_BODY_PARTS);
    for(unsigned int i=0;i<_N;i++)
      auxbest[i]=_xbest[i];
    auxbestf=_prev_best_fitness;
    for(unsigned int h=inilegs;h<_BODY_PARTS;h++)
	stepPart(_xbest,_xbest,&_prev_best_fitness,_ini_end_body_parts[h*2],_ini_end_body_parts[h*2+1],_params.stop_maxiter/_BODY_PARTS,_params.stop_maxeval/_BODY_PARTS);

    for(unsigned int i=0;i<_BODY_PARTS-inilegs;i++)
	stepPart(auxbest,auxbest,&auxbestf,_ini_end_body_parts[order[i]*2],_ini_end_body_parts[order[i]*2+1],_params.stop_maxiter/_BODY_PARTS,_params.stop_maxeval/_BODY_PARTS);
    
    if(auxbestf<_prev_best_fitness){
      _prev_best_fitness=auxbestf;
      for(unsigned int i=0;i<_N;i++)
	_xbest[i]=auxbest[i];
    }
    
    delete [] auxbest;
   
    //indicate to evaluators that frame analysis is finished
    SkeletonPose bestEstimation=getCurrentEstimation();
    _bodyModel.setPose(_auxBodyModel,bestEstimation);
    for(unsigned int i=0;i<evaluators.size();i++)
      evaluators[i]->frameFinished(&_auxBodyModel);
    
    _num_steps++;
}
*/
/*
*/

/**
*
*/

void CMAESHumanPartTracker::stepPart(double* x, double* xbest, double* best_fitness, unsigned int inibodypart, unsigned int endbodypart, double niterations, double nevals)throw(gu::Exception)
{
    double *xfinal;
    double *const*pop;
          
    //performs an iteration
    //it means sampling and propagation, and then evaluation by calling the  calculateParticleProb function

    arFunvals = cmaes_init(&evo, endbodypart-inibodypart+1, x+inibodypart, _stdini+inibodypart, 0, _params.lambda, niterations, nevals, NULL); 
    
    /* Iterate until stop criterion holds */
    while(!cmaes_TestForTermination(&evo))
    { 
      
      /* generate lambda new search points, sample population */
      pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

      /*
      if(checkfeasibility)
      for (int i = 0; i < cmaes_Get(&evo, "popsize"); ++i){ 
	for(unsigned int k=inibodypart,j=0;k<=endbodypart;k++,j++){
	    _xbest[k]=pop[i][j];
	}
	while (!isFeasible(_xbest)) 
	    cmaes_ReSampleSingle(&evo, i); 
      }
      */
      int Lambda=cmaes_Get(&evo, "lambda");
      /* evaluate the new search points using fitfun from above */ 
      #pragma omp parallel for
      for (int i = 0; i < Lambda; ++i) {
	 int threadId=omp_get_thread_num();
	 SkeletonPose pose;
	 fromPartVector2Pose(pop[i],inibodypart,endbodypart,pose);
	 //_bodyModel.setPose(auxBodyModel,pose);
	 //arFunvals[i]=_fitness_function->evaluate(&auxBodyModel,_TheViewSet);
	 _bodyModel.setPose(_auxBodyModelV[threadId],pose );
         arFunvals[i]=evaluators[threadId]->evaluate(&_auxBodyModelV[threadId],_TheViewSet);
	 
      }

      /* update the search distribution used for cmaes_SampleDistribution() */
      cmaes_UpdateDistribution(&evo, arFunvals);  
    }
    
    /* get best estimator for the optimum, xmean */
    xfinal = cmaes_GetNew(&evo, "xbestever"); /* "xbestever" might be used as well */
    double xfinal_fitness = cmaes_Get(&evo, "fbestever");
    
    if(xfinal_fitness < (*best_fitness)){
      *best_fitness=xfinal_fitness;
      for(unsigned int i=inibodypart,j=0;i<=endbodypart;i++,j++){
	xbest[i]=xfinal[j];
      }
    }
    /*
    else{
            *best_fitness=_prev_best_fitness;
      for(unsigned int i=inibodypart;i<=endbodypart;i++){
	xbest[i]=_xprev_best[i];
      }
    
    }*/
    /**/
    
    
    free(xfinal);
    
    cmaes_exit(&evo);	/* releases the memory */
    
}
/**
*
*/
SkeletonPose CMAESHumanPartTracker::getCurrentEstimation()
{
    SkeletonPose pose;
    
    fromPartVector2Pose(_xbest,0,_N-1,pose);
    
    return pose;

}

}
}

