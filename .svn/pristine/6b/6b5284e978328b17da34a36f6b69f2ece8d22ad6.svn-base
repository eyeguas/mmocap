#include "macmachainshumantracker.h"
#include <omp.h>
namespace gummocap
{
namespace algorithms
{
/**
*/
MACMAChainsHumanTracker::MACMAChainsHumanTracker ()
{
    _areParams=false;
    _TheViewSet=NULL;
}
/**
*/
MACMAChainsHumanTracker::~MACMAChainsHumanTracker ()
{
    if (_areParams) {
        
      delete _xbest;
      /*
      for (	unsigned int i=0;i<evaluators.size();i++)
            delete evaluators[i];
      */
    }
}

/**Sets the required params and initializes
*/
void MACMAChainsHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask)
{

    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    _fitness_function=fitness_func;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose,DOFmask);
    _N=(_reducedPose.size())*3+3;
    
    _xbest=new tChromosomeReal(_N);
    for (unsigned int i=0;i<_N;i++) {
        (*_xbest)[i]=0.0;
    }
    
    //Create independent evaluators , one for each thread
    //assert(evaluators.size()==0);
    
    for (int i=0;i<omp_get_max_threads();i++) { 
        _evaluators.push_back( _fitness_function->makeCopy());
        _auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
    }
    _num_steps=0;
    _areParams=true;
}


void MACMAChainsHumanTracker::step(guscene::view::ViewSet *VSet)throw(GUException)
{

    _TheViewSet=VSet;
    //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<_evaluators.size();i++)
        _evaluators[i]->startFrame(VSet,_num_steps);

    
    int threadId=omp_get_thread_num();
    // Algorithm initialization
    Random random(new SRandom(time(NULL)));
    HumanTrackerProblemInstance htp(&random, _N);
    htp.init(_auxBodyModelV[threadId],_TheViewSet, _evaluators[threadId]);
    ProblemPtr phtp = htp.get(_params.max_eval);
    for (unsigned i = 0; i < _N; ++i) {
	    phtp->setDomainValues(i, (*_xbest)[i]-_params.min_domain, (*_xbest)[i]+_params.max_domain, true);
	}
     
    SSGA *ssga = new SSGA(&random);
     ssga->setCross(new CrossBLX(0.5));
     ssga->setMutation(new MutationBGA());
     ssga->setSelect(new SelectNAM(3));
     ssga->setReplacement(new ReplaceWorst());

     IEA *ea = ssga;
     CMAESHansen *cmaes = new CMAESHansen("cmaesinit.par");
     cmaes->searchNeighborhood(0.5);
     ILocalSearch *ls = cmaes;

     MALSChains *ma = new MALSChains(ea, ls);
    // ma->setDebug();
     ma->setRestart(new RestartBest());
     Hybrid *hybrid = ma;
     hybrid->setEffortRatio(_params.lg_ratio); // L/G Ratio (local/global search ratio)
     hybrid->setIntensity(_params.ls_intensity);  //Number of fitness function evaluations required by the LS algorithm during its operation
    	
     EA alg(hybrid, phtp);
    
    
    //tFitness sum=0;
    tFitness fitness=0;
    
	for (unsigned int num = 1; num <= _params.stop_maxiter; num++) {
	    alg.apply((*_xbest), &fitness);
	    //printf("fitness:%Le\n", fitness);
	    //sum += fitness;
	}
	//printf("mean:%Le\n", sum/_params.stop_maxiter);
/*
	for (unsigned int i=0;i<_xbest->size();i++){
	  cout << (*_xbest)[i] << " ";
	}
	cout << endl;
*/	
    //performs an iteration
    
    /* Iterate until stop criterion holds */

    //indicate to evaluators that frame analysis is finished
    
    SkeletonPose bestEstimation=getCurrentEstimation();
    _bodyModel.setPose(_auxBodyModel,bestEstimation);
    for(unsigned int i=0;i<_evaluators.size();i++)
      _evaluators[i]->frameFinished(&_auxBodyModel);
    
    _num_steps++;
}
/**
*
*/
SkeletonPose MACMAChainsHumanTracker::getCurrentEstimation()
{
    SkeletonPose pose;

    fromVector2Pose((*_xbest),pose);

    return pose;

}

}
}

