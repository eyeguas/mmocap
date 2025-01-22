#include <gu/gurandom.h>
#include <gumocap/gupainter.h>
#include <gu/gutimemark.h>
#include "macmachainshumanparttracker.h"

namespace gummocap
{
namespace algorithms
{

/**
*/
MACMAChainsHumanPartTracker::MACMAChainsHumanPartTracker ()
{
}
/**
*/

MACMAChainsHumanPartTracker::~MACMAChainsHumanPartTracker ()
{
  if(_areParams){
  }
}

/**Sets the required params and initializes
*/
void MACMAChainsHumanPartTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask,unsigned int nBodyParts)
{
    
    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose,DOFmask);
    _N=(_reducedPose.size())*3+3;
    assert((unsigned int)_N==3*(1+_params.nbody_joints+_params.nlarm_joints+_params.nrarm_joints+_params.nlleg_joints+_params.nrleg_joints)+3);
    
    
    HumanPartTracker::setParams(nBodyParts,_params.nbody_joints,_params.nlarm_joints,_params.nrarm_joints,_params.nlleg_joints,_params.nrleg_joints);
    
    _xbest=new tChromosomeReal(_N);
    
    for(unsigned int i=0;i<_N;i++){
      (*_xbest)[i]=0.0;
    }

    _prev_best_fitness=INFINITE;
    _fitness_function=fitness_func;
    for (int i=0;i<omp_get_max_threads();i++) { 
        _evaluators.push_back( _fitness_function->makeCopy());
        _auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
    }
    
    _num_steps=0;
    _areParams=true;
}
/**
*
*/
void MACMAChainsHumanPartTracker::step(guscene::view::ViewSet *VSet)throw(gu::Exception){
  _TheViewSet=VSet;
  
  //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<_evaluators.size();i++)
        _evaluators[i]->startFrame(VSet,_num_steps);
  
    _prev_best_fitness=INFINITE;
    
    for(unsigned int h=0;h<_BODY_PARTS;h++)
	stepPart(_xbest,_xbest,&_prev_best_fitness,_ini_end_body_parts[h*2],_ini_end_body_parts[h*2+1],_params.stop_maxiter,_params.max_eval/_BODY_PARTS);

  
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

void MACMAChainsHumanPartTracker::stepPart(tChromosomeReal* x, tChromosomeReal* xbest, double* best_fitness, unsigned int inibodypart, unsigned int endbodypart, double niterations, double maxeval)throw(gu::Exception)
{
    // Algorithm initialization
    
    Random random(new SRandom(time(NULL)));
    
    HumanPartTrackerProblemInstance htp(&random, endbodypart-inibodypart+1);
    
    int threadId=omp_get_thread_num();
        
    htp.init(x,inibodypart,endbodypart,_auxBodyModelV[threadId],_TheViewSet, _evaluators[threadId]);
    
    ProblemPtr phtp = htp.get(maxeval);
    for (unsigned i = inibodypart, j=0; i <=endbodypart; ++i,j++) {
	    phtp->setDomainValues(j, (*x)[i]-_params.min_domain, (*x)[i]+_params.max_domain, true);
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
     //ma->setDebug();
     ma->setRestart(new RestartBest());
     Hybrid *hybrid = ma;
     hybrid->setEffortRatio(_params.lg_ratio); // L/G Ratio (local/global search ratio)
     hybrid->setIntensity(_params.ls_intensity);  //Number of fitness function evaluations required by the LS algorithm during its operation
    	
     EA alg(hybrid, phtp);
    
    tFitness fitness=0;
    tChromosomeReal sol(endbodypart-inibodypart+1);
    for(unsigned int i=inibodypart,j=0;i<=endbodypart;i++,j++){
	  sol[j]=(*x)[i];
	}
	for (unsigned int num = 1; num <= niterations; num++) {
	    alg.apply(sol, &fitness);
	}
	
     if(fitness < _prev_best_fitness){
      *best_fitness=fitness;
      for(unsigned int i=inibodypart,j=0;i<=endbodypart;i++,j++){
	(*xbest)[i]=sol[j];
      }
     }
     /*
     else{
       *best_fitness=_prev_best_fitness;
       for(unsigned int i=inibodypart;i<=endbodypart;i++){
	(*xbest)[i]=(*_xprev_best)[i];
       }
     }
     */
     
     
}
/**
*
*/
SkeletonPose MACMAChainsHumanPartTracker::getCurrentEstimation()
{
    SkeletonPose pose;
    
    fromPartVector2Pose((*_xbest),0,_N-1,pose);
    
    return pose;

}

}
}

/*
void MACMAChainsHumanPartTracker::stepParallel(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
    
    
    _prev_best_fitness=INFINITE;
    _xprev_best=_xbest;
    
  //   cout << "TRANS" << endl;
    stepPart(_xbest,_xbest,&_prev_best_fitness,_initrans,_endtrans,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
   //  cout << "BODY" << endl;
    stepPart(_xbest,_xbest,&_prev_best_fitness,_inibody,_endbody,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
  //  cout << "RARM" << endl;
    _parts_best_fitness[BODY]=_prev_best_fitness;
    _xprev_best=&_parts_best[BODY];
    for(unsigned int i=0;i<_N;i++)
      _parts_best[BODY][i]=(*_xbest)[i];

    stepPart(_xbest,_xbest,&_parts_best_fitness[RARM],_inirarm,_endrarm,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
  //  cout << "LARM" << endl;
    stepPart(_xbest,_xbest,&_parts_best_fitness[RARM],_inilarm,_endlarm,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
   //     cout << "RLEG" << endl;
    for(unsigned int i=0;i<_N;i++){
      _parts_best[RARM][i]=(*_xbest)[i];
      (*_xbest)[i]=_parts_best[BODY][i];
    }

    stepPart(_xbest,_xbest,&_parts_best_fitness[LARM],_inilarm,_endlarm,_params.stop_maxiter,_params.max_eval/10);
    stepPart(_xbest,_xbest,&_parts_best_fitness[LARM],_inirarm,_endrarm,_params.stop_maxiter,_params.max_eval/10);

    unsigned int best_arms_search=LARM;
    if(_parts_best_fitness[RARM]<_parts_best_fitness[LARM]){
	best_arms_search=RARM;
//	cout << "RL-ARMs better than LR-ARMs" << endl;
    }
//    else
//	cout << "LR-ARMs better than RL-ARMs" << endl;  
    for(unsigned int i=0;i<_N;i++){
      _parts_best[LARM][i]=(*_xbest)[i];
      if(best_arms_search==RARM)
	(*_xbest)[i]=_parts_best[RARM][i];
    }
    _prev_best_fitness=_parts_best_fitness[best_arms_search];
    _xprev_best=&_parts_best[best_arms_search];

    stepPart(_xbest,_xbest,&_parts_best_fitness[RLEG],_inirleg,_endrleg,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
 //   cout << "LLEG" << endl;
    stepPart(_xbest,_xbest,&_parts_best_fitness[RLEG],_inilleg,_endlleg,_params.stop_maxiter,_params.max_eval/10);
//    for(int i=0;i<_N;i++){
//      cout << _xbest[i] << " " ;
//    }
//    cout << endl;
    for(unsigned int i=0;i<_N;i++){
      _parts_best[RLEG][i]=(*_xbest)[i];
      (*_xbest)[i]=_parts_best[best_arms_search][i];
    }

    stepPart(_xbest,_xbest,&_parts_best_fitness[LLEG],_inilleg,_endlleg,_params.stop_maxiter,_params.max_eval/10);
    stepPart(_xbest,_xbest,&_parts_best_fitness[LLEG],_inirleg,_endrleg,_params.stop_maxiter,_params.max_eval/10);

    unsigned int best_legs_search=LLEG;
    if(_parts_best_fitness[RLEG]<_parts_best_fitness[LLEG]){
	best_legs_search=RLEG;
//	cout << "RL-LEGs better than LR-LEGs" << endl;
    }
//    else
//        cout << "LR-LEGs better than RL-LEGs" << endl;
    for(unsigned int i=0;i<_N;i++){
      _parts_best[LLEG][i]=(*_xbest)[i];
      if(best_legs_search==RLEG)
	(*_xbest)[i]=_parts_best[RLEG][i];
    }

    _prev_best_fitness=_parts_best_fitness[best_legs_search];
    _xprev_best=&_parts_best[best_legs_search];
}
*/

/**
*
*/
/*
void MACMAChainsHumanPartTracker::stepLimbs(unsigned int first_limb, unsigned int second_limb, tChromosomeReal xlimbs_in, tChromosomeReal xlimbs_out, unsigned int ini_first_limb, unsigned int end_first_limb, unsigned int ini_second_limb, unsigned int end_second_limb, double niterations, double numevals)throw(gu::Exception)
{
      stepPart(&xlimbs_in,&xlimbs_out,&_parts_best_fitness[first_limb],ini_first_limb,end_first_limb,niterations,numevals/2);
      stepPart(&xlimbs_out,&xlimbs_out,&_parts_best_fitness[first_limb],ini_second_limb,end_second_limb,niterations,numevals/2);
}
*/
/**
*
*/
/*
void MACMAChainsHumanPartTracker::stepParallelLimbs(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
    
    _prev_best_fitness=INFINITE;
    
    stepPart(_xbest,&(_parts_best[BODY]),&_parts_best_fitness[BODY],_initrans,_endtrans,_params.stop_maxiter,_params.max_eval/11);
    stepPart(&(_parts_best[BODY]),&(_parts_best[BODY]),&_parts_best_fitness[BODY],_inibody,_endbody,_params.stop_maxiter,_params.max_eval/11);

    _prev_best_fitness=_parts_best_fitness[BODY];
    _xprev_best=&_parts_best[BODY];
    for(unsigned int i=0;i<_N;i++){
      if((i>=_initrans)&&(i<=_endbody))
	(*_xbest)[i]=_parts_best[RLEG][i]=_parts_best[LLEG][i]=_parts_best[RARM][i]=_parts_best[LARM][i]=_parts_best[BODY][i];
      else
	_parts_best[RLEG][i]=_parts_best[LLEG][i]=_parts_best[RARM][i]=_parts_best[LARM][i]=(*_xbest)[i];
    }
    
    int next_limb[4];
    next_limb[0]=next_limb[2]=1;
    next_limb[1]=next_limb[3]=-1;
    int ini_next_limb[4];
    ini_next_limb[0]=ini_next_limb[2]=2;
    ini_next_limb[1]=ini_next_limb[3]=-2;
    #pragma omp parallel for
    for(int k=0;k<4;k++){
      unsigned int first_limb=LARM+k;
      unsigned int ini_first_limb=INILARM+(k*2);
      unsigned int end_first_limb=ENDLARM+(k*2);
      stepLimbs(first_limb,first_limb+next_limb[k],_parts_best[BODY],_parts_best[first_limb],_ini_end_body_parts[ini_first_limb],_ini_end_body_parts[end_first_limb],_ini_end_body_parts[ini_first_limb+ini_next_limb[k]],_ini_end_body_parts[end_first_limb+ini_next_limb[k]],_params.stop_maxiter,2*_params.max_eval/11);      
    }
    
    unsigned int best_arms_search=LARM;
    if(_parts_best_fitness[RARM]<_parts_best_fitness[LARM]){
	best_arms_search=RARM;
	//cout << "RL-ARMs better than LR-ARMs" << endl;
    }
    //else
	//cout << "LR-ARMs better than RL-ARMs" << endl;  
    for(unsigned int i=_inilarm;i<=_endrarm;i++){
      (*_xbest)[i]=_parts_best[best_arms_search][i];
    }
    

    unsigned int best_legs_search=LLEG;
    if(_parts_best_fitness[RLEG]<_parts_best_fitness[LLEG]){
	best_legs_search=RLEG;
	//cout << "RL-LEGs better than LR-LEGs" << endl;
    }
    //else
        //cout << "LR-LEGs better than RL-LEGs" << endl;
    for(unsigned int i=_inilleg;i<=_endrleg;i++){
      (*_xbest)[i]=_parts_best[best_legs_search][i];
    }

    stepPart(_xbest,_xbest,&_prev_best_fitness,_initrans,_endrleg,_params.stop_maxiter,_params.max_eval/11);
      
      _xprev_best=_xbest;
    
}
*/
