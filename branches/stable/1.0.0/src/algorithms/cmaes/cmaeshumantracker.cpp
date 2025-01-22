#include <gu/gurandom.h>
#include <gumocap/gupainter.h>
#include <gu/gutimemark.h>
#include "cmaeshumantracker.h"

namespace gummocap
{
namespace algorithms
{

/**
*/
CMAESHumanTracker::CMAESHumanTracker ()
{
    _areParams=false;
    _TheViewSet=NULL;
}
/**
*/
CMAESHumanTracker::~CMAESHumanTracker ()
{
    if (_areParams) {
        delete[] _xbest;
        delete[] _stdini;
	if (evo.state!=-1)
            cmaes_exit(&evo);	/* releases the memory */
        for (	unsigned int i=0;i<evaluators.size();i++)
            delete evaluators[i];

    }
}

/**Sets the required params and initializes
*/
void CMAESHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask)
{

    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose, DOFmask);
    _N=(_reducedPose.size())*3+3;
    _xbest=new double[_N];
    _stdini=new double[_N];
    for (unsigned int i=0;i<_N;i++) {
        _xbest[i]=0.0;
        _stdini[i]=p.std_ini;
    }
    //initXPos(_reducedPose); // Not necessary: all X are 0.
    arFunvals = cmaes_init(&evo, _N, _xbest, _stdini, 0, _params.lambda, _params.stop_maxiter, _params.stop_maxeval, NULL);
 
    _fitness_function=fitness_func;
    //Create independent evaluators , one for each thread
    assert(evaluators.size()==0);
    for (int i=0;i<omp_get_max_threads();i++) { 
        evaluators.push_back( _fitness_function->makeCopy());
        _auxBodyModelV.push_back( BodyModel<HumanSkeleton>(_bodyModel) );
    }
    _num_steps=0;
    _areParams=true;
}
/**
*
*
*/
void CMAESHumanTracker::step(guscene::view::ViewSet *VSet)throw(GUException)
{
    double *xfinal;
    double *const*pop;

    _TheViewSet=VSet;
    //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<evaluators.size();i++)
        evaluators[i]->startFrame(VSet,_num_steps);

    /* Initialize everything into the struct evo, 0 means default */
    //arFunvals = cmaes_init(&evo, 0, xbest, NULL, 0, 0, "initials.par");
    //printf("%s\n", cmaes_SayHello(&evo));
    //cmaes_ReadSignals(&evo, "signals.par");  /* write header and initial values */

    //performs an iteration
    //it means sampling and propagation, and then evaluation by calling the  calculateParticleProb function


    /* Iterate until stop criterion holds */
    while (!cmaes_TestForTermination(&evo))
    {
        /* generate lambda new search points, sample population */
        pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

        /* Here you may resample each solution point pop[i] until it
        becomes feasible, e.g. for box constraints (variable
        boundaries). function is_feasible(...) needs to be
        user-defined.
        Assumptions: the feasible domain is convex, the optimum is
        not on (or very close to) the domain boundary, initialX is
        feasible and initialStandardDeviations are sufficiently small
        to prevent quasi-infinite looping.
        */
        /*
        for (int i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
         while (!isFeasible(pop[i]))
           cmaes_ReSampleSingle(&evo, i);
        */
        /* evaluate the new search points using fitfun from above */
	int nIterations= cmaes_Get(&evo, "lambda");
	#pragma omp parallel for
        for (int i = 0; i <nIterations; ++i) {
            int threadId=omp_get_thread_num();
	    SkeletonPose pose;
            fromVector2Pose(pop[i],pose );
            _bodyModel.setPose(_auxBodyModelV[threadId],pose );
            arFunvals[i]=evaluators[threadId]->evaluate(&_auxBodyModelV[threadId],_TheViewSet);
        }

        /* update the search distribution used for cmaes_SampleDistribution() */
        cmaes_UpdateDistribution(&evo, arFunvals);

        /* read instructions for printing output or changing termination conditions */
        //cmaes_ReadSignals(&evo, "signals.par");
        //fflush(stdout); /* useful in MinGW */
    }

    //printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
    //cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

    /* get best estimator for the optimum, xmean */
    xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */

    for (unsigned int i=0;i<_N;i++) {
        _xbest[i]=xfinal[i];
    }

    free(xfinal);

    //cmaes_exit(&evo);
    arFunvals = cmaes_init_step(&evo, _N, _xbest, _stdini, 0, _params.lambda, NULL);
    
    //indicate to evaluators that frame analysis is finished
    SkeletonPose bestEstimation=getCurrentEstimation();
    _bodyModel.setPose(_auxBodyModelV[0],bestEstimation);
    for(unsigned int i=0;i<evaluators.size();i++)
      evaluators[i]->frameFinished(&_auxBodyModelV[0]);
    _num_steps++;

}
/**
*
*/
SkeletonPose CMAESHumanTracker::getCurrentEstimation()
{
    SkeletonPose pose;

    fromVector2Pose(_xbest,pose);

    return pose;

}
/**
*
*/
bool CMAESHumanTracker::isFeasible(double const *x) {
assert(false);

    return true;
}

/*
*
*
*/
bool CMAESHumanTracker::checkFeasibility(BodyModel<HumanSkeleton>& bm_sol,BodyModel<HumanSkeleton>& bm_org) {

    Point3Df spine_bottom=bm_sol.getJoint(HumanSkeleton::Spine_bottom).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Spine_bottom ).xyz;
    Point3Df spine_up=bm_sol.getJoint(HumanSkeleton::Spine_up).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Spine_up ).xyz;
    float body_radius=bm_org.getMesh(HumanSkeleton::Spine_bottom).width_radius*2;

    // DISTANCE BODY-ARMS (MUST BE GREATER THAN ARM RADIUS + BODY RADIUS)

    float arm_radius=bm_org.getMesh(HumanSkeleton::Lshoulder).width_radius*2;
    Point3Df lelbow=bm_sol.getJoint(HumanSkeleton::Lelbow).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Lelbow ).xyz;
    Point3Df midspine=(spine_up+spine_bottom)/2.;
    float larm_distance=midspine.distance(lelbow);
    if (larm_distance<(arm_radius+body_radius)) {
        //  cout << "LARM" << endl;
        return false;
    }
    //cout << "LARM DISTANCE:" << larm_distance << "THRESHOLD:" << arm_radius+spine_radius << endl;

    Point3Df relbow=bm_sol.getJoint(HumanSkeleton::Relbow).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Relbow ).xyz;
    float rarm_distance=midspine.distance(relbow);
    if (rarm_distance<(body_radius+arm_radius)) {
        //  cout << "RARM" << endl;
        return false;
    }

    // cout << "RARM DISTANCE:" << rarm_distance << "THRESHOLD:" << spine_radius+arm_radius << endl;
    Point3Df lwrist=bm_sol.getJoint(HumanSkeleton::Lwrist).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Lwrist ).xyz;
    float lwrist_distance=midspine.distance(lwrist);
    if (lwrist_distance<(body_radius+arm_radius)) {
        //  cout << "LWRIST" << endl;
        return false;
    }

    Point3Df rwrist=bm_sol.getJoint(HumanSkeleton::Rwrist).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Rwrist ).xyz;
    float rwrist_distance=midspine.distance(rwrist);
    if (rwrist_distance<(body_radius+arm_radius)) {
        //  cout << "RWRIST" << endl;
        return false;
    }

    // ARMS DISTANCE

    float arms_distance=relbow.distance(lelbow);
    if (arms_distance<(arm_radius*2)) {
        //  cout << "ARMS" << endl;
        return false;
    }

    //cout << "ARMS DISTANCE:" << arms_distance << "THRESHOLD:" << 2*arm_radius << endl;


    //KNEES DISTANCE

    Point3Df lknee=bm_sol.getJoint(HumanSkeleton::Lknee).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Lknee ).xyz;
    Point3Df rknee=bm_sol.getJoint(HumanSkeleton::Rknee).pose.Tmatrix*bm_sol.getJoint( HumanSkeleton::Rknee ).xyz;

    float leg_radius=bm_org.getMesh(HumanSkeleton::Lhip).width_radius;
    float knees_distance=rknee.distance(lknee);
    if (knees_distance<(2*leg_radius))
        return false;
    //cout << "KNEES DISTANCE:" << knees_distance << "THRESHOLD:" << 2*leg_radius << endl;

    return true;
}
/**
 * 
 */
void CMAESHumanTracker::setValidViews(vector<bool> vviews)throw(GUException)
{
  for (unsigned int i=0;i<evaluators.size();i++)
        evaluators[i]->setValidViews(vviews);
  
}

};
};

