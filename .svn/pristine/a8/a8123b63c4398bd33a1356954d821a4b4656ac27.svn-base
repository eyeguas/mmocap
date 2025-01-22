#ifndef _gummocap_MACMACHAINSHUMANTRACKER_
#define _gummocap_MACMACHAINSHUMANTRACKER_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <opencv/highgui.h>
#include <omp.h>
#include "humantracker.h"
#include "poseevaluators.h"

#include <time.h>
#include <gurealea/realea/problemcec2005.h>
#include <gurealea/realea/ilocalsearch.h>
#include <gurealea/realea/cross.h>
#include <gurealea/realea/ea.h>
#include <gurealea/macmachains/ssga.h>
#include <gurealea/macmachains/cmaeshan.h>
#include <gurealea/realea/srandom.h>
#include <gurealea/macmachains/malschains.h>
#include "humantrackerprobleminstance.h"

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;
using namespace realea;
using namespace gummocap::humantracker;
using namespace std;
using std::auto_ptr;


namespace gummocap
{
namespace algorithms {



class MACMAChainsHumanTracker: public HumanTracker
            /**\brief This class represents a human tracker using the MA-CMA-Chains algorithm.
            */
{

public:
    struct Params {
        Params(double maxiter, double max_ev, double min_dom, double max_dom, double ls_int, double lg_rat) {
            stop_maxiter=maxiter;
	    min_domain=min_dom;
	    max_domain=max_dom;
	    ls_intensity=ls_int;
	    lg_ratio=lg_rat;
	    max_eval=max_ev;
        }
        Params() {
            stop_maxiter=1;
	    min_domain=0;
	    max_domain=0;
	    ls_intensity=500;
	    lg_ratio=0.5;
	    max_eval=10000;
        }
        double stop_maxiter;	// Maximum number of iterations for each step
        double min_domain;	// Search domain: [xbest-min_domain,xbest+max_domain]
	double max_domain;
	double max_eval;	// Maximum number of evaluations per iteration = max_eval * dimensionality of the problem
	double ls_intensity;	// Local search intensity in number of evaluations
	double lg_ratio;	// Local search / Global search ratio
    };

    /**
    */
    MACMAChainsHumanTracker ();
    /**Releases the memory
    */
    ~MACMAChainsHumanTracker ();
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p,vector<bool> *DOFmask);
    /**Performs an step of the algorithm: a MA-CMA-Chains algorithm search.
    */
    void step(guscene::view::ViewSet *VSet)throw(GUException);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();
    /**
     */
    void enableDebug(bool enable){}
    

protected:
    Params _params;
    BodyModel<HumanSkeleton> _bodyModel,_auxBodyModel;
    vector<BodyModel<HumanSkeleton> > _auxBodyModelV;
    bool _areParams;
    gumocap::skeleton::SkeletonPoseWoT _reducedPose;
    guscene::view::ViewSet *_TheViewSet;

    tChromosomeReal* _xbest;
    
    unsigned int _N;
    PoseEvaluator* _fitness_function;
    vector<PoseEvaluator*> _evaluators;
    unsigned int _num_steps;

    void fromVector2PoseWoT(const tChromosomeReal &x, SkeletonPoseWoT &pose) {
        assert(_areParams);
        pose=_reducedPose;

        pose.Translation[0]=x[0];
        pose.Translation[1]=x[1];
        pose.Translation[2]=x[2];
        int index=3;
        //int index=0;
        // pose[0]=root
        for (unsigned int i=0;i<pose.size();i++) {
            pose[i].x=x[index++];
            pose[i].y=x[index++];
            pose[i].z=x[index++];
        }
    }

    /** converts a vector in a full skeleton pose
    */

    void fromVector2Pose(const tChromosomeReal &x, SkeletonPose &pose)
    {
        SkeletonPoseWoT redPose;
        fromVector2PoseWoT(x,redPose);
        SkeletonPoseWoT::fromWoT2Complete(redPose,pose);

    }


};
}
}
#endif

