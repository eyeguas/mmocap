#ifndef _gummocap_CMAESHUMANTRACKER_
#define _gummocap_CMAESHUMANTRACKER_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <opencv/highgui.h>
#include <omp.h>
#include "humantracker.h"
extern "C" {
#include "cmaes_interface.h"
}
#include "poseevaluators.h"

namespace gummocap
{
namespace algorithms {
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;

class CMAESHumanTracker: public HumanTracker
            /**\brief This class represents a human tracker using the CMAES algorithm.
            */
{

public:
    struct Params {
        Params(int lamb, double stddev_ini, double maxiter, double maxeval) {
            std_ini=stddev_ini;
            stop_maxiter=maxiter;
	    stop_maxeval=maxeval;
            lambda=lamb;
        }
        Params() {
            std_ini=0.0;
            stop_maxiter=1;
	    stop_maxeval=1;
            lambda=0;
        }
        double std_ini;		// Standard deviation for each point (parameter: rotation or translation) in the search
        double stop_maxiter;	// Maximum number of iterations for each step
        double stop_maxeval;	// Maximum number of evaluations for each step
        int lambda;
    };

    /**
    */
    CMAESHumanTracker ();
    /**Releases the memory
    */
    ~CMAESHumanTracker ();
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func, Params  p, vector<bool> *DOFmask);
    /**Performs an step of the algorithm: a CMAES search.
    */
    void step(guscene::view::ViewSet *VSet)throw(GUException);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();
    /**
     */
    void enableDebug(bool enable){}
    /**Checks if the body model solution is feasible.
    * bm_sol: body model solution
    * bm_org: body model original
    */
    bool checkFeasibility(BodyModel<HumanSkeleton>& bm_sol,BodyModel<HumanSkeleton>& bm_org);
    bool isFeasible(double const *x);

protected:

    /** converts a vector in a reduced skeleton pose
        */

    void fromVector2PoseWoT(double const *x, SkeletonPoseWoT &pose) {
        assert(_areParams);
        pose=_reducedPose;

        pose.Translation[0]=x[0];
        pose.Translation[1]=x[1];
        pose.Translation[2]=x[2];
        int index=3;
        // pose[0]=root
        for (unsigned int i=0;i<pose.size();i++) {
            pose[i].x=x[index++];
            pose[i].y=x[index++];
            pose[i].z=x[index++];
        }
    }

    /** converts a vector in a full skeleton pose
    */

    void fromVector2Pose(double const *x, SkeletonPose &pose)
    {
        SkeletonPoseWoT redPose;
        fromVector2PoseWoT(x,redPose);
        SkeletonPoseWoT::fromWoT2Complete(redPose,pose);

    }

    /** computes the feasibility distance
    * Returns true if the feasibility distance is greater or equal to the threshold.
    * Returns false otherwise.
    */

    bool distanceFeasible(Point3Df& p1, Point3Df& s1, Point3Df& s2, float threshold) {
        int t;
        float d=p1.distancePerpendicularPointSegment(s1,s2,t);
        //cout << p1.xyzh[0] <<" "<<p1.xyzh[1]<<" "<<p1.xyzh[2]<<" "<<t << endl;
        //cout << d << " "<< threshold << endl;
        if (t==0) {
            //cout << d << " "<< threshold << endl;
            if (d<(threshold)) {
                return false;
            }
        }
        return true;
    }
    /**Set a mask to select these views to be employed
      */
    void setValidViews(vector<bool> vviews)throw(GUException);
protected:
    Params _params;
    BodyModel<HumanSkeleton> _bodyModel,_auxBodyModel;
    vector<BodyModel<HumanSkeleton> > _auxBodyModelV;
    bool _areParams;
    gumocap::skeleton::SkeletonPoseWoT _reducedPose;
    guscene::view::ViewSet *_TheViewSet;


    unsigned int _num_steps;
    cmaes_t evo;
    double* arFunvals;
    double* _xbest;
    double* _stdini;
    unsigned int _N;
    PoseEvaluator* _fitness_function;
    vector<PoseEvaluator*> evaluators;



};
};
};
#endif

