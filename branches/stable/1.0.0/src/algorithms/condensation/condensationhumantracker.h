#ifndef _gummocap_CONDENSATIONHUMANTRACKER_
#define _gummocap_CONDENSATIONHUMANTRACKER_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <guai/condensationparallel.h>
#include "poseevaluators.h"
#include "humantracker.h"
namespace gummocap
{
namespace algorithms
{
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace guai::condensation;
using namespace gummocap::poseevaluators;
class CondensationHumanTracker  : public guai::condensation::CondensationParallel, public HumanTracker
{

public:
    struct Params {
        Params(int nPart=500) {
            nParticles=nPart;
        }
        //

        int nParticles;
	
    };

    /**
    */
    CondensationHumanTracker ();
    /**
    */
    ~ CondensationHumanTracker();
    /**Sets the required params and initializes
    */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet,PoseEvaluator* poseEval, Params  p,vector<bool> *DOFmask=NULL);
    /**Performs an step of the algorithm
    */
    void step(guscene::view::ViewSet *VSet)throw(GUException);
    /**draws the indicates particle in the viewset
    */
    void drawParticle(unsigned int i,guimage::ImageSet *ImgSet);


    /**
     */
    SkeletonPose getCurrentEstimation(){return getMeanPose();}
    /**
     */
    void enableDebug(bool enable){}

  void saveProbHist(string file);
    
protected:
    /**draws the indicates particle in the viewset
    */
    SkeletonPose getMeanPose();
    /**Reimplemeted from CONDENSATION. Creates a  particle for this problem
    */
    Particle  *createNewParticle();
    /** Reimplemeted from CONDENSATION.
    * Updates the particle. according to the prior probability
    */
    void variateParticle(Particle &src,Particle &dst) ;
    /** Reimplemeted from CONDENSATION. Calculates the prob of a particle
    */
    float calculateParticleProb(Particle &p,unsigned int ThreadId) ;


    /** converts a particle in a full skeleton pose
    */
    void fromParticle2PoseWoT(Particle &particle,SkeletonPoseWoT &pose);
    /** converts a particle in a full skeleton pose
    */
    void fromParticle2Pose(Particle &particle,SkeletonPose &pose);


   
private:

    int _nStepsDone;
    PoseEvaluator* _poseEvaluator;
    vector<PoseEvaluator*> _evaluators;
    Params _params;
    BodyModel<HumanSkeleton> _bodyModel;
    vector<BodyModel<HumanSkeleton> > _auxBodyModel;
    bool _areParams;
    gumocap::skeleton::SkeletonPoseWoT _reducedPose;
    guscene::view::ViewSet *_TheViewSet;
};

};
};
#endif

