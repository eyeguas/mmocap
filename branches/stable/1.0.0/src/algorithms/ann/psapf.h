#ifndef _gummocap_PSAPAF_H
#define _gummocap_PSAPAF_H
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include "annealedpfparallel.h"
#include "humantracker.h"
#include "poseevaluators.h"
#include <guscene/viewset2imageset.h>
#include <guimage/imagecomposer.h>
namespace gummocap
{
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton; 
using namespace gummocap::poseevaluators;
namespace algorithms {
/**\brief Class that implements the partitioned sample, annealed particle filter for human tracking
 */ 
class PSAPF: public AnnealedPFParallel, public HumanTracker
{
public:
    struct Params: public AnnealedPFParallel::Parameters
    {
        Params(int NParticles=300,int NLayers=5, double DesiredAlpha=0.5,double  StartAngleNoise=(2*M_PI)/180.,double StartTranslationNoise=0.01)
        {
            nParticles=NParticles;
            nLayers=NLayers;
            desiredAlpha=DesiredAlpha;
            startAngleNoise=StartAngleNoise;
            startTranslationNoise=StartTranslationNoise; 
        }
        void setNEvaluations(int nE,int NLayers=-1)
	{
	  if (NLayers==-1) nLayers=10;
	  else nLayers=NLayers;
	  nParticles=float(nE)/float(6*nLayers);
	}
	double startAngleNoise;//initial noise applied to angles
        double startTranslationNoise;//initial noise applied to angles
    };
    /**
     */
    PSAPF();
    /**
     */
    ~PSAPF();
    /**Sets the required params with the body model in the initial position
     */
    void setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet,PoseEvaluator* poseEval, Params params);
    /**Performs an step of the algorithm
        */
    void step(guscene::view::ViewSet *VSet)throw(GUException);
    /**Gets the current best solution
    */
    SkeletonPose getCurrentEstimation();
    
    //enabled debug
    void enableDebug(bool debug){_debug=debug;}
private:
    /**Calculates the probability of the particle passed
      * @param p particle
      * @param threadId id of the current thread performing the operation.
     */
    double calculateParticleProb(Particle &p,unsigned int threadId);
    /** converts a particle in a full skeleton pose
    */
    void fromParticle2PoseWoT(Particle &particle,SkeletonPoseWoT &pose);
    /** converts a particle in a full skeleton pose
    */
    void fromParticle2Pose(Particle &particle,SkeletonPose &pose);
    /** converts a particle in a full skeleton pose
    */
    void fromPoseWoT2Particle(SkeletonPoseWoT &pose,Particle &particle);
    /**
     */
    enum Masks{GLOBAL=0x0001,TORSO=0x0002,HEAD=0x0004,LSHOULDER=0x0008,LELBOW=0x0010,RSHOULDER=0x0020,RELBOW=0x0040,LHIP=0x0080,RHIP=0x0100,LKNEE=0x0200,RKNEE=0x0400};
    
    void getHierarchicalMask(int type,vector<bool> &mask);
    /**
     */
    void fromPose2Particle(SkeletonPose  &FullPose,Particle &particle);
    
    /**
    */
    void doSeach()throw(GUException);
    
    int _nStepsDone;
    bool _areParams;
    Params _params;
    PoseEvaluator* _poseEvaluator;
    vector<PoseEvaluator*> _evaluators;
    BodyModel<HumanSkeleton> _bodyModel;
    vector<BodyModel<HumanSkeleton> > _auxBodyModel;
    gumocap::skeleton::SkeletonPoseWoT _reducedPose;
    guscene::view::ViewSet *_TheViewSet;
    //debuggin
    bool _debug;
    guimage::ImageSet _TheViewSet_copy;
    guscene::utils::ViewSet2ImageSet VS2IS;
    guimage::ImageComposer ISComposer;
    void printParticles(ImageSet &IS);
};

};
};

#endif
