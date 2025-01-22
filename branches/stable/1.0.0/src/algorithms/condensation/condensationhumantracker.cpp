#include "condensationhumantracker.h"
#include <gu/gurandom.h>
#include <gumocap/gupainter.h>
#include <gu/gutimemark.h>
#include <guimage/imageutils.h>
#include <opencv/highgui.h>
#include <guscene/viewset2imageset.h>
#include <guimage/imagecomposer.h>
namespace gummocap
{
namespace algorithms {
/**
*/
CondensationHumanTracker::CondensationHumanTracker ()
{
    _areParams=false;
    _TheViewSet=NULL;


}
/**
*/
CondensationHumanTracker::~CondensationHumanTracker ()
{
    for (unsigned int i=0;i<_evaluators.size();i++)
        delete _evaluators[i];

}
/**Sets the required params and initializes
*/
void CondensationHumanTracker::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet,PoseEvaluator* poseEval, Params  p,vector<bool> *DOFmask )
{
    _params=p;
    _bodyModel=sk;
    _TheViewSet=VSet;
    _poseEvaluator=poseEval;
    //let us create as many evaluators as threads
    assert(_evaluators.size()==0);
    int nThreads=getNThreads();
    for (unsigned int i=0;i<nThreads;i++)_evaluators.push_back(_poseEvaluator->makeCopy());
    _auxBodyModel.resize(nThreads);
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose,DOFmask);
    //initialize creating number of particles
    init(_params.nParticles);
    _areParams=true;

    _nStepsDone=0;
    //debug
}

void CondensationHumanTracker::step(guscene::view::ViewSet *VSet)throw(GUException)
{

    _TheViewSet=VSet;
    //performs an iteration
       for (unsigned int i=0;i<_evaluators.size();i++)
          _evaluators[i]->startFrame(VSet,_nStepsDone);

    //it means sampling and propagation, and then evaluation by calling the  calculateParticleProb function
    iterate  ();
    //indicate to evaluators that frame analysis finish
    SkeletonPose bestEstimation=getCurrentEstimation();
     _bodyModel.setPose(_auxBodyModel[0],bestEstimation);
    for(unsigned int i=0;i<_evaluators.size();i++)
      _evaluators[i]->frameFinished(&_auxBodyModel[0]);
    
    
    
    _nStepsDone++;

}
/**Reimplemeted from CONDENSATION. Creates a  particle for this problem
   */
Particle  *CondensationHumanTracker::createNewParticle() {
    Particle *P=new Particle();
    //deterimne the
    //_reducedPose is initialized in setParams first
    P->resize(3+_reducedPose.size()*3); //Translation (3) +  x,y,z rotation of non terminal nodes
    (*P)[0]=_reducedPose.Translation[0];
    (*P)[1]=_reducedPose.Translation[1];
    (*P)[2]=_reducedPose.Translation[2];
    int index=3;
    for (unsigned int i=0;i<_reducedPose.size();i++) {
        (*P)[index++]=_reducedPose[i].x;
        (*P)[index++]=_reducedPose[i].y;
        (*P)[index++]=_reducedPose[i].z;
    }
//     variateParticle(*P,*P);
    return P;
}
/**
*
*/
void CondensationHumanTracker::variateParticle(Particle &src,Particle &dst)
{
    dst.resize(src.size());
    //dumyy variation just for fast checking
    //translation

    dst[0]=src[0]+ GURandom::nextGaussian()*0.01;
    dst[1]=src[1]+ GURandom::nextGaussian()*0.01;
    dst[2]=src[2]+ GURandom::nextGaussian()*0.01;
    float devGrad=(1.*M_PI)/180.;
    //variate only the shoulder and elbow

    for (unsigned int i=3;i<src.size();i++)
        dst[i]=src[i]+ GURandom::nextGaussian()*devGrad; 
}
/**
*
*/
void CondensationHumanTracker::drawParticle(unsigned int i,guimage::ImageSet *ImgSet)
{
    //get the real pose from the particle
    SkeletonPose pose;
    fromParticle2Pose( (*(*currParticles)[i]),pose);
    _bodyModel.setPose(_auxBodyModel[0],pose);
    //now, draw
    gumocap::utils::GuPainter::draw( * _auxBodyModel[0].getSkptr(),ImgSet);
}

/**
* HEre is the key, the evaluation of the particle
*/
float CondensationHumanTracker::calculateParticleProb(Particle &p,unsigned int ThreadId)
{
    SkeletonPose pose;
    fromParticle2Pose(p,pose);
    _bodyModel.setPose(_auxBodyModel[ThreadId],pose);
    //now, project all the vertices of the model in all images
    float res=1-_evaluators[ThreadId]->evaluate( &_auxBodyModel[ThreadId], _TheViewSet);
    return res;
}
/**
* Returns the best estimation
*/
SkeletonPose CondensationHumanTracker::getMeanPose()
{
    Particle p=getUnimodalEstimate  () ;
    SkeletonPose pose;
    fromParticle2Pose(p,pose);
    return pose;

}
/*
SkeletonPose CondensationHumanTracker::getBestPose()
{

    int bestIndex=0;
    float maxProb=(*currParticles)[bestIndex]->prob;
    for (unsigned int i=1;i<currParticles->size();i++)
        if (maxProb<(*currParticles)[i]->prob) {
            bestIndex=i;
            maxProb=(*currParticles)[i]->prob;
        }

    SkeletonPose pose;
    fromParticle2Pose(*(*currParticles)[bestIndex],pose);
    return pose;

}*/
void CondensationHumanTracker::saveProbHist(string file)
{

    int div=100;
    vector<float> hist(div);
    for (unsigned int i=0;i<div;i++) hist[i]=0;
    float maxProb=-1;
    for (unsigned int i=0;i<currParticles->size();i++)
        if (maxProb<(*currParticles)[i]->prob)
            maxProb=(*currParticles)[i]->prob;

    cout<<"maxProb="<<maxProb<<endl;

    for (unsigned int i=0;i<currParticles->size();i++)
    {
        float prob=(*currParticles)[i]->prob/maxProb;
        //get the corresponding bin
        int bin=int( prob*float(div));
//       cout<<"prob="<<prob<<endl;
        if (bin==div) bin=div-1;
        hist[bin]++;
    }
    ofstream fileOut(file.c_str());
    for (int i=0;i<div;i++)
        fileOut<<i<<" "<<hist[i]<<endl;
}
/** converts a particle in a full skeleton pose
   */
void CondensationHumanTracker::fromParticle2PoseWoT(Particle &particle,SkeletonPoseWoT &pose) {
    assert(_areParams);
    pose=_reducedPose;
    //now, copy data
    pose.Translation[0]=particle[0];
    pose.Translation[1]=particle[1];
    pose.Translation[2]=particle[2];
    int index=3;

    //TTTTTT
    for (unsigned int i=0;i<pose.size();i++) {
        pose[i].x=particle[index++];
        pose[i].y=particle[index++];
        pose[i].z=particle[index++];
    }
}
/** converts a particle in a full skeleton pose
*/
void CondensationHumanTracker::fromParticle2Pose(Particle &particle,SkeletonPose &pose)
{
    SkeletonPoseWoT redPose;
    fromParticle2PoseWoT(particle,redPose);//create a reduced pose first from the particle
    SkeletonPoseWoT::fromWoT2Complete(redPose,pose);//converts the reduced pose to full pose
}

};
};

