#include "psapf.h"
#include <algorithm>
#include <guimage/opencv/cvwrapper.h>
#include <opencv/highgui.h>
#include <opencv/cv.h>
#include <gumocap/gupainter.h>

namespace gummocap {
namespace algorithms
{
/**
 */
PSAPF::PSAPF()
{
    _areParams=false;
    _TheViewSet=NULL;
    _debug=false;
}
/**
 */
PSAPF::~PSAPF()
{
    for (unsigned int i=0;i<_evaluators.size();i++)
        delete _evaluators[i];
}
/**Sets the required params with the body model in the initial position
 */
void PSAPF::setParams(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet,PoseEvaluator* poseEval,  Params params)
{
    _bodyModel=sk;
    _TheViewSet=VSet;
    _poseEvaluator=poseEval;
    _params=params;
    assert(_evaluators.size()==0);
    int nThreads=getNThreads();
    for (unsigned int i=0;i<nThreads;i++)
        _evaluators.push_back(poseEval->makeCopy());

    _auxBodyModel.resize(nThreads);
    SkeletonPose  fullPose=sk.getSkptr()->getRestPose();
    gumocap::skeleton::SkeletonPoseWoT::fromComplete2WoT(fullPose,*sk.getSkptr(), _reducedPose);

    //Create the initial particle
    Particle initialParticle;
    fromPoseWoT2Particle(_reducedPose,initialParticle );
    cout<<"RP="<<_reducedPose.size()<<" "<<initialParticle.size()<<endl;

    //Create the the noise vector
    vector<double> noise(initialParticle.size());
    noise[0]=noise[1]=noise[2]=_params.startTranslationNoise;
    for (unsigned int i=3;i<noise.size();i++) noise[i]=_params.startAngleNoise;
    AnnealedPFParallel::setParams(initialParticle,noise,_params );

    //do a dummy evaluation
//     float res=_evaluators[0]->evaluate( &_bodyModel , _TheViewSet);
//     cout<<"res="<<res<<endl;
//     exit(0);
    _nStepsDone=0;
    _areParams=true;

}


/**Performs an step of the algorithm
    */
void PSAPF::step(guscene::view::ViewSet *VSet)throw(gu::Exception)
{
    _TheViewSet=VSet;
    //indicate to evaluators that frame analysis starts
    for (unsigned int i=0;i<_evaluators.size();i++)
        _evaluators[i]->startFrame(VSet,_nStepsDone);
    //performs an iteration
 
    //set all to one except these elements that will not be processed
    vector<bool> mask;
    getHierarchicalMask( GLOBAL|TORSO|HEAD,mask);
    setMask(mask);
    doSeach();
    //now, 
    getHierarchicalMask( LSHOULDER|LELBOW ,mask);
    setMask(mask);
    doSeach();
    //now, 
    getHierarchicalMask(RSHOULDER|RELBOW,mask);
    setMask(mask);
    doSeach();
    
    //now, 
    getHierarchicalMask(LSHOULDER|LELBOW,mask);
    setMask(mask);
    doSeach();
    
    //now, 
    getHierarchicalMask(RKNEE|RKNEE ,mask);
    setMask(mask);
    doSeach();

    //now, 
    getHierarchicalMask(LKNEE|LKNEE ,mask);
    setMask(mask);
    doSeach();

    
    //indicate to evaluators that frame analysis finish
    SkeletonPose bestEstimation=getCurrentEstimation();
     _bodyModel.setPose(_auxBodyModel[0],bestEstimation);
    for(unsigned int i=0;i<_evaluators.size();i++)
      _evaluators[i]->frameFinished(&_auxBodyModel[0]);

    _nStepsDone++;
}
void PSAPF::doSeach()throw(gu::Exception)
{
    //it means sampling and propagation, and then evaluation by calling the  calculateParticleProb function
    for ( AnnealedPFParallel::start();!AnnealedPFParallel::end(); ) {
        AnnealedPFParallel::nextStep();
        if (_debug) {
            cout<<"Current Exp="<<_currentExp<<endl;
            VS2IS(_TheViewSet,VS2IS.COLORINPUT )->copy(&_TheViewSet_copy);
            if (!ISComposer.areParamsSet()) ISComposer.setParams(&_TheViewSet_copy);
            printParticles(_TheViewSet_copy);
            cvNamedWindow("AnnealedPF");
            cvShowImage("AnnealedPF", guimage::opencv::getTmp( ISComposer(&_TheViewSet_copy) ));
            cvWaitKey(300);
        }
    }     
}

/**Calculates the probability of the particle passed
 * @param p particle
 * @param threadId id of the current thread performing the operation.
*/
double PSAPF::calculateParticleProb(Particle &p,unsigned int threadId)
{
    SkeletonPose pose;
    fromParticle2Pose(p,pose);
    _bodyModel.setPose(_auxBodyModel[threadId],pose);
    //now, project all the vertices of the model in all images
    float res=1-_evaluators[threadId]->evaluate( &_auxBodyModel[threadId], _TheViewSet);
    return res;
}
/**
 *
 */
SkeletonPose PSAPF::getCurrentEstimation()
{
    Particle *bestParticle=(*currParticles)[0];
    //calculates the mean position ponderated by the weight
    Particle mean((*currParticles)[0]->size());
    std::fill(mean.begin(),mean.end(),0);
    //now, iterate
    double probSum=0;
    for (unsigned int i=0;i<currParticles->size();i++) {
        if (bestParticle->prob <  ((*currParticles)[i])->prob )  bestParticle=((*currParticles)[i]);
        double prob=((*currParticles)[i])->prob;
        probSum+=((*currParticles)[i])->prob;
        for (unsigned int j=0;j<(*currParticles)[0]->size();j++)
            mean[j]+= (*(*currParticles)[i])[j]*prob ;
    }

    //now divide to gen the mean
    for (unsigned int i=0;i<mean.size();i++)      mean[i]/=probSum;
    //now, obtain the corresponding pose
    SkeletonPose retPose;
    fromParticle2Pose(mean,retPose );

//     fromParticle2Pose(*bestParticle,retPose );

    return retPose;
}
/** converts a particle in a full skeleton pose
   */
void PSAPF::fromParticle2PoseWoT(Particle &particle,SkeletonPoseWoT &pose)
{
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
void PSAPF::fromParticle2Pose(Particle &particle,SkeletonPose &pose)
{
    SkeletonPoseWoT redPose;
    fromParticle2PoseWoT(particle,redPose);//create a reduced pose first from the particle
    SkeletonPoseWoT::fromWoT2Complete(redPose,pose);//converts the reduced pose to full pose
}
/** converts a particle in a full skeleton pose
*/
void PSAPF::fromPoseWoT2Particle(SkeletonPoseWoT &pose,Particle &particle)
{
    particle.resize(3+ pose.size()*3);
    particle[0]=pose.Translation.xyzh[0];
    particle[1]=pose.Translation.xyzh[1];
    particle[2]=pose.Translation.xyzh[2];

    int index=0;
    for (int i=3 ;index<pose.size();index++) {
        particle[i++]=pose[index].x;
        particle[i++]=pose[index].y;
        particle[i++]=pose[index].z;
    }

}

/**
*/
void PSAPF::printParticles(ImageSet &IS)
{
    //take all particles and sort them according to prob in ascending order
    vector<Particle*> sorted(*currParticles);
    Particle*aux;
    float maxProb=0,minProb=9999;
    for (unsigned int i=0;i<sorted.size();i++) {
        if (maxProb<sorted[i]->prob ) maxProb=sorted[i]->prob ;
        if (minProb>sorted[i]->prob ) minProb=sorted[i]->prob ;
        for (unsigned int j=i+1;j<sorted.size();j++) {
            if (sorted[i]->prob >sorted[j]->prob ) {
                aux= sorted[i];
                sorted[i]=sorted[j];
                sorted[j]=aux;
            }
        }
    }
    //now, print them according to their prob with a color
    gu::math::LinearScale ColorScale(minProb,maxProb,0,1);
    int nParticlesToVisualize=sorted.size()*0.8;
    for (unsigned int i=sorted.size()-nParticlesToVisualize;i<sorted.size();i++) {
// 	cout << "(" <<i<<","<<(sorted[i]->prob)<<")"<<flush;
        SkeletonPose pose;
        fromParticle2Pose(*sorted[i],pose);
        _bodyModel.setPose(_auxBodyModel[0],pose);
        Point3Df color(1-ColorScale(sorted[i]->prob),0,ColorScale(sorted[i]->prob));
        gumocap::utils::GuPainter::draw(_auxBodyModel[0],&IS,color );
// 	cout<<"prob="<<prob<<endl;
    }
//     cout<<endl;
}

/** converts a particle in a full skeleton pose
*/
void PSAPF::fromPose2Particle(SkeletonPose  &FullPose, Particle &particle)
{
    SkeletonPoseWoT pose;
    SkeletonPoseWoT::fromComplete2WoT( FullPose,*_bodyModel.getSkptr(), pose);

    particle.resize(3+ pose.size()*3);
    particle[0]=pose.Translation.xyzh[0];
    particle[1]=pose.Translation.xyzh[1];
    particle[2]=pose.Translation.xyzh[2];
    int index=0;
    for (int i=3 ;index<pose.size();index++) {
        particle[i++]=pose[index].x;
        particle[i++]=pose[index].y;
        particle[i++]=pose[index].z;
    }

}
/**
*/

void PSAPF::getHierarchicalMask(int type,vector<bool> &mask)
{
    //set all to one except these elements that will not be processed
    gumocap::skeleton::SkeletonPose HierarchyMaks=_bodyModel.getSkptr()->getRestPose();
 
    if ((type & GLOBAL) !=0){
      HierarchyMaks.Translation=1;
    }
    
    if ((type & TORSO) !=0){
      for(unsigned int i=0;i<HumanSkeleton::Neck;i++)
	HierarchyMaks[i]=1;
    }
    if ((type & HEAD) !=0){      
      for(unsigned int i=HumanSkeleton::Neck;i<=HumanSkeleton::Head_tip;i++)
	HierarchyMaks[i]=1;
    }
    if ((type & LSHOULDER) !=0)
	HierarchyMaks[HumanSkeleton::Lshoulder]=1;   
    if ((type & LELBOW) !=0)
	HierarchyMaks[HumanSkeleton::Lelbow]=1;
    if ((type & RSHOULDER) !=0)
	HierarchyMaks[HumanSkeleton::Rshoulder]=1;   
    if ((type & RELBOW) !=0)
	HierarchyMaks[HumanSkeleton::Relbow]=1;
    if ((type & LHIP) !=0)
	HierarchyMaks[HumanSkeleton::Lhip]=1;
    if ((type & RHIP) !=0)
	HierarchyMaks[HumanSkeleton::Rhip]=1;
    if ((type & LKNEE) !=0)
	HierarchyMaks[HumanSkeleton::Lknee]=1;
    if ((type & RKNEE) !=0)
	HierarchyMaks[HumanSkeleton::Rknee]=1;
    
/*    
    if ((type & RARM) !=0){
      cout<<"RARM"<<endl;
      for(unsigned int i=HumanSkeleton::Rshoulder;i<=HumanSkeleton::Rhand_tip;i++)
	HierarchyMaks[i]=1;
    }
    if ((type & LLEG) !=0){
      cout<<"LLEG"<<endl;
      for(unsigned int i=HumanSkeleton::Lhip;i<=HumanSkeleton::Lfoot_tip;i++)
	HierarchyMaks[i]=1;
    }
    if ((type & RLEG) !=0){
      cout<<"RLEG"<<endl;
      for(unsigned int i=HumanSkeleton::Rhip;i<=HumanSkeleton::Rfoot_tip;i++)
	HierarchyMaks[i]=1;
    }
    */
    Particle pmask;
    fromPose2Particle(HierarchyMaks,pmask);

    //now, convert into a boolean mask
    mask.resize(pmask.size());
    for (unsigned int i=0;i< mask.size();i++){
        mask[i]=( pmask[i]==1);
    }
}
}
};
