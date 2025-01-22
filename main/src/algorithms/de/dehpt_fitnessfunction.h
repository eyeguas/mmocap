
#ifndef __DEHPT_FITNESSFUNCTION__
#define __DEHPT_FITNESSFUNCTION__

#include "deht_fitnessfunction.h"

namespace gummocap{
  
  namespace algorithms{

    class HPTFitnessFunction : public HTFitnessFunction
    /**
    * \brief This class implements the concept of a Fitness Function that obtains the fitness value 
    * \brief of a given solution for the hierarchical human tracker problem.
    */
    {

      public:
	/**
	*/
	HPTFitnessFunction(double* xbest, unsigned int N, unsigned int inibodypart, unsigned int endbodypart, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask);
	/**
	*/
	HPTFitnessFunction(double* xbest, unsigned int N, unsigned int inibodypart, unsigned int endbodypart, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask,vector<bool> *DOFXYZmask);
	/**Releases the memory
	*/
	~HPTFitnessFunction();
	/** Evaluates a solution
	*/
	long double fitness(const double *vars);
	/** Evaluates a solution
	*/
	long double fitness(const double *vars, int dimension);
	/** Returns the dimension of the solutions that is able to evaluate.
	*/
	int getDim();
	

      private:
	
	/** converts a vector in a reduced skeleton pose
    */
/*
    void fromPartVector2PoseWoT(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPoseWoT &pose){
      //assert(_areParams);
      pose=_reducedPose;
      unsigned int index=3, index_x=0;
      
      if(inibodypart==0){
	pose.Translation[0]=x[0];
	pose.Translation[1]=x[1];
	pose.Translation[2]=x[2];
	index_x=3;
      }
      else{
	pose.Translation[0]=_xbest[0];
	pose.Translation[1]=_xbest[1];
	pose.Translation[2]=_xbest[2];
      }
      
      // pose[0]=root
      for(unsigned int i=0;i<pose.size();i++){
	  if((index<inibodypart)||(index>endbodypart)){
	    pose[i].x=_xbest[index++];
	    pose[i].y=_xbest[index++];
	    pose[i].z=_xbest[index++];
	  }
	  else{
	    pose[i].x=x[index_x++];
	    pose[i].y=x[index_x++];
	    pose[i].z=x[index_x++];
	    index+=3;
	  }
      }
    }
    */

  void fromPartVector2PoseWoT(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPoseWoT &pose){
      //assert(_areParams);
      pose=_reducedPose;
      unsigned int index=3, index_x=0;
      
      if(inibodypart==0){
	pose.Translation[0]=x[0];
	pose.Translation[1]=x[1];
	pose.Translation[2]=x[2];
	index_x=3;
      }
      else{
	pose.Translation[0]=_xbest[0];
	pose.Translation[1]=_xbest[1];
	pose.Translation[2]=_xbest[2];
      }
      
      // pose[0]=root
      if(_reducedDOFXYZmask!=NULL){
	for(unsigned int i=0,j=0;i<pose.size();i++){
	  if((index<inibodypart)||(index>endbodypart)){
	     
	      pose[i].x=_xbest[index++];
	      pose[i].y=_xbest[index++];
	      pose[i].z=_xbest[index++];
	  }
	  else{
	    if((*_reducedDOFXYZmask)[j++]){
	      pose[i].x=x[index_x++];
	      index++;
	    }
	    else
	      pose[i].x=_xbest[index++];
	  
	    if((*_reducedDOFXYZmask)[j++]){
	      pose[i].y=x[index_x++];
	      index++;
	    }
	    else
	      pose[i].y=_xbest[index++];
	    if((*_reducedDOFXYZmask)[j++]){
	      pose[i].z=x[index_x++];
	      index++;
	    }
	    else
	      pose[i].z=_xbest[index++];
	    
	  }
	}
      }
      else{
      for(unsigned int i=0;i<pose.size();i++){
	  if((index<inibodypart)||(index>endbodypart)){
	    pose[i].x=_xbest[index++];
	    pose[i].y=_xbest[index++];
	    pose[i].z=_xbest[index++];
	  }
	  else{
	    pose[i].x=x[index_x++];
	    pose[i].y=x[index_x++];
	    pose[i].z=x[index_x++];
	    index+=3;
	  }
      }
      }
      
    }


    /** converts a vector in a full skeleton pose
    */
    
    void fromPartVector2Pose(double const *x, unsigned int inibodypart, unsigned int endbodypart, SkeletonPose &pose)
    {
      SkeletonPoseWoT redPose;
      fromPartVector2PoseWoT(x,inibodypart, endbodypart, redPose);
      SkeletonPoseWoT::fromWoT2Complete(redPose,pose);
      
    }

	
	unsigned int _inibodypart;
	unsigned int _endbodypart;
	double*      _xbest;
    };
    
  }

}

#endif
