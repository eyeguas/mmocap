#ifndef _gummocap_MetaEvaluator_H
#define _gummocap_MetaEvaluator_H
#include "poseevaluators.h"
namespace gummocap{
  
  namespace poseevaluators{
  
      /**Evaluates the performance of an evaluator
       */
    class MetaEvaluator{
      
    public:
      
      struct CameraConfidence{
	unsigned int camIndex;//index of the camera in the set
	float conf;

	
	CameraConfidence(){}
	CameraConfidence(unsigned int index,float confidence){camIndex=index;conf=confidence;}
	CameraConfidence(const CameraConfidence &C){ camIndex=C.camIndex; conf=C.conf;}
	
	
	struct Greater{
	  Greater(){}
	 bool operator()(const CameraConfidence &C1,const CameraConfidence &C2){return C1.conf>C2.conf;}  
	};
	
      };
      /**
       */
      MetaEvaluator();

      /**Evaluates the performace of the evaluator and returns the  list of cameras with their confidence
       * in descending order
       */
      vector<CameraConfidence> evaluate (PoseEvaluator *PE,BodyModel<HumanSkeleton> *body,ViewSet *vs);
      
      /**Returns a boolean vector with these views that have at least the confidence level indicated
       */
      vector<bool> getViewMaskConf(vector<CameraConfidence> &conf,float confLevel);
      /**Returns a boolean vector mask indicating the n best cameras if n>0. 
       * If n<0, the vector mask indicate the n worst cameras
       */
      vector<bool> getViewMasknCams(vector<CameraConfidence> &conf,int n);
    private:
      //confidence of cameras based on their position
      vector<float> positionConfidence;
      
      float getSin(Point3Df cam1, Point3Df cam2, Point3Df point);
    };
    
  }
};

#endif

