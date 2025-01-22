
#ifndef __DEHT_FITNESSFUNCTION__
#define __DEHT_FITNESSFUNCTION__

#include "de_fitnessfunction.h"
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include "poseevaluators.h"

namespace gummocap{
  
  namespace algorithms{
    using namespace gumocap;
    using namespace gumocap::models;
    using namespace gumocap::skeleton;
    using namespace gummocap::poseevaluators;
 
    
    class HTFitnessFunction : public FitnessFunction
    /**
    * \brief This class implements the concept of a Fitness Function that obtains the fitness value 
    * \brief of a given solution for the human tracker problem.
    */
    {
    public:
	/**
	*/
	HTFitnessFunction(unsigned int N, PoseEvaluator* fitness_function, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet* TheViewSet,vector<bool> *DOFmask);
	/**Releases the memory
	*/
	~HTFitnessFunction();
	/** Evaluates a solution
	*/
	long double fitness(const double *vars);
	/** Evaluates a solution
	*/
	long double fitness(const double *vars, int dimension);
	/** Returns the dimension of the solutions that is able to evaluate.
	*/
	int getDim();
	/** Sets the ViewSet.
	*/
	void setViewSet(ViewSet* VS);
	/**
	  * This function defines which is the best fitness value given two as parameters.
	  * It is useful to design maximising and minimisin problems.
	*/
	long double compare(long double f1, long double f2);
	ostringstream* getName();
	/** indicate to evaluator that frame analysis starts
	*/
	void startFrame(ViewSet* VS,unsigned int idFrame);
	/** indicate to evaluator that frame analysis ends
	*/
	void frameFinished(BodyModel<HumanSkeleton>* bodymodel);
	
	 /** converts a vector in a full skeleton pose
    */

    void fromVector2Pose(double const *x, SkeletonPose &pose)
    {
        SkeletonPoseWoT redPose;
        fromVector2PoseWoT(x,redPose);
        SkeletonPoseWoT::fromWoT2Complete(redPose,pose);

    }


	
      protected:

	/** converts a vector in a reduced skeleton pose
        */

    void fromVector2PoseWoT(double const *x, SkeletonPoseWoT &pose) {
        //assert(_areParams);
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

   	
      
	gumocap::skeleton::SkeletonPoseWoT _reducedPose;
	PoseEvaluator* _fitness_function;
	BodyModel<HumanSkeleton> _bodyModel; 
	BodyModel<HumanSkeleton> _auxBodyModel;
	guscene::view::ViewSet *_TheViewSet;
	unsigned int _N;
    };
    
  }
  
}

#endif
