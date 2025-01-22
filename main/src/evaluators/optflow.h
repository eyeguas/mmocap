#ifndef _gummocap_OptFlow_POSEEVALUATORS_
#define _gummocap_OptFlow_POSEEVALUATORS_

#include "poseevaluators.h"
#include "optflowpm.h"
#include <omp.h>
#include <vector>
#include <utility>
using namespace std;
namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class OptFlow: public PoseEvaluator
            /**\brief This class represents the OpticalFlow fitness function
            */
{

public:
    /**
     */
    OptFlow();
    /**
     */
    OptFlow(const  OptFlow&X);

    /**
     */
    ~OptFlow();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param parallel whether multiples threads must be employed
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs )throw(gu::Exception);
    /**
     */
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;
    /**
     */
    PoseEvaluator * makeCopy();
    /**
     */
    string getId() {
        return "optflow";
    }
    /**
    */
    void startFrame(ViewSet *vs,int frameId);
    /**Call when a frame has being processed passing the best estimation obtained
     */
    void frameFinished(BodyModel<HumanSkeleton> *body );
    /**
     */
    unsigned int numberOfEvaluations() {
        return _nEvals;
    }
private:
    unsigned int _nEvals;
    guscene::view::ViewSet *_TheViewSet;
    BodyModel<HumanSkeleton>* _bodyModel;
    bool _areParams;
    
    float getLinePointDist ( Point3D<float> uLstart,Point3D<float> uLend,Point3D<float> upoint) ;
    void find3dCorrespondeces(vector< vector<CvPoint2D32f> > &features,BodyModel<HumanSkeleton> *model);
    static vector< vector<CvPoint2D32f> > _currFeatures,_prevFeatures;
    static vector< vector< pair<int,int> > > feature3DCorrespondence;//indicates for each feature, to which point it belongs to
    static vector<vector<svector<bool> > > _Visibility;

    static bool _isFrameStarted,_isFrameFinished; 
    static int _currFrameId;
    static OptFlowPm _OptFlowInt;
};


}
}

#endif
