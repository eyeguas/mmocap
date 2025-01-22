#ifndef _gummocap_Surf_POSEEVALUATORS_
#define _gummocap_Surf_POSEEVALUATORS_

#include "poseevaluators.h"
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

class Surf: public PoseEvaluator
            /**\brief This class represents the Surf fitness function 
            */
{
 
public:
    /**
     */
    Surf();
    /**
     */
    Surf(const  Surf&X);

    /**
     */
    ~Surf();
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
        return "surf";
    }
    /**
    */
    void startFrame(ViewSet *vs,int frameId);
    /**Call when a frame has being processed passing the best estimation obtained
     */
    void frameFinished(BodyModel<HumanSkeleton> *body );
    /**
     */
    unsigned int numberOfEvaluations(){return _nEvals;}
private: 
    unsigned int _nEvals;
    guscene::view::ViewSet *_TheViewSet;
    BodyModel<HumanSkeleton>* _bodyModel;
    bool _areParams;
    bool _isFrameStarted;

    class SurfInternal {
    public:
        SurfInternal() {
            _areParams=false;
        }
        ~SurfInternal();
        void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs);
        void  startFrame(ViewSet *vs,int frameId);
	void frameFinished(BodyModel<HumanSkeleton> *body,ViewSet *vs );
	void saveVRMLModels(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs,Point3Df color,string baseName);
	
	
        float getLinePointDist ( Point3D<float> uLstart,Point3D<float> uLend,Point3D<float> upoint);
  
	
	
        void extractSurf(ViewSet *_Vs);
        void extractAndAssignSURF(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs);
// 	void getLinePointDist ( Point3D<float> uLstart,Point3D<float> uLend,Point3D<float> upoint,float &distCamInt,float &distInterPoint);

        void  findPairs(  CvSeq* objectKeypoints,  CvSeq* objectDescriptors,
                          CvSeq* imageKeypoints,  CvSeq* imageDescriptors,vector< pair<int,int> >& ptpairs );

        int naiveNearestNeighbor( const float* vec, int laplacian,
                                  const CvSeq* model_keypoints,
                                  const CvSeq* model_descriptors ,CvPoint2D32f point);
        double compareSURFDescriptors( const float* d1, const float* d2, double best, int length );
	void savePairingImages(string baseName);
	void assignSURF_3DPoints(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs);
        vector<CvSeq*> keypoints1,objectDescriptors1;
        vector<CvSeq*> keypoints2,objectDescriptors2;
        vector<CvSeq*>* _curr_keypoints,*_curr_objectDescriptors;
        vector<CvSeq*>* _prev_keypoints,*_prev_objectDescriptors;
        vector<vector<pair<int,int> > > _SurfTo3DMeshAssignments;
        vector<vector<pair<int,int> > > _SurfInterFrameAssignments; 
        vector<CvMemStorage*>  _storage;
        vector<IplImage*> _auxGrey;
	vector<vector<svector<bool> > > _Visibility;
        bool _areParams;
        int _currentFrameId;
	bool _frameFinished;
        gu::TimeMark Timer;
	
	//debug
	vector<IplImage *>_ImagesCp1,_ImagesCp2;
	vector<IplImage *> *_prevImgv,*_curImgv;
	float _maxDistPix;
    };

    static SurfInternal _SurfInt;
};


}
}

#endif
