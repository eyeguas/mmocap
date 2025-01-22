#include "optflow.h"
#include <gu/gutimestatistics.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>
#include <omp.h>
namespace gummocap {
namespace poseevaluators {
OptFlowPm OptFlow::_OptFlowInt;
vector< vector< pair<int,int> > > OptFlow::feature3DCorrespondence;
vector<vector<svector<bool> > > OptFlow::_Visibility;
vector< vector<CvPoint2D32f> > OptFlow::_currFeatures,OptFlow::_prevFeatures;
int OptFlow::_currFrameId=-1;
bool  OptFlow::_isFrameStarted=false;
bool  OptFlow::_isFrameFinished=false;
/*
*
*
*/
OptFlow::OptFlow() {
    _areParams=false;
}

/*
*
*
*/
OptFlow::~OptFlow() {


}
/*
*
*
*/
OptFlow::OptFlow(const  OptFlow&X) {
    _areParams=X._areParams;
    _nEvals=X._nEvals;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _currFeatures=X._currFeatures;
    _prevFeatures=X._prevFeatures;

}
/*
*
*
*/

PoseEvaluator * OptFlow::makeCopy() {
    return new OptFlow(*this);
}
/*
*
*
*/
void OptFlow::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs ) throw(gu::Exception) {
    _nEvals=0;
    _TheViewSet=vs;
    _bodyModel=body;
    _isFrameStarted=false;
    _currFrameId=-1;
    //allocation of memory
    //check all images are of same sime and allocate for only one
    IplImage *mask0=_TheViewSet->get(0)->getMaskImage();
    for (unsigned int v=1;v<_TheViewSet->size();v++) {
        IplImage *mask=_TheViewSet->get(v)->getMaskImage();
        if (mask->width!=mask0->width || mask->height!=mask0->height)
            throw gu::Exception("OptFlow::setParams images of different sizes not supported yet");
    }
    _OptFlowInt.setParams( vs);
#pragma omp  critical
    {
        _currFrameId=-1;
        _isFrameStarted=false;
	_isFrameFinished=false;
        _OptFlowInt.findFeatures(vs,_prevFeatures);
        find3dCorrespondeces(_prevFeatures,body);
    }
    _areParams=true;
}

/*
*
*
*/
void OptFlow::startFrame(ViewSet *vs,int frameId) {
    
//     cout<<"Thread id="<<omp_get_thread_num()<<" "<<omp_get_num_threads()<<" max="<<omp_get_max_threads()<<endl;
#pragma omp  critical
    if (_currFrameId!=frameId) {
        _currFrameId=frameId;
        _isFrameStarted=true;
	_isFrameFinished=false;
//         cout<<"Thread id IN="<<omp_get_thread_num()<<endl;
        _OptFlowInt.findMatches(vs,_currFeatures);
    }
    _TheViewSet=vs;
}

/*
*
*
*/
double OptFlow::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) {
 
    _nEvals++;
/*#pragma omp critical
      cout<<"Eval id IN="<<omp_get_thread_num()<<" ne="<<_nEvals<<endl;*/
    if (!_isFrameStarted) {
        cerr<<"OptFlow::evaluate startFrame must be called first:"<<__FILE__<<"  in line "<<__LINE__<<endl;
        exit(0);
    }
// cout<<"eval"<<endl;
    float distSum=0;
    int nPoints=0;
    //evaluation consists in deterining the distance from the points seen in the previous frame to
    //the current position in this frame.WE do this for all points matches in the two frames
    for (unsigned int v=0;v< vs->size();v++) {
        Point3Df camCenter;
        CameraParameters *cp=vs->get(v)->getCameraParameters();
        cp->getCameraLocation(camCenter[0],camCenter[1],camCenter[2]);
        guimage::Matrix3D TM;
        TM.createFromExtrinsicsParamsInv(cp->Extrinsics);
        for (unsigned int i=0;i<feature3DCorrespondence[v].size();i++) {
            //determine the point in the mesh
            if ( feature3DCorrespondence[v][i].first!=-1) {
                Point3Df point3dMesh=body->getMesh( feature3DCorrespondence[v][i].first ).getVertices()[feature3DCorrespondence[v][i].second];
                //current proyection
                CvPoint2D32f center=_currFeatures[v][i];
                //compute the distance betwen the line from the center and the point passing through
                //now, find a point 1 meter away from the camera projecting in the position indicated
                Point3Df _pointP;
                _pointP.xyzh[0]= (center.x-cp->Intrinsics._cx)/cp->Intrinsics._fx;
                _pointP.xyzh[1]= (cp->Intrinsics._cy-center.y)/cp->Intrinsics._fy;
                _pointP.xyzh[2]=1;
                //transofrm the point to be in world coordinates
                TM.transform(_pointP.xyzh[0],_pointP.xyzh[1],_pointP.xyzh[2]);
                float dist= getLinePointDist (camCenter,_pointP,point3dMesh);
                dist=dist*dist;
// 	      cout<<"d="<<dist<<","<<flush;
                distSum+=dist;
                nPoints++;
            }
        }
    }
    //compute the distance
    float meanDist=distSum/float(nPoints);

//     cout<<endl<<"meanDist="<<meanDist<<endl;

    return 1- exp( -meanDist);
}

/*
*
*
*/
void OptFlow::frameFinished(BodyModel<HumanSkeleton> *body ) {
//     cout<<"frameFinished Thread id="<<omp_get_thread_num()<<endl;
#pragma omp critical
    if (!_isFrameFinished)
    {
        _isFrameFinished=true;
        _isFrameStarted=false;
        _OptFlowInt.findFeatures(_TheViewSet, _prevFeatures); 
        find3dCorrespondeces(_prevFeatures,body);
    }
}

/**
*/
void  OptFlow:: find3dCorrespondeces(vector< vector<CvPoint2D32f> > &features,BodyModel<HumanSkeleton> *_body)
{
    feature3DCorrespondence.resize(_TheViewSet->size());
    for (unsigned int v=0;v<_TheViewSet->size();v++) {
        feature3DCorrespondence[v].resize(features[v].size());
    }
    _Visibility.resize(_TheViewSet->size());

    //now, find the assignment
#pragma omp parallel for
    for (unsigned int v=0;v<_TheViewSet->size();v++) {

        Point3Df camCenter;
        CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters();
        cp->getCameraLocation(camCenter.xyzh[0],camCenter.xyzh[1],camCenter.xyzh[2]);
        _body->computeVisibility_withintersections(camCenter,&_Visibility[v] );
        //_body->computeVisibilityLoc(camCenter,_Visibility[v] );
        guimage::Matrix3D TM;
        TM.createFromExtrinsicsParamsInv(cp->Extrinsics);
        for (unsigned int kp = 0; kp < features[v].size(); kp++ ) {
            CvPoint2D32f center=features[v][kp];
            //now, find a point 1 meter away from the camera projecting in the position indicated
            Point3Df _pointP;
            _pointP.xyzh[0]= (center.x-cp->Intrinsics._cx)/cp->Intrinsics._fx;
            _pointP.xyzh[1]= (cp->Intrinsics._cy-center.y)/cp->Intrinsics._fy;
            _pointP.xyzh[2]=1;
            //transofrm the point to be in world coordinates
            TM.transform(_pointP.xyzh[0],_pointP.xyzh[1],_pointP.xyzh[2]);
            //The line is camCenter-_pointP
            //find the nearest 3d point of the mesh to that line
            float minDist=1e10;
            feature3DCorrespondence[v][kp]=pair<int,int>(-1,-1);
            for (unsigned int part=0;part<_body->size();part++) {
                svector<Point3Df> &vertices=_body->getMesh(part).getVertices();

                for (unsigned int vert=0;vert<vertices.size();vert++) {
                    float dist=getLinePointDist (camCenter,_pointP,vertices[vert]);
                    if (_Visibility[v][part][vert]) {
                        if (dist<minDist && dist<0.01) {//assign if is near enough and the nearest
                            minDist=dist;
                            feature3DCorrespondence[v][kp].first=part;
                            feature3DCorrespondence[v][kp].second=vert;
                        }//
                    }
                }
            }
        }
    }
//     cout<<"done"<<endl;
}

/**
*/
float  OptFlow:: getLinePointDist ( Point3D<float> uLstart,Point3D<float> uLend,Point3D<float> upoint) {
    Point3D<float> inter;
    float Lstart[3]={uLstart[0],uLstart[1],uLstart[2]};
    float Lend[3]={uLend[0],uLend[1],uLend[2]};
    float point[3]={upoint[0],upoint[1],upoint[2]};
    float dp1[3],dp2[3];
    float num,den;
    dp1[0]= ( point[0] )- ( Lstart[0] );
    dp1[1]= ( point[1] )- ( Lstart[1] );
    dp1[2]= ( point[2] )- ( Lstart[2] );
    dp2[0]= ( Lstart[0] )- ( Lend[0] );
    dp2[1]= ( Lstart[1] )- ( Lend[1] );
    dp2[2]= ( Lstart[2] )- ( Lend[2] );
    num=dp1[0]*dp2[0]+dp1[1]*dp2[1]+dp1[2]*dp2[2];
    den= ( dp2[0]*dp2[0] ) + ( dp2[1]*dp2[1] ) + ( dp2[2]*dp2[2] );
    float t= ( -num/ ( den ) );
    inter[0]= float ( ( Lstart[0] ) + ( Lend[0]-Lstart[0] ) *t );
    inter[1]= float ( ( Lstart[1] ) + ( Lend[1]-Lstart[1] ) *t );
    inter[2]= float ( ( Lstart[2] ) + ( Lend[2]-Lstart[2] ) *t );
    return  inter.distance(upoint);
}



}

}
