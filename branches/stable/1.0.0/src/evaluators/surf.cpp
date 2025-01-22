#include "surf.h"
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

Surf::SurfInternal Surf::_SurfInt;
/*
*
*
*/
Surf::Surf() {
    _isFrameStarted=_areParams=false;

}

/*
*
*
*/
Surf::~Surf() {


}
/*
*
*
*/
Surf::Surf(const  Surf&X) {
    _areParams=X._areParams;
    _nEvals=X._nEvals;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _isFrameStarted=X._isFrameStarted;

}
/*
*
*
*/

PoseEvaluator * Surf::makeCopy() {
    return new Surf(*this);
}
/*
*
*
*/
void Surf::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs ) throw(GUException) {
    _nEvals=0;
    _TheViewSet=vs;
    _bodyModel=body;
    _isFrameStarted=false;
    //allocation of memory
    //check all images are of same sime and allocate for only one
    IplImage *mask0=_TheViewSet->get(0)->getMaskImage();
    for (unsigned int v=1;v<_TheViewSet->size();v++) {
        IplImage *mask=_TheViewSet->get(v)->getMaskImage();
        if (mask->width!=mask0->width || mask->height!=mask0->height)
            throw GUException("Surf::setParams images of different sizes not supported yet");
    }
    _SurfInt.setParams(body,vs);
//      _SurfInt.saveVRMLModels(body,vs,Point3Df(0,0,1),"init");
    _areParams=true;

}

/*
*
*
*/
void Surf::startFrame(ViewSet *vs,int frameId) {
    //extract new surf elements
    _TheViewSet=vs;
    _SurfInt.startFrame(vs,frameId);
    //
    _isFrameStarted=true;




}

/*
*
*
*/
double Surf::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) {
    _nEvals++;
    if (!_isFrameStarted) {
        cerr<<"Surf::evaluate startFrame must be called first:"<<__FILE__<<"  in line "<<__LINE__<<endl;
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
        for (unsigned int i=0;i<_SurfInt._SurfInterFrameAssignments[v].size();i++) {
            //determine the point in the mesh
            int prev=_SurfInt._SurfInterFrameAssignments[v][i].first;
            if (_SurfInt._SurfTo3DMeshAssignments[v][prev].first!=-1) {
                Point3Df point3dMesh=body->getMesh( _SurfInt._SurfTo3DMeshAssignments[v][prev].first ).getVertices()[_SurfInt._SurfTo3DMeshAssignments[v][prev].second];
                //current proyection
                int curr=_SurfInt._SurfInterFrameAssignments[v][i].second;
                CvSURFPoint* r = (CvSURFPoint*)cvGetSeqElem( (*_SurfInt._curr_keypoints)[v], curr );
                CvPoint center;
                int radius;
                center.x = cvRound(r->pt.x);
                center.y = cvRound(r->pt.y);
                radius = cvRound(r->size*1.2/9.*2);
                //compute the distance betwen the line from the center and the point passing through
                //now, find a point 1 meter away from the camera projecting in the position indicated
                Point3Df _pointP;
                _pointP.xyzh[0]= (center.x-cp->Intrinsics._cx)/cp->Intrinsics._fx;
                _pointP.xyzh[1]= (cp->Intrinsics._cy-center.y)/cp->Intrinsics._fy;
                _pointP.xyzh[2]=1;
                //transofrm the point to be in world coordinates
                TM.transform(_pointP.xyzh[0],_pointP.xyzh[1],_pointP.xyzh[2]);
                float dist=_SurfInt.getLinePointDist (camCenter,_pointP,point3dMesh);
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
void Surf::frameFinished(BodyModel<HumanSkeleton> *body ) {

//     _SurfInt.saveVRMLModels(body,_TheViewSet,Point3Df(1,0,1), "end");
    _SurfInt.frameFinished(body,_TheViewSet);
    //lets do a final check
    cout<<"Final val="<<evaluate(body,_TheViewSet)<<endl;

//     exit(0);
}

/*
*
*
*/
Surf::SurfInternal::~SurfInternal() {
    for (unsigned int i=0;i<_auxGrey.size();i++)
        cvReleaseImage(&_auxGrey[i]);
    for (unsigned int i=0;i< _storage.size();i++) {
        if (keypoints1[i]!=NULL) cvClearSeq(keypoints1[i]);
        if (objectDescriptors1[i]!=NULL)   cvClearSeq(objectDescriptors1[i]);
        if (keypoints2[i]!=NULL) cvClearSeq(keypoints2[i]);
        if (objectDescriptors2[i]!=NULL)   cvClearSeq(objectDescriptors2[i]);
    }
    for (unsigned int i=0;i<_storage.size();i++)
        cvReleaseMemStorage(&_storage[i]);

    //
    for (unsigned int i=0;i<_ImagesCp2.size();i++)
        cvReleaseImage(&_ImagesCp2[i]);
    for (unsigned int i=0;i<_ImagesCp1.size();i++)
        cvReleaseImage(&_ImagesCp1[i]);
}
/*
*
*
*/
void Surf::SurfInternal::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs) {
//     cvSaveImage("image-init.jpg",vs->get(0)->getBGRImage());
    //done only the first time by one of the threads
    _maxDistPix=(15 *vs->get(0)->getWidth()/640); //Change this!!
    
#pragma omp critical
    if (_areParams==false) {
        //create grey images
        for (unsigned int v=0;v<vs->size();v++)
            _auxGrey.push_back( cvCreateImage(cvGetSize(vs->get(v)->getBGRImage()),8,1));

        //allocate if not yet
        _storage.resize(vs->size());
        keypoints1.resize(vs->size());
        objectDescriptors1.resize(vs->size());
        keypoints2.resize(vs->size());
        objectDescriptors2.resize(vs->size());
        for (unsigned int i=0;i< vs->size();i++) {
            keypoints1[i]=objectDescriptors1[i]=NULL;
            keypoints2[i]=objectDescriptors2[i]=NULL;
            _storage[i]=   cvCreateMemStorage(0);
        }
        _curr_keypoints=&keypoints1;
        _curr_objectDescriptors=&objectDescriptors1;

        //
        for (unsigned int i=0;i<vs->size();i++) {
            _ImagesCp2.push_back( cvCloneImage(vs->get(i)->getBGRImage()));
            _ImagesCp1.push_back( cvCloneImage(vs->get(i)->getBGRImage()));
        }
        _curImgv=&_ImagesCp1;

        extractAndAssignSURF(body,vs);


        //swap
        _prev_keypoints=_curr_keypoints;
        _prev_objectDescriptors=_curr_objectDescriptors;
        _curr_keypoints=&keypoints2;
        _curr_objectDescriptors=&objectDescriptors2;
        _prevImgv=_curImgv;
        _curImgv=&_ImagesCp2;

        _areParams=true;
        _currentFrameId=-1;
        _frameFinished=true;
    }
}

void Surf::SurfInternal::startFrame(ViewSet *vs,int frameId) {
//     cvSaveImage("image-sf.jpg",vs->get(0)->getBGRImage());
    //only the first threads goes into
#pragma omp critical
    if (frameId!=_currentFrameId) {
        if (!_frameFinished) {
            cerr<<"Surf::SurfInternal::startFrame did not call to frameFinished"<<endl;
            exit(0);
        }
        //extract new surf elements
        extractSurf(vs);
        _SurfInterFrameAssignments.resize(vs->size());


#pragma omp parallel for
        for (unsigned int v=0;v<vs->size();v++) {
            findPairs( (*_prev_keypoints)[v], (*_prev_objectDescriptors)[v],(*_curr_keypoints)[v], (*_curr_objectDescriptors)[v],_SurfInterFrameAssignments[v] );
        }
        int totalMatches=0,totalKeyPoints=0;
        for (unsigned int i=0;i<_SurfInterFrameAssignments.size();i++) {

            totalMatches+=_SurfInterFrameAssignments[i].size();
            totalKeyPoints+=(*_prev_keypoints)[i]->total;
        }

        cout<<"matches  :"<<float(totalMatches)/float(totalKeyPoints)<<" "<<totalMatches<<endl;
        //Create a visible image of the pairing
//          savePairingImages("pairing-");

        _currentFrameId=frameId;
        _frameFinished=false;

    }
}

void Surf::SurfInternal::frameFinished(BodyModel<HumanSkeleton> *body,ViewSet *vs ) {
#pragma omp critical
    if (!_frameFinished) {
        //now, assign the SURF of the current frame to the mest body model obtained
        //for the next iteration
        assignSURF_3DPoints(body,vs);
        //finally, swap buffers
        swap( _curr_keypoints,_prev_keypoints);
        swap(_curr_objectDescriptors,_prev_objectDescriptors);

        swap(_prevImgv,_curImgv);
        _frameFinished=true;

    }
}

/*
*
*
*/

void Surf::SurfInternal::extractSurf(ViewSet *_Vs) {

    //deallocate previous data
    for (unsigned int i=0;i< _Vs->size();i++) {
        if ((*_curr_keypoints)[i]!=NULL)
            cvClearSeq((*_curr_keypoints)[i]);
        if ((*_curr_objectDescriptors)[i]!=NULL)
            cvClearSeq((*_curr_objectDescriptors)[i]);
    }

    //extract surf in each camera considering only the mask points
#pragma omp parallel for
    for (unsigned int i=0;i< _Vs->size();i++) {
        cvCvtColor( _Vs->get(i)->getBGRImage(),_auxGrey[i],CV_BGR2GRAY );
        CvSURFParams params = cvSURFParams(300, 1);
        cvExtractSURF( _auxGrey[i], _Vs->get(i)->getMaskImage() , & (*_curr_keypoints)[i], &(*_curr_objectDescriptors)[i], _storage[i], params );
    }

    for (unsigned int i=0;i< _Vs->size();i++)
        cvCopyImage( _Vs->get(i)->getBGRImage(),(*_curImgv)[i]);

}

void Surf::SurfInternal::assignSURF_3DPoints(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs) {
    _SurfTo3DMeshAssignments.resize(_Vs->size());
    for (unsigned int v=0;v<_Vs->size();v++) {
        _SurfTo3DMeshAssignments[v].resize((*_curr_keypoints)[v]->total);
    }
    _Visibility.resize(_Vs->size());

    //now, find the assignment
#pragma omp parallel for
    for (unsigned int v=0;v<_Vs->size();v++) {

        Point3Df camCenter;
        CameraParameters *cp=_Vs->get(v)->getCameraParameters();
        cp->getCameraLocation(camCenter.xyzh[0],camCenter.xyzh[1],camCenter.xyzh[2]);
        _body->computeVisibility_withintersections(camCenter,&_Visibility[v] );
        //_body->computeVisibilityLoc(camCenter,_Visibility[v] );
        guimage::Matrix3D TM;
        TM.createFromExtrinsicsParamsInv(cp->Extrinsics);
        for (int kp = 0; kp < (*_curr_keypoints)[v]->total; kp++ ) {
            CvSURFPoint* r = (CvSURFPoint*)cvGetSeqElem( (*_curr_keypoints)[v], kp );
            CvPoint center;
            int radius;
            center.x = cvRound(r->pt.x);
            center.y = cvRound(r->pt.y);
            radius = cvRound(r->size*1.2/9.*2);
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
            _SurfTo3DMeshAssignments[v][kp]=pair<int,int>(-1,-1);
            for (unsigned int part=0;part<_body->size();part++) {
                svector<Point3Df> &vertices=_body->getMesh(part).getVertices();

                for (unsigned int vert=0;vert<vertices.size();vert++) {
                    float dist=getLinePointDist (camCenter,_pointP,vertices[vert]);
                    if (_Visibility[v][part][vert]) {
                        if (dist<minDist && dist<0.01) {//assign if is near enough and the nearest
                            minDist=dist;
                            _SurfTo3DMeshAssignments[v][kp].first=part;
                            _SurfTo3DMeshAssignments[v][kp].second=vert;
                        }//
                    }
                }
            }
        }
    }
    cout<<"done"<<endl;

}
/*
*
*
*/

void Surf::SurfInternal::extractAndAssignSURF(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs) {
    extractSurf(_Vs);
    assignSURF_3DPoints(_body,_Vs);
}

/*
*
*
*/
void  Surf::SurfInternal::saveVRMLModels(BodyModel<HumanSkeleton> *_body,ViewSet *_Vs,Point3Df color,string baseName) {
//      return ;
    for (unsigned int v=0;v<_Vs->size();v++) {
        char name[100];
        sprintf(name,"%smodel-%d.wrl",baseName.c_str(),v);
        ofstream fileVRML(name);
        fileVRML<< "#VRML V2.0 utf8"<<endl;
        _body->saveToVRML(fileVRML,color);
        //now, add the selected points
        for (unsigned int i=0;i<_SurfTo3DMeshAssignments[v].size();i++) {
            if (_SurfTo3DMeshAssignments[v][i].first!=-1) {
                Point3Df point= _body->getMesh(_SurfTo3DMeshAssignments[v][i].first).getVertices()[_SurfTo3DMeshAssignments[v][i].second];
                fileVRML<<" Transform {"<<endl;
                fileVRML<<" translation " <<point[0]<<" "<<point[1]<<" "<<point[2]<<endl;
                fileVRML<<"children ["<<endl;
                fileVRML<<"Shape {"<<endl;
                fileVRML<<"appearance Appearance { material Material {diffuseColor 0.0 1 0.0  } }"<<endl;
                fileVRML<<"geometry Sphere {radius  0.01}	} ]}"<<endl;

            }
        }
    }
}
/*
*
*
*/
void Surf::SurfInternal::savePairingImages(string baseName ) {

//         for (unsigned int v=0;v<_Vs->size();v++) {
//             CameraParameters *cp=_Vs->get(v)->getCameraParameters();
//             for (unsigned int i=0;i<_SurfTo3DMeshAssignments[v].size();i++) {
//                 if (_SurfTo3DMeshAssignments[v][i].first!=-1) {
//                     Point3Df point= _body->getMesh(_SurfTo3DMeshAssignments[v][i].first).getVertices()[_SurfTo3DMeshAssignments[v][i].second];
//                     int xo,yo;
//                     cp->project3DPoint(point.xyzh[0],point.xyzh[1],point.xyzh[2],xo,yo);
//                     //draw
//                     cvRectangle( _Vs->get(v)->getBGRImage(),cvPoint(xo-2,yo-2),cvPoint(xo+2,yo+2),cvScalar(255,0,0),-1);
//                 }
//             }
//         }

    for (unsigned int v=0;v<(*_curr_keypoints).size();v++) {
        for (int i = 0; i < (*_curr_keypoints)[v]->total; i++ ) {
            CvSURFPoint* r = (CvSURFPoint*)cvGetSeqElem( (*_curr_keypoints)[v], i );
            CvPoint center;
            int radius;
            center.x = cvRound(r->pt.x);
            center.y = cvRound(r->pt.y);
            radius = cvRound(r->size*1.2/9.*2);
            cvRectangle( (*_curImgv)[v],cvPoint(center.x-2,center.y-2),cvPoint(center.x+2,center.y+2),cvScalar(255,0,0),-1);
        }
        for (int i = 0; i < (*_prev_keypoints)[v]->total; i++ ) {
            CvSURFPoint* r = (CvSURFPoint*)cvGetSeqElem( (*_prev_keypoints)[v], i );
            CvPoint center;
            int radius;
            center.x = cvRound(r->pt.x);
            center.y = cvRound(r->pt.y);
            radius = cvRound(r->size*1.2/9.*2);
            cvRectangle( (*_prevImgv)[v],cvPoint(center.x-2,center.y-2),cvPoint(center.x+2,center.y+2),cvScalar(255,0,0),-1);
        }
        //now, draw the matching
        //create a composed image
        IplImage *cvComposed=cvCreateImage( cvSize((*_prevImgv)[v]->width*2,(*_prevImgv)[v]->height),8,3);
        cvSetImageROI(cvComposed,cvRect(0,0,(*_prevImgv)[v]->width,(*_prevImgv)[v]->height));
        cvCopy((*_prevImgv)[v],cvComposed);
        cvSetImageROI(cvComposed,cvRect((*_prevImgv)[v]->width,0,(*_prevImgv)[v]->width,(*_prevImgv)[v]->height));
        cvCopy((*_curImgv)[v],cvComposed);
        cvResetImageROI(cvComposed);
        //now, write lines between the matches
        int width=(*_prevImgv)[v]->width;
        for (unsigned int a=0;a<_SurfInterFrameAssignments[v].size();a++) {
            CvSURFPoint* rprev = (CvSURFPoint*)cvGetSeqElem( (*_prev_keypoints)[v], _SurfInterFrameAssignments[v][a].first );
            CvSURFPoint* rcur = (CvSURFPoint*)cvGetSeqElem( (*_curr_keypoints)[v], _SurfInterFrameAssignments[v][a].second );
            CvPoint centerPrev;
            centerPrev.x = cvRound(rprev->pt.x);
            centerPrev.y = cvRound(rprev->pt.y);
            CvPoint centerCur;
            centerCur.x = cvRound(rcur->pt.x);
            centerCur.y = cvRound(rcur->pt.y);
            cvLine( cvComposed, cvPoint(centerPrev.x,centerPrev.y),cvPoint( width+centerCur.x,centerCur.y),cvScalar(255,255,255));

        }

        char name[100];
        sprintf(name,"%s-%d.jpg",baseName.c_str(),v);
        cvSaveImage(name,cvComposed);
        cvReleaseImage(&cvComposed);
    }
}


/**
*/
float  Surf::SurfInternal::getLinePointDist ( Point3D<float> uLstart,Point3D<float> uLend,Point3D<float> upoint) {
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

/**
*/
void Surf::SurfInternal::findPairs(  CvSeq* objectKeypoints,  CvSeq* objectDescriptors,
                                     CvSeq* imageKeypoints,  CvSeq* imageDescriptors, vector< pair<int,int> >& ptpairs ) {
    int i;
    CvSeqReader reader, kreader;
    cvStartReadSeq( objectKeypoints, &kreader );
    cvStartReadSeq( objectDescriptors, &reader );
    ptpairs.clear();

    for ( i = 0; i < objectDescriptors->total; i++ ) {
        const CvSURFPoint* kp = (const CvSURFPoint*)kreader.ptr;
        const float* descriptor = (const float*)reader.ptr;
        CV_NEXT_SEQ_ELEM( kreader.seq->elem_size, kreader );
        CV_NEXT_SEQ_ELEM( reader.seq->elem_size, reader );
        int nearest_neighbor = naiveNearestNeighbor( descriptor, kp->laplacian, imageKeypoints, imageDescriptors, kp->pt );
        if ( nearest_neighbor >= 0 ) {
            ptpairs.push_back (pair<int,int>(i,nearest_neighbor));
        }
    }
}

/**
*/
int
Surf::SurfInternal::naiveNearestNeighbor( const float* vec, int laplacian,
        const CvSeq* model_keypoints,
        const CvSeq* model_descriptors,CvPoint2D32f point ) {
    int length = (int)(model_descriptors->elem_size/sizeof(float));
    int i, neighbor = -1;
    double d, dist1 = 1e6, dist2 = 1e6;
    CvSeqReader reader, kreader;
    float pixDistance;
    cvStartReadSeq( model_keypoints, &kreader, 0 );
    cvStartReadSeq( model_descriptors, &reader, 0 );

    for ( i = 0; i < model_descriptors->total; i++ ) {
        const CvSURFPoint* kp = (const CvSURFPoint*)kreader.ptr;
        const float* mvec = (const float*)reader.ptr;
        CV_NEXT_SEQ_ELEM( kreader.seq->elem_size, kreader );
        CV_NEXT_SEQ_ELEM( reader.seq->elem_size, reader );
        if ( laplacian != kp->laplacian )
            continue;
        d = compareSURFDescriptors( vec, mvec, dist2, length );
        if ( d < dist1 ) {
            dist2 = dist1;
            dist1 = d;
            neighbor = i;
	    pixDistance=sqrt( (point.x-kp->pt.x)*(point.x-kp->pt.x)  +  (point.y-kp->pt.y)*(point.y-kp->pt.y));
        } else if ( d < dist2 )
            dist2 = d;
    }
    if ( dist1 < 0.6*dist2 && pixDistance<_maxDistPix )
        return neighbor;
    return -1;
}


/**
*/
double
Surf::SurfInternal::compareSURFDescriptors( const float* d1, const float* d2, double best, int length ) {
    double total_cost = 0;
    assert( length % 4 == 0 );
    for ( int i = 0; i < length; i += 4 ) {
        double t0 = d1[i] - d2[i];
        double t1 = d1[i+1] - d2[i+1];
        double t2 = d1[i+2] - d2[i+2];
        double t3 = d1[i+3] - d2[i+3];
        total_cost += t0*t0 + t1*t1 + t2*t2 + t3*t3;
        if ( total_cost > best )
            break;
    }
    return total_cost;
}


}

}
