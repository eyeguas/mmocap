#include "optflowpm.h"
#include <omp.h>
namespace gummocap {
namespace poseevaluators {
#define MAXPOINTS 500
OptFlowPm::OptFlowPm()
{
  _useSubPix=false;

}
/**
 */
OptFlowPm::~OptFlowPm()
{
    for (unsigned int i=0;i<_points1.size();i++) {
        delete _points1[i];
        delete _points2[i];
        delete _status[i];
        cvReleaseImage(&_grey1[i]);
        cvReleaseImage(&_grey2[i]);
        cvReleaseImage(&_eig[i]);
        cvReleaseImage(&_tmp[i]);
        cvReleaseImage(&_pyramid1[i]);
        cvReleaseImage(&_pyramid2[i]);
    }
}
/**
 */
void OptFlowPm::setParams(ViewSet *vs,bool useSubPixelRefinement)
{
    _useSubPix=useSubPixelRefinement;
    _flags=0;
    _points1.resize(vs->size());
    _points2.resize(vs->size());
    _nPoints.resize(vs->size());
    _grey1.resize(vs->size());
    _grey2.resize(vs->size());
    _pyramid1.resize(vs->size());
    _pyramid2.resize(vs->size());
    _eig.resize(vs->size());
    _tmp.resize(vs->size());
    _status.resize(vs->size());
    for (unsigned int i=0;i<vs->size();i++) {
        _points1[i]=new CvPoint2D32f[MAXPOINTS];
        _points2[i]=new CvPoint2D32f[MAXPOINTS];
        _status[i]=new char[MAXPOINTS];
        _grey1[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),8,1);
        _grey2[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),8,1);
        _eig[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),IPL_DEPTH_32F,1);
        _tmp[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),IPL_DEPTH_32F,1);
        _pyramid1[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),8,1);
        _pyramid2[i]=cvCreateImage(cvGetSize(vs->get(i)->getBGRImage()),8,1);
    }
    _currGrey=&_grey1;
    _prevGrey=&_grey2;
    _curPyramid=&_pyramid1;
    _prevPyramid=&_pyramid2;
    _curPoints=&_points1;
    _prevPoints=&_points2;

    _findFeaturesCall=false;

}
/**
  */
void OptFlowPm::findFeatures(ViewSet *vs,vector< vector<CvPoint2D32f> > &features)throw(GUException)
{ 
    double quality = 0.01;
    double min_distance = 10;
    int win_size = 5;
    features.resize(vs->size());
 #pragma omp parallel for
    for (unsigned int v=0;v<vs->size();v++) {
        _nPoints[v]=MAXPOINTS;
        cvCvtColor(vs->get(v)->getBGRImage(),(*_currGrey)[v],CV_BGR2GRAY);
        cvGoodFeaturesToTrack( (*_currGrey)[v], _eig[v], _tmp[v], (*_curPoints)[v], &_nPoints[v],
                               quality, min_distance, vs->get(v)->getMaskImage(), 3, 0, 0.04 );
	//Do not know why gives seg fault in version 2.1!
	if(_nPoints[v]>0 &&_useSubPix){//
        cvFindCornerSubPix( (*_currGrey)[v], (*_curPoints)[v], _nPoints[v],
                            cvSize(win_size,win_size), cvSize(-1,-1),
                            cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,20,0.03));
	
	}

        //now copy to the out vector
        features[v].resize( _nPoints[v]);
        for (  int p=0;p<_nPoints[v];p++)
            features[v][p]=(*_curPoints)[v][p];
    }
    _findFeaturesCall=true; 

}
/**
 */
void OptFlowPm::findMatches(ViewSet *vs,vector< vector<CvPoint2D32f> >&newLocations )throw(GUException)
{
    if (!_findFeaturesCall) 
      throw GUException("OptFlowPm::findMatches You need to call OptFlowPm::findFeatures");
    swap(_currGrey,_prevGrey);
    swap(_curPyramid,_prevPyramid);
    swap(_prevPoints,_curPoints);
    int win_size = 5; 
    CvPoint2D32f invalid=cvPoint2D32f(-1,-1);
    newLocations.resize(vs->size());
#pragma omp parallel for
    for (unsigned int v=0;v<vs->size();v++) {
        cvCvtColor(vs->get(v)->getBGRImage(),(*_currGrey)[v],CV_BGR2GRAY);
        cvCalcOpticalFlowPyrLK( (*_prevGrey)[v], (*_currGrey)[v], (*_prevPyramid)[v],(*_curPyramid)[v],
                                (*_prevPoints)[v], (*_curPoints)[v], _nPoints[v], cvSize(win_size,win_size), 3, _status[v], 0,
                                cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,20,0.03), _flags );
        //set to -1,-1 the points not detected
	newLocations[v].resize(_nPoints[v]);
        for (int i=0;i<_nPoints[v];i++) {
            if ( _status[v][i]==0) newLocations[v][i]=invalid;
            else newLocations[v][i]=(*_curPoints)[v][i];
        }

    }

    _flags |= CV_LKFLOW_PYR_A_READY;
 
}

}
};

