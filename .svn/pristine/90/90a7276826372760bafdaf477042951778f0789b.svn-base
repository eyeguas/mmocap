#include "gpu/gpu_xor.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>


// void __compute();

//#define _SUBREGION
namespace gummocap
{

 
namespace poseevaluators
{


/*
*
*
*/
Gpu_Xor::Gpu_Xor()
{
    _areParams=false;
}

/*
*
*
*/
Gpu_Xor::Gpu_Xor ( const  Gpu_Xor&X )
{
    _areParams=X._areParams;
    _nEvals=X._nEvals;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _proj = X._proj.clone();
    _meshAvgSize=X._meshAvgSize;
    _viewEval=X._viewEval;
    _validViews=X._validViews;
}
/*
*
*
*/

PoseEvaluator * Gpu_Xor::makeCopy()
{
    return new Gpu_Xor ( *this );
}
/*
*
*
*/
Gpu_Xor::~Gpu_Xor()
{

}
/*
*
*
*/
void Gpu_Xor::setParams ( BodyModel<HumanSkeleton> *body,ViewSet *vs )
{
    _TheViewSet=vs;
    _bodyModel=body;
    _nEvals=0;
    _viewEval.resize ( vs->size() );
    _validViews.resize ( vs->size() );
    for ( unsigned int i=0;i<_validViews.size();i++ ) _validViews[i]=true;
    //Determine the maximum number of points in a bodypart and reserve memory for fast proyection
    unsigned int maxP=0;
    for ( unsigned int i=0;i<body->size();i++ )
    {
        if ( maxP<body->getMesh ( i ).getVertices().size() )
            maxP=body->getMesh ( i ).getVertices().size();
    }
    ThePProj.reserve ( maxP );



    //check all images are of same sime and allocate for only one
    cv::Mat mask0=_TheViewSet->get ( 0 )->getMaskImage();
    for ( unsigned int v=1;v<_TheViewSet->size();v++ )
    {
        cv::Mat mask=_TheViewSet->get ( v )->getMaskImage();
        if ( mask.size()!=mask0.size())
        {
            cerr<<"Gpu_Xor::setParams images of different sizes not supported yet : "<<__FILE__<<":"<<__LINE__<<endl;
            exit ( 0 );
        }

        //allocate only for one image
        _proj.create( mask0.size(),mask0.type());
	_proj.setTo(cv::Scalar::all(0));
    }

    //get the average size of  polygon meshs in eachbody part
    _meshAvgSize.resize ( body->size() );
    for ( unsigned int part=0;part<body->size();part++ )
    {
        double meanSum=0;
        int NS=0;
        svector<PolygonMesh::MeshFace> &faces=body->getMesh ( part ).getFaces();
        svector<Point3Df> &vertices=body->getMesh ( part ).getVertices();
        for ( unsigned int f=0;f<faces.size();f++ )
        {
            //estimate only by the first two vertices
            Point3Df p1= vertices[faces[f].v[0]];
            Point3Df p2= vertices[faces[f].v[1]];
            meanSum+=p1.distance ( p2 );
            NS++;
        }
        if ( NS!=0 ) _meshAvgSize[part]=meanSum/double ( NS );
        else _meshAvgSize[part]=0;
//       cout<<"VS="<<_meshAvgSize[part]<<endl;
    }
    _areParams=true;

}

/*
*
*
*/
double Gpu_Xor::evaluate ( BodyModel<HumanSkeleton> *body,ViewSet *vs )
{
  _nEvals++;
//   __compute();

    //    double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
    double sumAvrg=0;
    _TheViewSet=vs;
    _bodyModel=body;
    int width=_TheViewSet->get(0)->getMaskImage()->width;
    int height=_TheViewSet->get(0)->getMaskImage()->height;
    //parallel version of the algorithm
    int totalErrors=0,totalPoints=0;
    //Here the main loop, the point analysis
    for (unsigned int v=0;v<_TheViewSet->size();v++) {
	_proj.setTo(cv::Scalar::all(0));//set syntetic image to zero
        CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters  ();
	float TMatrix[16];
	cp->getProjectionMatrix(TMatrix);
        //project the body part
        for (unsigned part=0;part<_bodyModel->size();part++) {
            gumocap::core::svector<Point3D<float> > &vertices=_bodyModel->getMesh(part).getVertices();
	    //project points
	    for(size_t  p=0;p<vertices.size();p++){
		Point3D<float> res;
		//translate the point to camera coordinates
		res.xyzh[0]=vertices[p].xyzh[0]* TMatrix[0]+ vertices[p].xyzh[1]* TMatrix[1]+ vertices[p].xyzh[2]* TMatrix[2] +TMatrix[3];
		res.xyzh[1]=vertices[p].xyzh[0]* TMatrix[4]+ vertices[p].xyzh[1]* TMatrix[5]+ vertices[p].xyzh[2]* TMatrix[6] +TMatrix[7];
		res.xyzh[2]=vertices[p].xyzh[0]* TMatrix[8]+ vertices[p].xyzh[1]* TMatrix[9]+ vertices[p].xyzh[2]* TMatrix[10] +TMatrix[11];
		//now, project
		int xo=( res.xyzh[0] / res.xyzh[2 ] ) *cp->Intrinsics._fx+cp->Intrinsics._cx ;
		int yo=( res.xyzh[1] / res.xyzh[2 ] ) *cp->Intrinsics._fy+cp->Intrinsics._cy ;
		//is into image limits , then draw
		if (xo>=0 && yo>=0 && xo<_proj.cols&& yo<_proj.rows){
		  _proj.at<uchar>(yo,xo)=255;
		}
	    }
        }


        cv::Mat mask=_TheViewSet->get(v)->getMaskImage();
 	cv::imshow("proj",_proj);
 	cv::waitKey(0);
	int count=xor_special(_proj,mask);
        totalErrors+=count;
        sumAvrg+=(float)count/(mask.cols*mask.rows);
    }
    float error=float(totalErrors)/float(totalPoints);

    return sumAvrg/float(_TheViewSet->size());

}
/**
*/
int Gpu_Xor::xor_special ( cv::Mat im1,cv::Mat im2)
{

  int count=0;
  assert(im1.type()==CV_8UC1);
  assert(im2.type()==CV_8UC1);
  for(int y=0;y<im1.rows;y++){
    uchar *im1_=im1.ptr<uchar>(y);
    uchar *im2_=im2.ptr<uchar>(y);
    for(int x=0;x<im1.cols;x++){
      if (im1_[x] xor im2_[x]) count++;
    }
  }
  
  return count;
}

}

}
