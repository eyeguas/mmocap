#include "xor.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>

//#define _SUBREGION
namespace gummocap
{


namespace poseevaluators
{


/*
*
*
*/
Xor::Xor()
{
    _areParams=false;
}

/*
*
*
*/
Xor::Xor ( const  Xor&X )
{
    _areParams=X._areParams;
    _nEvals=X._nEvals;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _proj =cvCloneImage ( X._proj );
    _meshAvgSize=X._meshAvgSize;
    _viewEval=X._viewEval;
    _validViews=X._validViews;
}
/*
*
*
*/

PoseEvaluator * Xor::makeCopy()
{
    return new Xor ( *this );
}
/*
*
*
*/
Xor::~Xor()
{
    if ( _proj!=NULL ) cvReleaseImage ( &_proj );

}
/*
*
*
*/
void Xor::setParams ( BodyModel<HumanSkeleton> *body,ViewSet *vs )
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
    IplImage *mask0=_TheViewSet->get ( 0 )->getMaskImage();
    for ( unsigned int v=1;v<_TheViewSet->size();v++ )
    {
        IplImage *mask=_TheViewSet->get ( v )->getMaskImage();
        if ( mask->width!=mask0->width || mask->height!=mask0->height )
        {
            cerr<<"Xor::setParams images of different sizes not supported yet : "<<__FILE__<<":"<<__LINE__<<endl;
            exit ( 0 );
        }

        //allocate only for one image
        _proj=cvCreateImage ( cvSize ( mask0->width,mask0->height ),mask0->depth,mask0->nChannels );
        cvSetZero ( _proj );

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
double Xor::evaluate ( BodyModel<HumanSkeleton> *body,ViewSet *vs )
{
  _nEvals++;

    //    double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
    double sumAvrg=0;
    _TheViewSet=vs;
    _bodyModel=body;
    int width=_TheViewSet->get(0)->getMaskImage()->width;
    int height=_TheViewSet->get(0)->getMaskImage()->height;
    //parallel version of the algorithm
    vector<int> patchSizes(_bodyModel->size());
    int totalErrors=0,totalPoints=0;
    //Here the main loop, the point analysis
    for (unsigned int v=0;v<_TheViewSet->size();v++) {
//          Timer.init();
	ThePProj.resetProjectedLimits();//reset for calculating the region where the model projects
        CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters  ();
        //compute the size of the projected patches  according to the distance to the camera
	for(int bp=0;bp<_meshAvgSize.size();bp++) {
	    Point3Df bodyPos= body->getJoint(bp).xyz;
	    float dist=cp->distanceToCamera (bodyPos[0],bodyPos[1],bodyPos[2]);
	    patchSizes[bp]=int(  ((_meshAvgSize[bp] * cp->Intrinsics._fx)/dist ) +0.5);
	    if (patchSizes[bp]%2==0) patchSizes[bp]++;//make it even
	    if (patchSizes[bp]>9) patchSizes[bp]=9;//set to maximum
      }
// cout<<"PS="<<patchSize<<endl;
        for (unsigned part=0;part<_bodyModel->size();part++) {
            gumocap::core::svector<Point3D<float> > &vertices=_bodyModel->getMesh(part).getVertices();
	    if (vertices.size()==0) continue; //avoid more computation in regions with no points
	    int patchSize=patchSizes[part];
            //project the body part
	    ThePProj.setROI(patchSize,patchSize,_proj->width -patchSize,_proj->height-patchSize); 	   
	    ThePProj.project(cp,(float*) vertices._data,vertices._size);
            //THIS IS THE LOOP THAT REPEATS MOST TIMES !
            //NOTE: memset is much faster that setting the points manually(specially for large regions >3)!
            //ok, now for each point
            int index=0;
            for (unsigned point=0;point<vertices.size() ;point++,index+=2) {
                if ( ThePProj.getInvalidV()[index]==0 /*_invalid_v[index]==0 */ ) { //the points has projected
                    //note : _projections[v][index] is the x projection and _projections[v][index+1] the y projection in the image
                    int yo=ThePProj.getProyections()[index+1];
                    int xo=ThePProj.getProyections()[index];
                    switch (patchSize) {
                    case 1:
                        (*(_proj->imageData+ _proj->widthStep * (yo)+xo))=255;
                        break;
                    case 3:
                            totalPoints+=9;
                            memset(_proj->imageData+ _proj->widthStep * (yo-1)+(xo-1),255,3);
                            memset(_proj->imageData+ _proj->widthStep * (yo)+(xo-1),255,3);
                            memset(_proj->imageData+ _proj->widthStep * (yo+1)+(xo-1),255,3);
                        break;
                    case 5:
                            totalPoints+=25;
                            memset(_proj->imageData+ _proj->widthStep * (yo-2)+(xo-2),255,5);
                            memset(_proj->imageData+ _proj->widthStep * (yo-1)+(xo-2),255,5);
                            memset(_proj->imageData+ _proj->widthStep * (yo)+(xo-2),255,5);
                            memset(_proj->imageData+ _proj->widthStep * (yo+1)+(xo-2),255,5);
                            memset(_proj->imageData+ _proj->widthStep * (yo+2)+(xo-2),255,5);
                        break;
                    case 7:
                            totalPoints+=49;
                            memset(_proj->imageData+ _proj->widthStep * (yo-3)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo-2)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo-1)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo+1)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo+2)+(xo-3),255,7);
                            memset(_proj->imageData+ _proj->widthStep * (yo+3)+(xo-3),255,7);
                        break;
		    case 9:
                            totalPoints+=64;
                            memset(_proj->imageData+ _proj->widthStep * (yo-4)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo-3)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo-2)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo-1)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo+1)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo+2)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo+3)+(xo-4),255,9);
                            memset(_proj->imageData+ _proj->widthStep * (yo+4)+(xo-4),255,9);
                        break;
                    }
                }
            }
        }


        IplImage *mask=_TheViewSet->get(v)->getMaskImage();
#ifndef SUBREGION_XOR
// 	cvNamedWindow("proj");
// 	cvShowImage("proj",_proj);
// 	cvWaitKey(0);
	int count=xor_special(_proj,mask);

#else

        //enlarge the region a bit (10%)
	float minX,minY,maxX,maxY;
	ThePProj.getProjectedLimits(minX,minY,maxX,maxY);
// 	cvRectangle(_proj,cvPoint(minX,minY),cvPoint(maxX,maxY),cvScalar(255,255,255),1);	
        int Ax=0.2*(maxX-minX);
        int Ay=0.1*(maxY-minY);
        minX-= Ax;minY-= Ay;
        maxX+= Ax;maxY+= Ay;
// 	cvRectangle(_proj,cvPoint(minX,minY),cvPoint(maxX,maxY),cvScalar(255,255,255),1);	
// 	cvNamedWindow("proj");
// 	cvShowImage("proj",_proj);
// 	cvWaitKey(0);
        int count=xor_special(_proj,mask ,minX,minY,maxX-1,maxY-1);
#endif
        totalErrors+=count;
        sumAvrg+=(float)count/(mask->width*mask->height);
//         Timer.end();
//         cout<<"ProjTime="<<Timer.avrg()<<endl;
    }
    float error=float(totalErrors)/float(totalPoints);

    return sumAvrg/float(_TheViewSet->size());

}
/**
*/
int Xor::xor_special ( IplImage *im1,IplImage *im2 )
{


    __m128i *end= ( __m128i * ) ( im1->imageData+im1->widthStep*im1->height );
    __m128i *ptr1= ( __m128i * ) im1->imageData;
    __m128i *ptr2= ( __m128i * ) im2->imageData;
    __m128i zeros=_mm_set_epi8 ( 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 );
    __m128i ones=_mm_set_epi8 ( 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 );

    unsigned char _aux[16] __attribute__ ( ( aligned ( 16 ) ) );
    short int *_siaux= ( short int * ) &_aux;
    __m128i *iaux= ( __m128i * ) &_aux;
    int sum=0;

    for ( ; ptr1<end;ptr1++,ptr2++ )
    {
        //first load

        __m128i r1= _mm_load_si128 ( ptr1 );
        __m128i r2= _mm_load_si128 ( ptr2 );
        __m128i raux;
        //does xor
        raux=_mm_xor_si128 ( r1,r2 );
        //afte the two following instructions, the result is raux[i]==0? 0:1
        raux=_mm_cmpeq_epi8 ( zeros,raux ); // raux[i]==0? 0xfff:0x000

        raux=_mm_andnot_si128 ( raux,ones ); //
        //use this instruction to add the value of the 8 first and the second 8 uchar and convert them into 16 bits
        raux= _mm_sad_epu8 ( raux,zeros );
        _mm_store_si128 ( iaux,raux );//store and sum
        sum=sum+_siaux[0]+_siaux[4];
        _mm_store_si128 ( ptr1,zeros );//Set zero in the im1
    }
    return sum;

}


/**
*/
int Xor::xor_special ( IplImage *im1,IplImage *im2,int minX,int minY,int maxX,int maxY )
{

//     cout<<"jjj"<<minX<<" "<<minY<<" - "<<maxX<<" "<<maxY<<endl;
    if ( minX<0 ) minX=0;
    if ( minY<0 ) minY=0;
    if ( maxX>=im1->width ) maxX=im1->width-1;
    if ( maxY>=im1->height ) maxY=im1->height-1;

    int XOffSet=minX;
    if ( XOffSet%16!=0 ) XOffSet-=XOffSet%16;
    int nXIterations=1+ ( ( maxX-minX ) /16 );
    if ( nXIterations*16+minX>=im1->width ) nXIterations--;
    __m128i *ptr1,*ptr2;
    __m128i zeros=_mm_set_epi8 ( 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 );
    __m128i ones=_mm_set_epi8 ( 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 );

    unsigned char _aux[16] __attribute__ ( ( aligned ( 16 ) ) );
    short int *_siaux= ( short int * ) &_aux;
    __m128i *iaux= ( __m128i * ) &_aux;
    int sum=0;
    for ( int y=minY ;y<=maxY;y++ )
    {
        ptr1= ( __m128i * ) ( im1->imageData+y*im1->widthStep+ XOffSet );
        ptr2= ( __m128i * ) ( im2->imageData+y*im2->widthStep+ XOffSet );
        for ( int x=0;x<nXIterations;x++,ptr1++,ptr2++ )
        {
            //first load
            __m128i r1= _mm_load_si128 ( ptr1 );
            __m128i r2= _mm_load_si128 ( ptr2 );
            __m128i raux;
            //does xor
            raux=_mm_xor_si128 ( r1,r2 );
            //afte the two following instructions, the result is raux[i]==0? 0:1
            raux=_mm_cmpeq_epi8 ( zeros,raux ); // raux[i]==0? 0xfff:0x000

            raux=_mm_andnot_si128 ( raux,ones ); //
            //use this instruction to add the value of the 8 first and the second 8 uchar and convert them into 16 bits
            raux= _mm_sad_epu8 ( raux,zeros );
            _mm_store_si128 ( iaux,raux );//store and sum
            sum=sum+_siaux[0]+_siaux[4];
            _mm_store_si128 ( ptr1,zeros );//Set zero in the im1
        }
    }
    return sum;

}

}

}
