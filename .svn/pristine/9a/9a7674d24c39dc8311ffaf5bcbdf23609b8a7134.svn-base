#include "stxorcolor.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>
#include <gumocap/gupainter.h>
#include <gustereo3/stereoimageutilities.h>
namespace gummocap
{

namespace poseevaluators {

/*
*
*
*/
StXorColor::StXorColor() {
    _areParams=false;
}

/*
*
*
*/
StXorColor::StXorColor(const  StXorColor&X)
{
    _areParams=X._areParams;
    _nEvals=X._nEvals;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _invalid_v=X._invalid_v;
    _projections=X._projections;
    _zvertices=X._zvertices;
    //Create new projetion images
    _proj =cvCloneImage(X._proj );
    _xorimg=cvCloneImage(X._xorimg  );
    //Color
    _ColorModel=X. _ColorModel;
    _auxColorModel=X._auxColorModel;
    _currColorModel=NULL;
    if (X._currColorModel!=NULL) {
        if (X._currColorModel== &X. _ColorModel) _currColorModel=&_ColorModel;
        else  _currColorModel=&_auxColorModel;
    }
}
/*
*
*
*/

PoseEvaluator * StXorColor::makeCopy()
{
    return new StXorColor(*this);
}
/*
*
*
*/
StXorColor::~StXorColor() {
    if (_areParams) {
        cvReleaseImage(&_proj);
        cvReleaseImage(&_xorimg);
    }
}
/*
*
*
*/
void StXorColor::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs) throw (GUException)
{
  
    //it is required to be stereo information in the images passed
    for (unsigned int i=0;i<vs->size();i++) {
        if ( vs->get(i)->getStereoImage()==NULL) throw GUException("StXorColor::setParams views must have stereo information");
        if ( vs->get(i)->getStereoImage()->hasStereo()==false) throw GUException("StXorColor::setParams stereo images must have stereo information calculated");
    }
    _TheViewSet=vs;
    _bodyModel=body;
   _nEvals=0;
    //Determine the maximum number of points in a bodypart and free previous possible memory
    _projections.clear();
    _invalid_v.clear();
    unsigned int maxP=0;
    for (unsigned int i=0;i<body->size();i++) {
        if (maxP<body->getMesh(i).getVertices().size())
            maxP=body->getMesh(i).getVertices().size();
    }
    //allocation of memory

    ///parallel version allocates for each view
    //check all images are of same sime and allocate for only one
    IplImage *mask0=_TheViewSet->get(0)->getMaskImage();
    for (unsigned int v=0;v<_TheViewSet->size();v++) {
        IplImage *mask=_TheViewSet->get(v)->getMaskImage();
        if (mask->width!=mask0->width || mask->height!=mask0->height)
            throw GUException("StXorColor::setParams images of different sizes not supported yet ");;
    }

    //allocate only for one image
    _proj=cvCreateImage(cvSize(mask0->width,mask0->height),mask0->depth,mask0->nChannels);
    _xorimg=cvCreateImage(cvSize(_proj->width,_proj->height),_proj->depth,_proj->nChannels);
    //allocate auxiliar vectors with the projection of the points
    _invalid_v.resize(maxP*2);
    _projections.resize(maxP*2);
    _zvertices.resize(maxP);
    cvSetZero(_proj);
    _areParams=true;


    //COlor
    _ColorModel.resize(body->size() );
    _auxColorModel.resize(body->size() );
    //now, evaluate writting in the color model so that it is created the first time
    _currColorModel=&_ColorModel;
    evaluate(body,vs);
    //finally, set to write in the auxModel
    _currColorModel=&_auxColorModel;
    //in the subsequent evaluations, it will be written in the aux model so as to compare

}

/*
*
*
*/


double StXorColor::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) {
    Point3Df camLoc;
    _nEvals++;
    Timer.init();

    double StereoError=0,XorError=0,ColorError=0;
    _TheViewSet=vs;
    _bodyModel=body;
 
    //parallel version of the algorithm
    //Here the main loop, the point analysis
    for (unsigned int v=0;v<_TheViewSet->size();v++) {
        CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters  ();
        float *stXYZ=_TheViewSet->get(v)->getStereoImage()->getDepthImage();
        unsigned char *stVal=_TheViewSet->get(v)->getStereoImage()->getValidityDepthMap();
        IplImage *_mask=_TheViewSet->get(v)->getMaskImage();
        cp->getCameraLocation ( camLoc.xyzh[0],camLoc.xyzh[1],camLoc.xyzh[2]);
        _bodyModel->computeVisibilityLoc (  camLoc );
        IplImage *bgr=_TheViewSet->get(v)->getBGRImage();	
        double sumStError=0;
        int nStereoPointsComputed=0;
	
	
        for (unsigned part=0;part<_bodyModel->size();part++) {
            gumocap::core::svector<Point3D<float> > &vertices=_bodyModel->getMesh(part).getVertices();
            //project the body part
            cp->project3DPoints_4bytes_nodist((float*) vertices._data ,_projections._data,vertices._size,_invalid_v._data );
            //now, determine the point location relative to the camera position. Only z component is required
            getZverticesRelativeToCamera(vertices,_zvertices,cp);
            svector<bool> &visible=_bodyModel->getMesh(part).getVisibility();
            //ok, now for each point calcuale the value
            int index=0;
            for (unsigned point=0;point<vertices.size() ;point++,index+=2) {
                if (_invalid_v[index]==0 && visible[point]  ) { //the points has projected
                    //note : _projections[v][index] is the x projection and _projections[v][index+1] the y projection in the image
                    int yo=_projections[index+1];
                    int xo=_projections[index];
                    //project points around the model
                    if (xo>=1 && xo+1<bgr->width && yo>=1 && yo+1<bgr->height) {
                        unsigned char *pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * yo;
                        pix[xo]=255;
                        pix[xo-1]=255;
                        pix[xo+1]=255;
                        pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo-1);
                        pix[xo]=255;
                        pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo+1);
                        pix[xo]=255;
                        ///calcuale  stereo error
                        unsigned char *maskpix=(unsigned char *)  _mask->imageData+ _mask->widthStep * yo;
                        int offSet= _proj->width * yo +xo;
                        if (maskpix[xo] &&  stVal[offSet]) {
                            //now, calculate the dsitance to the expected point
                            float dist=stXYZ[offSet*3+2]-_zvertices[point];
                            if (dist<0) dist=-dist;
                            sumStError+=1- exp ( -(dist*dist)/ (2*0.15*0.15) );
                            nStereoPointsComputed++;
                        }

                        //add to Color histogram
                        unsigned char *brgp=(unsigned char *)(bgr->imageData+bgr->widthStep*yo+(xo-1)*3);
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[2],brgp[1],brgp[0],1);
                    }
                }
            }
        }

        XorError+=(float)(xor_special(_proj,_mask))/float((_mask->width*_mask->height));
        StereoError+=sumStError/float(nStereoPointsComputed);
    }
    //color comparison
    float ColorSim=0;
    for (unsigned int i=0;i<_currColorModel->size();i++) {
        (*_currColorModel)[i].normalize();
        ColorSim+=(*_currColorModel)[i].compare( _ColorModel[i]/*, RGColorHistogram::Bhattacharyya */);
    }
    ColorError=1- (ColorSim/float(_currColorModel->size()));

    XorError/=float(_TheViewSet->size());
    StereoError/=float(_TheViewSet->size());
//     double joint=(1-XorError)*(1-StereoError)*(1-ColorError);
    return (ColorError+XorError+StereoError)/3;


}


/**
*/
int StXorColor::xor_special(IplImage *im1,IplImage *im2)
{


    __m128i *end=(__m128i *)(im1->imageData+im1->widthStep*im1->height);
    __m128i *ptr1=(__m128i *)im1->imageData;
    __m128i *ptr2=(__m128i *)im2->imageData;
    __m128i zeros=_mm_set_epi8( 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    __m128i ones=_mm_set_epi8( 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);

    unsigned char _aux[16] __attribute__((aligned(16)));
    short int *_siaux=(short int *)&_aux;
    __m128i *iaux=(__m128i *)&_aux;
    int sum=0;

    for ( ; ptr1<end;ptr1++,ptr2++) {
        //first load

        __m128i r1= _mm_load_si128(ptr1);
        __m128i r2= _mm_load_si128(ptr2);
        __m128i raux;
        //does xor
        raux=_mm_xor_si128(r1,r2);
        //afte the two following instructions, the result is raux[i]==0? 0:1
        raux=_mm_cmpeq_epi8(zeros,raux); // raux[i]==0? 0xfff:0x000

        raux=_mm_andnot_si128(raux,ones); //
        //use this instruction to add the value of the 8 first and the second 8 uchar and convert them into 16 bits
        raux= _mm_sad_epu8(raux,zeros);
        _mm_store_si128(iaux,raux);//store and sum
        sum=sum+_siaux[0]+_siaux[4];
        _mm_store_si128(ptr1,zeros);//Set zero in the im1
    }
    return sum;

}


/**
*/
void StXorColor::getZverticesRelativeToCamera(gumocap::core::svector<Point3D<float> > &vertices,svector<float> &vz,CameraParameters *cp)
{
//get the matrix that transform from the real world to camera ref system
    guimage::Matrix3D T;
    T. createFromExtrinsicsParams (cp->Extrinsics);
//we are only interested in the Z component
    float ML[4];
    ML[0]=T(2,0);
    ML[1]=T(2,1);
    ML[2]=T(2,2);
    ML[3]=T(2,3);
//transform the Z component to the camera ref

    for (int i=0;i<vertices.size();i++)
        vz._data[i]=  ML[0]* vertices._data[i].xyzh[0]+
                      ML[1]* vertices._data[i].xyzh[1]+
                      ML[2]* vertices._data[i].xyzh[2]+
                      ML[3];
}

}

}

