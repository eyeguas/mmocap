#include "xorcolor.h"
#include <gu/gutimestatistics.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>
namespace gummocap
{

namespace poseevaluators {

/*
*
*
*/
XorColor::XorColor() {
    _proj=NULL;
    _areParams=false;
}

/*
*
*
*/
XorColor::XorColor(const  XorColor&X)
{
    _areParams=X._areParams;
    _TheViewSet=X._TheViewSet;
    _bodyModel=X._bodyModel;
    _invalid_v=X._invalid_v;
    _projections=X._projections;
    _ColorModel=X. _ColorModel;
    _auxColorModel=X._auxColorModel;
    _nEvals=X._nEvals;
    _currColorModel=NULL;
    if (X._currColorModel!=NULL) {
        if (X._currColorModel== &X. _ColorModel) _currColorModel=&_ColorModel;
        else  _currColorModel=&_auxColorModel;
    }
    if (X._proj!=NULL)
        _proj=cvCloneImage(X._proj );
    else
        _proj=NULL;

}
/*
*
*
*/

PoseEvaluator * XorColor::makeCopy()
{
    return new XorColor(*this);
}
/*
*
*
*/
XorColor::~XorColor() {
    if (_proj!=NULL)
        cvReleaseImage(&(_proj));

}
/*
*
*
*/
void XorColor::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs ) throw(gu::Exception)
{
    _TheViewSet=vs;
    _bodyModel=body;
   _nEvals=0;
    //Determine the maximum number of points in a bodypart and free previous possible memory
    unsigned int maxP=0;
    for (unsigned int i=0;i<body->size();i++) {
        if (maxP<body->getMesh(i).getVertices().size())
            maxP=body->getMesh(i).getVertices().size();
    }
    //allocation of memory
    //check all images are of same sime and allocate for only one
    IplImage *mask0=_TheViewSet->get(0)->getMaskImage();
    for (unsigned int v=1;v<_TheViewSet->size();v++) {
        IplImage *mask=_TheViewSet->get(v)->getMaskImage();
        if (mask->width!=mask0->width || mask->height!=mask0->height)
            throw gu::Exception("XorColor::setParams images of different sizes not supported yet");
    }
    //allocate only for one image
    _proj=cvCreateImage(cvSize(mask0->width,mask0->height),mask0->depth,mask0->nChannels);

    //allocate auxiliar vectors with the projection of the points
    _invalid_v.resize(maxP*2);
    _projections.resize(maxP*2);
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
double XorColor::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) {

    assert(_areParams);
    _nEvals++;
    //    double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
    double sumAvrgXor=0;
    _TheViewSet=vs;
    _bodyModel=body;
    Point3Df camLoc;
    int width=_TheViewSet->get(0)->getMaskImage()->width;
    int height=_TheViewSet->get(0)->getMaskImage()->height;
    //vector with the index of the body parts in the correct order
    vector<pair<unsigned int,float> > partOrder;
    //reset color histogram
    for (unsigned int i=0;i<_currColorModel->size();i++) (*_currColorModel)[i].reset();


    vector<int> nPointsPart(_bodyModel->size());//number of points in each part
    for(unsigned int i=0;i<nPointsPart.size();i++) nPointsPart[i]=0;
    gu::TimeMark Timer;
    //Here the main loop, the point analysis
    for (unsigned int v=0;v<_TheViewSet->size();v++) {

        CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters  ();
        IplImage *bgr=_TheViewSet->get(v)->getBGRImage();
        //part must sorted according to distance to the camera
        sortBodyPart(partOrder,cp);
        //calculate the visibile points
        cp->getCameraLocation( camLoc.xyzh[0],camLoc.xyzh[1],camLoc.xyzh[2]);

//   	  _bodyModel->computeVisibilityLoc (  camLoc );
        _bodyModel->computeVisibility_withintersections (  camLoc );
        for (unsigned ipart=0;ipart<_bodyModel->size();ipart++) {
            unsigned int part=partOrder[ipart].first; //index of the next part to be evaluated 
	  
            gumocap::core::svector<Point3D<float> > &vertices=_bodyModel->getMesh(part).getVertices();
            //project the body part
            cp->project3DPoints_4bytes_nodist((float*) vertices._data ,_projections._data,vertices._size,_invalid_v._data );

            svector<bool> &visible=_bodyModel->getMesh(part).getVisibility();
            //ok, now for each point calcuale the value
            int index=0;
            for (unsigned point=0;point<vertices.size() ;point++,index+=2) {
                if (_invalid_v[index]==0 &&  visible[point]  ) { //the point is  projected and visible
                    //note : _projections[v][index] is the x projection and _projections[v][index+1] the y projection in the image
                    int yo=_projections[index+1];
                    int xo=_projections[index];
                    //do only if not in a border
                    if (xo>=1 && xo+1<width && yo>=1 && yo+1<height) {
                        unsigned char *pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * yo+ xo-1;
                        unsigned char *brgp=(unsigned char *)(bgr->imageData+bgr->widthStep*yo+(xo-1)*3);
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[2],brgp[1],brgp[0],255- *pix);
                        *pix=255;
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[5],brgp[4],brgp[3],255-*(pix+1));
                        *(pix+1)=255;
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[8],brgp[7],brgp[6],255-*(pix+2));
                        *(pix+2)=255;
                        pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo-1)+xo;
                        brgp=(unsigned char *)(bgr->imageData+bgr->widthStep*(yo-1)+(xo*3));
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[2],brgp[1],brgp[0],255-*pix)
                        *pix=255;
                        pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo+1)+xo;
                        brgp=(unsigned char *)(bgr->imageData+bgr->widthStep*(yo+1)+(xo*3));
                        _def_RG_ADD_RGB((*_currColorModel)[part],brgp[2],brgp[1],brgp[0],255-*pix)
                        *pix=255;
			 nPointsPart[part]++;
                    }
                }
            }
        }
        /*
         cvNamedWindow("proj");
         cvShowImage("proj", _proj);
         cvWaitKey(0);*/
        int count=xor_special(_proj,_TheViewSet->get(v)->getMaskImage());
        sumAvrgXor+=float(count)/float( width* height);
    }

    int totalPoints=0;
    for(unsigned int i=0;i<nPointsPart.size();i++)
	totalPoints+=nPointsPart[i];
//     cout<<"To="<<totalPoints<<endl;
    float ColorSim=0;
    for (unsigned int i=0;i<_currColorModel->size();i++) {
        (*_currColorModel)[i].normalize();
	float fact=float(nPointsPart[i])/float(totalPoints);
        ColorSim+=fact*(*_currColorModel)[i].compare( _ColorModel[i],RGColorHistogram::Bhattacharyya);
    }
//     ColorSim/=float(_currColorModel->size());//1 total equivalence, 0 no equivalence
//   cout<<" eval="<<sumAvrgXor/float(_TheViewSet->size())<<" "<<(1- ColorSim) <<endl;
    return 0.7*sumAvrgXor/float(_TheViewSet->size()) +   0.3*ColorSim ;

}

/**
*/
void XorColor::sortBodyPart(vector<pair<unsigned int,float> > &indexDist,CameraParameters *cp)
{
//lets calcualte the distance to the mean point of each skeleton's bone
    Point3Df camPos;
    cp->getCameraLocation(camPos.xyzh[0],camPos.xyzh[1],camPos.xyzh[2]);



//for each body part, determine the distance to camera
//The root is an special case bcause it has not parent. Solution: root always last.
    indexDist.resize( _bodyModel->size());
    indexDist[0].first=0;
    indexDist[0].second=9999999;//infinite distance

    for (unsigned int i=1;i<_bodyModel->size();i++) {
        indexDist[i].first=i;
        //determine first and last points of segment
        unsigned int parent=_bodyModel->getSkptr()->getNode(i).getParent();
        assert(parent!=-1);
        Point3Df start=_bodyModel->getJoint(parent).xyz;
        Point3Df end=_bodyModel->getJoint(i).xyz;
//         Point3Df mean=(start+end)/2.;
        indexDist[i].first=i;
        indexDist[i].second= min( start.distance(camPos),end.distance(camPos));
    }
//now, sort in ascending order
    for (unsigned int i=0;i<indexDist.size()-1;i++)
        for (unsigned int j=i+1;j<indexDist.size();j++)
            if (indexDist[i].second> indexDist[j].second) swap(indexDist[i],indexDist[j]);

}

/**
*/
int XorColor::xor_special(IplImage *im1,IplImage *im2)
{
//    if( !gu_isaligned( im1) ||!gu_isaligned( im2) )
//    {
//     cerr<<" XorColor::xor_special not aligned images"<<__FILE__<<":"<<__LINE__<<endl;
//     exit(0);
//    }
//    //yet, only for some sizes
//    if( im1->width%16!=0)
//    {
//     cerr<<" XorColor::xor_special not appropriated size for images"<<__FILE__<<":"<<__LINE__<<endl;
//     exit(0);
//    }


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
}

}
