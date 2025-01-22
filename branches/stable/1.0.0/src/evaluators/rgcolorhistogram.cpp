#include "rgcolorhistogram.h"
#include <guimage/alloc.h>
#include <cstdlib>
#include <cstring>
using namespace std;
namespace gummocap
{

namespace poseevaluators {
/**Empty constructor
*/
RGColorHistogram::RGColorHistogram()
{
    reset();
}


/**Copy constructor
*/
RGColorHistogram::RGColorHistogram ( const RGColorHistogram &CH )
{
    (*this)=CH;
}
/**
*/
RGColorHistogram & RGColorHistogram::operator= ( const   RGColorHistogram &CH )
{
    memcpy(_data,CH._data,size()*sizeof(float));
    return *this;
}
/**
*/
RGColorHistogram::~RGColorHistogram()
{

}


/**Resets the information in the object
*/
void RGColorHistogram::reset()
{
    memset(_data,0,size()*sizeof(float));
}
/**Adds a hsv color to the histogram
 */
void RGColorHistogram::addRGB ( int R,int G,int B,float val )
{
    //find the index
    if (R+G+B!=0) {
        float nr=float(R)/255.;
        float ng=float(G)/255.;
        float nb=float(B)/255.;
        float Y=0.299*nr+0.587*ng+0.114*nb;
        float Cr=(nr-Y)+0.5;
        float Cb=(nb-Y)+0.5;
        int indexR=(  Cr*float(_gummocap_RGCH_NDIV-1)  ) +0.5;
        int indexB=(  Cb*float(_gummocap_RGCH_NDIV-1)  )+0.5;
        int indexF=indexR*_gummocap_RGCH_NDIV+indexB;
        _data[indexF]+=val;
    }
    else _data[0]+=val;

}
/**Adds a hsv color to the histogram
 */
void RGColorHistogram::addRGB ( unsigned char *rgb,int ncolors,float val)
{
    for (int i=0;i<ncolors;i++,rgb+=3) {
        int R=rgb[0];
        int G=rgb[1];
        int B=rgb[2];

        if (R+G+B!=0) {
            float nr=float(R)/255.;
            float ng=float(G)/255.;
            float nb=float(B)/255.;
            float Y=0.299*nr+0.587*ng+0.114*nb;
            float Cr=(nr-Y)+0.5;
            float Cb=(nb-Y)+0.5;
            int indexR=(  Cr*float(_gummocap_RGCH_NDIV-1)  ) +0.5;
            int indexB=(  Cb*float(_gummocap_RGCH_NDIV-1)  )+0.5;
            int indexF=indexR*_gummocap_RGCH_NDIV+indexB;
            _data[indexF]+=val;
        }
        else _data[0]+=val;
    }
}

/**
*/
double RGColorHistogram::compare (   RGColorHistogram &Model,ComparisonType type )
{
    switch ( type )
    {
    case Intersection:
    {
        double sum=0;
        for ( int i=0;i<size();i++ )
            sum+= _data[i]< Model._data[i] ? _data[i]:Model._data[i];
        return sum;
    }
    break;
    case Bhattacharyya:
    {
        Model.normalize();
        normalize();

        double sum=0;
        for (   int i=0;i<size();i++ )
            sum+= sqrt (_data[i]*Model._data[i]);
        return sqrt(1-sum);
    }break;
    default:
        cerr<<" RGColorHistogram::compare invalid comparison method"<<endl;
    };
    return -1;
}

/**
*/
void RGColorHistogram::normalize()
{
    float _votesSum=0;
    for ( int i=0;i<size();i++ )
        _votesSum+=_data[i];
    if ( _votesSum==0 )return;
    float invNVotes=1./float ( _votesSum );
    for ( int i=0;i<size();i++ )
        _data[i]*=invNVotes;
}

/**
*/
void RGColorHistogram::add(const RGColorHistogram &hist)
{
    for (int i=0;i<size();i++)
        _data[i]+=hist._data[i];
}

/**Updates this color histogram using the one passed as this=this*(1-innovation)+innovation*CH
 *@param CH the new color histogram
 *@param alpha update ratio
 */
void RGColorHistogram::update ( RGColorHistogram &CH,double innovation )
{
    assert ( size() ==CH.size() );
    double innovation_1=1-innovation;
    double sum=0;
    normalize();
    CH.normalize();

    for ( int i=0;i<size();i++ )
    {
        _data[i]= _data[i]*innovation_1+innovation*CH._data[i];
        sum+= _data [i];
    }
    //ahora normalizamos
    if (sum!=0 )
    {
        double invSum=1./sum;
        for ( int i=0;i<size();i++ )
            _data [i]*=invSum;
    }
}
};
};

