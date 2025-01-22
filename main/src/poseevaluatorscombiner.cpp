#include <poseevaluatorscombiner.h>


namespace gummocap
{

namespace poseevaluators {
/**
 *
 */
PoseEvaluatorCombiner::PoseEvaluatorCombiner()
{
    _areParams=false;
    _destroyMem=false;
}


/**
 */
PoseEvaluatorCombiner::PoseEvaluatorCombiner( const PoseEvaluatorCombiner &X)
{
    for (unsigned int i=0;i<X._eval.size();i++)
        _eval.push_back( X._eval[i]->makeCopy());
    _areParams=X._areParams;
    _destroyMem=true;
    _mode=X._mode;
    _nEvals=X._mode;

}
/**
 *
 */
PoseEvaluatorCombiner::~PoseEvaluatorCombiner()
{
    if (_destroyMem)
        for (unsigned int i=0;i<_eval.size();i++)
            delete  _eval[i];
}

/**
 *
 */
void PoseEvaluatorCombiner::setParams(CombMode mode)
{
    _mode=mode;
    _areParams=true;
}

/**
 *
 */
void PoseEvaluatorCombiner::add(PoseEvaluator *pe)
{
    _eval.push_back(pe);
}

/**
 *
 */
double PoseEvaluatorCombiner::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs)
{ 
    _nEvals++;
    double  vals[_eval.size()];
    for (unsigned int i=0;i<_eval.size();i++)
        vals[i]=_eval[i]->evaluate(body,vs);
    //now aggregate
    double res=0;
    switch (_mode) {

    case  ADD: {
        res=0;
        for (unsigned int i=0;i<_eval.size();i++)
            res+=vals[i];
    }break;
    case  AVRG: {
        res=0;
        for (unsigned int i=0;i<_eval.size();i++)
            res+=vals[i];
        res/=float(_eval.size());
    }break;
    case MULT:
        res=1;
        for (unsigned int i=0;i<_eval.size();i++)
            res=res* (1.-vals[i]);
        res=1.-res;
        break;
    }; 
    return res;
}


/**
 *
 */
PoseEvaluator * PoseEvaluatorCombiner::makeCopy()
{
    return new PoseEvaluatorCombiner(*this);
}
/**
 *
 */
string PoseEvaluatorCombiner::getId()
{
    string combined="comb://"+_eval[0]->getId();
    for (unsigned int i=1;i<_eval.size();i++)
        combined+="+"+_eval[i]->getId();
    return combined;
}

/**
 *
 */
void PoseEvaluatorCombiner::startFrame(ViewSet *vs,int frameId)
{
    for (unsigned int i=0;i<_eval.size();i++)
        _eval[i]->startFrame(vs,frameId);
}
/**
 *
 */
void PoseEvaluatorCombiner::frameFinished(BodyModel<HumanSkeleton> *body )
{
    for (unsigned int i=0;i<_eval.size();i++)
        _eval[i]->frameFinished(body);
}


}
}

