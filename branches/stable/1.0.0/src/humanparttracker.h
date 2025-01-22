#ifndef __HUMANPARTTRACKER__
#define __HUMANPARTTRACKER__

#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include <omp.h>
#include "poseevaluators.h"
#include <vector>

namespace gummocap
{
namespace algorithms
{
  using namespace std;  
  using namespace gumocap;
  using namespace gumocap::models;
  using namespace gumocap::skeleton;
  using namespace gummocap::poseevaluators;
  /**
 * \brief Base class for all human part trackers
 */
class HumanPartTracker
{
  public:
    
    /**
    */
    HumanPartTracker();
    /**
    */
    ~HumanPartTracker();
    /**Sets the required params and initializes
    */
    void setParams(unsigned int nBodyParts, unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints);
  
  protected:
    
    void setIniEndBodyParts6();
    void setIniEndBodyParts12();
    
    unsigned int _initrans;
    unsigned int _endtrans;
    unsigned int _inirot;
    unsigned int _endrot;
    unsigned int _inibody;
    unsigned int _endbody;
    unsigned int _inilarm;
    unsigned int _endlarm;
    unsigned int _inirarm;
    unsigned int _endrarm;
    unsigned int _inilleg;
    unsigned int _endlleg;
    unsigned int _inirleg;
    unsigned int _endrleg;
    vector<unsigned int> _ini_end_body_parts;
    
    
    bool _enabledParams;
    
    double _prev_best_fitness;
    
    static const double INFINITE=100000000.0;
    unsigned int _BODY_PARTS;  
    
    enum BodyPart6{BODYO=0,TORSO,LARM,RARM,LLEG,RLEG};
    enum BodyPart12{BODY_P=0,BODY_O,B_TORSO,LU_ARM,LL_ARM,RU_ARM,RL_ARM,HEAD,LU_LEG,LL_LEG,RU_LEG,RL_LEG};
    enum IniEndBodyParts{INITRANS=0, ENDROT, INIBODY, ENDBODY, INILARM, ENDLARM, INIRARM, ENDRARM, INILLEG, ENDLLEG, INIRLEG, ENDRLEG};
};

}

}

#endif
