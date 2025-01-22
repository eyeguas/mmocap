#include "algorithms/algorithms.h"
#include <iostream>
#include <cmath>
#include <gu/guexception.h>
#include <guimage/imagesetvideoreader.h>
#include <guscene/bckviewsetproducer.h>
#include <guscene/stbckviewsetproducer.h>
#include <guscene/scenedisplay.h>
#include <guscene/sceneparams.h>
#include <guscene/viewset2imageset.h>
#include <guscene/bckviewsetproducervis.h>
#include <gumocap/gupainter.h>
#include <gumocap/markererrorevaluator.h>
#include <gumocap/humanskeleton.h>
#include <gustereo3/stereosetvideoreader.h>
#include <guavi/aviwriter.h>
#include <getopt.h>
#include <evaluators/xor.h>
#include <evaluators/stev.h>
#include <evaluators/countforeg.h>
#include <evaluators/fogvh.h>
#include <evaluators/surf.h>
#include <evaluators/optflow.h>
#include <evaluators/metaevaluator.h>
#include <gustereo3/stereoimageutilities.h>
#include <poseevaluatorscombiner.h>
using namespace gummocap;
using namespace gummocap::poseevaluators;
using namespace gummocap::algorithms;
using namespace gumocap;
using namespace gumocap::utils;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace guimage;
using namespace guscene::view::monocular;
using namespace guscene::view;
using namespace guscene;
using namespace guavi;
using namespace std;


BodyModel<HumanSkeleton> ThePersonModel;
SkeletonPose TheCurrenBestPose;
GUTimeMark Timer;
///Video I/O
string TheInputVideoFile,TheInputModelFile,TheBckModelFile,TheOutAviFile;
int inputType=0;//0 vis, 1 svs
int TheVideo_fps=30;
ViewSetProducer *TheViewSetReader;
guscene::ViewSet *TheViewSet;
guimage::ImageSet TheISCopy;
guimage::ImageComposer ImgComposer;
guscene::utils::ViewSet2ImageSet VS2IS; //converts form viewset to imageset so that imageset routines can be employed
guscene::utils::ViewSet2ImageSet::ImageType VS2IS_TypeOfView=guscene::utils::ViewSet2ImageSet::COLORINPUT;
guavi::AviWriter NewVideo;
bool videoEnabled=false;
string TheBckVisFile;
vector<pair<int,float> > ThresHolds;

///Algoritm's parameters and evaluators
HumanTracker *TheHumanTracker=NULL;
string TheFitnessFunction;
float TheBckGlobalThresHold=3;
int TheNBckFrames=40;
gummocap::poseevaluators::PoseEvaluator* FF;
string strTheHumanTrackerAlgorithm="con";//default algorithm
string hierarchicalStrategy="6-step";
bool HS=false; // default no hierarchical strategy
unsigned int numHsteps=6;
gummocap::algorithms::CondensationHumanTracker::Params con_params;
gummocap::algorithms::AnnealedHumanTracker::Params ann_params;
gummocap::algorithms::PSAPF::Params psapf_params;
////////// CMAES ////////////////
double CMAESStdDev=0.01;
double CMAESHSStdDev=0.01;
unsigned int CMAESLambda=0;
////////// MA-CMA-Chains ///////
double MACMASearchRadius=0.05;
double MACMALocalIntensity=100;
double MACMALGRatio=0.9;
////////// PSO ////////////////
unsigned int PSONoOfAgents=10;
double PSOMaximumVelocity=0.025;
double PSOSearchRadius=0.05;
double PSOInitialWeight=2.0;
////////// DE ////////////////
double DESearchRadius=0.05;
unsigned int DEPopulationSize=10;
unsigned int DEnLPC=5;
double DE_CR=0.9;
double DE_F=0.5;
double NumFunctionFitnessEvaluations=0.0;
double NumIterations=0.0;
vector<bool> *_DOFMask=NULL; //mask for limiting the number of joints
unsigned int nJointsRightArm=HumanSkeleton::NJOINTS_RIGHT_ARM;
unsigned int nJointsLeftArm=HumanSkeleton::NJOINTS_LEFT_ARM;
unsigned int nJointsRightLeg=HumanSkeleton::NJOINTS_RIGHT_LEG;
unsigned int nJointsLeftLeg=HumanSkeleton::NJOINTS_LEFT_LEG;
unsigned int nJointsBodyandHead=HumanSkeleton::NJOINTS_BODY_AND_HEAD;
vector<bool> * createReducedModelMask();
///Error Analisys
string TheErrorFile,TheErrorFileOut;
bool TheDoAnalyzeError=false;
gumocap::markers::MarkerErrorEvaluator<HumanSkeleton> TheErrorAnalyzer;
///Misc
bool TheDebugFlag=false;
bool TheFinishFlag=false;
int timeWaitKey=0;
bool TheNoXflag=false;

///Manual mode
bool TheManualModeFlag=false;
void printControls();
void processKey ( char key );
void manualModelManipulation ( SkeletonPose &pose );

///STEREO debug
int  writeVRML ( string path,View *TV,BodyModel<HumanSkeleton> *bodyModel );
///POSE ESTIMATION
int PoseEstimation=-1;

///Meta Evaluation
MetaEvaluator MEvaluator;
float maxConf=-1,nCamsConf=0;
/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
void usage()
{
    cout<<"This program tests the human tracking approach\n\n";
    cout<<" --in-vis in.vis | --in-svs in.svs : specifies the input"<<endl;
    cout<<" -m model.hbm"<<endl;
    cout<<" -b backgrond.bck "<<endl;
    cout<<" -v video.avi"<<endl;
    cout<<" --video_fps <value>: number of frames per second employed for the video"<<endl;

    cout<<" -a algorithm (con: condensation,ann: annealed pf,psapf,cmaes,ma-cma,pso,de)"<<endl;
    cout<<"\t Parameters Condensation:"<<endl;
    cout<<"\t\t --con_np <val> : number of particles"<<endl;
    cout<<"\t Parameters Annealed ParticleFilter and PSAPF:"<<endl;
    cout<<"\t\t --ann_np <val> : number of particles"<<endl;
    cout<<"\t\t --ann_nl <val> : number of layers"<<endl;
    cout<<"\t Parameters CMAES:"<<endl;
    cout<<"\t\t --cms_dev <val> : standard deviation"<<endl;
    cout<<"\t\t --cms_n <val> : population size"<<endl;
    cout<<"\t Parameters MA-CMA:"<<endl;
    cout<<"\t\t --mcm_sr <val> : search radius"<<endl;
    cout<<"\t\t --mcm_lsi <val> : local search intensity"<<endl;
    cout<<"\t\t --mcm_lgr <val> : local / global search ratio"<<endl;
    cout<<"\t Parameters PSO:"<<endl;
    cout<<"\t\t --pso_na <val> : number of agents"<<endl;
    cout<<"\t\t --pso_mv <val> : maximum velocity"<<endl;
    cout<<"\t\t --pso_sr <val> : search radius"<<endl;
    cout<<"\t Parameters DE:"<<endl;
    cout<<"\t\t --de_ps <val> : population size"<<endl;
    cout<<"\t\t --de_lpc <val> : number of leading, placing and correcting solutions"<<endl;
    cout<<"\t\t --de_sr <val> : search radius"<<endl;
    cout<<" -f fitness_function (xor,fog,fogvh,st,stxor,stxorcolor,xorcolor,xorgrad,glxor,surf,optflow)"<<endl;
    cout<<"\t\t  Combinations can be used example:  xor+surf combines both"<<endl;
    cout<<" -e number of function fitness evaluations (stop criterion)"<<endl;
    cout<<" -t background thresHold"<<endl;
    cout<<" -s Stream:Value: specifies a paticular theshold for a video stream."<<endl;
    cout<<"\t Example: -i 0:1.3 -i 1:3  . Sets the theshold 1.3 for stream 0- and 3 for stream 1"<<endl;
    cout<<" --nbck val: indicates the number of frames employed for background generation is a bckground model file is not provided"<<endl;
    cout<<" -d enabled debug information"<<endl;
    cout<<" --showfrg: show foreground instead of color information"<<endl;
    cout<<" --nowait: do not stop between frames"<<endl;
    cout<<" --sequential: disables parallel processing"<<endl;
    cout<<" --noX: does not show results in windows"<<endl;
    cout<<" --bckVis <file.vis>: uses a precalcualted background video"<<endl;
    cout<<" --errorFile <file.csv>: do error calculation"<<endl;
    cout<<" --errorFileOut <file.csv>: saves the error in a separte file"<<endl;
    cout<<" --manualMode : swtiches to manual mode"<<endl;
    cout<<" --reduceModel: uses a reduced human model with less DOF "<<endl;
    cout<<" --hs <strategy> (6-step, 12-step) "<<endl;
    cout<<"\t\t Hierarchical strategy:"<<endl;
    cout<<"\t\t 6-step: body, torso and neck, left arm, right arm, left leg, right leg."<<endl;
    cout<<"\t\t 12-step: body, torso, left upper arm, left lower arm, right upper arm, "<<endl;
    cout<<"\t\t\t\t right lower arm, neck, left upper leg, left lower leg, right upper leg, right lower leg."<<endl;
    cout<<endl<<endl;
    cout<<" --pose_estimation <frame>: use for testing only the pose estimation in a particular frame"<<endl;
    cout<<" --maxConf <value> : estimates the confidence  of each view and selects these proviindg at least the confidence indicated."<<endl;
    cout<<"\t value is in the range [0,1]. 0 low confidence and 1 total confidence (all cameras)"<<endl;
    cout<<" --nCamsConf <value> :Uses the number of cameras specified."<<endl;
}

/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
static const char short_options [] = "h:m:b:v:e:f:t:s:a:d";

static const struct option
            long_options [] =
{
    { "help",           no_argument,   NULL,                 'h' },
    { "model",     required_argument,   NULL,           'm' },
    { "backgrond",     required_argument,   NULL,           'b' },
    { "video",     required_argument,   NULL,           'v' },
    { "fitnessFunction",     required_argument,   NULL,           'f' },
    { "evals",     required_argument,   NULL,           'e' },
    { "thresHold",     required_argument,   NULL,           't' },
    { "streamThreshold",   required_argument,   NULL,           's' },
    { "algorithm",   required_argument,   NULL,           'a' },
    { "debug",   required_argument,   NULL,           'd' },
    { "nbck",   required_argument,   NULL,           300 },
    { "con_np",   required_argument,   NULL,           301 },
    { "ann_np",   required_argument,   NULL,           302 },
    { "ann_nl",   required_argument,   NULL,           303 },
    { "cms_dev",   required_argument,   NULL,           315 },
    { "cms_n",   required_argument,   NULL,           316 },
    { "mcm_sr",   required_argument,   NULL,           317 },
    { "mcm_lsi",   required_argument,   NULL,           318 },
    { "mcm_lgr",   required_argument,   NULL,           319 },
    { "pso_na",   required_argument,   NULL,           320 },
    { "pso_mv",   required_argument,   NULL,           321 },
    { "pso_sr",   required_argument,   NULL,           322 },
    { "de_ps",   required_argument,   NULL,           323 },
    { "de_lpc",   required_argument,   NULL,           324 },
    { "de_sr",   required_argument,   NULL,           325 },
    { "in-vis",     required_argument,   NULL,           304 },
    { "in-svs",     required_argument,   NULL,           305 },
    { "showfrg",     no_argument,   NULL,           306 },
    { "nowait",     no_argument,   NULL,           307 },
    { "sequential",     no_argument,   NULL,           308 },
    { "noX",     no_argument,   NULL,           309 },
    { "bckVis",     required_argument,   NULL,           310 },
    { "errorFile",     required_argument,   NULL,           311 },
    { "errorFileOut",     required_argument,   NULL,           500 },

    { "manualMode",     no_argument,   NULL,           312 },
    { "reduceModel",     no_argument,   NULL,           313 },
    { "hs",     required_argument,   NULL,           314 },
    { "pose_estimation",     required_argument,   NULL,           600 },
    { "maxConf",     required_argument,   NULL,           601 },
    { "nCamsConf",     required_argument,   NULL,           602 },
    { "video_fps",     required_argument,   NULL,           603 },

    { 0, 0, 0, 0 }
};


/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
void readArguments ( int argc,char **argv )
{
    for ( ;; )
    {
        int index;
        int c;
        c = getopt_long ( argc, argv,
                          short_options, long_options,
                          &index );

        if ( -1 == c )
            break;
        switch ( c )
        {
        case 0:
            break;
        case 'h':
            usage ();
            exit ( EXIT_SUCCESS );
            break;
        case 'm':
            TheInputModelFile=optarg;
            break;
        case 'b':
            TheBckModelFile=optarg;
            break;
        case 'e':
            NumFunctionFitnessEvaluations= ( double ) atoi ( optarg );
            break;
        case 'f':
            TheFitnessFunction=optarg;
            break;
        case 'v':
            TheOutAviFile=optarg;
            break;
        case 't':
            TheBckGlobalThresHold=atof ( optarg );
            break;
        case 's':
        {
            pair<int,float> val;
            int n=sscanf ( optarg,"%d:%f",&val.first,&val.second );
            if ( n!=2 )
            {
                cerr<<"-s argument must be in the way stream:value"<<endl;
                usage ();
                exit ( EXIT_FAILURE );
            }
            ThresHolds.push_back ( val );
        }
        break;
        case 'a':
            strTheHumanTrackerAlgorithm=optarg;
            break;
        case 'd':
            TheDebugFlag=true;
            break;
        case 300:
            TheNBckFrames=atoi ( optarg );
            break;
        case 301:
            con_params.nParticles=atoi ( optarg );
            break;
        case 302:
            ann_params.nParticles=atoi ( optarg );
            psapf_params.nParticles=atoi ( optarg );
            break;
        case 303:
            ann_params.nLayers=atoi ( optarg );
            psapf_params.nLayers=atoi ( optarg );
            break;
        case 304:
            TheInputVideoFile=optarg;
            inputType=0;
            break;
        case 305:
            TheInputVideoFile=optarg;
            inputType=1;
            break;
        case 306:
            VS2IS_TypeOfView=guscene::utils::ViewSet2ImageSet::FOREGROUND;
            break;
        case 307:
            timeWaitKey=100;
            break;
        case 308:
            omp_set_num_threads ( 1 );
            break;
        case 309:
            TheNoXflag=true;
            break;
        case 310:
            TheBckVisFile=optarg;
            break;
        case 311:
            TheErrorFile=optarg;
            TheDoAnalyzeError=true;
            break;
        case 312:
            TheManualModeFlag=true;
            break;
        case 313:
            _DOFMask=createReducedModelMask();
            break;
        case 314:
            HS=true;
            if ( !strcmp ( optarg,"6-step" ) )
                numHsteps=6;
            else if ( !strcmp ( optarg,"12-step" ) )
                numHsteps=12;
            else
                cout << "Assuming a 6-step hierarchical strategy..." << endl;
            break;
        case 315:
            CMAESStdDev=atof ( optarg );
            CMAESHSStdDev=CMAESStdDev;
            break;
        case 316:
            CMAESLambda=atoi ( optarg );
            break;
        case 317:
            MACMASearchRadius=atof ( optarg );
            break;
        case 318:
            MACMALocalIntensity=atof ( optarg );
            break;
        case 319:
            MACMALGRatio=atof ( optarg );
            break;
        case 320:
            PSONoOfAgents=atoi ( optarg );
            break;
        case 321:
            PSOMaximumVelocity=atof ( optarg );
            break;
        case 322:
            PSOSearchRadius=atof ( optarg );
            break;
        case 323:
            DEPopulationSize=atoi ( optarg );
            break;
        case 324:
            DEnLPC=atoi ( optarg );
            break;
        case 325:
            DESearchRadius=atof ( optarg );
            break;
        case 500:
            TheErrorFileOut=optarg;
            break;
        case 600:
            PoseEstimation=atoi ( optarg );
            break;
        case 601:
            maxConf=atof ( optarg );
            break;
	case 602:
	    nCamsConf=atoi ( optarg );
	    break;
	case 603:
	      TheVideo_fps=atoi ( optarg );
	    break;
        default:
            usage ();
            exit ( EXIT_FAILURE );
        };
    }
}


/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
gummocap::poseevaluators::PoseEvaluator* getIndividualFitnessEvaluator ( string type ) throw ( GUException )
{

    if ( type=="xor" )
    {
        gummocap::poseevaluators::Xor *FF=new Xor();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }

    else if ( type=="fog" )
    {
        gummocap::poseevaluators::CountForeground *FF=new CountForeground();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }
    else if ( type=="fogvh" )
    {
        gummocap::poseevaluators::FogVH *FF=new FogVH();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }
    else if ( type=="st" )
    {
        gummocap::poseevaluators::StEv *FF=new gummocap::poseevaluators::StEv();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }
    else if ( type=="surf" )
    {
        gummocap::poseevaluators::Surf *FF=new gummocap::poseevaluators::Surf();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }
    else if ( type=="optflow" )
    {
        gummocap::poseevaluators::OptFlow *FF=new gummocap::poseevaluators::OptFlow();
        FF->setParams ( &ThePersonModel,TheViewSet );
        return FF;
    }
    else
    {
        throw GUException ( "getFitnessEvaluator invalid evaluator descriptor :"+type );
    }
}



/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
gummocap::poseevaluators::PoseEvaluator* getFitnessEvaluator ( string type ) throw ( GUException )
{
//detect if several evaluator are indicated
    bool foundPLus=false;
    bool foundMult=false;

    if ( type.find ( "+" ) !=string::npos ) foundPLus=true;
    if ( type.find ( "*" ) !=string::npos ) foundMult=true;
    if ( ( foundMult|foundPLus ) ==false )  //none found
    {
        return getIndividualFitnessEvaluator ( type );
    }
    else
    {
        if ( foundPLus&foundMult ) throw GUException ( "Error combining evaluators use only + or *" );
        //several must be employed
        PoseEvaluatorCombiner *PEC =new PoseEvaluatorCombiner();

        if ( foundPLus )
            PEC->setParams ( PoseEvaluatorCombiner::AVRG );
        else if ( foundMult )
            PEC->setParams ( PoseEvaluatorCombiner::MULT );

        int prevPos=0;
        unsigned int pos=min ( type.find_first_of ( "+",0 ),type.find_first_of ( "*",0 ) );
        cout<<"POS="<<pos<<endl;
        while ( pos!=string::npos )
        {
            string stype=type.substr ( prevPos,pos-prevPos );
            cout<<"ADDING ="<<stype<<endl;
            PEC->add ( getIndividualFitnessEvaluator ( stype ) );
            prevPos=pos+1;
            pos=min ( type.find_first_of ( "+",prevPos ),type.find_first_of ( "*",prevPos ) );
        }
        //if not, take the last element
        string stype=type.substr ( prevPos,type.size()-1 );
        cout<<"ADDING  ="<<stype<<endl;
        PEC->add ( getIndividualFitnessEvaluator ( stype ) );
        return PEC;
    }
}
/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
gummocap::algorithms::HumanTracker* getHumanTracker ( string type,BodyModel<HumanSkeleton> &bodyModel,ViewSet *VS,gummocap::poseevaluators::PoseEvaluator* posEv ) throw ( GUException )
{


    if ( type=="con" )
    {
        gummocap::algorithms::CondensationHumanTracker *Alg=new gummocap::algorithms::CondensationHumanTracker();
        if ( NumFunctionFitnessEvaluations>0 )
            con_params.nParticles=NumFunctionFitnessEvaluations;
        Alg->setParams ( bodyModel, VS,posEv, con_params,_DOFMask );
        return Alg;
    }
    else if ( type=="ann" )
    {
        gummocap::algorithms::AnnealedHumanTracker *Alg=new gummocap::algorithms::AnnealedHumanTracker();
        if ( NumFunctionFitnessEvaluations>0 ) ann_params.setNEvaluations ( NumFunctionFitnessEvaluations );
        Alg->setParams ( bodyModel, VS,posEv, ann_params,_DOFMask );
        return Alg;
    }
    else if ( type=="cmaes" )
    {
        if ( NumFunctionFitnessEvaluations>0.0 )
            NumIterations=1e299;
        else
        {
            NumIterations=100;
            NumFunctionFitnessEvaluations=1e299;
        }
        if ( HS )
        {
            gummocap::algorithms::CMAESHumanPartTracker *Alg=new gummocap::algorithms::CMAESHumanPartTracker();
            Alg->setParams ( bodyModel, VS,posEv, CMAESHumanPartTracker::Params ( CMAESLambda,CMAESHSStdDev,NumIterations,NumFunctionFitnessEvaluations,nJointsBodyandHead,nJointsLeftArm,nJointsRightArm,nJointsLeftLeg,nJointsRightLeg,CMAESHSStdDev,CMAESHSStdDev,CMAESHSStdDev,CMAESHSStdDev,CMAESHSStdDev ),_DOFMask,numHsteps );
            return Alg;
        }
        else
        {
            gummocap::algorithms::CMAESHumanTracker *Alg=new gummocap::algorithms::CMAESHumanTracker();
            Alg->setParams ( bodyModel, VS,posEv, CMAESHumanTracker::Params ( CMAESLambda,CMAESStdDev,NumIterations,NumFunctionFitnessEvaluations ),_DOFMask );
            return Alg;
        }

    }
    else if ( type=="ma-cma" )
    {
        if ( NumFunctionFitnessEvaluations<=0.0 )
            NumFunctionFitnessEvaluations=20*48;
        NumIterations=1;

        if ( HS )
        {
            gummocap::algorithms::MACMAChainsHumanPartTracker *Alg=new gummocap::algorithms::MACMAChainsHumanPartTracker();
            Alg->setParams ( bodyModel, VS,posEv, MACMAChainsHumanPartTracker::Params ( NumIterations,NumFunctionFitnessEvaluations,MACMASearchRadius,MACMASearchRadius,MACMALocalIntensity,MACMALGRatio,nJointsBodyandHead,nJointsLeftArm,nJointsRightArm,nJointsLeftLeg,nJointsRightLeg ),_DOFMask,numHsteps );
            return Alg;
        }
        else
        {
            gummocap::algorithms::MACMAChainsHumanTracker *Alg=new gummocap::algorithms::MACMAChainsHumanTracker();
            Alg->setParams ( bodyModel, VS,posEv, MACMAChainsHumanTracker::Params ( NumIterations,NumFunctionFitnessEvaluations,MACMASearchRadius,MACMASearchRadius,MACMALocalIntensity,MACMALGRatio ),_DOFMask );
            return Alg;
        }

    }
    else if ( type=="pso" )
    {
        if ( NumFunctionFitnessEvaluations<=0.0 )
            NumFunctionFitnessEvaluations=5000;
        NumIterations=1;
        if ( HS )
        {
            gummocap::algorithms::PSOHumanPartTracker *Alg=new gummocap::algorithms::PSOHumanPartTracker();
            Alg->setParams ( bodyModel, VS,posEv, PSOHumanPartTracker::Params ( PSONoOfAgents,PSOMaximumVelocity,PSOSearchRadius,PSOInitialWeight,NumFunctionFitnessEvaluations,nJointsBodyandHead,nJointsLeftArm,nJointsRightArm,nJointsLeftLeg,nJointsRightLeg ),_DOFMask,numHsteps );
            return Alg;
        }
        else
        {
            gummocap::algorithms::PSOHumanTracker *Alg=new gummocap::algorithms::PSOHumanTracker();
            Alg->setParams ( bodyModel, VS,posEv, PSOHumanTracker::Params ( PSONoOfAgents,PSOMaximumVelocity,PSOSearchRadius,PSOInitialWeight,NumFunctionFitnessEvaluations ),_DOFMask );
            return Alg;
        }
    }
    else if ( type=="de" )
    {
        if ( NumFunctionFitnessEvaluations<=0.0 )
            NumFunctionFitnessEvaluations=5000;
        NumIterations=1;
        if ( HS )
        {
            gummocap::algorithms::DEHumanPartTracker *Alg=new gummocap::algorithms::DEHumanPartTracker();
            Alg->setParams ( bodyModel, VS,posEv, DEHumanPartTracker::Params ( nJointsBodyandHead,nJointsLeftArm,nJointsRightArm,nJointsLeftLeg,nJointsRightLeg,DEPopulationSize,DEnLPC,DESearchRadius,DESearchRadius,DE_CR,DE_F,NumFunctionFitnessEvaluations ),_DOFMask,numHsteps );
            return Alg;
        }
        else
        {
            gummocap::algorithms::DEHumanTracker *Alg=new gummocap::algorithms::DEHumanTracker();
            Alg->setParams ( bodyModel, VS,posEv, DEHumanTracker::Params ( DEPopulationSize,DEnLPC,DESearchRadius,DESearchRadius,DE_CR,DE_F,NumFunctionFitnessEvaluations ),_DOFMask );
            return Alg;
        }
    }
    else if ( type=="psapf" )
    {
        gummocap::algorithms::PSAPF *Alg=new gummocap::algorithms::PSAPF();
        if ( NumFunctionFitnessEvaluations>0.0 )
            psapf_params.setNEvaluations ( NumFunctionFitnessEvaluations );
        Alg->setParams ( bodyModel, VS,posEv,psapf_params );
        return Alg;
    }

    else
    {
        throw GUException ( "getHumanTracker invalid algorithm descriptor :"+type );
    }
    return NULL;
}


/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
ViewSetProducer * getViewSetProducer() throw ( GUException )
{
    switch ( inputType )
    {
    case 0:
    {
        //open vis file and do initial background learning
        ImageSetVideoReader *VisVideoReader=new ImageSetVideoReader;
        VisVideoReader->connect ( TheInputVideoFile );
        if ( TheBckVisFile=="" )
        {

            BckViewSetProducer *TheVSReader=new BckViewSetProducer();
            //Prepare background
            BckViewSetProducer::Params p;
            if ( TheBckModelFile!="" )
                TheVSReader->loadBckModel ( TheBckModelFile );//loads background from file
            else
                p.setNBckFrames ( TheNBckFrames );
            p.thresHold=TheBckGlobalThresHold;
            if ( ThresHolds.size() !=0 )
                p.thresHolds=ThresHolds;
            //now, go to initial frame

            TheVSReader->setParams ( p,VisVideoReader );
            return TheVSReader;
        }
        else
        {
            BckViewSetProducerVis *TheVSReader=new BckViewSetProducerVis();
            //Prepare background
            BckViewSetProducerVis::Params p;
            //now, go to initial frame
            TheVSReader->setParams ( p,VisVideoReader,TheBckVisFile );
            return TheVSReader;
        }

        //use the extended wrapper


//          return TheVSReader;
    }
    break;

    case 1:
    {
        //open vis file and do initial background learning
        gustereo3::StereoSetVideoReader *StReader=new gustereo3::StereoSetVideoReader();
        StReader->connect ( TheInputVideoFile );
        guscene::view::stereo::StBckViewSetProducer  *TheVSReader=new guscene::view::stereo::StBckViewSetProducer();

        //Prepare background
        guscene::view::stereo::StBckViewSetProducer::Params p;
        if ( TheBckModelFile!="" )
            TheVSReader->loadBckModelMono ( TheBckModelFile );//loads background from file
        else
            p.setNBckFrames ( TheNBckFrames );

        p.thresHoldColor =TheBckGlobalThresHold;
        p.doExtrinsiscTransformation=false;
        //now, go to initial frame
        TheVSReader->setParams ( p,StReader );
        return TheVSReader;
    }
    break;

    default:
        throw GUException ( "getViewSetProducer() invalid input type" );

    }
};
/////////////////////////////////////////////
//
//
/////////////////////////////////////////////
int main ( int argc,char **argv )
{
    try
    {
//         if ((argc<4)||(argc>5)) throw GUException("Usage: human_bodymodel.hbm  video.vis background.bck [out.avi]");
        readArguments ( argc,argv );
        if ( TheInputVideoFile=="" ) throw GUException ( "input required" );
        if ( TheInputModelFile=="" ) throw GUException ( "-m required" );
        if ( TheFitnessFunction=="" ) throw GUException ( "-f required" );
        if ( TheOutAviFile!="" )	  videoEnabled=true;

        //read person model
        ThePersonModel.readFromFile ( TheInputModelFile );
        //prepare the view set reader
        TheViewSetReader=getViewSetProducer();
        TheViewSet=TheViewSetReader->produceNext();
        //chech that external parameters are valid
        for ( unsigned int i=0;i<TheViewSet->size();i++ )
        {
            if ( TheViewSet->get ( i )->getCameraParameters()->Extrinsics.isValid() ==false )
                throw GUException ( "Invalid extrinsics parameters in any of the streams given" );
        }
        if ( ( FF=getFitnessEvaluator ( TheFitnessFunction ) ) ==NULL )
            throw GUException ( "Invalid fitness function" );


        char key;

        //if want to estimate the pose only in a particular frame, then
        //go to the frame and initializa the pose to some nearby one
        if ( PoseEstimation!=-1 )
        {

            if ( !TheDoAnalyzeError ) throw GUException ( "In PoseEstimation mode, the --errorFile option must be provided" );
            cout<<"Skipping to frame "<<PoseEstimation<<endl;
            while ( PoseEstimation!= TheViewSetReader->currentFrameIndex() && TheViewSetReader->eof() ==false )
            {
                TheViewSet=TheViewSetReader->produceNext();
                cout<<"."<<flush;
            }
            if ( TheViewSetReader->eof() ) throw GUException ( "The frame indicated for pose estimation is not correct. The video is not that long" );
            //now, calculate a plausible translation so that the model is near the real position
            //to avoid very difficult situations
            //for that purpose, get the location of the first marker which is at the pelvis and
            //use it as global translation for the model
            gumocap::markers::SeqMarkerSet _SqMarkers;
            _SqMarkers.readFromCSVFile ( TheErrorFile );
            gumocap::markers::MarkerSet *MS=_SqMarkers.getByFrame ( PoseEstimation );
            if ( MS==NULL ) throw GUException ( "There is not groundtruth info for the frame indicated" );
            if ( MS->getNValidMarkers() !=MS->size() ) throw GUException ( "Groundtruth for the frame indicated is not valid" );
            //
            BodyModel<HumanSkeleton> auxModel=ThePersonModel;
            TheCurrenBestPose=ThePersonModel.getSkptr()->getRestPose();
            TheCurrenBestPose.Translation= ( *MS ) [0]-ThePersonModel.getSkptr()->getNodeInfo ( 0 ).xyz;//TheCurrenBestPose.Translation-(*MS)[0];
            auxModel.setPose ( ThePersonModel,TheCurrenBestPose );
        }


//Create the tracker

        TheHumanTracker=getHumanTracker ( strTheHumanTrackerAlgorithm,ThePersonModel,TheViewSet,FF );
        if ( TheDebugFlag )
            TheHumanTracker->enableDebug ( TheDebugFlag );

        if ( TheDoAnalyzeError )
            TheErrorAnalyzer.setParams ( TheErrorFile,*ThePersonModel.getSkptr(),TheErrorFileOut );

        //initializes image composer
        VS2IS ( TheViewSet,VS2IS_TypeOfView )->copy ( &TheISCopy );
//created a single image by joining all of them
        ImgComposer.setParams ( &TheISCopy,320,240 );


        int index=0;
        TheCurrenBestPose=ThePersonModel.getSkptr()->getRestPose();


        GUTimeMark TimerComplete;
        do
        {
            TimerComplete.init();
//ok,lets show the scene
//first, make a copy of the viewset. We'll copy in a ImageSet because   the routines that paints
//the model are based on ImageSets and not ViewSets
            VS2IS ( TheViewSet,VS2IS_TypeOfView )->copy ( &TheISCopy );


            if ( !TheManualModeFlag )
            {
                //now, iterate the algorithm
                cout<<"start"<<endl;
                Timer.init();
                TheHumanTracker->step ( TheViewSet );
                cout<<"Time="<<Timer.end() <<endl;
                TheCurrenBestPose=TheHumanTracker->getCurrentEstimation();
            }
            else
            {
                manualModelManipulation ( TheCurrenBestPose );
            }




            BodyModel<HumanSkeleton> auxModel=ThePersonModel;

            ThePersonModel.setPose ( auxModel,TheCurrenBestPose );
            gumocap::utils::GuPainter::draw ( auxModel,&TheISCopy );
            //meta evaluation
            if ( maxConf!=-1 ||  nCamsConf!=0)
            {
                vector<MetaEvaluator::CameraConfidence> cc=MEvaluator.evaluate ( FF,&auxModel,TheViewSet );
		vector<bool> validViews;
		if ( maxConf!=-1)
		  validViews=MEvaluator.getViewMaskConf(cc,maxConf);
		else
		  validViews=MEvaluator.getViewMasknCams(cc,nCamsConf);
		int nValidViews=0;
		for(unsigned int i=0;i<validViews.size();i++) 
		  if(validViews[i]) 
		    nValidViews++;
 		cout<<"<conf> frame: "<<TheViewSetReader->currentFrameIndex() <<", nviews : "<< nValidViews<<" </conf>"<<endl;
                TheHumanTracker->setValidViews ( validViews );
            }
            //Create a composed image to show
            guimage::Image *composed=ImgComposer ( &TheISCopy );

            //save to avi
            if ( videoEnabled )
            {
                if ( ! NewVideo.isConnected() )
                {
                    NewVideo.setParams ( composed->getWidth(),composed->getHeight(),ImageFormats::BGR,TheVideo_fps );
                    NewVideo.connect ( TheOutAviFile );
                }
                NewVideo.consume ( composed->getImageData() );
            }
            //analyze error
            if ( TheDoAnalyzeError )
            {
                float error=TheErrorAnalyzer.evaluateFrame ( TheCurrenBestPose,TheViewSetReader->currentFrameIndex() );
                cout<<"error in frame ="<<TheViewSetReader->currentFrameIndex() <<" "<<error<<endl;
            }

            //show image
            if ( !TheNoXflag )
            {
                cvNamedWindow ( "image" );
                cvShowImage ( "image",guimage::opencv::getTmp ( composed ) );
                key=cvWaitKey ( timeWaitKey );
                processKey ( key );
            }

//next step
            TheViewSet=TheViewSetReader->produceNext();
            index++;

            cout<<"TOTAL TIME="<<TimerComplete.end() <<endl;
        }
        while ( !TheFinishFlag && !TheViewSetReader->eof() && PoseEstimation==-1 );
        if ( videoEnabled )
            NewVideo.disconnect();

    }
    catch ( std::exception &ex )
    {
        std::cout<<ex.what() <<endl;
        usage();
        return EXIT_FAILURE;
    }
}


///////////////////////////////////////
//
//
///////////////////////////////////////
void printControls()
{
    cout<<"Controls:"<<endl;
    cout<<" h: this info"<<endl;
    cout<<" c: toggles view between color and foreground"<<endl;

}

///////////////////////////////////////
//
//
///////////////////////////////////////
void processKey ( char key )
{
    switch ( key )
    {

    case 'h':
        printControls();
        break;
    case 'c':
        if ( VS2IS_TypeOfView==guscene::utils::ViewSet2ImageSet::FOREGROUND )
            VS2IS_TypeOfView=guscene::utils::ViewSet2ImageSet::COLORINPUT;
        else  VS2IS_TypeOfView=guscene::utils::ViewSet2ImageSet::FOREGROUND;
        break;
    case 'q':
    case 27:
        TheFinishFlag=true;
        break;
    }
};




/**
*/
int  writeVRML ( string path,View *TV,BodyModel<HumanSkeleton> *bodyModel )
{

    ofstream file ( path.c_str() );
    file<< "#VRML V2.0 utf8"<<endl;
    //bool drawAxes=true;
    file<< "Background {   skyColor [ 0.5 0.5 0.5 ] }"<<endl;

    for ( unsigned int part=0;part<bodyModel->size();part++ )
        bodyModel->getMesh ( part ).saveToVRML ( file,false,Point3Df ( 1,0,0 ) );

    file.close();
    //now, add the stereo info
    StereoImage *st=TV->getStereoImage();
    StereoImage copy;
    st->copy ( &copy );
    copy.doExtrinsic3DTransformation();
    StereoImageUtilities::createVRMLModel ( "out.wrl",&copy , ( unsigned char* ) TV->getMaskImage()->imageData );
    char cmd[100];
    sprintf ( cmd,"cat out.wrl >> %s",path.c_str() );
    int res=system ( cmd );

    return res;
}
/**
 */
unsigned int TheCurrentJoint=0;
float TheAngleFactor=M_PI/100.;
float TheMovFactor=0.05;
void processKeyManualManipulation ( char key )
{

    switch ( key )
    {
    case '2':
        TheCurrentJoint++;
        if ( TheCurrentJoint>=TheCurrenBestPose.size() ) TheCurrentJoint=TheCurrenBestPose.size()-1;
        cout<<"Current Joint="<<ThePersonModel.getJoint ( TheCurrentJoint ).name<<endl;
        break;
    case '1':
        TheCurrentJoint--;
        if ( TheCurrentJoint<0 ) TheCurrentJoint=0;
        cout<<"Current Joint="<<ThePersonModel.getJoint ( TheCurrentJoint ).name<<endl;
        break;
    case '+':
        TheMovFactor+=0.005;
        cout<<" TheAngleFactor="<<TheAngleFactor<< " TheMovFactor="<<TheMovFactor<<endl;
        break;
    case '-':
        TheMovFactor-=0.005;
        if ( TheMovFactor<=0 ) TheMovFactor=0.005;
        cout<<" TheAngleFactor="<<TheAngleFactor<< " TheMovFactor="<<TheMovFactor<<endl;
        break;
    case 81:
        TheCurrenBestPose[TheCurrentJoint].y+=TheAngleFactor;
        break;
    case 82:
        TheCurrenBestPose[TheCurrentJoint].x-=TheAngleFactor;
        break;
    case 83:
        TheCurrenBestPose[TheCurrentJoint].y-=TheAngleFactor;
        break;
    case 84:
        TheCurrenBestPose[TheCurrentJoint].x+=TheAngleFactor;
        break;

    case 'q':
        TheCurrenBestPose[TheCurrentJoint].z+=TheAngleFactor;
        break;
    case 'w':
        TheCurrenBestPose[TheCurrentJoint].z-=TheAngleFactor;
        break;
    case 'a':
        TheCurrenBestPose.Translation[0]+=TheMovFactor;
        break;
    case 's':
        TheCurrenBestPose.Translation[0]-=TheMovFactor;
        break;
    case 'd':
        TheCurrenBestPose.Translation[1]+=TheMovFactor;
        break;
    case 'f':
        TheCurrenBestPose.Translation[1]-=TheMovFactor;
        break;
    case 'g':
        TheCurrenBestPose.Translation[2]+=TheMovFactor;
        break;
    case 'h':
        TheCurrenBestPose.Translation[2]-=TheMovFactor;
        break;
    };
}


/**
 */
void manualModelManipulation ( SkeletonPose &pose )
{
//process in manual mode
    bool endOfManualProcess=false;
    TheCurrentJoint=0;
    FF->startFrame ( TheViewSet,TheViewSetReader->currentFrameIndex() );
    while ( !endOfManualProcess )
    {

        VS2IS ( TheViewSet,VS2IS_TypeOfView )->copy ( &TheISCopy );
        BodyModel<HumanSkeleton> auxModel=ThePersonModel;
        ThePersonModel.setPose ( auxModel,TheCurrenBestPose );
        //evaluate the model
        float res=FF->evaluate ( &auxModel, TheViewSet );
        cout<<"Evaluation rsult="<<res<<endl;
        //draw the model and show
        gumocap::utils::GuPainter::draw ( auxModel,&TheISCopy );
        //Create a composed image to show
        guimage::Image *composed=ImgComposer ( &TheISCopy );
        cvNamedWindow ( "image" );
        cvShowImage ( "image",guimage::opencv::getTmp ( composed ) );
        char key=cvWaitKey ( 0 );
        if ( key==32 )  endOfManualProcess=true;
        else processKeyManualManipulation ( key );
//          cout<<"key="<<(int)key<<endl;
    }
    cvSetMouseCallback ( "image",NULL );
}


vector<bool> * createReducedModelMask()
{
    unsigned  int size=ThePersonModel.getSkptr()->getRestPose().size();
    vector<bool> *mask=new   vector<bool>();
    mask->resize ( size );
    //set all to true but these we do not want
    for ( unsigned int i=0;i<mask->size();i++ ) ( *mask ) [i]=true;
    //
    ( *mask ) [HumanSkeleton::Spine_up]=false;
    nJointsBodyandHead--;
    ( *mask ) [HumanSkeleton::Lwrist]=false;
    nJointsLeftArm--;
    ( *mask ) [HumanSkeleton::Rwrist]=false;
    nJointsRightArm--;
    ( *mask ) [HumanSkeleton::Lankle]=false;
    nJointsLeftLeg--;
    ( *mask ) [HumanSkeleton::Rankle]=false;
    nJointsRightLeg--;
    return mask;
}
