#ifndef __LRVTRACK_ALL_HPP__
#define __LRVTRACK_ALL_HPP__
#include <map>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#if CV_VERSION_MAJOR>=3
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/videoio/videoio_c.h>
#include <opencv2/core/utility.hpp>
#endif

#include <opencv2/ml/ml.hpp>
#ifdef LRVTRACK_WITH_CUDA
#include <opencv2/gpu/gpu.hpp>
#elif defined(LRVTRACK_WITH_OPENCL)
#include <opencv2/core/ocl.hpp>
#endif
#include <opencv2/core/ocl.hpp>
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <fstream>
#include "larvaObject.hpp"
#include "lrvTrackFit.hpp"
#include <boost/timer/timer.hpp>
#ifdef _WIN32
#include <winsock.h>
//#include "lrvTrackConfig.h"
#endif

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

//std::string VIDEO_TYPE=".avi";
std::string VIDEO_TYPE=".mp4";
//int VIDEO_CODEC=CV_FOURCC('Z','L','I','B');
//int VIDEO_CODEC=CV_FOURCC('I', '4', '2', '0');
//int VIDEO_CODEC=CV_FOURCC('Y', 'U', 'Y', '2');
//int VIDEO_CODEC=CV_FOURCC('Y', '4', '2', '2');
//int VIDEO_CODEC=CV_FOURCC('Y', 'V', '1', '2');
//int VIDEO_CODEC=CV_FOURCC('U', 'Y', 'V', 'Y');
//int VIDEO_CODEC=CV_FOURCC('I', 'Y', 'U', 'V');
//int VIDEO_CODEC=CV_FOURCC('Y', '8', '0', '0');
//int VIDEO_CODEC=CV_FOURCC('H', 'D', 'Y', 'C');
//int VIDEO_CODEC=CV_FOURCC('H','F','Y','U');
//int VIDEO_CODEC=CV_FOURCC('L','A','G','S');
//int VIDEO_CODEC=CV_FOURCC('I', 'Y', 'U', 'V');
//int VIDEO_CODEC=CV_FOURCC('I', '4', '2', '0');
//int VIDEO_CODEC=CV_FOURCC('J', 'P', 'G', 'L');
//int VIDEO_CODEC=CV_FOURCC('Y', 'U', 'V', '4');
//int VIDEO_CODEC=CV_FOURCC('H','F','Y','U');
//int VIDEO_CODEC=CV_FOURCC('F','F','V','1');
int VIDEO_CODEC=CV_FOURCC('H','2','6','4');
//int VIDEO_CODEC=CV_FOURCC('S','V','Q','3');
//int VIDEO_CODEC=CV_FOURCC('M','P','4','S');
//int VIDEO_CODEC=CV_FOURCC('M','J','P','G');
//int VIDEO_CODEC=CV_FOURCC('F','M','P','4');
//int VIDEO_CODEC=CV_FOURCC('L','M','P','4');
//int VIDEO_CODEC=CV_FOURCC('M','P','4','2');
//int VIDEO_CODEC=CV_FOURCC('A','V','C','1');
//int VIDEO_CODEC=CV_FOURCC('M','P','G','4');
//int VIDEO_CODEC=CV_FOURCC('D','A','V','C');

std::ofstream summary;
std::ofstream csvfile;
int START_FRAME=0;

// *** FLAGS ***
int LRVTRACK_DSTEP=15;
size_t LRVTRACK_PETRIDISH=90;
double LRVTRACK_WSTEP=0.02;
int LRVTRACK_VERBOSE_LEVEL=1;
std::string LRVTRACK_RESULTS_FOLDER;
std::string LRVTRACK_DATE;
std::string LRVTRACK_NAME;
std::string LRVTRACK_SAVE_VIDEO;
std::string LRVTRACK_SAVE_PROCESSED_VIDEO="";
std::string LRVTRACK_FILE_INPUT;
std::string LRVTRACK_INPUT_METADATA;
std::string LRVTRACK_ODOR_LR;
std::string LRVTRACK_ROI_INPUT;
// Larvae File Descriptors for output
std::map<size_t,std::ofstream> lfds;

int  LRVTRACK_CAMERA_INPUT;
int  LRVTRACK_ODOUR_CUPS=0;
int  LRVTRACK_CONTOUR_RESOLUTION=150;
size_t LRVTRACK_MIN_OBJ_SIZE=20;
size_t LRVTRACK_MAX_OBJ_SIZE=2000;
size_t LRVTRACK_THREADS=8;
size_t LRVTRACK_MIN_OUTPUT_DURATION=0;
size_t LRVTRACK_FRAME_WIDTH=0;
size_t LRVTRACK_FRAME_HEIGHT=0;
double LRVTRACK_MODEL_DURATION=3.5;
double LRVTRACK_GAMMA=1.0;
double LRVTRACK_BRIGHTNESS=0.0;
double LRVTRACK_CONTRAST=1.0;
double LRVTRACK_SMOOTHING=2.0;
double LRVTRACK_MPP=0.129;
bool LRVTRACK_USE_MODEL=false;
bool LRVTRACK_INVERT=false;
bool LRVTRACK_PARALLEL=false;
bool LRVTRACK_NORMALIZE=false;
bool LRVTRACK_CHOREOGRAPHY_OUTPUT=false;
bool LRVTRACK_EXTRACT_OFFLINEBG=false;
bool LRVTRACK_EXTRACT_OFFLINEBG_MIN=false;
bool LRVTRACK_CSV_OUTPUT=false;
bool LRVTRACK_SHOW_LIVE=false;
bool LRVTRACK_SHOW_SKELETON=false;
bool LRVTRACK_SHOW_CONTOUR=false;
bool LRVTRACK_SHOW_ORIENTATION=false;
bool LRVTRACK_SHOW_CENTROID=true;
bool LRVTRACK_SHOW_HEAD_TAIL=true;
bool LRVTRACK_SHOW_TAGS=true;
// *** FLAGS ***

//VALUES FOR VARIOUS COMPARISONS
double LARVA_SIZE_COMPARISON_FACTOR=1.3;
double LARVA_CENTRE_COMPARISON_FACTOR=1.05;
double LARVA_MAHALANOBIS_THRESHOLD=2.3;
double LARVA_OBJECT_LENGTH=40;
double IS_LARVA_THRESHOLD=253.2;
double COLLISION_DURATION_THRESHOLD=48;
size_t HISTORY_SIZE=10;

double VIDEO_FPS=16;
size_t CURRENT_FRAME=0;
static size_t TOTAL_FRAMES=0;
size_t LARVAE_COUNT=0;
size_t MAX_LARVAE_DETECTED=0;

size_t FRAME_COLS;
size_t FRAME_ROWS;

std::map<size_t, std::vector<size_t> > parent_blobs;
std::map<size_t, std::vector<size_t> > children_blobs;
std::map<size_t,larvaObject> detected_larvae;
std::map<size_t,std::vector<size_t> > reincarnations;
std::vector<collisionModel> larvaeModels;

std::map<size_t,vector<size_t> > lost_blobs;
std::map<size_t,vector<size_t> > lost_blobs_ring;
std::map<size_t,vector<size_t> > lost_blobs_odor_left;
std::map<size_t,vector<size_t> > lost_blobs_odor_right;
std::map<size_t,vector<size_t> > partial_lost_blobs;
std::vector<size_t> certain_blobs;

std::vector<cv::Vec3f> circles;
std::vector<cv::Vec3f> cups;
cvb::CvBlobs cupBlobs;
cvb::CvBlobs dishBlob;
cv::Mat cupContours;
size_t cupContoursWhitePix;

cvb::CvBlobs NEW;

typedef struct {
  size_t frame;
  bool converging;
  std::vector<size_t> from;
  std::vector<size_t> to;
} event;

std::vector<event> events;

//********* Used by the new tracking algorithm ************************
// map of assignments of the blobs in the previous frame.
// Assignments are such that:
//    [ IDP of Previous Frame, [ IDN1, ..., IDNN]] IDs of
//    new Frame assigned to IDP
std::map<size_t, std::vector<size_t> > assignedPrevious;
std::map<size_t, size_t> assignedPreMap;
std::vector<size_t> newClusters;
std::map<size_t,std::vector<size_t> > newDiverging;

// map of assignments of the blobs in the current frame.
// Assignments are such that:
//    [ IDP of Current Frame, [ID_NN, IDN1, ..., IDNN]] IDs of
//    [ IDP of Current Frame, [ID_NN, IDN1, ..., IDNN]] IDs of
//    Previous Frame assigned to IDP
//    ID NN
std::map<size_t, std::vector<size_t> > assignedNew;
std::vector<size_t> newInFrame;
//**********************************************************************

/*cpu_timer tS;
cpu_timer tP;
cpu_times FrameEllapsedTime;
cpu_times CurrentTime;*/

double PIXEL_SIZE_IN_MM=0.14450867052023;

const uint64_t factorial_vec[21]={1,
                            1,
                            2,
                            6,
                            24,
                            120,
                            720,
                            5040,
                            40320,
                            362880,
                            3628800,
                            39916800,
                            479001600,
                            6227020800,
                            87178291200,
                            1307674368000,
                            20922789888000,
                            355687428096000,
                            6402373705728000,
                            121645100408832000,
                            2432902008176640000};


cv::Mat processedFrame;
cv::Mat greyFrame;
cv::Mat unprocessedFrame;
cv::Mat colorFrame;
cv::Mat bgFrame;
cv::Mat previousFrame;
cv::Mat previousOrigFrame;
cv::Mat toRingTemplate;
IplImage *labelImg;
int DEBUG_INFO=0;

double Wlength=1.0;
double Wsize=0.2;

size_t bestCircle=0;

#endif
