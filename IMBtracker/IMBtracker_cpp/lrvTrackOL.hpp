#ifndef __LRVTRACKOL_HPP__
#define __LRVTRACKOL_HPP__
#include <string>
#include <utility>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <numeric>
#include <vector>
#include <string>
#include "lrvTrack.hpp"
#include "cvblob.h"
#include "lrvTrackBase.hpp"
#include "blobUtils.hpp"
#include "larvaDistanceMap.hpp"
#include "lrvTrackDebug.hpp"
#include "lrvTrackFit.hpp"
#include <fstream>
#ifndef _WIN32
#include <sys/time.h>
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace cv;
using namespace std;

void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<size_t, larvaObject> &NEW_LARVA);

double mh_dist(size_t N,size_t C);
double kn_dist(size_t N,size_t C);

/*
 * Class to perform the computations for each larva.
 * The class is needed to perform these operations in parallel.
 * The format is based on this answer:
 * http://answers.opencv.org/question/3730/how-to-use-parallel_for/
 */
class larvaeUpdateBody : public ParallelLoopBody
{
private:
  cvb::CvBlobs &In;
  cvb::CvBlobs &Prev;
  tbb::concurrent_hash_map<size_t, larvaObject> &NEW_LARVA;
  cvb::CvBlobs::iterator it;

public:
  larvaeUpdateBody(cvb::CvBlobs &IIn,
      cvb::CvBlobs &IPrev,
      tbb::concurrent_hash_map<size_t, larvaObject> &n
      ):
    In(IIn),
    Prev(IPrev),
    NEW_LARVA(n)
  {}
    void operator ()(const Range& range) const
    {
        for (int i = range.start; i < range.end; ++i)
        {
          cvb::CvBlobs::iterator it=In.begin();
          advance(it,i);
            updateOneLarva(In,Prev,it,NEW_LARVA);
        }
    }
#if CV_VERSION_MAJOR<3
    void operator ()(const BlockedRange& range) const
    {
        for (int i = range.begin(); i < range.end(); ++i)
        {
          cvb::CvBlobs::iterator it=In.begin();
          advance(it,i);
            updateOneLarva(In,Prev,it,NEW_LARVA);
        }
    }
#endif
};


/*
 * class to encode the possible mappings for the larvae before and after a collision.
 */
class lrvMapping {

    double dst;

  public:
    //The actual mapping
    // newID -> detectedLarvaeID
    pair<size_t,size_t> mapping;
    //vector<size_t> candidates;
    size_t nlrv;
    size_t plrv;
    //Constructor function sets up the mapping
    lrvMapping(size_t a,size_t b)
    {
      mapping=make_pair(a,b);
      nlrv=mapping.first;
      plrv=mapping.second;
      dst=-1;
    }
    void setDistance()
    {
       //dst=kn_dist(nlrv,plrv);
       dst=mh_dist(nlrv,plrv);
    }
    void print()
    {
      cerr << "M[" << nlrv << "->" << plrv << "]" << " #" << dst << endl;
    }
    double getDistance()
    {
      if(dst==-1)
        dst=mh_dist(nlrv,plrv);
        //dst=kn_dist(nlrv,plrv);

      return dst;
    }
};

// Vector of mappings: covers the possible assignments.
typedef vector<lrvMapping> ltAssignments;

void pair_powersets(vector<lrvMapping> &IN,
    vector<ltAssignments > &OUT);


#endif
