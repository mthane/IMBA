#ifndef __LRVTRACK_FIT_HPP
#define __LRVTRACK_FIT_HPP
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <vector>
#include "cvblob.h"
#include "larvaObject.hpp"

//We direct our approximation initially
//to the following spine points (p[0-11], 0:head 11:tail):
//  p[8]: main reference point
//  angleg: global angle of p[8] p[6]
//  angles: spine angles at p[2],p[4],p[6],p[8]
//
//  We need points p[0],p[2],p[4],p[6],p[8],p[11]
//  and their corresponding width points.
//
//  Angles are in degrees and anticlockwise, measured
//  from the old vector to the new.


using namespace cv;
using namespace cvb;
using namespace std;

extern size_t FRAME_COLS;
extern size_t FRAME_ROWS;
extern double LRVTRACK_WSTEP;
extern bool LRVTRACK_PARALLEL;
extern Mat greyFrame;
extern Mat greyFrame;
extern std::map<size_t,larvaObject> detected_larvae;

#define LRVFIT_ANGLES 4
#define IPOL 3
#define PADDING 4

class fitData
{
  public:
  Point2f midtail;
  double length;
  double width;
  vector<double> angles;
  double global_angle;

  fitData()
  {};
  fitData(Point2f mp,
      double l,
      double w,
      double a1,
      double a2,
      double a3,
      double a4,
      double ga):
    midtail(mp),
    length(l),
    width(w),
    global_angle(ga)
  {
    angles.push_back(a1);
    angles.push_back(a2);
    angles.push_back(a3);
    angles.push_back(a4);
  }
  fitData(Point2f mp,
      double l,
      double w,
      vector<double> a,
      double ga):
    midtail(mp),
    length(l),
    width(w),
    angles(a),
    global_angle(ga){}
};

class lrvFit
{
  private:
    size_t orig_width;
    size_t orig_length;
    std::vector<Point2f> intSpine;
    vector<double> ga;
    vector<double> a1;
    vector<double> a2;
    vector<double> a3;
    vector<double> a4;
    vector<double> wl;
    double degree_step=15*0.0174532925199432957692369076848;
    double wstep=LRVTRACK_WSTEP;
    size_t PAD=10;
  public:
    std::vector<Point2f> spine;
    std::vector<Point2f> contour;
    std::vector<Point> cpoints;
    size_t ID;
    size_t minx;
    size_t maxx;
    size_t miny;
    size_t maxy;
    fitData larvaFitData;
    void setloops(int ds,double ws);
    void setloops();
    void createContourFromFit();
    void setupSpine();
    void pointToPoint(Point2f &p1, double angle, double d, Point2f &p2);
    void createContourFromFit(Mat &bg);
    void calculateContourPoints(Point2f &a,
        Point2f &b,
        Point2f &c,
        double b_index,
        double width,
        Point &cl,
        Point &cr);

    void calculateContour2f();

    void calculateContourPoint2f(Point2f &a,
        Point2f &b,
        Point2f &c,
        double b_index,
        double width,
        Point2f &cl,
        Point2f &cr);

    void calculateContourPoints(Point2f &a,
        Point2f &b,
        Point2f &c,
        Point2f &bp,
        double b_index,
        double width,
        Point &cl,
        Point &cr);

    void csvLine(size_t CURRENT_FRAME, 
        size_t SFRAME,
        size_t EFRAME,
        size_t VIDEO_FPS, 
        cv::Point2f &cc, 
        double ppm, 
        std::string &csvline);

    void paintPoly(Mat &ROI, std::vector<Point> f,size_t fsize);
    void generate(std::vector<fitData> &l);
    void setup(larvaObject &l,
        size_t minx,
        size_t maxx,
        size_t miny,
        size_t maxy,
        size_t FRAME);
    void setup(larvaObject &l,
        size_t new_minx,
        size_t new_maxx,
        size_t new_miny,
        size_t new_maxy,
        int dstep, double wstep, size_t FRAME=0);
    //~lrvFit();
    lrvFit(){}
    lrvFit(larvaObject &l,int dstep,double wstep,
            size_t minx,size_t maxx, size_t miny, size_t maxy, size_t FRAME=0);
    lrvFit(larvaObject &l,
            size_t minx,size_t maxx, size_t miny, size_t maxy, size_t FRAME=0);
    double errorFunc(Mat &ref);
    double errorFunc(Mat &ref, Mat &bg);
    void showContour();
    double optimize(Mat &ref);
    double optimizeAndReturn(Mat &ref,Mat &ret);
    double optimize(Mat &ref, Mat &bg);
    double optimize(Mat &ref,Mat &bg, Mat &ret);
    void createMatfromFit(Mat &larvaFitContour,
                                Mat &bg,
                                bool verbose=false
                                );
    void createMatfromFit(Mat &larvaFitContour,
                          bool verbose=false
                          );
    void filterAngle(std::vector<double> &a,
        double &val,
        double lim,
        double add);
    void changeBaseCoords(size_t minx,
                          size_t maxx,
                          size_t miny,
                          size_t maxy);
};

ostream& operator<<(ostream& cout,lrvFit &obj);

class collisionModel {
  public:
  bool SUCCESS;
  size_t SFRAME;
  size_t EFRAME;
  void createBg(int exclude,size_t FRAME,Mat &bg);
  std::vector<std::vector<lrvFit> > larvae_models; //larvae_models[larva][frame]
  std::vector<Mat> optimizedFitImages;
  std::vector<double> frameError;
  collisionModel();
  collisionModel(std::vector<size_t> clarvae, cvb::CvBlob &b,
                 size_t START_FRAME,size_t END_FRAME);
  void updateModel(cvb::CvBlob &blob,size_t FRAME,Mat &out);
};

#endif
