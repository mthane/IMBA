#ifndef __LRVTRACK_BLOBUTILS_HPP__
#define __LRVTRACK_BLOBUTILS_HPP__
#include <opencv2/core/core.hpp>

#if CV_VERSION_MAJOR>=3
#include <opencv2/imgproc/imgproc.hpp>
#endif

#include "cvblob.h"
#define ROI_PADDING 0

unsigned long long getmsofday();
extern double LRVTRACK_SMOOTHING;
extern size_t LRVTRACK_FRAME_WIDTH;
extern size_t LRVTRACK_FRAME_HEIGHT;

void spline3(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di);

void spline4(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &cm,
           std::vector<float> &vcurv,
           double &cvariance);

void spline2(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &cm,
           std::vector<float> &vcurv);

void nospline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv,
           double &cvariance);

void sg_spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv,
           double &cvariance);

namespace lrvTrack{
  
inline  double p2fdist(const double x1,const double y1,const double x2, const double y2)
  {
    double xdiff=x1 - x2;
    double ydiff=y1 - y2;
    return sqrt(xdiff*xdiff + ydiff*ydiff);
  }

template <typename point_type_a, typename point_type_b >
  inline    double p2fdist(const point_type_a &a, const point_type_b &b)
    {
      double xdiff=a.x - b.x;
      double ydiff=a.y - b.y;
      return sqrt(xdiff*xdiff + ydiff*ydiff);
    }

template <typename point_type_a, typename point_type_b >
  inline    double p2fdist(point_type_a &&a, const point_type_b &b)
    {
      double xdiff=a.x - b.x;
      double ydiff=a.y - b.y;
      return sqrt(xdiff*xdiff + ydiff*ydiff);
    }

template <typename point_type_a, typename point_type_b >
  inline    double p2fdist(const point_type_a &a, point_type_b &&b)
    {
      double xdiff=a.x - b.x;
      double ydiff=a.y - b.y;
      return sqrt(xdiff*xdiff + ydiff*ydiff);
    }

template <typename point_type_a, typename point_type_b >
  inline    double p2fdist(point_type_a &&a, point_type_b &&b)
    {
      double xdiff=a.x - b.x;
      double ydiff=a.y - b.y;
      return sqrt(xdiff*xdiff + ydiff*ydiff);
    }
  
}

void blobToBlobOverlap(cvb::CvBlob &blob1,
                       cvb::CvBlob &blob2,
                       double &b1ratio,
                       double &b2ratio);

double diff(cv::Point2f &a, cv::Point2f &b);
void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point2f> &points,size_t PAD=0);

double vecVariance(std::vector<double> &c);

void pointsToContourVector(cvb::CvBlob &blob,
                         std::vector<cv::Point2f> &kpoints,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points);

bool createSimpleROI(cv::Mat &img,
                     size_t minx,
                     size_t miny,
                     size_t maxx,
                     size_t maxy,
                     size_t PADDING,
                     cv::Mat &res);

/*void blobToContourVector(cvb::CvBlob &blob,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points,
                         std::vector<size_t> &baseContourPointIdx);*/

void blobToContourVector(cvb::CvBlob &p,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points);

void tile2same(cv::Mat &a, cv::Mat &b, cv::Mat &r);

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob);

void createBlobContour(cv::Mat &ROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int PADDING=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=8,
                        cv::Scalar bg=cv::Scalar(0));

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        std::vector<std::pair<cv::Point2f,cv::Point2f> > &cpairs,
                        std::vector<cv::Point2f> &spine,
                        int type=CV_8UC1,
                        int padding=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=8,
                        cv::Scalar bg=cv::Scalar(0)
                        );

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int padding=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=8,
                        cv::Scalar bg=cv::Scalar(0)
                        );

void createLarvaContour_custom(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int minx=0,
                        int maxx=0,
                        int miny=0,
                        int maxy=0,
                        int PAD=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=8,
                        cv::Scalar bg=cv::Scalar(0),
			bool verbose=false
    );

/*
void createLarvaContourCV(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int PADDING=0,
                        bool FILL=true);

*/
cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &e, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d);
cv::Point2f px5Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d, cv::Point2f &e);

void smoothPoints(std::vector<cv::Point2f> &origPoints, std::vector<cv::Point2f> &smoothened,cv::Mat &frame, size_t rep);

void spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w, 
           int n, 
           float xp1, 
           float yp1, 
           float xpn, 
           float ypn, 
           std::vector<float> &x2, 
           std::vector<float> &y2);

void splint(std::vector<cv::Point2f> &cp,
            std::vector<float> &x2a, 
            std::vector<float> &y2a, 
            std::vector<float> &d, 
            int n, 
            float t, 
            int khi,
            int klo,
            cv::Point2f &np);

/*void sspline(std::vector<float> &cp,
           std::vector<float> &t,
           std::vector<float> &w,
           int n, 
           float s,
           std::vector<float> &ax, 
           std::vector<float> &bx, 
           std::vector<float> &cx, 
           std::vector<float> &dx);

void ssplint(std::vector<float> &d,
           float t,
           int klo,
           std::vector<float> &ax, 
           std::vector<float> &bx, 
           std::vector<float> &cx, 
           std::vector<float> &dx,
           std::vector<float> &ay, 
           std::vector<float> &by, 
           std::vector<float> &cy, 
           std::vector<float> &dy,
           cv::Point2f &np);
*/

void createLarvaContourPacked(cv::Point &first, 
                              size_t &size,
                              std::string &STR,
                              cvb::CvBlob &blob);

void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type=CV_8UC1,
                              int padding=0);

void lengthAreaPerimeter(double a,double b,double &length, double &width);

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame);

double getPerimeter(cvb::CvBlob &blob);

double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob,cv::Mat &grey_frame);
double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob,cv::Mat &grey_frame,cv::Mat &preFrame);
double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD=0);
double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
double angleD( cv::Point2f pt1, cv::Point2f pt0, cv::Point2f pt2 );
double angleC( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
cv::Point2f perp(cv::Point2f &p1, cv::Point2f &p2, 
                 cv::Point2f &p3, double mag=1.0, 
                 bool right=true);
bool intersection(cv::Point2f o1, cv::Point2f p1, 
                  cv::Point2f o2, cv::Point2f p2, cv::Point2f &r);

void derivVec(std::vector<cv::Point2f> &in,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d);

void derivVec(std::vector<cv::Point2f> &in,
    std::vector<cv::Point2f> &d1);

void deriv2Vec(std::vector<cv::Point2f> &d1,
    std::vector<cv::Point2f> &d2);

void deriv2Vec(
    std::vector<cv::Point2f> &d1,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d2);

void curvVec(std::vector<cv::Point2f> &in,
             std::vector<float> &p,
             std::vector<float> &c);

void curvVec(std::vector<cv::Point2f> &in,
    std::vector<float> &out);

void contourAngles(std::vector<cv::Point2f> &in,
                   std::vector<double> &angles,
                   size_t S=2);

void getBestCurvature(std::vector<float> &c,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax
                      );
void getBestCurvatureS(std::vector<float> &c,
                      std::map<float,size_t> &curvatures,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax,
                      double &variance
                      );
namespace std{
template<typename data>
  data getNeighboursAvg(std::vector<data> &in,int idx, int range, data initval)
  {
    data sum=initval;
    for(int i=(idx-range/2);i<=(idx+range/2);i++)
    {
      int id=i;
      if(i<0)
        id=in.size()+i;
      if(i>(int) in.size()-1)
        id=i-in.size();
      sum+=in[id];
    }
    data retval=sum*(1.0/range);
    return retval;
  }

template<typename data>
  void smoothVecMap(std::vector<data> &in,
                 std::map<data,size_t> &out,
                 std::vector<data> &out2,
                 int range, 
                 data initval)
  {
    for(size_t i=0;i<in.size();i++)
    {
      data ret;
      ret=getNeighboursAvg(in,i,range,initval);
      out[ret]=i;
      out2.push_back(ret);
    }
  }

template<typename data>
  void findLocalMaxima(std::vector<data> &m,std::vector<size_t> &loc,bool maxima=true)
  {
      if(m.size()<3)
        return;

    struct iv{
      data value;
      size_t index;
      iv(data v, size_t i)
      {
        value=v;
        index=i;
      }
    };

    std::vector<iv> d;
    if(maxima==true)
    {
      long idx=0;
      for(idx=0;idx<m.size();idx++)
      {
        //Check previous points
        long p=idx;
        long n=idx;
        while(m[idx]==m[p])
        {
          p--;
          if(p<0)
            p=m.size()+p;
        }
        while(m[idx]==m[n])
        {
          n++;
          if(n>=m.size())
            n=n-m.size();
        }
        if(m[idx]>m[n] && m[idx]>m[p])
        {
          d.emplace_back(iv(m[idx],idx));
        }
      }
      struct by_value_max{
        bool operator()(iv const &a, iv const &b) {
          return a.value > b.value;
        }
      };
      std::sort(d.begin(),d.end(),by_value_max());
    }
    else
    {
      long idx=0;
      for(idx=0;idx<m.size();idx++)
      {
        //Check previous points
        long p=idx;
        long n=idx;
        while(m[idx]==m[p])
        {
          p--;
          if(p<0)
            p=m.size()+p;
        }
        while(m[idx]==m[n])
        {
          n++;
          if(n>=m.size())
            n=n-m.size();
        }
        if(m[idx]<m[n] && m[idx]<m[p])
        {
          d.emplace_back(iv(m[idx],idx));
        }
      }
      struct by_value_min{
        bool operator()(iv const &a, iv const &b) {
          return a.value < b.value;
        }
      };
      std::sort(d.begin(),d.end(),by_value_min());
    }

    for(size_t i=0;i<d.size();i++)
      loc.push_back(d[i].index);
  }

template<typename data>
  data smoothVal(std::vector<data> &in,
                 int idx,int size)
  {
    int start,end;
    start=idx-((size-1)/2);
    end=idx+((size-1)/2);
    if(start<0)
      start=0;
    if(end>in.size()-1)
      end=in.size()-1;

    int d=end-start+1;
    data sum=in[start];
    for(int i=start+1;i<=end;i++)
    {
      sum+=in[i];
    }
    return ((double)(1.0/d))*sum;
  }
}

template<typename data>
  void smoothVec(std::vector<data> &in,
                 std::vector<data> &out, 
                 int range, 
                 data initval)
  {
    for(size_t i=0;i<in.size();i++)
    {
      data ret;
      ret=getNeighboursAvg(in,i,range,initval);
      out.push_back(ret);
    }
  }

#endif
