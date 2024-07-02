#include <libalglib/stdafx.h>
#include "blobUtils.hpp"
#include <boost/dynamic_bitset.hpp>
#include <libalglib/interpolation.h>
#include <sys/time.h>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "lrvTrackDebug.hpp"
#include <opencv2/highgui/highgui.hpp>
#ifdef LRV_TRACK_VISUAL_DEBUG
#include <opencv2/plot.hpp>
#endif

using namespace boost::accumulators;
using namespace lrvTrack;

double vecVariance(std::vector<double> &c)
{
  accumulator_set<double, stats<tag::variance> > acc;
  for_each(c.begin(), c.end(), boost::bind<void>(boost::ref(acc), _1));

  return sqrt(variance(acc));
}

unsigned long long getmsofday()
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (long long)tv.tv_sec*1000*1000 + tv.tv_usec;
}

bool isLeft(cv::Point2f &a, cv::Point2f &b, cv::Point2f &c){
  return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) > 0;
}
std::string printVector(alglib::real_1d_array &vec)
{
  if (vec.length()==0)
    return "";
  std::stringstream sstm;
  //bool const is_number= std::is_arithmetic<number>::value;
  //static_assert( is_number, "Provided type is not an arithmetic type");
  size_t i=0;
  sstm << "[" ;
  sstm << vec[i] ;
  ++i;
  for( ; i < vec.length(); ++i)
    sstm << ","<< vec[i];
  sstm << "]";
  return sstm.str();
}

/*double p2fdist(double x1,double y1, double x2, double y2)
{
  double xdiff=x1 - x2;
  double ydiff=y1 - y2;
  return sqrt(xdiff*xdiff + ydiff*ydiff);
}*/

/*double p2fdist(cv::Point2f a, cv::Point2f b)
{
  double xdiff=a.x - b.x;
  double ydiff=a.y - b.y;
  return sqrt(xdiff*xdiff + ydiff*ydiff);
}*/

//cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c)
cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &e, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d)
{
  cv::Vec3b va=f.at<cv::Vec3b>(a);
  cv::Vec3b vb=f.at<cv::Vec3b>(b);
  cv::Vec3b vc=f.at<cv::Vec3b>(c);
  /*cv::Vec3b va(1,1,1);
  cv::Vec3b vb(1,1,1);
  cv::Vec3b vc(1,1,1);*/
  
  
  float ua=(va[0]+va[1]+va[2])/3;
  float ub=(vb[0]+vb[1]+vb[2])/3;
  float uc=(vc[0]+vc[1]+vc[2])/3;

  float sum=ua+ub+uc;
  float wa=ua/sum;
  float wb=ub/sum;
  float wc=uc/sum;

  return wa*a+wb*b+wc*c;
  //return b;
}

void derivVec(std::vector<cv::Point2f> &in,
    std::vector<cv::Point2f> &d1)
{
  std::vector<double> idx;
  d1.push_back(in[0]-in.back());
  for(size_t i=1;i<in.size();i++)
  {
    d1.push_back(in[i]-in[i-1]);
  }
}

void derivVec(std::vector<cv::Point2f> &in,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d)
{
  std::vector<double> idx;
  d.push_back(cv::Point(0,0));
  for(size_t i=1;i<in.size();i++)
  {
    double R=1/(p[i]-p[i-1]);
    d.push_back((in[i]-in[i-1])*R);
  }
  d[0]=(in[0]-in.back())*(1/(p[p.size()-1]-p[p.size()-2]));

}

void deriv2Vec(std::vector<cv::Point2f> &d1,
    std::vector<cv::Point2f> &d2)
{
  std::vector<double> idx;
  d1.push_back(d1[0]-d1.back());
  for(size_t i=1;i<d1.size();i++)
  {
    d2.push_back((d1[i]-d1[i-1]));
  }
}

void deriv2Vec(
    std::vector<cv::Point2f> &d1,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d2)
{
  d2.push_back(cv::Point(0,0));
  for(size_t i=1;i<d1.size();i++)
  {
    double R=1/(p[i]-p[i-1]);
    d2.push_back((d1[i]-d1[i-1])*R);
  }
  d2[0]=(d1[0]-d1.back())*(1/(p[p.size()-1]-p[p.size()-2]));
}

void curvVec(std::vector<cv::Point2f> &in,
             std::vector<float> &out)
{

  std::vector<cv::Point2f> d1,d2;
  derivVec(in,d1);
  deriv2Vec(d1,d2);
  for(size_t i=0;i<in.size();i++)
  {
    float val=((float)(d1[i].x*d2[i].y-d1[i].y*d2[i].x))/
        pow(d1[i].x*d1[i].x+d1[i].y*d1[i].y,1.5);
    out.push_back(val);
  }
}

void curvVec(std::vector<cv::Point2f> &in,
             std::vector<float> &p,
             std::vector<float> &c)
{

  std::vector<cv::Point2f> d1,d2;
  derivVec(in,p,d1);
  deriv2Vec(d1,p,d2);
  for(size_t i=0;i<in.size();i++)
  {
    //std::cerr << d1[i].x << "," << d1[i].y << "," << d2[i].x << "," << d2[i].y << std::endl;
    float val=(d1[i].x*d2[i].y-d1[i].y*d2[i].x)/
        pow(d1[i].x*d1[i].x+d1[i].y*d1[i].y,1.5);
    c.push_back(val);
  }
}


cv::Point2f px5Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d, cv::Point2f &e)
{
  /*cv::Vec3b va=f.at<cv::Vec3b>(a);
  cv::Vec3b vb=f.at<cv::Vec3b>(b);
  cv::Vec3b vc=f.at<cv::Vec3b>(c);
  cv::Vec3b vd=f.at<cv::Vec3b>(d);
  cv::Vec3b ve=f.at<cv::Vec3b>(e);
  
  float ua=(va[0]+va[1]+va[2])/3;
  float ub=(vb[0]+vb[1]+vb[2])/3;
  float uc=(vc[0]+vc[1]+vc[2])/3;
  float ud=(vd[0]+vd[1]+vd[2])/3;
  float ue=(ve[0]+ve[1]+ve[2])/3;

  float sum=ua+ub+uc+ue+ud;
  float wa=ua/sum;
  float wb=ub/sum;
  float wc=uc/sum;
  float wd=ud/sum;
  float we=ue/sum;

  return wa*a+wb*b+wc*c+wd*d+we*e;*/
  cv::Point2f NP=0.2*(a+b+c+d+e);
  return NP;
}


void smoothPoints(std::vector<cv::Point2f> &origPoints, std::vector<cv::Point2f> &smoothened,cv::Mat &frame, size_t rep)
{
  std::vector<cv::Point2f> temp1=origPoints;
  for (size_t k=0;k<rep;k++)
  {
    cv::Point2f NP=px3Smooth(frame, temp1[temp1.size()-2], temp1.back(), temp1[0], temp1[1], temp1[2]); 
    smoothened.push_back(NP);
    for(size_t i=1;i<temp1.size();i++)
    {
      if(temp1[i-1]!=temp1[i])
      {
        cv::Point2f NP;
        if(i==1)
          NP=px3Smooth(frame,temp1.back(),temp1[i-1],temp1[i],temp1[i+1],temp1[i+2]);
        else if(i<temp1.size()-2)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[i+1],temp1[i+2]);
        else if (i==temp1.size()-1)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[0],temp1[1]);
        else if (i==temp1.size()-2)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[i+1],temp1[0]);

        if(NP!=smoothened.back())
        {
          smoothened.push_back(NP);
        }
      }
    }
    if(k!=rep-1)
    {
      temp1=smoothened;
      smoothened.clear();
    }
  }
}


double diff(cv::Point2f &a, cv::Point2f &b)
{
  return (fabs((double) a.x-b.x)+fabs((double)a.y-b.y));
}

void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point2f> &points, size_t PAD)
{
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&p.contour);
  cvb::CvContourPolygon::iterator a=cntPoly->begin();
  while(a!=cntPoly->end())
    {
      points.push_back((cv::Point2f)(*a++)-cv::Point2f(PAD,PAD));
    }
  delete cntPoly;
}

void pointsToContourVector(cvb::CvBlob &blob,
                         std::vector<cv::Point2f> &kpoints,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points)
{
  
  cv::Mat ROI;
  if(!createSimpleROI(
      frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      PAD,
      ROI))
    return;

  cv::Point2f PADpoint(blob.minx-PAD,blob.miny-PAD);
  for(size_t i=1;i<kpoints.size();i++)
  {
    cv::LineIterator it(ROI,kpoints[i-1]-PADpoint,kpoints[i]-PADpoint,8);
    for(int j = 1; j < it.count; j++)
    {
      cv::Point2f newPoint;
      newPoint.x=it.pos().x;
      newPoint.y=it.pos().y;
      points.push_back(newPoint+PADpoint);
      ++it;
    }
  }

  cv::LineIterator it(ROI,kpoints.back()-PADpoint,kpoints[0]-PADpoint,8);
  for(int i = 0; i < it.count; i++)
  {
    cv::Point2f newPoint;
    newPoint.x=it.pos().x;
    newPoint.y=it.pos().y;
    points.push_back(newPoint+PADpoint);
    ++it;
  }
}

void blobToContourVector(cvb::CvBlob &blob,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points)
{
  std::vector<cv::Point2f> kpoints;
  blobToPointVector(blob,kpoints);

  cv::Mat ROI;
  if(!createSimpleROI(
      frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      PAD,
      ROI))
    return;

  cv::Point2f PADpoint(blob.minx-PAD,blob.miny-PAD);
  for(size_t i=1;i<kpoints.size();i++)
  {
    cv::LineIterator it(ROI,kpoints[i-1]-PADpoint,kpoints[i]-PADpoint,8);
    for(int j = 1; j < it.count; j++)
    {
      cv::Point2f newPoint;
      newPoint.x=it.pos().x;
      newPoint.y=it.pos().y;
      points.push_back(newPoint+PADpoint);
      ++it;
    }
  }

  cv::LineIterator it(ROI,kpoints.back()-PADpoint,kpoints[0]-PADpoint,4);
  for(int i = 0; i < it.count; i++)
  {
    cv::Point2f newPoint;
    newPoint.x=it.pos().x;
    newPoint.y=it.pos().y;
    points.push_back(newPoint+PADpoint);
    ++it;
  }
}

void blobToBlobOverlap(cvb::CvBlob &blob1,
                       cvb::CvBlob &blob2,
                       double &b1ratio,
                       double &b2ratio)
{
  size_t rminx = std::min(blob1.minx,blob2.minx);
  size_t rminy = std::min(blob1.miny,blob2.miny);
  size_t rmaxx = std::max(blob1.maxx,blob2.maxx);
  size_t rmaxy = std::max(blob1.maxy,blob2.maxy);
  cv::Mat lrvROI1,lrvROI2;
  createLarvaContour_custom(lrvROI1,
                            blob1,
                            CV_8UC1,
                            rminx,
                            rmaxx,
                            rminy,
                            rmaxy);
  size_t b1pixels=countNonZero(lrvROI1);
  createLarvaContour_custom(lrvROI2,
                            blob2,
                            CV_8UC1,
                            rminx,
                            rmaxx,
                            rminy,
                            rmaxy);
  size_t b2pixels=countNonZero(lrvROI2);

  if(lrvROI1.cols != lrvROI2.cols || 
     lrvROI1.rows != lrvROI2.rows)
  {	
  	b1ratio=0;
	b2ratio=0;
	return;
  }
  cv::Mat lrvROI=lrvROI1 & lrvROI2;

  size_t commonPixels=countNonZero(lrvROI);
  b1ratio=(double)commonPixels/(double)b1pixels;
  b2ratio=(double)commonPixels/(double)b2pixels;

}

void lengthAreaPerimeter(double a,double b,double &length, double &width)
{
  double x1=-0.0459441*sqrt(a*a+25.1327*b)+0.0649747*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
  double y1=(0.0530516*(0.028715*a*a*sqrt(a*a+25.1327*b)-(0.866025*a*a*b)/sqrt(a*a+25.1327*b)+0.00574301*pow((a*a+25.1327*b),1.5)-1.01036*b*sqrt(a*a+25.1327*b)-(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)+0.0324874*a*a*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0281349*a*sqrt(a*a+25.1327*b)*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0162437*pow((-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)-2.24537*b*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

  double x2=-0.0459441*sqrt(a*a+25.1327*b)-0.0649747*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
  double y2=(0.0530516*(0.028715*a*a*sqrt(a*a+25.1327*b)-(0.866025*a*a*b)/sqrt(a*a+25.1327*b)+0.00574301*pow((a*a+25.1327*b),1.5)-1.01036*b*sqrt(a*a+25.1327*b)-(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)-0.0324874*a*a*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0281349*a*sqrt(a*a+25.1327*b)*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0162437*pow((-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)+2.24537*b*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

  double x3=0.0459441*sqrt(a*a+25.1327*b)+0.0649747*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
  double y3=(0.0530516*(-0.028715*a*a*sqrt(a*a+25.1327*b)+(0.866025*a*a*b)/sqrt(a*a+25.1327*b)-0.00574301*pow((a*a+25.1327*b),1.5)+1.01036*b*sqrt(a*a+25.1327*b)+(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)+0.0324874*a*a*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0281349*a*sqrt(a*a+25.1327*b)*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0162437*pow(((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)-2.24537*b*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

  double x4=0.0459441*sqrt(a*a+25.1327*b)-0.0649747*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
  double y4=(0.0530516*(-0.028715*a*a*sqrt(a*a+25.1327*b)+(0.866025*a*a*b)/sqrt(a*a+25.1327*b)-0.00574301*pow((a*a+25.1327*b),1.5)+1.01036*b*sqrt(a*a+25.1327*b)+(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)-0.0324874*a*a*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0281349*a*sqrt(a*a+25.1327*b)*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0162437*pow(((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)+2.24537*b*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

  if(x1>0 && y1>0)
  {
    length = std::max(x1,y1);
    width = std::min(x1,y1);
  }
  if(x2>0 && y2>0)
  {
    length = std::max(x2,y2);
    width = std::min(x2,y2);
  }
  if(x3>0 && y3>0)
  {
    length = std::max(x3,y3);
    width = std::min(x3,y3);
  }
  if(x4>0 && y4>0)
  {
    length = std::max(x4,y4);
    width = std::min(x4,y4);
  }

}

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob)
{

  if(!createSimpleROI(
      frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      ROI_PADDING,
      ROI))
    return;

  //cv::normalize(ROI,ROI,0,255,cv::NORM_MINMAX);
}

/* BROKEN!! IF FIXED MAY BE FASTER
void createLarvaContourCV(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type,
                        int PADDING,
                        bool FILL)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);

  lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING)+(2*PADDING),
                        blob.maxx-blob.minx+(2*ROI_PADDING)+(2*PADDING),
                        type);

  cv::Point *ContourPoints[1];
  ContourPoints[0]=(cv::Point*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point)
                   );

  std::vector<cv::Point2f> cntVec;

  for (size_t i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
    ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
    cntVec.push_back(ContourPoints[0][i]);
  }

  sizes[0]=static_cast<int> (cntPoly->size());
  std::vector<std::vector<cv::Point2f> > cvec;
  cvec.push_back(cntVec);

  cv::Scalar color;
  if(type==CV_8UC1)
    color=cv::Scalar(255);
  else
    color=cv::Scalar(255,255,255);

  drawContours(lrvROI,cvec,0,color,1,4,cv::noArray(),INT_MAX,cv::Point(PADDING,PADDING));
}
*/

bool createSimpleROI(cv::Mat &img,
                     size_t minx,
                     size_t miny,
                     size_t maxx,
                     size_t maxy,
                     size_t PADDING,
                     cv::Mat &res)
{
  try{
    if (minx-PADDING>0 &&
        miny-PADDING>0 &&
        maxx+1+PADDING>LRVTRACK_FRAME_WIDTH &&
        maxy+1+PADDING>LRVTRACK_FRAME_HEIGHT)
    {
      res=img(cv::Rect(minx-PADDING,
            miny-PADDING,
            maxx-minx+1+2*PADDING,
            maxy-miny+1+2*PADDING));
    }
    else
    {
      res=img(
          cv::Rect(minx,
            miny,
            maxx-minx+1,
            maxy-miny+1));

          copyMakeBorder(res,
            res,
            PADDING,
            PADDING,
            PADDING,
            PADDING,
            cv::BORDER_CONSTANT,
            cv::Scalar(0));
          }
          }
          catch(...)
          {
          std::cerr << "Error creating ROI" << std::endl;
          return false;
          }
          return true;
}

void tile2same(cv::Mat &a, cv::Mat &b, cv::Mat &r)
{
  r=cv::Mat(a.rows,a.cols+b.cols,a.type());
  cv::Rect r1(0,0,a.cols,a.rows);
  cv::Rect r2(a.cols,0,b.cols,a.rows);
  cv::Mat e1=r(r1);
  cv::Mat e2=r(r2);
  a.copyTo(e1);
  b.copyTo(e2);
}

void createLarvaContour_custom(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type,
                        int minx,
                        int maxx,
                        int miny,
                        int maxy,
                        int PAD,
                        bool FILL,
                        cv::Scalar color,
                        int connectivity,
                        cv::Scalar bg,
			bool verbose)
{
  cv::Mat lrvROIlocal;
  int dxu,dyu,dxb,dyb;
  int MAXBORDER=blob.area;
  /*if(verbose)
  {
	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : blob area: " << blob.area << std::endl;
	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : blob minx: " << blob.minx << std::endl;
	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : blob maxx: " << blob.maxx << std::endl;
	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : blob miny: " << blob.miny << std::endl;
	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : blob maxx: " << blob.maxy << std::endl;
  }*/
  createLarvaContour(lrvROIlocal,blob,type,0,FILL,color,connectivity,bg);
  //if(verbose)
  //{
  //    	  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : lrvROIlocal: " << lrvROIlocal.size() << std::endl;
  //}
  dxu=(int) blob.minx-(int) minx;
  dyu=(int) blob.miny-(int) miny;
  dxb=(int) maxx-(int) blob.maxx;
  dyb=(int) maxy-(int) blob.maxy;
  if(! ( dxu+PAD<0 || dyu+PAD<0 || dxb+PAD<0 || dyb+PAD<0))
  {
    if(dxu>MAXBORDER)
      dxu=0;
    if(dyu>MAXBORDER)
      dyu=0;
    if(dxb>MAXBORDER)
      dxb=0;
    if(dyb>MAXBORDER)
      dyb=0;
    //if(verbose)
    //{
    //  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : cmb: " << dyu << "," << dyb <<"," << dxu <<"," << dxb << " PAD: " << PAD << std::endl;
    //}
    cv::copyMakeBorder(lrvROIlocal,lrvROI,dyu+PAD,dyb+PAD,dxu+PAD,dxb+PAD,cv::BORDER_CONSTANT,cv::Scalar(0));

  }
  else
  {
    if(dxu<0)
    {
      //remove dxu cols from the left of lrvROIlocal
      cv::Mat tmpROI = lrvROIlocal.colRange(-1*dxu,lrvROIlocal.cols);
      lrvROIlocal=tmpROI;
      dxu=0;
    }
    if(dxb<0)
    {
      //remove dxb cols from the right of lrvROIlocal
      cv::Mat tmpROI= lrvROIlocal.colRange(0,lrvROIlocal.cols+dxb);
      lrvROIlocal=tmpROI;
      dxb=0;
    }
    if(dyu<0)
    {
      //remove dyu rows from the top of lrvROIlocal
      cv::Mat tmpROI = lrvROIlocal.rowRange(-1*dyu,lrvROIlocal.rows);
      lrvROIlocal=tmpROI;
      dyu=0;
    }
    if(dyb<0)
    {
      //remove dyb rows from the bottom of lrvROIlocal
      cv::Mat tmpROI = lrvROIlocal.rowRange(0,lrvROIlocal.rows+dyb);
      lrvROIlocal=tmpROI;
      dyb=0;
    }
    if(dxu>MAXBORDER)
      dxu=0;
    if(dyu>MAXBORDER)
      dyu=0;
    if(dxb>MAXBORDER)
      dxb=0;
    if(dyb>MAXBORDER)
      dyb=0;
    //if(verbose)
    //{
    //  std::cerr << "DBG: CREATELARVACONTOUR_CUSTOM : cmb: " << dyu << "," << dyb <<"," << dxu <<"," << dxb << " PAD: " << PAD << std::endl;
    //}
    cv::copyMakeBorder(lrvROIlocal,lrvROI,dyu+PAD,dyb+PAD,dxu+PAD,dxb+PAD,cv::BORDER_CONSTANT,cv::Scalar(0));

  }
}

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        std::vector<std::pair<cv::Point2f,cv::Point2f> > &cpairs,
                        std::vector<cv::Point2f> &spine,
                        int type,
                        int PADDING,
                        bool FILL,
                        cv::Scalar color,
                        int connectivity,
                        cv::Scalar bg)
{
  int sizes[1];

  if (lrvROI.empty())
    lrvROI=cv::Mat(blob.maxy-blob.miny+1+(2*PADDING),
                        blob.maxx-blob.minx+1+(2*PADDING),
                        type,bg);

  cv::Point2f bp(blob.minx-PADDING,blob.miny-PADDING);
  cv::Point *ContourPoints[1];
  size_t ContourSz=2*cpairs.size()+2;
  ContourPoints[0]=(cv::Point*) malloc(
                     ContourSz*sizeof(cv::Point)
                   );

  size_t i=1;
  ContourPoints[0][0]=spine[0]-bp;
  ContourPoints[0][ContourSz/2]=spine.back()-bp;
  cv::Point pl=ContourPoints[0][0];
  cv::Point pr=ContourPoints[0][0];
  for(auto &p: cpairs)
  {
    ContourPoints[0][i]=p.first-bp;
    ContourPoints[0][ContourSz-i]=p.second-bp;
    if(!FILL)
    {
      cv::line(lrvROI,pl,ContourPoints[0][i],color,1,connectivity);
      cv::line(lrvROI,pr,ContourPoints[0][ContourSz-i],color,1,connectivity);
      pl=ContourPoints[0][i];
      pr=ContourPoints[0][ContourSz-i];
    }
    i++;
  }
  if(!FILL)
  {
    cv::line(lrvROI,ContourPoints[0][ContourSz/2],ContourPoints[0][ContourSz/2-1],color,1,connectivity);
    cv::line(lrvROI,ContourPoints[0][ContourSz/2],ContourPoints[0][ContourSz/2+1],color,1,connectivity);
  }


  if(FILL)
  {
    cv::fillPoly(lrvROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        color,
        connectivity
        );
  }
  free(ContourPoints[0]);
}
void createBlobContour(cv::Mat &ROI,
                        cvb::CvBlob &blob,
                        int type,
                        int PADDING,
                        bool FILL,
                        cv::Scalar color,
                        int connectivity,
                        cv::Scalar bg)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);

  if (ROI.empty())
    ROI=cv::Mat(blob.maxy-blob.miny+1+(2*PADDING),
                        blob.maxx-blob.minx+1+(2*PADDING),
                        type,bg);

  cv::Point *ContourPoints[1];
  ContourPoints[0]=(cv::Point*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point)
                   );
  sizes[0]=static_cast<int> (cntPoly->size());

  for (size_t i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x;
    ContourPoints[0][i].y=(*cntPoly)[i].y;
    if(!FILL && i>0)
    {
      cv::line(ROI,ContourPoints[0][i-1],ContourPoints[0][i],color,1,connectivity);
    }
  }
  if(!FILL)
    cv::line(ROI,ContourPoints[0][cntPoly->size()-1],ContourPoints[0][0],color,1,connectivity);


  if(FILL)
  {
    cv::fillPoly(ROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        color,
        connectivity
        );
  }
  free(ContourPoints[0]);
  delete(cntPoly);

}

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type,
                        int PADDING,
                        bool FILL,
                        cv::Scalar color,
                        int connectivity,
                        cv::Scalar bg)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);

  if (lrvROI.empty())
    lrvROI=cv::Mat(blob.maxy-blob.miny+1+(2*PADDING),
                        blob.maxx-blob.minx+1+(2*PADDING),
                        type,bg);

  cv::Point *ContourPoints[1];
  ContourPoints[0]=(cv::Point*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point)
                   );
  sizes[0]=static_cast<int> (cntPoly->size());

  for (size_t i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+PADDING;
    ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+PADDING;
    if(!FILL && i>0)
    {
      cv::line(lrvROI,ContourPoints[0][i-1],ContourPoints[0][i],color,1,connectivity);
    }
  }
  if(!FILL)
    cv::line(lrvROI,ContourPoints[0][cntPoly->size()-1],ContourPoints[0][0],color,1,connectivity);


  if(FILL)
  {
    cv::fillPoly(lrvROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        color,
        connectivity
        );
  }
  free(ContourPoints[0]);
  delete(cntPoly);

}


// Finds the intersection of two lines, or returns false.
// The lines are defined by (o1, p1) and (o2, p2).
bool intersection(cv::Point2f o1, cv::Point2f p1, cv::Point2f o2, cv::Point2f p2,
    cv::Point2f &r)
{
  cv::Point2f x = o2 - o1;
  cv::Point2f d1 = p1 - o1;
  cv::Point2f d2 = p2 - o2;

  float cross = d1.x*d2.y - d1.y*d2.x;
  if (fabs(cross) < /*EPS*/1e-8)
    return false;

  double t1 = (x.x * d2.y - x.y * d2.x)/cross;
  r = o1 + d1 * t1;
  if(r.x<std::min(o1.x,p1.x) || r.x<std::min(o2.x,p2.x)
      || r.x>std::max(o1.x,p1.x) || r.x>std::max(o2.x,p2.x)
      || r.y<std::min(o1.y,p1.y) || r.y<std::min(o2.y,p2.y)
      || r.y>std::max(o1.y,p1.y) || r.y>std::max(o2.y,p2.y))
    return false;
  return true;
}

cv::Point2f perp(cv::Point2f &p1, cv::Point2f &p2, cv::Point2f &p3, double len, bool right)
{
  cv::Point2f b(0,0);
  cv::Point2f v1=p1-p2;
  cv::Point2f v2=p3-p2;
  double l1=sqrt(v1.x*v1.x+v1.y*v1.y);
  double l2=sqrt(v2.x*v2.x+v2.y*v2.y);
  cv::Point2f nv1=v1*(1/l1);
  cv::Point2f nv2=v2*(1/l2);
  cv::Point2f s=nv1+nv2;
  double d=sqrt(s.x*s.x+s.y*s.y);
  if(d<1e-8 || d<1e-8)
  {
    if(right)
    {
      if(nv2.x*nv2.y>0)
        s.y=-nv2.y;
      else
        s.x=-nv2.x;
    }
    else
    {
      if(nv2.x*nv2.y>0)
        s.x=-nv2.x;
      else
        s.y=-nv2.y;
    }
    return p2+len*s;
  }
  else 
  {
    s=s*(1/d);
    if(right)
      return p2+len*s;
    else
      return p2+len*(-1*s);
  }
}


void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type,
                              int PADDING)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);
  lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING)+(2*PADDING),
                        blob.maxx-blob.minx+(2*ROI_PADDING)+(2*PADDING),
                        type);
  cv::Point2f *ContourPoints[1];
  ContourPoints[0]=(cv::Point2f*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point2f)
                   );
  cv::Scalar color;
  if(type==CV_8UC1)
    color=cv::Scalar(255);
  else
    color=cv::Scalar(255,255,255);

  sizes[0]=static_cast<int> (cntPoly->size());
  for (size_t i=0; i<cntPoly->size(); ++i)
    {
      ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
      ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
      cv::circle(lrvROI,
                 ContourPoints[0][i], // circle centre
                 0,       // circle radius
                 color, // color
                 -1);              // thickness
    }
  free(ContourPoints[0]);
  delete(cntPoly);
}

void createLarvaContourPacked(cv::Point &first, 
                              size_t &size,
                              std::string &STR,
                              cvb::CvBlob &blob)
{
  std::vector<cv::Point2f> contourPoints;
  blobToPointVector(blob, contourPoints);
  first.x=(int) contourPoints.back().x;
  first.y=(int) contourPoints.back().y;
  std::stringstream outline;
  size_t acc=0;
  size_t cntSize=0;
  cv::Mat lrvROI;
  createLarvaContour(lrvROI,blob,CV_8UC1,0,false);
  int cnt3=0;

  std::vector<cv::Point2f>::reverse_iterator P=contourPoints.rbegin()+1;
  std::vector<cv::Point2f>::reverse_iterator pP=contourPoints.rbegin();
  for(;pP!=contourPoints.rend();)
  {
    cv::Point a;
    a.x= (int) P->x - blob.minx;
    a.y= (int) P->y - blob.miny;
    cv::Point pre;
    pre.x= (int) pP->x - blob.minx;
    pre.y= (int) pP->y - blob.miny;

    cv::LineIterator it(lrvROI, pre, a, 4);
    ++it;
    std::vector<uchar> buf(it.count);
    cv::Point cur;
    for(int i = 1; i < it.count; i++, ++it)
    {
      acc <<= 2;
      //contourSize++;
      cur.x=(int) it.pos().x;
      cur.y=(int) it.pos().y;
      cv::Point d=cur-pre;
      if(d.x==1 && d.y==0)
      {
        acc |= 0x1;
      }
      else if(d.y==1 && d.x==0)
      {
        acc |= 0x3;
      }
      else if(d.x==-1 && d.y==0)
      {
        acc |= 0x0;
      }
      else if(d.y==-1 && d.x==0)
      {
        acc |= 0x2;
      }
      else
      {
        std::cerr << "Error creating 4-connected bitset!!!" << std::endl;
      }
      cnt3++;
      cntSize++;
      lrvROI.at<uchar>(cur)=100;
      lrvROI.at<uchar>(a)=50;
      if (cnt3==3)
      {
        cnt3=0;
        outline << (uchar) (acc+48);
        acc=0;
      }
      pre.x=cur.x;
      pre.y=cur.y;
    }
    ++pP;
    ++P;
    if(P==contourPoints.rend())
      P=contourPoints.rbegin();
  }

  if (cnt3==1)
  {
    acc <<= 2;
    acc <<= 2;
    acc |= 0x1;
    cnt3=0;
    outline << (char) (acc+48);
  }
  if (cnt3==2)
  {
    acc <<= 2;
    cnt3=0;
    outline << (char) (acc+48);
  }
  STR=outline.str();
  size=cntSize;
}


double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  double dx1 = pt1.x - pt0.x;
  double dy1 = pt1.y - pt0.y;
  double dx2 = pt2.x - pt0.x;
  double dy2 = pt2.y - pt0.y;
  return acos((dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10));
}

// Angle within 0-2*Pi clockwise from last point to first point around middle point
// Y-axes reversed to match opencv Mat axes
double angleD( cv::Point2f pt2, cv::Point2f pt0, cv::Point2f pt1 )
{
  cv::Point2f newPT1=pt1-pt0;
  double arc1=atan2(newPT1.y,newPT1.x);
  double arc2=atan2(pt2.y,pt2.x);
  cv::Point2f newPT2=pt2-pt0;
  cv::Point2f rt;
  rt.x=newPT2.x * cos(-arc1) - newPT2.y*sin(-arc1);
  rt.y=newPT2.x * sin(-arc1) + newPT2.y*cos(-arc1);
  arc2=atan2(rt.y,rt.x);
  if(arc2<0)
    arc2=2*CV_PI+arc2;
  return arc2;
}

double angleC( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  //Absolute Tail/Head Angle
  double MULVAL=180.0/CV_PI;
  cv::Point2f newPT1=pt1-pt0;
  cv::Point2f newPT2=pt2-pt0;
  double d1=newPT1.x*newPT1.x+newPT1.y*newPT1.y;
  double d2=newPT2.x*newPT2.x+newPT2.y*newPT2.y;
  double arc1=MULVAL*acos(newPT1.dot(newPT2)/sqrt(d1*d2));
  return arc1;
}

double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD)
{
  double angle = cvb::cvAngle(blob);

  double x1,y1,x2,y2;
  double cx,cy;
  double lengthLine = MAX(blob->maxx-blob->minx, blob->maxy-blob->miny)/2.;

  cx=blob->centroid.x-blob->minx+PAD;
  cy=blob->centroid.y-blob->miny+PAD;

  x1=cx-lengthLine*cos(angle);
  y1=cy-lengthLine*sin(angle);
  x2=cx+lengthLine*cos(angle);
  y2=cy+lengthLine*sin(angle);
  cv::line(ROIimg,
           cv::Point2f(int(x1),int(y1)),
           cv::Point2f(int(x2),int(y2)),
           cv::Scalar(0,255,0));
  return angle;
}

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame)
{
  cv::Mat ROI;
  //TODO: Fix when the Padding exceeds the image size!!!

  cv::Mat ROIcopy;
  if(!createSimpleROI(grey_frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      ROI_PADDING,
      ROIcopy))
    return 0;

  ROIcopy.copyTo(ROI);
  ROI=ROI&larvaROI;
  //lrvTrackNormalize(ROI, ROI, 0, 255, CV_MINMAX );
  double nz=cv::norm(ROI,cv::NORM_L1);
  return nz;
}

inline bool baseXsmaller(cv::Point2f a, cv::Point2f b)
{
  return(a.x<b.x);
}

inline bool baseYsmaller(cv::Point2f a, cv::Point2f b)
{
  return(a.y<b.y);
}

void sortPoints(std::vector<cv::Point2f> &cp,
                bool sortX)
{
  if(sortX)
  {
    std::sort(cp.begin(),cp.end(),baseXsmaller);
  }
  else
  {
    std::sort(cp.begin(),cp.end(),baseYsmaller);
  }
}

//Angles of the contour points with step 2
void contourAngles(std::vector<cv::Point2f> &in,
                   std::vector<double> &angles,
                   size_t S)
{
  for(size_t i=0; i<in.size();i++)
  {
    long p=i-S;
    long n=i+S;
    if(p<0)
      p=in.size()+p;
    if(n>in.size()-1)
      n=n-in.size();
    angles.push_back(angleD(in[p],in[i],in[n]));
  }
}

void getBestCurvaturePoints(std::vector<cv::Point2f> &np)
{
  double dmin,dmax,amax;
  for(int pidx=0 ; pidx<np.size() ; pidx++)
  {
    int lidx=pidx+1;
    if(lidx>=np.size()) lidx=0;
  }
}

void getBestCurvatureS(std::vector<float> &curv,
                      std::map<float,size_t> &curvatures,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax,
                      double &variance
                      )
{
  //Curvatures doubly smoothened
  std::vector<float> c;
  //Temporary var to store smoothened curvatures
  std::vector<float> c_tmp1, c_tmp2, c_tmp3;// c_tmp4;
  //Size of filter
  //std::cout << printVector(curv) << std::endl; 
  int sf=(curv.size()*0.05);
  //std::cout << printVector(curv) << std::endl;
  smoothVec(curv,c_tmp1,sf,(float)0.0);
  smoothVec(c_tmp1,c_tmp2,sf,(float)0.0);
  smoothVec(c_tmp2,c_tmp3,sf,(float)0.0);
  //smoothVecMap(c_tmp,curvatures,c,sf,(float)0.0);
  
  smoothVecMap(c_tmp3,curvatures,c,sf,(float)0.0);
  
  //map storing maxima: curvature value, index
  std::map<float,std::vector<size_t> > maxima;
  
  //initialize starting from 0
  bool prePos=false;
  //for position 0: if previous is seq then prePos==true
  //(i.e. we're going up...)
  if(c[0]-c.back()>=0)
    prePos=true;

  //set max and min to current position
  double max=curv[0];
  double min=curv[0];
  
  for(size_t i=1;i<c.size();i++)
  {
    //Check for min, max
    if(curv[i]>max)
      max=curv[i];
    if(curv[i]<min)
      min=curv[i];

    //check if we're going up or stable
    if(c[i]-c[i-1]>=0 )
    {
      //if yes, set the prePos to true
      prePos=true;
      continue;
    }
    //if we're not going up
    if(c[i]-c[i-1]<0 )
    {
      //If we have been going up, it means we have a local maximum
      //save it in our maxima
      if(prePos)
        maxima[c[i-1]].push_back(i-1);
      //Set prePos to false (we're not going up anymore)
      prePos=false;
    }
    //mean+=curv[i];
  }
  variance=max-min;
  //mean=mean/curv.size();
  
  /*double s2,s3=0.0;
  for(size_t i=0;i<curv.size();i++)
  {
    s2=s2+(curv[i]-mean)*(curv[i]-mean); 
    s3=s3+(curv[i]-mean); 
  }
  variance = (s2-s3*s3/c.size())/(c.size()-1);*/

  //Check for position 0 again?
  if(c[0]-c.back()<0 )
  {
    if(prePos)
      maxima[c.back()].push_back(c.size()-1);
  }
  //If we didn't find enough maxima quit!!
  if(maxima.size()<2)
  {
    //std::cerr << " NOT ENOUGHT MAXIMA! " << std::endl;
    return;
  }

  std::map<float,std::vector<size_t> >::reverse_iterator f=maxima.rbegin();
  float m1=10.0;
  int i1;
  size_t sz1;
  //Look for largest maxima bigger than the threshold
  //Since map is sorted we don't really need to do this no?
  float minCurv=3.0;
  while(m1>minCurv)
  {
    if(f==maxima.rend())
      return;
    m1=f->first;
    sz1=f->second.size();
    i1=f->second[0];
    f++;
  }
  int i2=i1;
  int sz2;
  float m2;
  i2=i1;
  //Trying to find largest next maxima 
  //that is of distance more than a fourth of the total points
  int minDist = (double)curv.size()/6.0;
  while(abs(i2-i1)<minDist || (c.size()-abs(i2-i1))<minDist)
  {
    if(f==maxima.rend())
      return;
    m2=f->first;
    sz2=f->second.size();
    i2=f->second[0];
    f++;
  }
  double P=0.8;

  int p1l=INT_MIN;
  int p1r=INT_MIN;

  int p2l=INT_MIN;
  int p2r=INT_MIN;

  //Go downwards from i1 looking for a minimum
  for(int i=i1;i>i1-50;i--)
  {
    int j=i;
    if(j<0)
    {
      j=c.size()+i;
    }
    if(c[j]<P*m1)
    {
        p1l=i;
        break;
    }
  }

  //Go upwards from i1 looking for a minimum
  for(int i=i1;i<i1+50;i++)
  {
    int j=i;
    if(j>(int) c.size()-1)
    {
      j=i-c.size();
    }
    if(c[j]<P*m1)
    {
      p1r=i;
      break;
    }
  }
  //from the two minima get middle
  if(p1r!=INT_MIN && p1l!=INT_MIN)
  {
    i1=(p1r+p1l)/2;
    if(i1>=(int) c.size())
      i1=i1-c.size();
    if(i1<0)
      i1=c.size()+i1;
  }

  //same for the other local maximum
  for(int i=i2;i>i2-50;i--)
  {
    int j=i;
    if(j<0)
    {
      j=c.size()+i;
    }
    if(c[j]<P*m2)
    {
      p2l=i;
      break;
    }
  }

  for(int i=i2;i<i2+50;i++)
  {
    int j=i;
    if(i>(int) c.size()-1)
    {
      j=i-c.size();
    }
    if(c[j]<P*m2)
    {
      p2r=i;
      break;
    }
  }

  if(p2r!=INT_MIN && p2l!=INT_MIN)
  {
    i2=(p2r+p2l)/2;
    if(i2>=(int) c.size())
      i2=i2-c.size();
    if(i2<0)
      i2=c.size()+i2;
  }
  //return our two best options
  di[0]=i1;
  dmax[0]=m1;

  di[1]=i2;
  dmax[1]=m2;

}

void getBestCurvature(std::vector<float> &c,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax
                      )
{
  bool positive=true;
  for(size_t i=1;i<c.size();i++)
  {
    if(c[i]-c[i-1]>=0 )
    {
      positive=true;
      continue;
    }
    else if(positive)
    {
      positive=false;
      std::vector<size_t>::iterator dj=di.begin();
      for(std::vector<float>::iterator j=dmax.begin();j!=dmax.end();j++)
      {
        if(c[i-1]>(*j))
        {
          dmax.insert(j,c[i-1]);
          di.insert(dj,i);
          dmax.pop_back();
          di.pop_back();
          dj++;
          break;
        }
        dj++;
      }
    }
  }
}

void extractCentripetal(
    alglib::real_1d_array &x,
    alglib::real_1d_array &y,
    alglib::real_1d_array &ad,
    std::vector<float> &d)
{
  size_t step=1;
  float csqrt=0.0;
  for(size_t i=step;i<x.length();i+=step)
  {
      csqrt+=sqrt(p2fdist(x[i],y[i],x[i-1],y[i-1]));
  }
  csqrt+=sqrt(p2fdist(x[x.length()-1],y[y.length()-1],x[0],y[0]));
  ad[0]=0;
  d.push_back(0);
  size_t i=1;
  for(i=1;i<x.length()-1;i++)
  {
    float newD = ad[i-1] + sqrt(p2fdist(x[i-1],y[i-1],x[i],y[i]))/csqrt;
    ad[i]=newD;
    d.push_back(ad[i]);
  }
  i=x.length()-1;
  ad[i]=1.0;
  d.push_back(1.0);
}

void splint(std::vector<cv::Point2f> &cp,
            std::vector<float> &x2a, 
            std::vector<float> &y2a, 
            std::vector<float> &d,
            float t, 
            int khi,
            int klo,
            cv::Point2f &np)
{
  float b,a;
  
  double h=d[khi]-d[klo];
  a=(d[khi]-t)/h;
  b=(t-d[klo])/h;
  np.x=a*cp[klo].x+b*cp[khi].x+((a*a*a-a)*x2a[klo]+(b*b*b-b)*x2a[khi])*(h*h)/6.0;
  np.y=a*cp[klo].y+b*cp[khi].y+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
void spline3(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           /*int d0,
           int dm,*/
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures)
{
  std::vector<float> x2;
  std::vector<float> y2;
  /*std::cerr << cp.size() << std::endl;
  std::cerr << printVector(cp) << std::endl;
  std::cerr << d.size() << std::endl;
  std::cerr << printVector(d) << std::endl;
  std::cerr << w.size() << std::endl;
  std::cerr << printVector(w) << std::endl;*/

  spline(cp,d,w,n,FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX,x2,y2);

  double t0=0;
  double t;
  double step=(1.0-t0)/RES;
  std::vector<float> dmax(di.size(),0);
  std::vector<float> curvature;
  std::vector<float> pvals;
  std::vector<cv::Point2f> NP;
  size_t didx=0;
  for (int i=0;i<RES;i++)
  {
    t=t0+i*step;
    //double dx,dy;
    //double d2x,d2y;
    //pspline2calc(p,t,x,y);
    cv::Point2f N;
    if(d[didx]<t)
      didx++;
    splint(cp,x2,y2,d,t,didx+1,didx,N);
    NP.push_back(N);
    pvals.push_back(t);
    //alglib::pspline2diff2( p, t, x, dx, d2x, y, dy, d2y);

    //dp.push_back(cv::Point2f(d2x,d2y));
    //double val=(dx*d2y-dy*d2x)/pow(dx*dx+dy*dy,1.5);
    //curvature.push_back(val);
  }
  pvals.push_back(t+step);
  smoothVec(NP,np,13,cv::Point2f(0,0));

  std::vector<float> scurvature;
  curvVec(np,pvals,scurvature);
  //getBestCurvatureS(scurvature,curvatures,di,dmax,cvariance);
}
// Savitzky Golay filtered points
void sg_spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv,
           double &cvariance)
{
  double sigma=LRVTRACK_SMOOTHING;
  cv::Mat G;
  int order=3;
  int sz=cp.size()*0.05*2+1;
  int N = cp.size();
  int W = cp.size() / 16 * 2 + 1;
  std::vector<int> order_range = {0, 1, 2, 3};
  int half_window = (W-1)/2;

  std::vector<int> window_range;
  for(int i=0;i<W+1;i++)
	window_range.push_back(-half_window+i);
  // Order 3, order_range size = 4
  cv::Mat b(W+1,order+1,CV_32F);
  for (int i=0;i<W+1;i++)
  {
	  for(int j=0;j<order+1;j++)
		  b.at<float>(i,j) = pow(window_range[i],order_range[j]);
  }
  /*std::cout << " SAVITZKY GOLAY M(): " << std::endl; 
  std::cout << b << std::endl;
  */	
  cv::Mat inv;
  cv::invert(b,inv,cv::DECOMP_SVD);
  cv::Mat m = inv.row(0);
  //std::cout << " m.cols " << m.cols << std::endl;
  //Mat fm;
  //cv::flip(m,fm,0);
  //std::cout << m << std::endl;
  std::vector<cv::Point2f> tcp;
  std::vector<cv::Point2f> tnp;
  tcp.insert(tcp.end(),cp.begin(),cp.end());
  tcp.insert(tcp.end(),cp.begin(),cp.end());
  tcp.insert(tcp.end(),cp.begin(),cp.end());

  cv::Point anchor(m.cols - 2, -1);
  cv::filter2D(tcp, tnp, -1, m, anchor);
  /*std::cout << "TNP: " << std::endl;*/


  //cv::transpose(cv::getGaussianKernel(sz,sigma,CV_64FC1),G);
  //cv::flip(G,G,0);
  //cv::Point anchor(G.cols - sz -1, G.rows -0 -1);
  //cv::filter2D(tcp,tnp,-1,G,anchor);
  //cv::filter2D(tcp,tnp,-1,G);
  np.insert(np.begin(),tnp.begin()+cp.size(),tnp.begin()+2*cp.size());
  //std::cout << "CP SIZE: " << cp.size() << std::endl;
  //std::cout << "NP SIZE: " << np.size() << std::endl;
  //std::cout << printVector(cp) << std::endl;
  //std::cout << printVector(np) << std::endl;

  /*for (int i=0;i<np.size();i++)
  {
	std::cout << np[i] << " ";
  }
  std::cout << std::endl;*/
  
  curvVec(np,vcurv);
  std::vector<float> dmax(di.size(),0);
  getBestCurvatureS(vcurv,curvatures,di,dmax,cvariance);

}

void nospline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv,
           double &cvariance)
{
  double sigma=LRVTRACK_SMOOTHING;
  cv::Mat G;
  int sz=cp.size()*0.05*2+1;
  std::vector<cv::Point2f> tcp;
  std::vector<cv::Point2f> tnp;
  tcp.insert(tcp.end(),cp.begin(),cp.end());
  tcp.insert(tcp.end(),cp.begin(),cp.end());
  tcp.insert(tcp.end(),cp.begin(),cp.end());
  cv::transpose(cv::getGaussianKernel(sz,sigma,CV_64FC1),G);
  //cv::flip(G,G,0);
  //cv::Point anchor(G.cols - sz -1, G.rows -0 -1);
  //cv::filter2D(tcp,tnp,-1,G,anchor);
  cv::filter2D(tcp,tnp,-1,G);
  np.insert(np.begin(),tnp.begin()+cp.size(),tnp.begin()+2*cp.size());
  //std::cout << "CP SIZE: " << cp.size() << std::endl;
  //std::cout << "NP SIZE: " << np.size() << std::endl;
  //std::cout << printVector(cp) << std::endl;
  //std::cout << printVector(np) << std::endl;
  curvVec(np,vcurv);
  std::vector<float> dmax(di.size(),0);
  getBestCurvatureS(vcurv,curvatures,di,dmax,cvariance);

}

void spline4(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv,
           double &cvariance)
{
  alglib::real_1d_array x,y;                 
  alglib::real_1d_array ad,aw;               
  alglib::ae_int_t info;                     
  alglib::ae_int_t m=n;
  alglib::real_1d_array nu;
  alglib::integer_1d_array inu;
  inu.setlength(0);
  nu.setlength(0);
  //alglib::ae_int_t ink=0;                     
  size_t extra=5;
  x.setlength(m+2*extra);                          
  y.setlength(m+2*extra);                          
  ad.setlength(m+2*extra);                         
  aw.setlength(m+2*extra);                         
  alglib::spline1dinterpolant sx;            
  alglib::spline1dinterpolant sy;            
  alglib::spline1dfitreport rep;             
  double rho=LRVTRACK_SMOOTHING/10;
  alglib::pspline2interpolant p;
  std::vector<cv::Point2f> NP;
  for(int i=0;i<extra;i++)
  {
    x[i]=cp[cp.size()-extra+i].x;
    y[i]=cp[cp.size()-extra+i].y;
    //aw[i]=w[cp.size()-extra+i];
    aw[i]=1e-15;
  }
  for(int i=0;i<m;i++)
  {
    x[i+extra]=cp[i].x;
    y[i+extra]=cp[i].y;
    aw[i+extra]=w[i];
    //aw[i]=1e-15;
  }
  for(int i=0;i<extra;i++)
  {
    x[extra+m+i]=cp[i].x;
    y[extra+m+i]=cp[i].y;
    aw[i+m+extra]=w[i];
    aw[i]=1e-15;
  }
  extractCentripetal(x,y,ad,d);
  try{
    alglib::spline1dfitpenalized(ad,x,m+2*extra,rho,info,sx,rep);
    alglib::spline1dfitpenalized(ad,y,m+2*extra,rho,info,sy,rep);
    //alglib::spline1dfitpenalizedw(ad,x,aw,m+2*extra,rho,info,sx,rep);
    //alglib::spline1dfitpenalizedw(ad,y,aw,m+2*extra,rho,info,sy,rep);
    //alglib::spline1dfitcubicwc(ad,x,aw,m+extra,nu,nu,inu,ink,m+extra+2,info,sx,rep);
    //alglib::spline1dfitcubicwc(ad,y,aw,m+extra,nu,nu,inu,ink,m+extra+2,info,sy,rep);
  }
  catch(alglib::ap_error &e)
  {
    std::cerr << "Spline Error:" << e.msg << std::endl;
    std::cerr << printVector(x) << std::endl;
    std::cerr << printVector(y) << std::endl;
    std::cerr << printVector(ad) << std::endl;
    return;
  }
  //alglib::spline1dfithermite(ad,x,m+extra,info,sx,rep);
  //alglib::spline1dfithermite(ad,y,m+extra,info,sy,rep);
  
  double t0=ad[extra];
  double tn=ad[extra+m];
  double t;
  double step=(tn-t0)/(RES);
  std::vector<float> dmax(di.size(),0);
  std::vector<float> curvature;
  std::vector<float> pvals;
  std::vector<double> adfactor;
  adfactor.push_back(step);
  size_t j=0;
  for (int i=0;i<RES;i++)
  {
    //step=adfactor[j];
    t=t0+i*step;
    double vx,vy;
    vx = spline1dcalc(sx, t);
    vy = spline1dcalc(sy, t);
    cv::Point2f N(vx,vy);
    np.push_back(N);
    pvals.push_back(t);
    while(t<ad[j])
      j++;
  }
  pvals.push_back(t+step);
  curvVec(np,pvals,vcurv);
  getBestCurvatureS(vcurv,curvatures,di,dmax,cvariance);
}

void spline2(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &curvatures,
           std::vector<float> &vcurv)
{
  alglib::real_2d_array xy;
  xy.setlength(n+2,2);
  //int rot=10;
  //xyrot.setlength(n+2,2);
  alglib::pspline2interpolant p;
  std::vector<cv::Point2f> NP;
  for(int i=0;i<n;i++)
  {
    xy[i][0]=cp[i].x;
    xy[i][1]=cp[i].y;
  }
  xy[n][0]=xy[0][0];
  xy[n][1]=xy[0][1];
  xy[n+1][0]=xy[1][0];
  xy[n+1][1]=xy[1][1];

  alglib::pspline2buildperiodic(xy,n,1,2,p);
  
  double t0=0;
  double t;
  double step=(1.0-t0)/(RES);
  std::vector<float> dmax(di.size(),0);
  std::vector<float> curvature;
  std::vector<float> pvals;
  for (int i=0;i<RES;i++)
  {
    t=t0+i*step;
    double x,y;
    alglib::pspline2calc(p,t,x,y);
    cv::Point2f N(x,y);
    NP.push_back(N);
    pvals.push_back(t);
  }
  pvals.push_back(t+step);
  smoothVec(NP,np,5,cv::Point2f(0,0));

  curvVec(np,pvals,vcurv);
  double variance;
  getBestCurvatureS(vcurv,curvatures,di,dmax,variance);
}

void spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n, 
           float xp1, 
           float yp1, 
           float xpn, 
           float ypn, 
           std::vector<float> &x2, 
           std::vector<float> &y2)
{
  int i,k;
  float px,sig,qnx,unx;
  float py,qny,uny;
  
  std::vector<float> yu;
  std::vector<float> xu;
  yu.reserve(n);
  xu.reserve(n);
  x2.reserve(n);
  y2.reserve(n);
  if(xp1 > 0.99e30)
  {
    x2[0]=0.0;
    xu[0]=0.0;
  }
  else
  {
    x2[0]=-0.5;
    xu[0]=(3.0*(d[1]-d[0]))*(w[0]*(cp[1].x-cp[0].x)/(d[1]-d[0])-xp1);
  }

  if(yp1 > 0.99e30)
  {
    y2[0]=0.0;
    yu[0]=0.0;
  }
  else
  {
    y2[0]=-0.5;
    yu[0]=(3.0*(d[1]-d[0]))*(w[0]*(cp[1].y-cp[0].y)/(d[1]-d[0])-yp1);
  }

  for(i=1;i<=n-2;i++)
  {
    sig=(d[i]-d[i-1])/(d[i+1]-d[i-1]);
    px=sig*x2[i-1]+2.0;
    py=sig*y2[i-1]+2.0;
    x2[i]=((sig-1.0)/px);
    y2[i]=((sig-1.0)/py);

    xu[i]=w[i]*(cp[i+1].x-cp[i].x)/(d[i+1]-d[i]) - w[i]*(cp[i].x - cp[i-1].x)/(d[i]-d[i-1]);
    yu[i]=w[i]*(cp[i+1].y-cp[i].y)/(d[i+1]-d[i]) - w[i]*(cp[i].y - cp[i-1].y)/(d[i]-d[i-1]);

    xu[i]=(6.0*xu[i]/(d[i+1]-d[i-1])-sig*xu[i-1])/px;
    yu[i]=(6.0*yu[i]/(d[i+1]-d[i-1])-sig*yu[i-1])/py;
  }
  if (xpn>0.99e30)
    qnx=unx=0.0;
  else
  {
    qnx=0.5;
    unx=(3.0/(d[n-1]-d[n-2]))*(xpn-(w[n-1]*cp[n-1].x-cp[n-2].x)/(d[n-1]-d[n-2]));
  }
  
  if (ypn>0.99e30)
    qny=uny=0.0;
  else
  {
    qny=0.5;
    uny=(3.0/(d[n-1]-d[n-2]))*(ypn-(w[n-1]*cp[n-1].y-cp[n-2].y)/(d[n-1]-d[n-2]));
  }

  x2[n-1]=(unx-qnx*xu[n-2])/(qnx*x2[n-1]+1.0);
  y2[n-1]=(uny-qny*yu[n-2])/(qny*y2[n-1]+1.0);
  for(k=n-2;k>=0;k--)
  {
    x2[k]=x2[k]*x2[k+1]+xu[k];
    y2[k]=y2[k]*y2[k+1]+yu[k];
  }
}

void fin(int n, 
         int m2,
         float h,
         float p,
         std::vector<float> &t,
         std::vector<float> &a,
         std::vector<float> &b,
         std::vector<float> &c,
         std::vector<float> &d,
         std::vector<float> &u,
         std::vector<float> &v,
         std::vector<float> &cp)
{
  for (int i=1;i<=n;i++)
  {
    a[i]=cp[i]-p*v[i];
    c[i]=u[i];
  }
  for (int i=1;i<=m2;i++)
  {
    h=t[i+1]-t[i];
    d[i]=(c[i+1]-c[i])/(3*h);
    b[i]=(a[i+1]-a[i])/h-(h*d[i]+c[i])*h;
  }
}

/*void sspline(std::vector<cv::Point2f> &cp,
           std::vector<float> &t,
           std::vector<float> &w,
           int n, 
           float s,
           newPoints np)

{
  float* xpointsx;
  float* xpointsy;
  xpointsx=(float *)malloc(n*sizeof(float));
  xpointsy=(float *)malloc(n*sizeof(float));

  for(size_t i=0;i<xpoints.size();i++)
  {
    xpointsx[i]=xpoints[i].x;
    xpointsy[i]=xpoints[i].y;
  }
  float **v;
  float **a;
  v=(float **)malloc(n*sizeof(float *));
  a=(float **)malloc(n*sizeof(float *));
  for(size_t i=0;i<xpoints.size();i++)
  {
    v[i]=(float *)malloc(7*sizeof(float));
    a[i]=(float *)malloc(4*sizeof(float));
  }

  smooth_(&x,&t[0],&w[0],&n,&s,&v,&a);
}*/

double getPerimeter(cvb::CvBlob &blob)
{
  std::vector<cv::Point2f> cntPoints;
  blobToPointVector(blob,cntPoints);
  return arcLength(cntPoints, true);
}

double getSurroundingSize(cv::Point2f &point, 
    cvb::CvBlob &blob, 
    cv::Mat &grey_frame)
{
  cv::Mat larvaImg,lrvROI;
  size_t PADDING=4;
  cv::Mat rROI,ROI;
  if(!createSimpleROI(grey_frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      PADDING,
      rROI))
    return 0;
  rROI.copyTo(ROI);
  //cvtColor(rROI,ROI,CV_BGR2GRAY);

  createLarvaContour(lrvROI, blob,CV_8UC1,PADDING);
  //createLarvaContour(contourBlack, blob,CV_8UC1,PADDING);

  cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
  //cv::imshow("Original",lrvROI);
  //cv::waitKey(-1);
  cv::dilate(lrvROI,lrvROI,element);
  //cv::imshow("Eroded",lrvROI);
  //cv::waitKey(-1);
  cv::bitwise_and(ROI,lrvROI,ROI);
  //cv::equalizeHist(ROI, ROI);
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  size_t xd=blob.maxx-blob.minx;
  size_t yd=blob.maxy-blob.miny;
  size_t sd = sqrt(xd*xd+yd*yd)/4;
  cv::circle(cROI, cv::Point2f(point.x+PADDING,point.y+PADDING),sd,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
  //cv::imshow("headtail",area);
  //cv::waitKey(-1);
  //cv::equalizeHist( area, area);
  //lrvTrackNormalize(area,area,0,255,cv::NORM_MINMAX);
  double nz=cv::norm(area,cv::NORM_L1);
  double nc=cv::countNonZero(area);
  return nz/nc;
}

double getSurroundingSize(cv::Point2f &point, 
    cvb::CvBlob &blob, 
    cv::Mat &grey_frame,
    cv::Mat &prev_frame)
{
  cv::Mat larvaImg,lrvROI;
  size_t PADDING=4;
  cv::Mat rROI,ROI;
  cv::Mat preROI;
  if(!createSimpleROI(grey_frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      PADDING,
      rROI))
    return 0;
  cvtColor(rROI,ROI,CV_BGR2GRAY);

  if(!prev_frame.empty())
  { 
    if(!createSimpleROI(prev_frame,
      blob.minx,
      blob.miny,
      blob.maxx,
      blob.maxy,
      0,
      preROI))
    return 0;

    cv::copyMakeBorder(preROI,
        preROI,
        PADDING,
        PADDING,
        PADDING,
        PADDING,
        cv::BORDER_CONSTANT,cv::Scalar(255));
  }
  createLarvaContour(lrvROI, blob,CV_8UC1,PADDING);
  //createLarvaContour(contourBlack, blob,CV_8UC1,PADDING);
  cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
  cv::erode(lrvROI,lrvROI,element);
  //cv::erode(preROI,preROI,element);
  //lrvTrackNormalize(lrvROI, lrvROI, 0, 255, CV_MINMAX );
  cv::bitwise_and(ROI,lrvROI,ROI);
  //cv::bitwise_and(preROI,lrvROI,preROI);
  /*
  cv::Mat dbg;
  ROI.copyTo(dbg);
  cv::circle(dbg,
      cv::Point2d(point.x+PADDING,point.y+PADDING),
      0,
      cv::Scalar(255),
      -1);

  cv::resize(dbg,dbg,cv::Size(),8,8,cv::INTER_NEAREST);
  cv::imshow("headtail",dbg);
  cv::waitKey(1);
  */
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  cv::circle(cROI, cv::Point2f(point.x+PADDING,point.y+PADDING),12,cv::Scalar(255),-1);
  //cv::Mat test=ROI&cROI;
  cv::Mat area=ROI&cROI;
  //cv::Mat prearea=preROI&cROI;
  cv::equalizeHist( area, area);
  //cv::equalizeHist( prearea, prearea);
  //lrvTrackNormalize(area,area,0,255,cv::NORM_MINMAX);
  //lrvTrackNormalize(prearea,prearea,0,255,cv::NORM_MINMAX);
  //cv::adaptiveBilateralFilter(test,area,cv::Size(5,5),10);
  //cv::bilateralFilter(test,area,4,8,2);
  /*
  cv::resize(area,dbg,cv::Size(),8,8,cv::INTER_NEAREST);
  cv::imshow("headtail",dbg);
  cv::waitKey(1);
  */
  double nz=cv::norm(area,cv::NORM_L1);
  double pnz=0; //cv::norm(prearea,cv::NORM_L1);
  double nc=cv::countNonZero(area);
  double pnc=0;//cv::countNonZero(prearea);
  return (nz+pnz)/(nc+pnc);
  //return nz;
}
