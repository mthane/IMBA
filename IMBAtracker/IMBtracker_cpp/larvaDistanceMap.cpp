#include "larvaDistanceMap.hpp"
#include "lrvTrackDebug.hpp"
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <utility>
#include <opencv2/highgui/highgui.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <opencv2/video/tracking.hpp>
//#include "gnuplot_i.hpp"
#include <libalglib/ap.h>
#include <libalglib/stdafx.h>
#include <libalglib/interpolation.h>
//#define LRV_TRACK_VISUAL_DEBUG
#ifdef LRV_TRACK_VISUAL_DEBUG
#include <opencv2/plot.hpp>
#endif

using namespace lrvTrack;

void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
  std::cout << std::endl << "Press any key to continue..." << std::endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    std::cout << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}


void lBFS(int p1,
          std::vector<cv::Point2f> &Points ,
          larvaDistanceMap &Distances
         )
{
  std::queue<int> Q;
  Q.push(p1);
  double MAX=Distances.MaxDist;
  while (!Q.empty())
    {
      int cur=Q.front();
      Q.pop();
      for(size_t i=0; i<Points.size() ; ++i)
        {
          //PointPair cur_pi = PointPair(Points[cur],Points[i]);
          PointPair p1_pi= PointPair(Points[p1],Points[i]);
          //PointPair p1_cur= PointPair(Points[p1],Points[cur]);
          if (Distances[cur][i]>0)
            {
              if (Distances[p1][i]<0)
                Q.push(i);
              float multres=Distances[cur][i]*Distances[p1][cur];
              float sqrtres;
              sqrtres=sqrt(multres);
              double newDst = Distances[cur][i]+Distances[p1][cur] +
                              2*sqrtres;
              if (Distances[p1][i]>newDst)
                {
                  Distances[p1][i]=newDst;
                  Distances[i][p1]=newDst;
                }
              if (Distances[p1][i]>MAX)
                {
                  MAX=Distances[p1][i];
                  Distances.MaxDistPoints=p1_pi;
                  Distances.MaxDist=MAX;
                }
            }
        }
    }
}

void pointWithDist(
                   cv::Point2f &a, 
                   cv::Point2f &b, 
                   cv::Point2f &ret,
                   double Dist)
{
  double abDist=p2fdist(a,b);
  double l=Dist/abDist;
  ret.x=a.x+(l*(b.x-a.x));
  ret.y=a.y+(l*(b.y-a.y));
}

void assignPoint(cv::Point2f &a,cv::Point2f &b)
{
  a.x=b.x;
  a.y=b.y;
}

bool distanceMatch(double dist1, double dist2)
{
  double E=1e-6;
  if(dist1>dist2-E && dist1<dist2+E)
    return true;
  else
    return false;
}

void showPoints(cvb::CvBlob blob,
                std::vector<cv::Point2f > points,
                std::vector<cv::Point2f > curSpine,
                size_t i,
                size_t fwdIdx,
                size_t bwdIdx,
                size_t pfwdIdx,
                size_t pbwdIdx,
                cv::Point2f bwdPoint,
                cv::Point2f fwdPoint,
                size_t PADDING
                )
{
  cv::Mat contour;
  cv::Mat lcontour;

  createLarvaContour(contour,blob,CV_8UC3,PADDING);
  cv::circle(contour,
      cv::Point2d(points[i].x-blob.minx+PADDING,points[i].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(points[bwdIdx].x-blob.minx+PADDING,points[bwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[pbwdIdx].x-blob.minx+PADDING,points[pbwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[fwdIdx].x-blob.minx+PADDING,points[fwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,250,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[pfwdIdx].x-blob.minx+PADDING,points[pfwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,250,0),
      -1);
  cv::circle(contour,
      cv::Point2d(bwdPoint.x-blob.minx+PADDING,bwdPoint.y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(fwdPoint.x-blob.minx+PADDING,fwdPoint.y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(curSpine.back().x-blob.minx+PADDING,curSpine.back().y-blob.miny+PADDING),
      0,
      cv::Scalar(0,0,255),
      -1);
  cv::resize(contour,lcontour,cv::Size(),8,8,cv::INTER_NEAREST);
  //cv::imshow("Spine3",lcontour);
  cv::waitKey(1);

  std::cerr << "Point[i] : " << points[i] << std::endl;
  std::cerr << "Point[bwdIdx] : " << points[bwdIdx] << std::endl;
  std::cerr << "Point[pbwdIdx] : " << points[pbwdIdx] << std::endl;
  std::cerr << "Point[fwdIdx] : " << points[fwdIdx] << std::endl;
  std::cerr << "Point[pfwdIdx] : " << points[pfwdIdx] << std::endl;
  std::cerr << "bwdPoint : " << bwdPoint << std::endl;
  std::cerr << "fwdPoint : " << fwdPoint << std::endl;
  std::cerr << "curSpine.back() : " << curSpine.back() << std::endl;
  std::cerr << "i: " << i << std::endl;
  std::cerr << "bwdIdx: " << bwdIdx << std::endl;
  std::cerr << "pbwdIdx: " << pbwdIdx << std::endl;
  std::cerr << "fwdIdx: " << fwdIdx << std::endl;
  std::cerr << "pfwdIdx: " << pfwdIdx << std::endl;
  std::cerr << "Dist fwdIdx->pfwdIdx: " << p2fdist(points[fwdIdx],points[pfwdIdx]) << std::endl;

  std::cerr << "Dist bwdIdx->pbwdIdx: " << p2fdist(points[bwdIdx],points[pbwdIdx]) << std::endl;
}

void getArcSimple(size_t A,
                size_t B,
                std::vector<cv::Point2f> &contour,
                std::vector<cv::Point2f> &Arc,
                bool left
                )
{
  size_t I=A;
  int step;
  //size_t contour_size=contour.size();
  if(left)
    step=-1;
  else
    step=1;

  if(left==true)
    if(I==0)
      I=contour.size()-1;
    else
      I+=step;
  else 
  {
    if(I==contour.size()-1)
      I=0;
    else
      I+=step;
  }
  
  while(I!=B)
  {
    Arc.push_back(contour[I]);
    if(left==true)
      if(I==0)
        I=contour.size()-1;
      else
        I+=step;
    else 
    {
      if(I==contour.size()-1)
        I=0;
      else
        I+=step;
    }
  }  
}

void getLeftArc(size_t A,
                size_t B,
                std::vector<cv::Point2f> &contour,
                std::vector<cv::Point2f> &leftArc,
                std::vector<float> &centripetal,
                std::vector<float> &leftArcCentripetal,
                std::vector<float> &leftArcDistanceMap
                )
{
  size_t I=A;
  
  while(I!=B)
  {
    if(I==A)
    {
      leftArcDistanceMap.push_back(0);
      leftArc.push_back(contour[A]);
      leftArcCentripetal.push_back(centripetal[A]);
      if(I==0)
        I=contour.size()-1;
      else
        I--;
      continue;
    }
    leftArc.push_back(contour[I]);
    leftArcCentripetal.push_back(centripetal[I]);
      
    if(I==0)
    {
      double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[I+1]);
      leftArcDistanceMap.push_back(newDst);
      I=contour.size()-1;
    }
    else
    {
      if(I==contour.size()-1)
      {
        double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[0]);
        leftArcDistanceMap.push_back(newDst);
      }
      else
      {
        double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[I+1]);
        leftArcDistanceMap.push_back(newDst);
      }
      --I;
    }
  }
  if(B!=contour.size()-1)
  {
    double newDst=leftArcDistanceMap.back()+p2fdist(contour[B],contour[B+1]);
    leftArcDistanceMap.push_back(newDst);
  }
  else
  {
    double newDst=leftArcDistanceMap.back()+p2fdist(contour[B],contour[0]);
    leftArcDistanceMap.push_back(newDst);
  }
  leftArc.push_back(contour[B]);
  leftArcCentripetal.push_back(centripetal[B]);
}

void getRightArc(size_t A,
                size_t B,
                std::vector<cv::Point2f> &contour,
                std::vector<cv::Point2f> &rightArc,
                std::vector<float> &centripetal,
                std::vector<float> &rightArcCentripetal,
                std::vector<float> &rightArcDistanceMap
                )
{
  size_t I=A;
  
  while(I!=B)
  {
    if(I==A)
    {
      rightArcDistanceMap.push_back(0);
      rightArc.push_back(contour[A]);
      rightArcCentripetal.push_back(centripetal[A]);
      if(I==contour.size()-1)
        I=0;
      else
        I++;
      continue;
    }

    rightArc.push_back(contour[I]);
    rightArcCentripetal.push_back(centripetal[I]);
    if(I==contour.size()-1)
    {
      I=0;
      rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                    p2fdist(contour[0],contour.back()));
    }
    else
    {
      if(I==0)
      {
        rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                     p2fdist(contour[0],contour.back()));
      }
      else
      {
        rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                      p2fdist(contour[I],contour[I-1]));
      }
      ++I;
    }
  }
  if(B!=0)
    rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                  p2fdist(contour[B],contour[B-1]));
  else
    rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                  p2fdist(contour[B],contour.back()));
  rightArc.push_back(contour[B]);
  rightArcCentripetal.push_back(centripetal[B]);


}

size_t getNextPointInArc(std::vector<cv::Point2f> &arc,
                       std::vector<float> &arcDistanceMap,
                       double distance,
                       cv::Point2f &newPoint
    )
{
  size_t i=0;
  while(distance>arcDistanceMap[i] && i<arcDistanceMap.size())
    ++i;

  if(i>0)
  {
    double d=distance-arcDistanceMap[i-1];
    double diip1=arcDistanceMap[i]-arcDistanceMap[i-1];
    double r=d/diip1;
    newPoint=arc[i-1]+r*(arc[i]-arc[i-1]);
    r=1;
  }
    return i;
}

double distanceBetweenSpines(larvaDistanceMap &a,larvaDistanceMap &b)
{
  double adist=0.0;
  double rdist=0.0;
  for(size_t i=0;i<a.Spine.size();i++)
  {
    adist+=p2fdist(a.Spine[i],b.Spine[i]);
    rdist+=p2fdist(a.Spine[i],b.Spine[a.Spine.size()-1-i]);
  }
  if(rdist>adist)
    return adist;
  else
    return rdist;
}

size_t findContourPointClosestTo(
                        cv::Point2f &p,
                        std::vector<cv::Point2f> &v
                        )
{
  double min=65535;
  size_t retval;
  size_t step=v.size()/10;
  for(size_t i=0;i<v.size();i+=step)
  {
    double cdist=p2fdist(p,v[i]);
    if(cdist<min)
    {
      min=cdist;
      retval=i;
    }
  }

  for(int i=retval-step;i<retval+step;++i)
  {
    size_t j=i;
    if(i<0)
      j=v.size()+i;
    if(i>v.size()-1)
      j=i-v.size();
    double cdist=p2fdist(v[j],p);
    if(cdist<min)
    {
      min=cdist;
      retval=j;
    }
  }
  return retval;
  
}

bool getPerpendicularPoints(cv::Point2f &p1, 
                            cv::Point2f &p2,
                            cv::Point2f &p3,
                            std::vector<cv::Point2f> &leftArc,
                            std::vector<cv::Point2f> &rightArc,
                            double maxDist,
                            cv::Point2f &L,
                            cv::Point2f &R
    )
{

  cv::Point2f lp=perp(p1,p2,p3,maxDist, 
      false);
  cv::Point2f rp=perp(p1,p2,p3,maxDist, 
      true);
  size_t li=1,ri=1;
  cv::Point2f lb,le,rb,re;

  lb=leftArc[li-1]; le=leftArc[li];
  rb=rightArc[ri-1]; re=rightArc[ri];

  double minlDist=DBL_MAX;
  double minrDist=DBL_MAX;
  while(li<leftArc.size())
  {
    double d;
    cv::Point2f intl;
    if(intersection(lb,le,lp,rp,intl))
    {
      d=p2fdist(intl,p2);
      if(d<minlDist)
      {
        minlDist=d;
        L=intl;
      }
    }
    lb=leftArc[li-1]; le=leftArc[li];
    ++li;
  }
  while(ri<rightArc.size())
  {
    double d;
    cv::Point2f intr;
    if(intersection(rb,re,lp,rp,intr))
    {
      d=p2fdist(intr,p2);
      if(d<minrDist)
      {
        minrDist=d;
        R=intr;
      }
    }
    rb=rightArc[ri-1]; re=rightArc[ri];
    ++ri;
  }
 return true;
}

/*bool getIntersectionsWithContour(
  size_t idx,
  cv::Point2f &p1,
  cv::Point2f &p2,
  cv::Point2f &p3,
  larvaDistanceMap &candidateMap,
  std::vector<cv::Point2f> &leftArc,
  std::vector<cv::Point2f> &rightArc,
  double maxDist,
  bool noNearest)
{
  cv::Point2f int_left;
  cv::Point2f int_right;
  std::vector<PointPair> &P=candidateMap.fSpinePairs;
  std::vector<double> &W=candidateMap.Widths;
  size_t cli;
  size_t cri;
  bool retval=true;
  maxDist*=2;

  if(idx>1 && idx<SPINE_SEGMENTS-2)
  {
    //Main body we choose the closest point
    cv::Point2f &cl=leftArc[findContourPointClosestTo(p2,leftArc)];
    cv::Point2f &cr=rightArc[findContourPointClosestTo(p2,rightArc)];
    //if(p2fdist(cl,p2)<p2fdist(cr,p2))
    if(rightArc.size()>leftArc.size())
    {
      //Left arc has the nearest point
      cv::Point2f xp=cl+3*(p2-cl);
      size_t ri=1;
      cv::Point2f rb=rightArc[ri-1];
      cv::Point2f re=rightArc[ri];
      double minrDist=DBL_MAX;
      int_left=cl;
      while(ri<rightArc.size())
      {
        double d;
        cv::Point2f intr;
        //Intersection between: 
        //previous rightArc point, current rightArc point
        //closest point on leftArc and extendedPoint
        if(intersection(rb,re,xp,cl,intr))
        {
          d=p2fdist(intr,p2);
          if(d<minrDist)
          {
            minrDist=d;
            int_right=intr;
          }
        }
        rb=rightArc[ri-1]; re=rightArc[ri];
        ++ri;
      }
    }
    else
    {
      //Right arc has the nearest point
      cv::Point2f xp=cr+3*(p2-cr);
      size_t li=1;
      cv::Point2f lb=leftArc[li-1];
      cv::Point2f le=leftArc[li];
      double minlDist=DBL_MAX;
      int_right=cr;
      while(li<leftArc.size())
      {
        double d;
        cv::Point2f intl;
        //Intersection between: 
        //previous leftArc point, current leftArc point
        //closest point on rightArc and extendedPoint
        if(intersection(lb,le,xp,cr,intl))
        {
          d=p2fdist(intl,p2);
          if(d<minlDist)
          {
            minlDist=d;
            int_left=intl;
          }
        }
        lb=leftArc[li-1]; le=leftArc[li];
        ++li;
      }
    }
  }
  else
  {
    cv::Point2f lp=perp(p1,p2,p3,maxDist, 
        false);
    cv::Point2f rp=perp(p1,p2,p3,maxDist, 
        true);
    size_t li=1,ri=1;
    cv::Point2f lb,le,rb,re;

    lb=leftArc[li-1]; le=leftArc[li];
    rb=rightArc[ri-1]; re=rightArc[ri];

    double minlDist=DBL_MAX;
    double minrDist=DBL_MAX;
    while(li<leftArc.size())
    {
      double d;
      cv::Point2f intl;
      if(intersection(lb,le,lp,rp,intl))
      {
        d=p2fdist(intl,p2);
        if(d<minlDist)
        {
          minlDist=d;
          int_left=intl;
        }
      }
      lb=leftArc[li-1]; le=leftArc[li];
      ++li;
    }
    while(ri<rightArc.size())
    {
      double d;
      cv::Point2f intr;
      if(intersection(rb,re,lp,rp,intr))
      {
        d=p2fdist(intr,p2);
        if(d<minrDist)
        {
          minrDist=d;
          int_right=intr;
        }
      }
      rb=rightArc[ri-1]; re=rightArc[ri];
      ++ri;
    }
  }
  P[idx-1].first=int_left;
  P[idx-1].second=int_right;
  W[idx-1]=p2fdist(int_left,int_right);
  if(W[idx-1]>candidateMap.WidthDist)
    candidateMap.WidthDist=W[idx-1];
  //p2=0.5*(int_left+int_right);

  return true;
}*/

/*void reSpine(
    std::vector<cv::Point2f> &leftArc,
    std::vector<cv::Point2f> &rightArc,
    std::vector<cv::Point2f> &contour,
    larvaDistanceMap &candidateMap,
    double maxDist,
    bool noNearest)
{
  bool haveFailedIntersections=false;
  std::vector<float> spineDistanceMap;
  cv::Point2f prePoint=1*candidateMap.Spine[0];
  float psum=0.0;

  for(size_t i=1;i<SPINE_SEGMENTS-1;++i)
  {
    if(i>=2)
    {
      cv::Point2f &p0=candidateMap.Spine[i-2];
      cv::Point2f &p1=candidateMap.Spine[i-1];
      cv::Point2f &p2=candidateMap.Spine[i];
      if(getIntersectionsWithContour(i-1,p0,p1,p2,candidateMap,
            leftArc,rightArc,maxDist,noNearest)==false)
        haveFailedIntersections=false;
      candidateMap.Angles[i-1]=
        angleD(p0,p1,p2);
    [>std::cerr << "Inserting to Angles of size: " 
      << candidateMap.Angles.size() 
      << " in position " << i-1 << std::endl;<]
      if(candidateMap.maxAngle<candidateMap.Angles.back())
      {
        candidateMap.maxAngle=candidateMap.Angles.back();
        candidateMap.maxAngleLocation=i-1;
      }
    }
  }
  cv::Point2f &p0=candidateMap.Spine[SPINE_SEGMENTS-3];
  cv::Point2f &p1=candidateMap.Spine[SPINE_SEGMENTS-2];
  cv::Point2f &p2=candidateMap.Spine[SPINE_SEGMENTS-1];
  if(getIntersectionsWithContour(SPINE_SEGMENTS-2,p0,p1,p2,candidateMap,
        leftArc,rightArc,maxDist,noNearest)==false)
    haveFailedIntersections=false;
  candidateMap.WidthAvg=candidateMap.WidthAvg/candidateMap.Widths.size();

  spineDistanceMap.clear();
  for(auto &sp: candidateMap.Spine)
  {
    spineDistanceMap.push_back(
        psum+p2fdist(sp,prePoint));
    prePoint=1*sp;
    psum=spineDistanceMap.back();
  }
  double dist=spineDistanceMap.back();

  double step=dist/(SPINE_SEGMENTS-1);
  size_t idx;
  for(size_t i=1;i<SPINE_SEGMENTS-1;++i)
  {
    cv::Point2f np;
    getNextPointInArc(candidateMap.Spine,
                      spineDistanceMap,
                      i*step,
                      np);
    [>
    std::cerr << "Inserting to Spine of size: " 
      << candidateMap.Spine.size() 
      << " in position " << i << std::endl;
      <]
    candidateMap.Spine[i]=1*np;
  }

  size_t sidx=1; 
  [>if(haveFailedIntersections)
  {
    for(auto &p: candidateMap.fSpinePairs)
    {
      candidateMap.Spine[sidx]=0.5*(p.first+p.second);
      sidx++;
    }
  }<]

  size_t Widx=0;
  for(auto &w: candidateMap.Widths)
  {
    if(Widx<2 && Widx>SPINE_SEGMENTS-2 )
      Widx++;
    else
    {
      double D=(w-candidateMap.WidthAvg);
      candidateMap.WidthSD+=D*D;
      Widx++;
    }

  }
  candidateMap.WidthSD=sqrt(candidateMap.WidthSD/candidateMap.Widths.size());

  candidateMap.MaxDist=cv::arcLength(candidateMap.Spine,false);
}*/

/*//We assume leftArc and rightArc start with A
//and end with B included.
void spineAndContour( 
    size_t A,
    size_t B,
    std::vector<cv::Point2f> &leftArc,
    std::vector<float> &leftArcCurvatures,
    std::vector<cv::Point2f> &rightArc,
    std::vector<float> &rightArcCurvatures,
    std::vector<cv::Point2f> &contour,
    std::vector<float> &leftArcDistanceMap,
    std::vector<float> &rightArcDistanceMap,
    larvaDistanceMap &candidateMap
  )
{
  double stepLeft=leftArcDistanceMap.back()/(SPINE_SEGMENTS-1);
  double stepRight=rightArcDistanceMap.back()/(SPINE_SEGMENTS-1);
  std::vector<cv::Point2f> wholeSpine;
  wholeSpine.resize(SPINE_SEGMENTS);
  candidateMap.Spine.resize(SPINE_SEGMENTS);
  candidateMap.spinePairs.resize(SPINE_SEGMENTS-2);
  candidateMap.Widths.resize(SPINE_SEGMENTS-2);
  candidateMap.fSpinePairs.resize(SPINE_SEGMENTS-2);
  candidateMap.Angles.resize(SPINE_SEGMENTS-2);
  size_t lpap=0; // previous_arc_point left
  size_t rpap=0; // previous_arc_point left
  size_t lpapr=0; // previous_arc_point left
  size_t rpapr=0; // previous_arc_point left
  candidateMap.WidthDist=0;
  double maxDist=0;

  [>double lMinVal=0.0;
    double rMinVal=0.0;
    cv::minMaxLoc(leftArcCurvatures,&lMinVal);
    cv::minMaxLoc(rightArcCurvatures,&rMinVal);
    double cMin=std::min(lMinVal,rMinVal);
    if(cMin<0)
    {
    cv::add(leftArcCurvatures,(-1.1)*cMin,leftArcCurvatures);
    cv::add(rightArcCurvatures,(-1.1)*cMin,rightArcCurvatures);
    }
    cv::resize(leftArcCurvatures,
    leftArcCurvatures,
    cv::Size(SPINE_SEGMENTS-2,1),
    0,
    0,
    cv::INTER_LANCZOS4);
    cv::resize(rightArcCurvatures,
    rightArcCurvatures,
    cv::Size(SPINE_SEGMENTS-2,1),
    0,
    0,
    cv::INTER_LANCZOS4);

    std::vector<double> cRatio;
    cv::divide(leftArcCurvatures,rightArcCurvatures,cRatio);
    double cRm=cv::mean(cRatio)[0];
    for(auto &v:cRatio)
    v=v/cRm;<]

  wholeSpine[0]=1*contour[A];
  wholeSpine[SPINE_SEGMENTS-1]=1*contour[B];
  for(size_t i=1;i<SPINE_SEGMENTS-1;++i)
  {
    cv::Point2f cPointLeft;
    cv::Point2f cPointRight;
    lpap=getNextPointInArc(leftArc,
        leftArcDistanceMap,
        i*stepLeft,/[>cRatio[i-1],
        cPointLeft);

    rpap=getNextPointInArc(rightArc,
        rightArcDistanceMap,
        i*stepRight,
        cPointRight);

    cv::Point2f newPoint=(cPointLeft+cPointRight)*0.5;
    candidateMap.spinePairs[i-1]=std::make_pair(cPointLeft,cPointRight);
    double dmax=p2fdist(cPointRight,cPointLeft);
    candidateMap.Widths[i-1]=dmax;
    if(dmax>maxDist)
      maxDist=dmax;
    wholeSpine[i]=1*newPoint;
  }
  candidateMap.Spine=wholeSpine;

  reSpine(leftArc,rightArc,contour,candidateMap,maxDist,true);
  //for(size_t i=0;i<4;i++)
  //{
  //  reSpine(leftArc,rightArc,contour,candidateMap,maxDist,false);
  //}
}*/

/*void getCandidateSpine2(cvb::CvBlob &blob,
                       size_t A,
                       size_t B, 
                       std::vector<cv::Point2f> &contour,
                       larvaDistanceMap &candidateMap,
                       std::map<float,size_t> &curvatures,
                       size_t RESOLUTION=SPINE_SEGMENTS,
                       std::vector<cv::Point2f> *heads=NULL,
                       std::vector<cv::Point2f> *tails=NULL,
                       std::vector<cvb::CvBlob> *blobs=NULL)
{
  std::vector<cv::Point2f> leftArc;
  std::vector<float> leftArcDistanceMap;
  std::vector<float> leftArcCurvatures;
  std::vector<cv::Point2f> rightArc;
  std::vector<float> rightArcDistanceMap;
  std::vector<float> rightArcCurvatures;
  std::vector<cv::Point2f> wholeSpine;

  getArcSimple(A,B,contour,leftArc,true);
  getArcSimple(A,B,contour,rightArc,false);

  candidateMap.Spine.clear();
  wholeSpine.clear();
  wholeSpine.push_back(contour[A]);
  candidateMap.MaxDist=0;
  candidateMap.WidthDist=0;
  candidateMap.maxAngle=0;
  candidateMap.WidthAvg=0;
  candidateMap.firstHalfWidthsSum=0;
  candidateMap.secondHalfWidthsSum=0;

  std::map<float,size_t>::reverse_iterator curvIT=curvatures.rbegin();
  size_t CURVTEST=1;
  int TOTAL=curvatures.size();
  candidateMap.curvatureBias=0.0;
  for(size_t i=0; i<CURVTEST;++i)
  {
    size_t dA=std::min(abs(curvIT->second-A),TOTAL-abs(curvIT->second-A));
    size_t dB=std::min(abs(curvIT->second-B),TOTAL-abs(curvIT->second-B));
    if(dA<dB)
    {
      candidateMap.curvatureBias+=1.0/CURVTEST;
    }
    curvIT++;
  }

  std::map<size_t, std::pair<size_t,float> > minDist;
  std::vector<cv::Point2f> *smallest_arc;
  std::vector<cv::Point2f> *longest_arc;
  if(leftArc.size()<rightArc.size())
  {
    smallest_arc=&leftArc;
    longest_arc=&rightArc;
  }
  else
  {
    smallest_arc=&rightArc;
    longest_arc=&leftArc;
  }

  size_t preli=0;
  size_t li=0;
  for(size_t si=0; si<smallest_arc->size();si++)
  {
    if((*smallest_arc)[si]==(*longest_arc)[li])
    {
      std::cerr << "lv==rv should not happen" << std::endl;
    }
    double dist=p2fdist((*smallest_arc)[si],(*longest_arc)[li]);
    double distAl=p2fdist((*longest_arc)[li],contour[A]);
    double distAs=p2fdist((*smallest_arc)[li],contour[A]);
    double distD=distAl-distAs;
    double newdist=dist-distD;
    bool entry=true;
    while(newdist<=dist || entry)
    {
      entry=false;
      li++;
      dist=newdist;
      distAl=p2fdist((*longest_arc)[li],contour[A]);
      distAs=p2fdist((*smallest_arc)[li],contour[A]);
      distD=distAl-distAs;
      newdist=p2fdist((*smallest_arc)[si],
          (*longest_arc)[li])-distD;
    }
    li--;
    if(li==preli)
      continue;
    preli=li;

    if((*longest_arc)[li]==contour[B])
      break;

    cv::Point2f &cPointLeft=(*smallest_arc)[si];
    cv::Point2f &cPointRight=(*longest_arc)[li];
    cv::Point2f newPoint=(cPointLeft+cPointRight)*0.5;
    candidateMap.spinePairs.push_back
      (std::make_pair(cPointLeft,cPointRight));
    candidateMap.Widths.push_back(p2fdist(cPointLeft,cPointRight));
    candidateMap.WidthAvg+=candidateMap.Widths.back();
    //cv::line(FOO,cPointLeft*2,cPointRight*2,cv::Scalar(0,0,200));
    if(candidateMap.Widths.back()>candidateMap.WidthDist)
    {
      candidateMap.WidthDist=candidateMap.Widths.back();
    }
    if(li>0.85*leftArc.size())
      candidateMap.firstHalfWidthsSum+=candidateMap.Widths.back();

    if(li<0.15*leftArc.size())
      candidateMap.secondHalfWidthsSum+=candidateMap.Widths.back();

    wholeSpine.push_back(newPoint);
  }

  candidateMap.WidthAvg=candidateMap.WidthAvg/candidateMap.Widths.size();
  for(size_t i=0;i<candidateMap.Widths.size();++i)
  {
    double D=(candidateMap.Widths[i]-candidateMap.WidthAvg);
    candidateMap.WidthSD+=D*D;
  }
  candidateMap.WidthSD=sqrt(candidateMap.WidthSD/candidateMap.Widths.size());

  size_t i=RESOLUTION-1;
  wholeSpine.push_back(contour[B]);
  if(wholeSpine.size()==0)
  {
    candidateMap.MaxDist=-1;
    return;
  }
  candidateMap.MaxDist=cv::arcLength(wholeSpine,false);
  std::vector<float> spineDistanceMap;

  spineDistanceMap.push_back(0);
  spineDistanceMap.push_back(p2fdist(wholeSpine[1],wholeSpine[0]));
  for (size_t i=2;i<wholeSpine.size()-1;++i)
    spineDistanceMap.push_back(spineDistanceMap.back()+
        p2fdist(wholeSpine[i],wholeSpine[i-1]));

  spineDistanceMap.push_back(candidateMap.MaxDist);


  candidateMap.Spine.push_back(wholeSpine[0]);
  spineAndContour(leftArc, rightArc,
      wholeSpine,
      candidateMap,
      spineDistanceMap);
  candidateMap.Angles.push_back(angleD(
        candidateMap.Spine[SPINE_SEGMENTS-3],
        candidateMap.Spine[SPINE_SEGMENTS-2],
        candidateMap.Spine.back()));
  if(candidateMap.maxAngle<candidateMap.Angles.back())
  {
    candidateMap.maxAngle=candidateMap.Angles.back();
    candidateMap.maxAngleLocation=i-1;
  }
}*/

void presize(std::vector<cv::Point2f> &a,
             std::vector<cv::Point2f> &o,
             size_t n)
{

  alglib::real_2d_array xy;
  xy.setlength(a.size(),2);
  alglib::pspline2interpolant p;
  std::vector<cv::Point2f> NP;
  for(int i=0;i<a.size();i++)
  {
    xy[i][0]=a[i].x;
    xy[i][1]=a[i].y;
  }

  alglib::pspline2build(xy,a.size(),1,2,p);
  
  double t0=0;
  double t;
  double step=1.0/(n-1);
  if(o==a)
    o.clear();
  for (int i=0;i<n;i++)
  {
    t=t0+i*step;
    double x,y;
    alglib::pspline2calc(p,t,x,y);
    cv::Point2f N(x,y);
    o.push_back(N);
  }
}

void getCandidateSpine(cvb::CvBlob &blob,
                       size_t A,
                       size_t B, 
                       std::vector<cv::Point2f> &contour,
                       larvaDistanceMap &candidateMap,
                       std::map<float,size_t> &curvatures,
                       std::vector<float> &curvatureMat,
                       size_t RESOLUTION=SPINE_SEGMENTS,
                       std::vector<cv::Point2f> *heads=NULL,
                       std::vector<cv::Point2f> *tails=NULL,
                       std::vector<cvb::CvBlob> *blobs=NULL)
{
  std::vector<cv::Point2f> leftArc;
  std::vector<cv::Point2f> leftArcRev;
  std::vector<float> leftArcDistanceMap;
  std::vector<float> leftArcCurvatures;
  std::vector<cv::Point2f> rightArc;
  std::vector<float> rightArcDistanceMap;
  std::vector<float> rightArcCurvatures;
  getLeftArc(A,B,contour,
             leftArc,curvatureMat,
             leftArcCurvatures,
             leftArcDistanceMap);
  getRightArc(A,B,contour,
              rightArc,curvatureMat,
              rightArcCurvatures,
              rightArcDistanceMap);

  size_t samples=std::max(leftArc.size(),rightArc.size());
  //if(samples!=leftArc.size())
  presize(leftArc,leftArc,samples);
  //else
  presize(rightArc,rightArc,samples);
  std::vector<cv::Point2f> cSpine;
  cv::addWeighted(leftArc,0.5,rightArc,0.5,0,cSpine);
  double spineLength=cv::arcLength(cSpine,false);
  double distanceStep=spineLength/(SPINE_SEGMENTS-1);
  double cdist=0.0;
  cv::Point2f pp=cSpine[0];
  cv::Point2f mp=cSpine[1];
  //size_t idx=0;
  double wMax=0.0;
  //double maxDist=100;
  candidateMap.Spine.push_back(cSpine[0]);

  //for(int idx=2;idx<cSpine.size()-1;idx++)
  for(int idx=1;idx<cSpine.size()-1;idx++)
  {
    cv::Point2f &p=cSpine[idx];
    //cv::Point2f &np=cSpine[idx+1];
    double ndist=cdist+p2fdist(pp,p);
    cdist+=p2fdist(pp,p);
    if(ndist>distanceStep)
    {
      cdist=ndist-distanceStep;

      candidateMap.Spine.push_back(pp);
      //cv::Point2f lp,rp;
      /*getPerpendicularPoints(pp, 
                            p,
                            np,
                            leftArc,
                            rightArc,
                            maxDist,
                            lp,
                            rp);*/
      candidateMap.spinePairs.push_back(std::make_pair(leftArc[idx],
                                                       rightArc[idx]));
      /*candidateMap.spinePairs.push_back(std::make_pair(lp,
                                                       rp));*/
      double cw=p2fdist(leftArc[idx],rightArc[idx]);
      //double cw=p2fdist(lp,rp);
      candidateMap.Widths.push_back(cw);
      if(cw>wMax)
        wMax=cw;
    }
    pp=p;
  }
  if(candidateMap.Spine.size()<SPINE_SEGMENTS)
    candidateMap.Spine.push_back(cSpine.back());

  candidateMap.Spine.back()=cSpine.back();
    
  std::vector<double> wmean;
  std::vector<double> std;
  cv::meanStdDev(candidateMap.Widths,wmean,std);
  candidateMap.WidthSD=std[0];
  candidateMap.WidthDist=wMax;
  candidateMap.MaxDist=spineLength;
}

/*void getCandidateSpine2(cvb::CvBlob &blob,
                       size_t A,
                       size_t B, 
                       std::vector<cv::Point2f> &contour,
                       larvaDistanceMap &candidateMap,
                       std::map<float,size_t> &curvatures,
                       std::vector<float> &curvatureMat,
                       size_t RESOLUTION=SPINE_SEGMENTS,
                       std::vector<cv::Point2f> *heads=NULL,
                       std::vector<cv::Point2f> *tails=NULL,
                       std::vector<cvb::CvBlob> *blobs=NULL)
{
  std::vector<cv::Point2f> leftArc;
  std::vector<cv::Point2f> leftArcRev;
  std::vector<float> leftArcDistanceMap;
  std::vector<float> leftArcCurvatures;
  std::vector<cv::Point2f> rightArc;
  std::vector<float> rightArcDistanceMap;
  std::vector<float> rightArcCurvatures;
  std::vector<cv::Point2f> wholeSpine;

  getLeftArc(A,B,contour,
             leftArc,curvatureMat,
             leftArcCurvatures,
             leftArcDistanceMap);
  getRightArc(A,B,contour,
              rightArc,curvatureMat,
              rightArcCurvatures,
              rightArcDistanceMap);

  float dleft=leftArcDistanceMap.back();
  float dright=rightArcDistanceMap.back();

  if(dleft > 4*dright || dright > 4*dleft) 
  {
    candidateMap.MaxDist=0;
    return;
  }
  
  double stepLeft=dleft/(SPINE_SEGMENTS-1);
  double stepRight=dright/(SPINE_SEGMENTS-1);
  candidateMap.Spine.clear();
  wholeSpine.clear();
  candidateMap.MaxDist=0;
  candidateMap.WidthDist=0;
  candidateMap.maxAngle=0;
  candidateMap.WidthAvg=0;
  candidateMap.firstHalfWidthsSum=0;
  candidateMap.secondHalfWidthsSum=0;

  spineAndContour( 
      A,B,
      leftArc,
      leftArcCurvatures,
      rightArc,
      rightArcCurvatures,
      contour,
      leftArcDistanceMap,
      rightArcDistanceMap,
      candidateMap
  );

  if(candidateMap.Spine.size()==0)
  {
    candidateMap.MaxDist=-1;
    return;
  }

  std::map<float,size_t>::reverse_iterator curvIT=curvatures.rbegin();
  size_t CURVTEST=1;
  int TOTAL=curvatures.size();
  candidateMap.curvatureBias=0.0;

}
*/
/*bool checkFinish(
                size_t &pointsCrossed,
                std::vector<cv::Point2f> &points,
                std::vector<cv::Point2f> &curSpine,
                size_t i,
                size_t &fwdIdx,
                size_t &bwdIdx,
                size_t &pfwdIdx,
                size_t &pbwdIdx,
                cv::Point2f &bwdPoint,
                cv::Point2f &fwdPoint,
                double &curDist
                )
{
      if(pointsCrossed==points.size()-1)
      {
        fwdIdx++;
        if(fwdIdx==points.size())
          fwdIdx=0;
        curSpine.push_back(points[fwdIdx]);
        curDist+=p2fdist(points[fwdIdx],curSpine.back());
        return true;
      }
      if(pointsCrossed==points.size()-2)
      {
        fwdIdx++;
        if(fwdIdx==points.size())
          fwdIdx=0;

        if(bwdIdx==0)
          bwdIdx=points.size();
        bwdIdx--;

        double fangle=angleD(curSpine[curSpine.size()-2],
               curSpine[curSpine.size()-1],
               points[fwdIdx]);

        double bangle=angleD(curSpine[curSpine.size()-2],
               curSpine[curSpine.size()-1],
               points[bwdIdx]);

        double MINANGLE=140;
        if (abs(fangle)<MINANGLE&& abs(bangle)<MINANGLE)
        {
          cv::Point2f np=(points[bwdIdx]+points[fwdIdx])*0.5;
          curDist+=p2fdist(np,curSpine.back());
          curSpine.push_back(np);
          return true;
        }
        else if(abs(fangle)<MINANGLE)
        {
          double distb=p2fdist(curSpine.back(),points[bwdIdx]);
          curSpine.push_back(points[bwdIdx]);
          curDist+=distb;
          return true;
        }
        else if(abs(bangle)<MINANGLE)
        {
          double distf=p2fdist(curSpine.back(),points[fwdIdx]);
          curSpine.push_back(points[fwdIdx]);
          curDist+=distf;
          return true;
        }
        else
        {
          double distf=p2fdist(curSpine.back(),points[fwdIdx]);
          double distb=p2fdist(curSpine.back(),points[bwdIdx]);

          if(distb>distf)
          {
            curSpine.push_back(points[bwdIdx]);
            curDist+=distb;
          }
          else
          {
            curSpine.push_back(points[fwdIdx]);
            curDist+=distf;
          }
          return true;
        }
      }
      return false;
}*/

void getPreviousPoint(std::vector<cv::Point2f> &a,
                  size_t idx,
                  cv::Point2f &n)
{
  if(idx==0)
    n=a[a.size()-1];
  else
    n=a[idx-1];
}

void expandcontour(cv::Mat &frame,
                   std::vector<cv::Point2f> &points,
                   std::vector<cv::Point2f> &outpoints)
{
  for(size_t i=0;i<points.size();++i)
  {
    //check points lrud
    cv::Point2f up(points[i].x,points[i].y-1);
    cv::Point2f down(points[i].x,points[i].y+1);
    cv::Point2f left(points[i].x-1,points[i].y);
    cv::Point2f right(points[i].x+1,points[i].y);

    cv::Vec3b uval= frame.at<cv::Vec3b>(cv::Point2f(up));
    cv::Vec3b dval= frame.at<cv::Vec3b>(cv::Point2f(down));
    cv::Vec3b lval= frame.at<cv::Vec3b>(cv::Point2f(left));
    cv::Vec3b rval= frame.at<cv::Vec3b>(cv::Point2f(right));
    
    if(uval[0]<dval[0] && uval[0] < lval[0] && uval[0] <rval[0])
    {
      outpoints.push_back(up);
    }
    else if(dval[0]<uval[0] && dval[0] < lval[0] && dval[0] < rval[0])
    {
      outpoints.push_back(down);
    }
    else if(lval[0]<uval[0] && lval[0] < dval[0] && lval[0] <rval[0])
    {
      outpoints.push_back(left);
    }
    else if(rval[0]<uval[0] && rval[0] < uval[0] && rval[0] <lval[0])
    {
      outpoints.push_back(right);
    }
    else
    {
      outpoints.push_back(points[i]);
    }
  }
}

void extractCentripetal(
                       std::vector<cv::Point2f> &ppoints,
                       std::vector<cv::Point2f> &xpoints,
                       std::vector<float> &d,
                       float &csqrt)
{
  size_t step=1;
  xpoints.push_back(ppoints[0]);
  for(size_t i=step;i<ppoints.size();i+=step)
  {
    cv::Point2f NP=ppoints[i];
    if(NP!=xpoints.back() && NP!=xpoints[0])
    {
      csqrt+=sqrt(p2fdist(NP,xpoints.back()));
      xpoints.push_back(NP);
    }
    else
    {
      if( i<ppoints.size()-1 && 
          ppoints[i+1]!=xpoints.back() && 
          ppoints[i+1] != xpoints[0] )
      {
        NP=ppoints[i+1];
        csqrt+=sqrt(p2fdist(NP,xpoints.back()));
        xpoints.push_back(NP);
      }
    }
  }
  csqrt+=sqrt(p2fdist(xpoints.back(),xpoints[0]));
  d.push_back(0);
  for(size_t i=1;i<xpoints.size();++i)
  {
    float newD = d.back() + sqrt(p2fdist(xpoints[i-1],xpoints[i]))/csqrt;
    d.push_back(newD);
  }
  d.push_back(d.back()+sqrt(p2fdist(xpoints.back(),xpoints[0]))/csqrt);
}

/*void extractCurvature(
                       std::vector<cv::Point2f> &ppoints,
                       std::vector<float> &curvature)
{
}*/

void smoothAndExtractCentripetal(cv::Mat &frame,
                                 std::vector<cv::Point2f> &ppoints,
                                 std::vector<cv::Point2f> &xpoints,
                                 std::vector<float> &d,
                                 float &csqrt)
{
  cv::Point2f NP=px5Smooth(frame, 
                           ppoints[ppoints.size()-2],
                           ppoints.back(), 
                           ppoints[0],
                           ppoints[1],
                           ppoints[2]); 
  xpoints.push_back(NP);
  for(size_t i=1;i<ppoints.size();++i)
  {
    if(ppoints[i-1]!=ppoints[i])
    {
      cv::Point2f NP;
      if(i==1)
        NP=px5Smooth(frame,ppoints.back(),ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if(i<ppoints.size()-2)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if (i==ppoints.size()-1)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[0],ppoints[1]);
      else if (i==ppoints.size()-2)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[0]);

      if(NP!=xpoints.back())
      {
        if(NP==xpoints[0])
          break;
        csqrt+=sqrt(p2fdist(NP,xpoints.back()));
        xpoints.push_back(NP);
      }
    }
  }
  csqrt+=sqrt(p2fdist(xpoints.back(),xpoints[0]));
  d.push_back(0);
  for(size_t i=1;i<xpoints.size();++i)
  {
    float newD = d.back() + sqrt(p2fdist(xpoints[i-1],xpoints[i]))/csqrt;
    d.push_back(newD);
  }
  d.push_back(d.back()+sqrt(p2fdist(xpoints.back(),xpoints[0]))/csqrt);
}

void getNextPoint(std::vector<cv::Point2f> &a,
                  size_t idx,
                  cv::Point2f &n)
{
  if(idx==a.size()-1)
    n=a[0];
  else
    n=a[idx+1];
}

void fixContour(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    size_t RES,
    cv::Mat &frame,
    cv::Mat &previousFrame,
    std::vector<cv::Point2f> *heads,
    std::vector<cv::Point2f> *tails,
    std::vector<cvb::CvBlob> *blobs)
{
  //int ID=blob.label;
  std::vector<cv::Point2f> cntPoints;
  std::vector<cv::Point2f> tmpPoints;
  std::vector<cv::Point2f> ppoints;
  std::vector<cv::Point2f> newPoints;
  std::vector<float> x2;
  std::vector<float> y2;
  std::vector<float> d;
  double cvariance;
  cv::Point2f bp(blob.minx,blob.miny);
  int PAD=0;
  //blobToPointVector(blob,ppoints);
  blobToContourVector(blob,frame,2,ppoints);
  bool RECONSTRUCT_SPINE=false;
  std::string rec="NOREC";
  if(ppoints.size()<0.1*RES)
  //if(ppoints.size()<0.8*RES)
  {
    ppoints.clear();
    RECONSTRUCT_SPINE=true;
    blobToPointVector(blob,ppoints);
    rec="REC";
  }

  std::vector<cv::Point2f> xpoints;
  //fixContourSimple(blob,frame);
  //cv::approxPolyDP(tmpPoints,ppoints,0.3,true);
  cv::Point2f nh(0,0),nt(0,0);
  if(heads==(std::vector<cv::Point2f> *) -14)
  {
    cvb::CvBlob &pblob=blobs->back();
    cv::Point2f ph=heads->back();
    cv::Point2f pt=tails->back();
    cv::Mat cROI=cv::Mat(frame,
        cv::Rect(pblob.minx-PAD,
          pblob.miny-PAD,
          pblob.maxx-pblob.minx+1+2*PAD,
          pblob.maxy-pblob.miny+1+2*PAD
          )
        );
  cv::Mat ROI;
  cvtColor(cROI,ROI,CV_BGR2GRAY);
    cv::Mat pROI=cv::Mat(previousFrame,
        cv::Rect(pblob.minx-PAD,
          pblob.miny-PAD,
          pblob.maxx-pblob.minx+1+2*PAD,
          pblob.maxy-pblob.miny+1+2*PAD
          )
        );
    cv::Mat d1,d2;
    //double diff1,diff2;
    //double total1, total2;
    d1=ROI-pROI;
    d2=pROI-ROI;
    cv::Mat pcROI;
    cvtColor(pROI,pcROI,CV_GRAY2BGR);
    cv::circle(pcROI, 
      (ph)+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);
    cv::circle(pcROI, 
      (pt)+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);
  }
  float csqrt=0.0;
  //smoothPoints(cntPoints,ppoints,frame,1);
  //smoothAndExtractCentripetal(frame,ppoints,xpoints,d,csqrt);
  extractCentripetal(ppoints,xpoints,d,csqrt);
  //std::cerr << "ppoints: " << ppoints.size() << std::endl;
  //std::cerr << "xpoints: " << xpoints.size() << std::endl;
  /*csqrt+=sqrt(p2fdist(xpoints[1],xpoints[0]));
  csqrt+=sqrt(p2fdist(xpoints[2],xpoints[1]));*/

  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[1],xpoints[0]))/csqrt);*/
  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[2],xpoints[1]))/csqrt);*/
  //xpoints.push_back(xpoints[0]);
  //xpoints.push_back(xpoints[1]);
  /*xpoints.push_back(xpoints[2]);*/
  std::vector<float> w;
  //float AvgGrVal=0.0;
  std::vector<float> dx;
  std::vector<float> dy;
  std::vector<float> d2x;
  std::vector<float> d2y;
  std::vector<float> cvt;
  d2x.push_back(0);
  d2y.push_back(0);
  cvt.push_back(0);
  std::map<float,size_t> curvatures;
  std::vector<size_t> dIdx(3,0);

  cv::Mat ROI=cv::Mat(frame,
              cv::Rect(blob.minx-PAD,
                       blob.miny-PAD,
                       blob.maxx-blob.minx+1+2*PAD,
                       blob.maxy-blob.miny+1+2*PAD
                      )
             );
  cv::Mat cROI,final;
  ROI.copyTo(cROI);
#ifdef LRV_TRACK_VISUAL_DEBUG
  normalize(cROI,cROI,0,255,CV_MINMAX);
#endif
  float wmin=FLT_MAX;
  float wmax=FLT_MIN;

  for(size_t i=0;i<xpoints.size();i++)
  {
    cv::Vec3b val=cROI.at<cv::Vec3b>((xpoints[i]-bp)+cv::Point2f(PAD,PAD));
    w.push_back((255-(val[0]+val[1]+val[2]))/255);
    if(w.back()<wmin)
      wmin=w.back();
    if(w.back()>wmax)
      wmax=w.back();
  }

  for(size_t i=0;i<xpoints.size();i++)
  {
    w[i]=(w[i]-wmin)/(wmax-wmin);
  }

#ifdef LRV_TRACK_VISUAL_DEBUG
  /*for(size_t i=0;i<xpoints.size();++i)
  {
    cv::Vec3b val= cROI.at<cv::Vec3b>((xpoints[i]-bp)+cv::Point2f(PAD,PAD));
    cv::Scalar sval(0,0,std::min(val[0]+val[1]+val[2]+50,255));
    cv::circle(cROI, 
      (xpoints[i]-bp)+cv::Point2f(PAD,PAD),
      0,
      sval,-1);
  }*/
  static int MULT=8;
  cv::resize(cROI,cROI,cv::Size(),MULT,MULT,cv::INTER_NEAREST);
#endif
  //cv::imshow("LRV1",cROI);

//if(heads==NULL)
  dIdx.resize(2);

std::vector<float> vcurv;
//if(xpoints.size()<1000*RES)
if(RECONSTRUCT_SPINE)
{
    //std::cout << "Y" << std::endl;
  try{
    spline4(xpoints,
        d,
        w,
        xpoints.size(),
        RES,
        newPoints,
        dIdx,
        curvatures,
        vcurv,
        cvariance);
  }
  catch(alglib::ap_error &e)
  {
    std::cerr << "Spline Error: [" << rec << "] "  << e.msg << std::endl;
    std::cerr << printVector(xpoints) << std::endl;
    std::cerr << printVector(d) << std::endl;
    return;
  }
  catch(...)
  {
    std::cerr << "Spline Error [" << blob.label << ", " << xpoints.size() << ", " << rec << "]" << std::endl;
    return;
  }
}
else
{
    //std::cout << "N" << std::endl;
    sg_spline(xpoints,
    //nospline(xpoints,
        d,
        w,
        xpoints.size(),
        RES,
        newPoints,
        dIdx,
        curvatures,
        vcurv,
        cvariance);
}


  bp=cv::Point2f(blob.minx-0.5,blob.miny-0.5);
#ifdef LRV_TRACK_VISUAL_DEBUG
  for(size_t i=0;i<newPoints.size();++i)
  {
    cv::circle(cROI, 
        MULT*(newPoints[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        1,
        cv::Scalar(0,255,0),-1);
  }
  cv::circle(cROI, 
      MULT*(newPoints[0]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      2,
      cv::Scalar(0,0,255),-1);
#endif
  double maxDist;
  larvaDistanceMap bestDist(newPoints);

  int k=0;
  k++;
  std::vector<cv::Point2f> spine;
  double maxLength=-1;
  if(!(nh.x==0 && nh.y==0 && nt.x==0 && nt.y == 0))
  {
    int hi,ti;
    cv::Point2f gnh=nh+bp;
    cv::Point2f gnt=nt+bp;
    hi=findContourPointClosestTo(gnh,newPoints);
    ti=findContourPointClosestTo(gnt,newPoints);
    dIdx[0]=hi;
    dIdx[1]=ti;
  }
  for (size_t i=0; i<dIdx.size()-1;++i)
  {
    for(size_t j=i+1; j<dIdx.size();j++)
    {
      std::vector<cv::Point2f> cSpine;
      std::vector<PointPair> pPairs;
      size_t di=dIdx[i];
      size_t dj=dIdx[j];
      if(blobs!=NULL && blobs->size()>0)
      {
        cv::Point2f bp((*blobs)[blobs->size()-1].minx,
            (*blobs)[blobs->size()-1].miny);
        cv::Point2f gH=heads->back()+bp;
        if(p2fdist(newPoints[dIdx[i]],gH)
            >
            p2fdist(newPoints[dIdx[j]],gH)
          )
        {
          di=dIdx[j];
          dj=dIdx[i];
        }
      }
      if(di == dj)
        continue;
      larvaDistanceMap candidateMap(newPoints);
      int sf=71;
      std::vector<float> c_tmp,c;
      vcurv[0]=(vcurv.back()+vcurv[1])/2;
      smoothVec(vcurv,c_tmp,sf,(float)0.0);
      smoothVec(c_tmp,c,sf,(float)0.0);

      /*double processNoiseC=0.00000001;
      double measurementNoiseC=0.000001;
      double errorCov=0.001;
      if(heads && heads->size()>2)
      {
        cv::KalmanFilter KF(6, 2, 0);
        cv::Mat_<float> state(6, 1); [> (x, y, Vx, Vy) <]
        cv::Mat processNoise(6, 1, CV_32F);
        cv::Mat_<float> measurement(2,1); measurement.setTo(cv::Scalar(0));

        cv::Point2f BP(blobs->at(blobs->size()-3).minx,blobs->at(blobs->size()-3).miny);
        cv::Point2f headP=heads->at(heads->size()-3);
        headP+=BP;
        
        KF.statePre.at<float>(0) = headP.x;
        KF.statePre.at<float>(1) = headP.y;
        KF.statePre.at<float>(2) = 0;
        KF.statePre.at<float>(3) = 0;
        KF.statePre.at<float>(4) = 0;
        KF.statePre.at<float>(5) = 0;
        KF.transitionMatrix = (cv::Mat_<float>(6, 6) << 1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1);
        KF.measurementMatrix = (cv::Mat_<float>(2, 6) << 1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
        setIdentity(KF.processNoiseCov, cv::Scalar::all(processNoiseC));
        setIdentity(KF.measurementNoiseCov, cv::Scalar::all(measurementNoiseC));
        setIdentity(KF.errorCovPost, cv::Scalar::all(errorCov));
        
       cv::Point2f statePt(state(0),state(1));
       cv::Mat prediction = KF.predict();

       BP=cv::Point2f(blobs->at(blobs->size()-2).minx,blobs->at(blobs->size()-2).miny);
       headP=heads->at(heads->size()-2);
       headP+=BP;
       measurement(0)=headP.x;
       measurement(1)=headP.y;
       cv::Mat estimated = KF.correct(measurement);

       prediction = KF.predict();

       BP=cv::Point2f(blobs->at(blobs->size()-1).minx,blobs->at(blobs->size()-1).miny);
       headP=heads->at(heads->size()-1);
       headP+=BP;
       measurement(0)=headP.x;
       measurement(1)=headP.y;
       estimated = KF.correct(measurement);

       size_t hi;
       if(p2fdist(newPoints[di],headP)<=p2fdist(newPoints[dj],headP))
         hi=di;
       else
         hi=dj;

       measurement(0)=newPoints[hi].x;
       measurement(1)=newPoints[hi].y;

       estimated = KF.correct(measurement);
       cv::Point2f kp(estimated.at<float>(0),estimated.at<float>(1));
       size_t nth=findContourPointClosestTo(kp,newPoints);
       if(hi==di)
         di=nth;
       else
         dj=nth;
      }
      if(tails && tails->size()>2)
      {
        cv::KalmanFilter KF(6, 2, 0);
        cv::Mat_<float> state(6, 1); [> (x, y, Vx, Vy) <]
        cv::Mat processNoise(6, 1, CV_32F);
        cv::Mat_<float> measurement(2,1); measurement.setTo(cv::Scalar(0));

        cv::Point2f BP(blobs->at(blobs->size()-3).minx,blobs->at(blobs->size()-3).miny);
        cv::Point2f tailP=tails->at(tails->size()-3);
        tailP+=BP;
        
        KF.statePre.at<float>(0) = tailP.x;
        KF.statePre.at<float>(1) = tailP.y;
        KF.statePre.at<float>(2) = 0;
        KF.statePre.at<float>(3) = 0;
        KF.statePre.at<float>(4) = 0;
        KF.statePre.at<float>(5) = 0;
        KF.transitionMatrix = (cv::Mat_<float>(6, 6) << 1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1);
        KF.measurementMatrix = (cv::Mat_<float>(2, 6) << 1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
        setIdentity(KF.processNoiseCov, cv::Scalar::all(processNoiseC));
        setIdentity(KF.measurementNoiseCov, cv::Scalar::all(measurementNoiseC));
        setIdentity(KF.errorCovPost, cv::Scalar::all(errorCov));
        
       cv::Point2f statePt(state(0),state(1));
       cv::Mat prediction = KF.predict();

       BP=cv::Point2f(blobs->at(blobs->size()-2).minx,blobs->at(blobs->size()-2).miny);
       tailP=tails->at(tails->size()-2);
       tailP+=BP;
       measurement(0)=tailP.x;
       measurement(1)=tailP.y;
       cv::Mat estimated = KF.correct(measurement);

       prediction = KF.predict();

       BP=cv::Point2f(blobs->at(blobs->size()-1).minx,blobs->at(blobs->size()-1).miny);
       tailP=tails->at(tails->size()-1);
       tailP+=BP;
       measurement(0)=tailP.x;
       measurement(1)=tailP.y;
       estimated = KF.correct(measurement);

       size_t hi;
       if(p2fdist(newPoints[di],tailP)<=p2fdist(newPoints[dj],tailP))
         hi=di;
       else
         hi=dj;

       measurement(0)=newPoints[hi].x;
       measurement(1)=newPoints[hi].y;

       estimated = KF.correct(measurement);
       cv::Point2f kp(estimated.at<float>(0),estimated.at<float>(1));
       size_t nth=findContourPointClosestTo(kp,newPoints);
       if(hi==di)
         di=nth;
       else
         dj=nth;
      }
      if(di==dj)
        continue;*/
      getCandidateSpine(blob,
          di,
          dj,
          newPoints,
          candidateMap,
          curvatures,
          c,
          50,
          heads,
          tails,
          blobs);

      if(maxLength<=candidateMap.MaxDist)
      {
        maxLength=candidateMap.MaxDist;
        bestDist=candidateMap;
      }
    }
  }

  Distances=bestDist;
  maxLength=Distances.MaxDist;

  maxDist=maxLength;

#ifdef LRV_TRACK_VISUAL_DEBUG
  std::cout << "DIDX SIZE: " << dIdx.size() << std::endl;
  std:: cout << "DIDX: " << printVector(dIdx) << std::endl;
    cv::circle(cROI, 
        MULT*(newPoints[dIdx[0]]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        7,
        cv::Scalar(0,255,255),-1);
    cv::circle(cROI, 
        MULT*(newPoints[dIdx[1]]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        7,
        cv::Scalar(0,255,255),-1);
  for(size_t i=1;i<Distances.Spine.size();++i)
  {
    cv::circle(cROI, 
        MULT*(Distances.Spine[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        3,
        cv::Scalar(0,0,255),-1);
    cv::line(cROI, 
        MULT*(Distances.Spine[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        MULT*(Distances.Spine[i-1]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        cv::Scalar(255,255,0),
        2
        );
  }
/*for(auto &pair: Distances.spinePairs)
{
  cv::line(cROI, 
      MULT*(pair.first-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      MULT*(pair.second-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      cv::Scalar(255,0,0),
      1
      );
}*/
  if(Distances.Spine.size() > 0)
    cv::circle(cROI, 
      MULT*(Distances.Spine.back()-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      4,
      cv::Scalar(255,0,255),-1);

  for(size_t i=0;i<dIdx.size();++i)
  {
    cv::circle(cROI, 
        MULT*(newPoints[dIdx[i]]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        3,
        cv::Scalar(0,255-(i*4),0),-1);
  }
  cv::circle(cROI, 
    MULT*(newPoints[dIdx.front()]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
    2,
    cv::Scalar(255,50,50),-1);
  
  /*cv::circle(cROI, 
      MULT*(newPoints[0]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      5,
      cv::Scalar(0,0,255),-1);*/

  //if((blob.label==0) && blobs!=NULL && CURRENT_FRAME>60 && CURRENT_FRAME<160 && blobs->size()>63)
  if(blob.label==2)
  {
    size_t i;
    cv::Mat plot_result;
    std::cout << "HERE START" << std::endl;
    std::cout << "=============================" << std::endl;
    std::vector<float> c_tmp1, c_tmp2, c_tmp3;// c_tmp4;
    int sf=(vcurv.size()*0.05);
    smoothVec(vcurv,c_tmp1,sf,(float)0.0);
    smoothVec(c_tmp1,c_tmp2,sf,(float)0.0);
    smoothVec(c_tmp2,c_tmp3,sf,(float)0.0);
    std::cout << printVector(c_tmp3) << std::endl;
    cv::Mat curvature_data = cv::Mat(c_tmp3).reshape(0,c_tmp3.size());
    curvature_data.convertTo(curvature_data,CV_64F);
    cv::Ptr<cv::plot::Plot2d> plot = cv::plot::createPlot2d(curvature_data);
    plot->setPlotLineColor(cv::Scalar( 50, 50, 255 ) );
    plot->render( plot_result );
    std::cerr << "HERE END" << std::endl;

    cv::Point p1 = cv::Point(
		    (double)dIdx[0]/c_tmp3.size() * plot_result.cols,
		    c_tmp3[dIdx[0]]/(
		      *std::max_element(c_tmp3.begin(),
		      c_tmp3.end()) - *std::min_element(
		      c_tmp3.begin(),c_tmp3.end())) * plot_result.rows);
    cv::Point p2 = cv::Point(
		    (double)dIdx[1]/c_tmp3.size() * plot_result.cols,
		    c_tmp3[dIdx[1]]/(
		      *std::max_element(c_tmp3.begin(),
		      c_tmp3.end()) - *std::min_element(
		      c_tmp3.begin(),c_tmp3.end())) * plot_result.rows);
    cv::circle(plot_result,p1,3,
		    cv::Scalar(255,0,255));
    cv::circle(plot_result,p2,3,
		    cv::Scalar(255,0,255));
    
    std::cout << "cROI: " << cROI.size() << std::endl;
    std::cout << "plot_result: " << plot_result.size() << std::endl;
    if(cROI.cols > plot_result.cols)
    {
	cv::Size ns=cv::Size(
			cROI.cols,
			cROI.cols/plot_result.cols * plot_result.rows);
        cv::resize(plot_result,plot_result,ns);
    }
    else if(cROI.cols < plot_result.cols)
    {
	cv::Size ns = cv::Size(plot_result.cols,
		size_t(((double)plot_result.cols/cROI.cols) * cROI.rows));
        cv::resize(cROI,cROI,ns);
    }
    cv::Mat C;
    std::cout << "cROI: " << cROI.size() << std::endl;
    std::cout << "plot_result: " << plot_result.size() << std::endl;
    cv::vconcat(cROI,plot_result,C);
    cv::imshow("C",C);
    /*for(i=0;i<Distances.Spine.size()-1;i++)
      std::cout << Distances.Spine[i] << ", " ;
    std::cout << Distances.Spine[i] << std::endl;

    for(i=0;i<Distances.Angles.size()-1;i++)
      std::cerr << Distances.Angles[i] << ", " ;
    std::cerr << Distances.Angles.back() << std::endl;*/
    cv::waitKey(-1);
    std::cout << CURRENT_FRAME << "," << blob.label << "," << Distances.WidthDist << "," <<
      Distances.WidthSD << std::endl;
    //g1.remove_tmpfiles();
    //std::cerr << "dIdx: " <<  dIdx[0] << " , " << dIdx[1] << std::endl;
    i=0;
    //std::vector<float> c;
    std::vector<float> c_tmp;
    //int sf=21;
    std::vector<double> idx;
    /*for(size_t i=0;i<vcurv.size();i++)
    {
      idx.push_back(i);
    }
    vcurv[0]=(vcurv.back()+vcurv[1])/2;
    Gnuplot g1("linespoints");
    Gnuplot g2("linespoints");
    smoothVec(vcurv,c_tmp,sf,(float)0.0);
    smoothVec(c_tmp,c,sf/2,(float)0.0);
    std::stringstream s1,s2;
    s1 << "curvature_smoothened_(" << CURRENT_FRAME << ")";
    s2 << "curvature_(" << CURRENT_FRAME << ")";
    //g1.set_yrange(-0.8,0.8);
    //g2.set_yrange(-1.5,1.5);
    g1.savetops(s1.str());
    g2.savetops(s2.str());
    g1.plot_xy(idx,c,s1.str());
    g2.plot_xy(idx,vcurv,s2.str());*/
    //cv::waitKey(100);
    //g1.remove_tmpfiles();
    //g2.remove_tmpfiles();
  }
#endif

if(Distances.Spine.size()==0)
{
  BOOST_LOG_TRIVIAL(debug) << "Error: generating spine: id: " << blob.label 
    << "Frame: " << CURRENT_FRAME; 
  Distances.Spine.resize(2);
  return;
}
  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=Distances.Spine[2].x;
  Distances.p20.y=Distances.Spine[2].y;
  Distances.p80.x=Distances.Spine[9].x;
  Distances.p80.y=Distances.Spine[9].y;
  Distances.MidPoint=0.5*(Distances.Spine[5]+Distances.Spine[6]);
  Distances.curvatureVariance=cvariance;
  return;
}

void computeSpine(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    cv::Mat &frame
    )
{
  cv::Mat contour,lcontour;
  //createLarvaContour(contour,blob,CV_8UC1,0,false);
  //cv::resize(contour,lcontour,cv::Size(),16,16,cv::INTER_CUBIC);
  //std::vector<cv::Point2f> &Dpoints=Distances.points;
  std::vector<cv::Point2f> points;
  std::vector<cv::Point2f> procPoints;
  std::vector<cv::Point2f> PlainPoints;
  std::vector<float> t;
  float csqrt;
  int PAD=2;
  blobToContourVector(blob,frame,PAD,points);
  //cv::approxPolyDP(Dpoints,points,0.4,true);
  //smoothVec(points,PlainPoints,3,cv::Point2f(0,0));
  smoothAndExtractCentripetal(frame,points,PlainPoints,t,csqrt);
  std::vector<float> scurvature;
  curvVec(PlainPoints,t,scurvature);
  /*for(size_t i=0;i<PlainPoints.size();i++)
  {
  cv::Vec3b val= frame.at<cv::Vec3b>((PlainPoints[i]));
  cv::Scalar sval(0,0,val[0]+val[1]+val[2]);
    cv::circle(frame, 
        PlainPoints[i],
        0,
        sval,-1);
  }
  cv::imshow("test",frame);*/
  std::vector<size_t> dIdx(10,0);
  std::vector<float> dmax(10,0);
  cv::Point2f bp(blob.minx,blob.miny);

  getBestCurvature(scurvature,
      dIdx,
      dmax
      );
  
  //cv::approxPolyDP(Dpoints,procPoints,0.,true);
  //cv::approxPolyDP(Dpoints,PlainPoints,0.8,true);

  //std::cerr << "computeSpine: points: " << points.size() << std::endl;
  //int PADDING=2;
  //double RES=0.5;
  //Distances.WidthDist=DBL_MAX;
  /*PointPair wpoints;
  double maxDist=0;
  std::vector<double> widths;

  cv::Mat foo;
  std::map<double,int> anglesPoints;
  //frame.copyTo(foo);
  for (size_t i=0; i<PlainPoints.size(); ++i)
  {
    double A;
    if(i==0)
    {
      A=(angleC(PlainPoints.back(),PlainPoints[i],PlainPoints[i+1]));
    }
    else if(i==PlainPoints.size()-1)
    {
      A=(angleC(PlainPoints[i-1],PlainPoints[i],PlainPoints.front()));
    }
    else
    {
      A=(angleC(PlainPoints[i-1],PlainPoints[i],PlainPoints[i+1]));
    }
    while(anglesPoints.find(A)!=anglesPoints.end())
    {
      A=A+0.000001;
    }
    anglesPoints[A]=i;*/
    /*cv::circle(foo, 
      PlainPoints[i],
      0,
      cv::Scalar(0,0,255),-1);*/
 // }

  std::vector<size_t> &candidatePoints=dIdx;
  /*std::map<double,int>::iterator f=anglesPoints.begin();
  //frame.copyTo(foo);
  size_t maxel=2<PlainPoints.size()?2:PlainPoints.size();
  for(size_t i=0;i<maxel;i++)
  {
    bool jump=false;
    for(size_t j=0;j<candidatePoints.size();j++)
    {
      if(abs(candidatePoints[j]-f->second)<2)
      {
        jump=true;
        break;
      }
    }
    if (jump)
    {
      i--;
      f++;
      continue;
    }
    candidatePoints.push_back(f->second);*/
  /*cv::circle(foo, 
      PlainPoints[candidatePoints.back()],
      0,
      cv::Scalar(0,250,0),-1);
    f++;
  }*/

  int k=0;
  k++;
  double maxDist;
  double MaxLength=-1;
  for (size_t i=0; i<candidatePoints.size()-1;i++)
  {
    for(size_t j=i+1; j<candidatePoints.size();j++)
    {
      larvaDistanceMap candidateMap(PlainPoints);
      std::vector<cv::Point2f> candidateSpine;
      /*getCandidateSpine(candidatePoints[i],
                        candidatePoints[j],
                        PlainPoints,
                        candidateMap.Spine,
                        candidateMap.Widths,
                        candidateMap.Angles,
                        candidateMap.spinePairs,
                        candidateMap.maxAngle,
                        candidateMap.maxAngleLocation,
                        candidateMap.firstHalfWidthsSum,
                        candidateMap.secondHalfWidthsSum,
                        candidateMap.WidthDist,
                        candidateMap.MaxDist,
                        30);
                        */
      if(MaxLength<candidateMap.MaxDist)
      {
        MaxLength=candidateMap.MaxDist;
        Distances=candidateMap;
      }
        
    }
  }
  maxDist=MaxLength;

  
  std::vector<PointPair>::iterator SPit;
  /*size_t i=0;
  for(SPit=Distances.spinePairs.begin();
      SPit!=Distances.spinePairs.end();
      SPit++)
  {
    cv::circle(foo, 
        SPit->first,
        0,
        cv::Scalar(250,20*i,0),-1);
    cv::circle(foo, 
        SPit->second,
        0,
        cv::Scalar(250,20*i,0),-1);
    i++;
  }

  cv::imshow("endpoints",foo);
  cv::waitKey(1);*/

  if(Distances.Spine.size()==0)
  {
    std::cerr << "Compute Spine: No spine constructed." << std::endl;
    return;
  }
  if(Distances.Spine.front()==Distances.Spine.back())
  {
    std::cerr << "Head and Tail are the same!!!" << std::endl;
  }

  //std::vector<cv::Point2f> spinepoints=spine;
  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=Distances.Spine[2].x;
  Distances.p20.y=Distances.Spine[2].y;
  Distances.p80.x=Distances.Spine[8].x;
  Distances.p80.y=Distances.Spine[8].y;
  Distances.MidPoint.x=Distances.Spine[5].x;
  Distances.MidPoint.y=Distances.Spine[5].y;
}

void computeInnerDistances(cvb::CvBlob &blob,
                           larvaDistanceMap &Distances,
                           cv::Point2f &MidPoint)
{
  cv::Mat contour;
  createLarvaContour(contour,blob);
  cv::Mat workingContour;
  contour.copyTo(workingContour);
  std::vector<cv::Point2f> &points=Distances.points;
  std::vector<cv::Point2f> SimplePoints;
  cv::approxPolyDP(points,SimplePoints,0.8,true);
  cv::Point2f MP(MidPoint.x+blob.minx,MidPoint.y+blob.miny);
  points=SimplePoints;
  double mWidth=0;
  Distances.WidthDist=DBL_MAX;

  int origNZ=countNonZero(contour);
  double MAX=0;
  for (size_t i=0; i<points.size(); ++i)
    {
      cv::Point2f p1=cv::Point2f(points[i].x,points[i].y);

      Distances[i][i]=0;
      for (size_t j=i+1; j<points.size(); ++j)
        {
          cv::Point2f p2(points[j].x,points[j].y);
          cv::line(workingContour,
                   p1,
                   p2,
                   cv::Scalar(255));
          if (countNonZero(workingContour)>origNZ)
            {
              Distances[i][j]=-1;
              Distances[j][i]=-1;
            }
          else
            {
              PointPair p1p2 = PointPair(p1,p2);
              double xdiff=points[i].x-points[j].x;
              double ydiff=points[i].y-points[j].y;
              double sqDst=xdiff*xdiff+ydiff*ydiff;
              Distances[i][j]=sqDst;
              Distances[j][i]=sqDst;

              if (MAX<Distances[i][j])
                {
                  MAX=Distances[i][j];
                  Distances.MaxDistPoints=p1p2;
                  Distances.MaxDist=MAX;
                }
              if(abs((int)(j-i))>1)
                {
                  std::vector<cv::Point2f> widthArc;
                  widthArc.push_back(p1);
                  widthArc.push_back(MP);
                  widthArc.push_back(p2);
                  mWidth=cv::arcLength(widthArc, false);
                  if (Distances.WidthDist > mWidth)
                    {
                      Distances.WidthDist=mWidth;
                      //Distances.WidthDistPoints=p1p2;
                    }
                }
            }
          contour.copyTo(workingContour);
        }
    }

  for (size_t i=0; i<points.size(); ++i)
    {
      cv::Point2f p1(points[i].x,points[i].y);
      lBFS(i,points,Distances);
    }
}

