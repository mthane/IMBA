#include "lrvTrackOL.hpp"
#ifdef LRVTRACK_CMAKE
#include "lrvTrackConfig.h"
#else
#define LRVTRACK_VERSION_MAJOR 0
#define LRVTRACK_VERSION_MINOR 0
#endif
#include "lrvTrackPartitionGenerator.hpp"
#include <boost/tokenizer.hpp>
#include <numeric>
#include <algorithm>
#include <functional>
#include "lpsolver.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <vector>

using namespace lrvTrack;

//For debugging to extract type:
std::string type2str(int type) {
  std::string r;

  uchar depth = type & CV_MAT_DEPTH_MASK;
  uchar chans = 1 + (type >> CV_CN_SHIFT);

  switch ( depth ) {
    case CV_8U:  r = "8U"; break;
    case CV_8S:  r = "8S"; break;
    case CV_16U: r = "16U"; break;
    case CV_16S: r = "16S"; break;
    case CV_32S: r = "32S"; break;
    case CV_32F: r = "32F"; break;
    case CV_64F: r = "64F"; break;
    default:     r = "User"; break;
  }

  r += "C";
  r += (chans+'0');

  return r;
}

int frameCount = 0;

//VideoWriter combOut;

Point2f filterMidpoint(vector<larvaDistanceMap> &in,
                    int idx, size_t size)
{
  int start,end;
  start=idx-((size-1)/2);
  end=idx+((size-1)/2);
  if(start<0)
    start=0;
  if(end>(int)in.size()-1)
    end=in.size()-1;

  int d=end-start+1;
  Point2f sum=in[start].MidPoint;
  for(int i=start+1;i<=end;i++)
  {
    sum+=in[i].MidPoint;
  }
  return ((double)(1.0/d))*sum;
}

Point2f filterCentroid(vector<CvBlob> &in,
                    int idx, size_t size)
{
  int start,end;
  start=idx-((size-1)/2);
  end=idx+((size-1)/2);
  if(start<0)
    start=0;
  if(end>(int)in.size()-1)
    end=in.size()-1;

  int d=end-start+1;
  Point2f sum=Point2f(in[start].centroid.x,in[start].centroid.y);
  for(int i=start+1;i<=end;i++)
  {
    sum+=Point2f(in[i].centroid.x,in[i].centroid.y);
  }
  return ((double)(1.0/d))*sum;
}

Point2f filterPoint(vector<Point2f> &in,
                    int idx, size_t size,
                    vector<CvBlob> &b)
{
  int start,end;
  start=idx-((size-1)/2);
  end=idx+((size-1)/2);
  if(start<0)
    start=0;
  if(end>(int)in.size()-1)
    end=in.size()-1;

  int d=end-start+1;
  Point2f sum=in[start]+Point2f(b[start].minx,b[start].miny);
  for(int i=start+1;i<=end;i++)
  {
    sum+=in[i]+Point2f(b[i].minx,b[i].miny);
  }
  return ((double)(1.0/d))*sum; // TS: mittlerer bewegungsrichtung je schritt
}

/*
 * Function to draw spine points of larvaObject lrv on frame img
 */
void drawSpinePoints(Mat &img, larvaObject &lrv,size_t idx)
{
  if(lrv.lrvDistances.size()==0)
    return;
  if(idx>lrv.lrvDistances.size())
    return;
  vector<Point2f> const &spine=lrv.lrvDistances[idx].Spine;
  if(spine.size()==0)
    return;
  vector<Point2f>::const_iterator it=spine.begin();
  for (;it!=spine.end();++it)
  {
    circle(img,
        *it,
        0,
        Scalar(255,0,255),
        -1);
  }
}

void drawSpinePoints(Mat &img, larvaObject &lrv)
{
  drawSpinePoints(img,lrv,lrv.lrvDistances.size()-1);
}

/*
 * Function to print brief contents of detected_larvae
 */
void dumpDetectedLarvae()
{
  BOOST_LOG_TRIVIAL(trace) << "Size: " << detected_larvae.size();
  std::map<size_t,larvaObject>::const_iterator dlit;
  dlit=detected_larvae.begin();
  stringstream d;
  d << "Contents: ";
  for(;dlit!=detected_larvae.end();++dlit)
  {
    d << dlit->first;
    d << " ";
    (dlit->second).dump();
  }
  BOOST_LOG_TRIVIAL(trace) << d.str() << endl;
}

void readIni()
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(LRVTRACK_INPUT_METADATA, pt);
  LRVTRACK_ODOR_LR = pt.get<std::string>("Trial Data.OdorA");
  if(LRVTRACK_ODOR_LR=="left" || LRVTRACK_ODOR_LR=="right")
  	cout << "Odor location from metadata: " << LRVTRACK_ODOR_LR << endl;
  else
  	LRVTRACK_ODOR_LR = pt.get<std::string>("Trial Data.OdorA_Side");
  	cout << "Odor location from metadata: " << LRVTRACK_ODOR_LR << endl;
}

void writeIni(cv::Point2f odorA_centroid)
{
  stringstream s;
  s << odorA_centroid.x << ", " << odorA_centroid.y;
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(LRVTRACK_INPUT_METADATA, pt);
  pt.put<std::string>("Trial Data.OdorA_Position",s.str());
  boost::property_tree::ini_parser::write_ini(LRVTRACK_INPUT_METADATA, pt);
}

void writeIni(double x, double y)
{
  stringstream s;
  s << x << ", " << y;
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(LRVTRACK_INPUT_METADATA, pt);
  pt.put<std::string>("Trial Data.OdorALocation",s.str());
  boost::property_tree::ini_parser::write_ini(LRVTRACK_INPUT_METADATA, pt);
}

void detectHeadTail()
{
  // Loop
}

/*void detectHeadTailLarva(larvaObject &l)
{
  double distance;
  for(size_t i=0;i<l.centroid_speed_x.size();i++)
  {
    Point2f distz=l.lrvDistances[i].Spine[0] - l.lrvDistances[i].MidPoint;
    Point2f distb=l.lrvDistances[i].Spine.back() - l.lrvDistances[i].MidPoint;
    Point2f mspeed(l.midpoint_speed_x,l.midpoint_speed_y);
  }
}*/

/*
 * Quick function for returning the average of a vector
 */
double avgVec(vector<double> &vec)
{
  return (accumulate(vec.begin(),vec.end(),0)/vec.size());
}

/*
 * Quick function for returning the average of the last N
 * values of a vector.
 */
namespace std{
  template<typename data>
    double avgNVec(vector<data> &vec,size_t N=HISTORY_SIZE)
    {
      //double SUM=0;
      size_t range;
      if (vec.size()>=N)
        range=N;
      else
        range=vec.size();

      return (accumulate(vec.rbegin(),vec.rbegin()+range,0)/range);
    }
}

/*
 * Simple factorial (we just need small numbers so this should be quick enough)
 */
uint64_t factorial(size_t n)
{
      return factorial_vec[n];
}

/*
 * Number of combinations k of n
 */
size_t kofn(size_t k, size_t n)
{
  return factorial(n)/(factorial(k)*factorial(n-k));
}

/*
 * Stirling number of the second kind
 * Set of N elements partitioned into k-blocks
 */

size_t stirling_2(size_t n, size_t k)
{
  int sum=0;
  for(size_t j=0;j<=k;++j)
    sum+=(int)(2*(((int) j+1)%2)-1) * (int)(kofn(j,k)*pow(k-j,n));
  return sum/factorial(k);
}

/* This is not exactly the powerset. It is used for the detection of clustering
 * or divergence based on the combination that matches the data best.
 * We do not need the sets of 1 element and the complete SET in there
 * Input:
 *  IN: Vector of int (representing larvae IDs) to create the powersets from.
 *  OUT: Vector of vectors with the powersets (it does not contain powersets
 *       of 1 nor the complete vector as a set).
 */
void powersets(vector<size_t> &IN, vector<vector<size_t> > &OUT){
  for (size_t i=2 ; i<IN.size();i++)
  {
    vector<size_t> pointers;
    for(size_t k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (size_t j=0 ; j<kofn(i,IN.size());j++)
    {
      vector<size_t> cvec;
      for(size_t idx=0;idx<i;idx++)
      {
        cvec.push_back(IN[pointers[idx]]);
      }
      OUT.push_back(cvec);
      for(size_t inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          size_t add=0;
          for(size_t res=inc;res<i;res++)
          {
            add++;
            pointers[res]=pointers[inc-1]+add;
          }
          break;
        }
      }
    }
  }
}

/*
 * Function to handle various keystrokes for each frame
 */
void handleKeys(bool &STEP,//bool &STOP
                bool &TRACK, bool &SHOWTAGS)
{
    int k;
    if (STEP==true)
      k=46;
    else
      k=waitKey(1);
    if (k>=0)
    {
      if (k==' ')
      {
        if(!TRACK)
          while (waitKey(1)!=' ')
          {
            //No OP
          }
      }
      if (k=='.')
      {
        STEP=true;
        while ((k=cv::waitKey(1))>=127 || k<0)
        {
        }
        if (k!='.')
          STEP=false;
      }
      if (k=='s')
      {
        SHOWTAGS=!SHOWTAGS;
      }
      if (k=='t')
      {
        TRACK=!TRACK;
      }
      if (k=='x')
      {
        exit(0);
      }
    }
}

inline double w(double x)
{
  double d=(2*x-1)*(2*x-1);
  if(x>=0.65)
    return sqrt(1-d*d);
  if(x<0.65)
    return sqrt(1-pow(1-(2*x/1.3),3));
  return 0;
}

void calculateContourPoints(Point2f &a,
                            Point2f &b,
                            Point2f &c,
                            Point2f &bp,
                            double b_index,
                            double width,
                            Point &cl,
                            Point &cr)
{
  //Make b our reference point
  Point2f r=a-b;
  Point2f f=c-b;
  Point2f sp=r+f;
  Point2f lPoint,rPoint;
  double nwidth;
  nwidth=w(b_index)*width;
  if(fabs(sp.x)<0.01 && fabs(sp.y)<0.01) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=(double)-r.x/(double)r.y;
      double xc1=sqrt((double)(nwidth*nwidth)/(double)(1+pslope*pslope));
      lPoint.x=xc1;
      lPoint.y=pslope*xc1;
      rPoint.x=-xc1;
      rPoint.y=-pslope*xc1;
    }
    else if(r.x==0 && r.y==0)
    {
      cerr << "Weird case" << endl;
      exit(0);
    }
    else if(r.x==0)
    {
      lPoint.x=-nwidth;
      lPoint.y=0;
      rPoint.x=+nwidth;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=-nwidth;
      rPoint.x=0;
      rPoint.y=+nwidth;
    }

  }
  else
  {
    double ratio=nwidth/sqrt(sp.x*sp.x+sp.y*sp.y);
    lPoint.x=ratio*sp.x;
    lPoint.y=ratio*sp.y;
    rPoint.x=-lPoint.x;
    rPoint.y=-lPoint.y;
  }

  double avals[4]={atan2f(r.y,r.x),
                    atan2f(-r.y,-r.x),
                    atan2f(rPoint.y,rPoint.x),
                    atan2f(lPoint.y,lPoint.x)};
  if(avals[0]<0)
    avals[0]=2*CV_PI+avals[0];
  if(avals[1]<0)
    avals[1]=2*CV_PI+avals[1];
  if(avals[2]<0)
    avals[2]=2*CV_PI+avals[2];
  if(avals[3]<0)
    avals[3]=2*CV_PI+avals[3];
  if((avals[2]>avals[0] || avals[2]<avals[1]) &&
     (avals[3]<avals[0] || avals[3]>avals[1]))
  {
    //points are correct lPoint is to the left, rPoint is to the right
    cl=lPoint+b-bp;
    cr=rPoint+b-bp;
  }
  else if((avals[2]<=avals[0] || avals[2]>=avals[1]) &&
          (avals[3]>=avals[0] || avals[3]<=avals[1]))
  {
    //points are reverse lPoint is to the left, rPoint is to the right
    cl=rPoint+b-bp;
    cr=lPoint+b-bp;
  }
  Point2f ppointa=a-bp;
  Point2f ppointb=b-bp;
  Point2f ppointc=c-bp;
  Point2f ppointl=lPoint+b-bp;
  Point2f ppointr=rPoint+b-bp;
}

void createMatfromSpine(Mat &larvaFitContour,
                        larvaObject &o, size_t idx)
{
  CvBlob &b=o.blobs[idx];
  larvaDistanceMap &d=o.lrvDistances[idx];
  Mat tmp=cv::Mat(b.maxy-
                  b.miny+1,
                  b.maxx-b.minx+1,
                  CV_8UC1,Scalar(0));

  Point2f bp(b.minx,b.miny);
  vector<Point2f> &spine = d.Spine;

  std::vector<Point> cpoints;
  cpoints.resize(2*spine.size()-2);
  cpoints[0]=spine[0]-bp;
  cpoints[cpoints.size()/2]=spine.back()-bp;
  for(int i=1;i<(int)spine.size()-1;i++)
  {
    Point2f a,b,c;
    a=spine[i-1]-bp;
    b=spine[i]-bp;
    c=spine[i+1]-bp;
    calculateContourPoints(spine[i-1],
                           spine[i],
                           spine[i+1],
                           bp,
                           (double) i/spine.size(),
                           (double) o.width_mean/2.0,
                           cpoints[i],
                           cpoints[cpoints.size()-i]);
    Point2f l=cpoints[i];
    Point2f r=cpoints[cpoints.size()-i];
    circle(tmp,cpoints[i],0,Scalar(100),-1);
    circle(tmp,cpoints[cpoints.size()-i],0,Scalar(100),-1);
    circle(tmp,b-bp,0,Scalar(150),-1);
    circle(tmp,a-bp,0,Scalar(150),-1);
    circle(tmp,c-bp,0,Scalar(150),-1);
  }

  size_t csz=cpoints.size();
  const Point *p=&cpoints[0];
  fillPoly(tmp,(const Point **) &p,(int *) &csz,1,Scalar(255));
  for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(50));
  }
  /*for(size_t i=0;i<cpoints.size();i++)
  {
    //line(tmp,cpoints[i-1],cpoints[i],Scalar(100));
    circle(tmp,cpoints[i],0,Scalar(100),-1);
  }*/
  tmp.copyTo(larvaFitContour);
}

void createMatfromSpine(Mat &larvaFitContour,
                        larvaObject &o)
{
  Mat tmp=cv::Mat(o.blobs.back().maxy-
                  o.blobs.back().miny+1,
        o.blobs.back().maxx-
        o.blobs.back().minx+1,
        CV_8UC1,Scalar(0));

  Point2f bp(o.blobs.back().minx,o.blobs.back().miny);
  vector<Point2f> &spine = o.lrvDistances.back().Spine;

  std::vector<Point> cpoints;
  cpoints.resize(2*spine.size()-2);
  cpoints[0]=spine[0]-bp;
  cpoints[cpoints.size()/2]=spine.back()-bp;
  cout << cpoints[0] << cpoints[cpoints.size()/2] << endl;
  for(int i=1;i<(int)spine.size()-1;i++)
  {
    circle(tmp,
           0.5*(spine[i-1]+spine[i])-bp,
           w((double)(i)/(double)spine.size())*(o.width_mean/2.0),
           Scalar(255),
           -1);
    circle(tmp,
           spine[i]-bp,
           w((double)i/(double)spine.size())*(o.width_mean/2.0),
           Scalar(255),-1);
  }

  //size_t csz=cpoints.size();
  //const Point *p=&cpoints[0];
  //cout << printVector(cpoints) << endl;
  //fillPoly(tmp,(const Point **) &p,(int *) &csz,1,Scalar(255));
  for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(100));
  }
  tmp.copyTo(larvaFitContour);
}

bool check_contour(larvaObject &o)
{
  Mat larvaFitContour;
  Mat larvaContour;
  if(o.lrvDistances.size()<1)
    return false;
  createMatfromSpine(larvaFitContour,o);
  createLarvaContour(larvaContour,
      o.blobs.back(),
      o.lrvDistances.back().spinePairs,
      o.lrvDistances.back().Spine,
      CV_8UC3,0,false,
      Scalar(255),8);

  imshow("MatFromSpine",larvaFitContour);
  imshow("Larva",larvaContour);
  waitKey(-1);
  return true;
}

bool check_roundness(size_t lID,
		     double area,
		     double length,
		     double length_mean,
		     double width,
		     double width_mean,
		     double perimeter,
		     double perimeter_mean,
                     double rndMean,
                     double curRnd,
                     //double RndRat,
                     double minHTDist_mean,
                     double distMin,
                     //double variance,
                     bool IN
                     )
{
    /*BOOST_LOG_TRIVIAL(debug) << "RND, " << CURRENT_FRAME
			<< ", " << lID
			<< ", " << area
			<< ", " << length
			<< ", " << length_mean
			<< ", " << width
			<< ", " << width_mean
			<< ", " << perimeter
			<< ", " << perimeter_mean
			<< ", " << curRnd;*/

	/*BOOST_LOG_TRIVIAL(debug) << "RND, " << CURRENT_FRAME
        << ", " << lID
	<< ", " << area
	<< ", " << length
	<< ", " << length_mean
	<< ", " << width
	<< ", " << width_mean
	<< ", " << perimeter
	<< ", " << perimeter_mean
        << ", " << curRnd;*/
      /*if(!IN)
      {
        if(curRnd<=rndMean*0.80 && distMin>minHTDist_mean*1.20)
        //if(curRnd<=2.9 && distMin>minHTDist_mean*1.20)
        {
          return true;
        }
        else
          return false;
      }
      else
      {
        if(curRnd>=rndMean*0.80 && distMin<=minHTDist_mean*1.20)
        //if(curRnd>=2.7 && distMin<=minHTDist_mean*1.20)
        {
          return true;
        }
        else
          return false;
      }*/
	//These numbers were received from one 9cm dish and one 15cm dish
	double ROUNDNESS_LIM=3.55;
	double LW_LIM=3.35;
	double LP_LIM=0.42;
	if(curRnd<ROUNDNESS_LIM &
	   length/width<LW_LIM &
	   length/width>0.001 & //Reasonable numbers...
	   length/perimeter<LP_LIM &
	   length/perimeter>0.001)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double is_larva(larvaObject &f)
{
  for(size_t i=0;i<f.blobs.size();i++)
  {
    if(!f.lrvDistances.empty() && !f.lrvDistances[i].Spine.empty())
    {
      Mat shouldLook;
      Mat actuallyLooks;
      createMatfromSpine(shouldLook,f,i);
      createLarvaContour(actuallyLooks,
                         f.blobs[i]);
    }
  }
  return 0.0;
}

void brightnessContrastGamma(//Mat &f,
    Mat &o, double b, double c, double g)
{
  if(c!=1.0 || b!=0.0)
  {
    addWeighted(o,0,o,c,b-255,o);
    //addWeighted(o,c,o,c,b,o);
  }
  if(g!=1.0)
  {
    o.convertTo(o, CV_32F);
    pow(o,g,o);
    o.convertTo(o,CV_8UC1);
  }
  /*Mat fgROI=Mat(o.rows , o.cols,o.depth(),Scalar(0));
  circle(fgROI,
      Point2f(o.rows/2.0,o.cols/2.0),
      o.rows/5,Scalar(255),-1);
  Mat f;
  f=o&fgROI;
  equalizeHist(f,f);
  imshow("eqhist", f);*/
}

/*#ifdef LRVTRACK_WITH_OPENCL
void brightnessContrastGammaOCL(ocl::oclMat &f,
                                ocl::oclMat &o,
                                double b,
                                double c,
                                double g)
{
  if(c!=1.0 || b!=0.0)
  {
    ocl::addWeighted(o,0,o,c,b,o);
  }
  if(g!=1.0)
  {
    o.convertTo(o, CV_32F);
    ocl::pow(o,g,o);
    o.convertTo(o,CV_8UC1);
  }
}
#endif*/

void showModels()
{
  cv::Point2f cc=Point2f(circles[bestCircle][0],
      circles[bestCircle][1]);
  double ppm;
  if(LRVTRACK_MPP==0.0)
    ppm=LRVTRACK_PETRIDISH/(2*circles[bestCircle][2]);
  else
    ppm=LRVTRACK_MPP;
  for(auto &m: larvaeModels)
  {
    if(!m.SUCCESS)
      continue;
    if(m.SFRAME<=CURRENT_FRAME && m.EFRAME>=CURRENT_FRAME)
    {
      for(auto &l: m.larvae_models)
      {
        string data;
        size_t c_index=CURRENT_FRAME-m.SFRAME;
        l[c_index].csvLine(CURRENT_FRAME,
                  m.SFRAME,
                  m.EFRAME,
                  VIDEO_FPS,
                  cc,
                  ppm,
                  data);
        if(data != "")
        {
          lfds[l[c_index].ID] << data;
        }

        size_t IDX=CURRENT_FRAME-m.SFRAME;
        int i=0;
        for (auto &p:l[IDX].spine)
        {
          circle(unprocessedFrame,
              p,
              0,
              Scalar(0,255,0),
              -1);
          if(p!=l[IDX].spine[0])
            line(unprocessedFrame,p,l[IDX].spine[i-1],Scalar(0,0,255));
          i++;
        }
        circle(unprocessedFrame,
            l[IDX].spine[0],
            1,
            Scalar(255,255,0),
            -1);
        circle(unprocessedFrame,
            l[IDX].spine.back(),
            1,
            Scalar(0,0,255),
            -1);
      }
    }
  }
}

void getAnglesFromSpine(vector<Point2f> &spine, vector<double> &angles, bool reverse)
{
  if(!reverse)
  {
    for(auto it=spine.begin()+1;it!=spine.end()-1;it++)
    {
      double a=angleD(*(it-1),*it,*(it+1));
      angles.push_back(a);
    }
  }
  else
  {
    for(auto it=spine.rbegin()+1;it!=spine.rend()-1;it++)
    {
      double a=angleD(*(it-1),*it,*(it+1));
      angles.push_back(a);
    }
  }
}

bool checkLarvaLength(size_t ID)
{
  larvaObject &o=detected_larvae[ID];
  size_t length=o.end_frame-o.start_frame+1;
  if(o.larva_ID == o.updated_ID)
  {
      return (length>LRVTRACK_MIN_OUTPUT_DURATION);
  }
  else
  {
    for(auto &l: detected_larvae)
    {
      larvaObject &p=l.second;
      length+=p.lifetimeWithStats;
      //length+=o.lifetimeWithStats;
    }
    for(auto &m: larvaeModels)
    {
      if(m.SUCCESS)
      {
        for(auto &lm: m.larvae_models)
        {
          if(lm[0].ID==o.updated_ID)
          {
            length+=(m.EFRAME-m.SFRAME);
            break;
          }
        }
      }
    }
    return (length>LRVTRACK_MIN_OUTPUT_DURATION);
  }
}

void make_Video(size_t id1, size_t id2, size_t id3)
{
  size_t cidx1=CURRENT_FRAME-detected_larvae[id1].start_frame;
  size_t cidx2=CURRENT_FRAME-detected_larvae[id2].start_frame;
  size_t cidx3=CURRENT_FRAME-detected_larvae[id3].start_frame;

  Mat l1,l2,l3;
  Mat lf1,lf2,lf3;

  Mat buzzF,buzz;
  greyFrame(cv::Rect(0,
        0,
        140,
        140)).copyTo(buzzF);

  createLarvaROI(greyFrame,l1,detected_larvae[id1].blobs[cidx1]);
  createLarvaROI(greyFrame,l2,detected_larvae[id2].blobs[cidx2]);
  createLarvaROI(greyFrame,l3,detected_larvae[id3].blobs[cidx3]);

  int Mr=max(l1.rows,l2.rows);
  Mr=max(Mr,l3.rows);
  int Mc=max(l1.cols,l2.cols);
  Mc=max(Mc,l3.cols);
  int M=max(Mc,Mr)+2;

  int QSize=80;

  resize(buzzF,buzz,Size(M,M));
  copyMakeBorder(l1,lf1,(M-l1.rows)/2,M-(l1.rows+(M-l1.rows)/2),
                        (M-l1.cols)/2,M-(l1.cols+(M-l1.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(l2,lf2,(M-l2.rows)/2,M-(l2.rows+(M-l2.rows)/2),
                        (M-l2.cols)/2,M-(l2.cols+(M-l2.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(l3,lf3,(M-l3.rows)/2,M-(l3.rows+(M-l3.rows)/2),
                        (M-l3.cols)/2,M-(l3.cols+(M-l3.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(lf1,lf1,(QSize-lf1.rows)/2,QSize-(lf1.rows+(QSize-lf1.rows)/2),
                        (QSize-lf1.rows)/2,QSize-(lf1.cols+(QSize-lf1.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(lf2,lf2,(QSize-lf2.rows)/2,QSize-(lf2.rows+(QSize-lf2.rows)/2),
                        (QSize-lf2.rows)/2,QSize-(lf2.cols+(QSize-lf2.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(lf3,lf3,(QSize-lf3.rows)/2,QSize-(lf3.rows+(QSize-lf3.rows)/2),
                        (QSize-lf3.rows)/2,QSize-(lf3.cols+(QSize-lf3.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(buzz,buzz,(QSize-buzz.rows)/2,QSize-(buzz.rows+(QSize-buzz.rows)/2),
                        (QSize-buzz.rows)/2,QSize-(buzz.cols+(QSize-buzz.cols)/2),
                        BORDER_CONSTANT,
                        Scalar(0));

  copyMakeBorder(lf1,lf1,2,2,2,2,
                        BORDER_CONSTANT,
                        Scalar(155));
  copyMakeBorder(lf2,lf2,2,2,2,2,
                        BORDER_CONSTANT,
                        Scalar(155));

  copyMakeBorder(lf3,lf3,2,2,2,2,
                        BORDER_CONSTANT,
                        Scalar(155));

  copyMakeBorder(buzz,buzz,2,2,2,2,
                        BORDER_CONSTANT,
                        Scalar(155));

  Mat f;
  Mat intF1,intF2;
  tile2same(lf1,lf2,intF1);
  tile2same(lf3,buzz,intF2);
  transpose(intF1,intF1);
  transpose(intF2,intF2);
  tile2same(intF1,intF2,f);
  Mat fout;
  cvtColor(f,fout,CV_GRAY2BGR);
  //combOut << fout ;
  imshow("Comb",f);
  waitKey(1);
}

void showTags2()
{
  stringstream frm;
  frm << CURRENT_FRAME;
  map<size_t,larvaObject>::iterator it=detected_larvae.begin();
  circle(unprocessedFrame,
      Point2f(circles[bestCircle][0],circles[bestCircle][1]),
      int(circles[bestCircle][2]),
      Scalar(0,255,0),1);
  int PAD=2;
  for(auto &b:cupBlobs)
  {
    Mat cupROI;
      createBlobContour(unprocessedFrame,
                         *b.second,
                         CV_8UC3,
                         PAD,
                         false,
                         Scalar(0,0,255),
                         8);
  }
  while (it!=detected_larvae.end())
  {
    if(it->second.start_frame>CURRENT_FRAME ||
        it->second.end_frame<CURRENT_FRAME)
    {
      it++;
      continue;
    }
    stringstream sstm;
    size_t c_index=CURRENT_FRAME-it->second.start_frame;
    cvb::CvBlob* blob=
      &(it->second.blobs[c_index]);
    if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
    {
      //printBlobFile(detected_larvae[it->first]);
    }
    sstm << (*it).second.updated_ID;
    if(!it->second.isCluster && it->second.round_flag[c_index])
    {
      sstm << "!";
    }
    putText(unprocessedFrame,
        sstm.str(),
        Point2f(blob->centroid.x+12,blob->centroid.y+12),
        FONT_HERSHEY_PLAIN,
        0.8,
        Scalar(255,255,255),
        1,
        CV_AA);

    Mat larvaROI;
    if(!createSimpleROI(unprocessedFrame,
          blob->minx,
          blob->miny,
          blob->maxx,
          blob->maxy,
          ROI_PADDING+PAD,
          larvaROI))
      break;
    Point2f cr(blob->centroid.x-blob->minx,blob->centroid.y-blob->miny);
    //larvaSkel testSkel(larvaROI,cr);
    //testSkel.drawSkeleton(larvaROI,Scalar(0,0,255));
    //imshow("ROI",larvaROI);
    //waitKey(1);
    //if(is_larva(blob)<IS_LARVA_THRESHOLD)
    //{
    /*
       circle(frame,
       Point2f(blob->centroid.x,blob->centroid.y),
       0,
       Scalar(255,0,0),
       -1);
       */
    //}
      cerr << "!!!! Save tracking-data1" << endl;
    if(it->second.isCluster==false)
    {
            cerr << "!!!! Save tracking-data2" << endl;
      circle(larvaROI,
          Point2f(
            it->second.heads[c_index].x+PAD,
            it->second.heads[c_index].y+PAD),
          2,
          Scalar(255,255,0),
          -1);

      circle(larvaROI,
          Point2f(
            it->second.tails[c_index].x+PAD,
            it->second.tails[c_index].y+PAD),
          2,
          Scalar(0,0,255),
          -1);

      createLarvaContour(larvaROI,
                         *blob,
                         it->second.lrvDistances[c_index].spinePairs,
                         it->second.lrvDistances[c_index].Spine,
                         CV_8UC3,PAD,false,
                         Scalar(0,255,0),8);
      drawSpinePoints(unprocessedFrame,it->second,c_index);
      cv::Point2f cc=Point2f(circles[bestCircle][0],
          circles[bestCircle][1]);
      double ppm;
      if(LRVTRACK_MPP==0.0)
        ppm=LRVTRACK_PETRIDISH/(2*circles[bestCircle][2]);
      else
        ppm=LRVTRACK_MPP;

      cerr << "!!!! Save tracking-data3" << endl;
      
      if(LRVTRACK_USE_MODEL)
      {
        fs::path data_folder(LRVTRACK_DATE+"-data");
        if(!fs::exists(data_folder))
        {
          fs::create_directory(data_folder);
        }
      }
      else
      {
        fs::path data_folder(LRVTRACK_DATE+"-data");
        if(!fs::exists(data_folder))
        {
          fs::create_directory(data_folder);
        }
      }
      //Here we must check whether the complete added track length is
      //long enough.
      //if((it->second.end_frame-it->second.start_frame)>LRVTRACK_MIN_OUTPUT_DURATION)
      if (checkLarvaLength(it->second.larva_ID))
      {
        larvaObject &l=it->second;
        stringstream f;
        if(LRVTRACK_USE_MODEL)
        {
          f << LRVTRACK_DATE+"-data/" << l.updated_ID << ".csv";
        }
        else
          f << LRVTRACK_DATE+"-data/" << l.updated_ID << ".csv";
        fs::path filename(f.str());
        //stringstream data;

        if(lfds[l.updated_ID].is_open()==false)
          lfds[l.updated_ID].open(filename.string(),fstream::out | fstream::app);
            string data;
        l.csvLine(CURRENT_FRAME,
            VIDEO_FPS,
            cc,
            ppm,
            data);
        if(data != "")
          lfds[l.updated_ID] << data;
	else
	{
	  BOOST_LOG_TRIVIAL(debug) << "SKIP: "
		  << c_index << ", " << l.updated_ID << endl;
	}
        /*if(l.blobs.size()-1==c_index)
        {
          lfds[l.updated_ID].close();
        }*/
      }
    }
    it++;
  }
  if(LRVTRACK_USE_MODEL)
    showModels();
  if(CURRENT_FRAME==TOTAL_FRAMES)
  {
    for(auto &f:lfds)
    {
      f.second.close();
    }
  }
}

void showTags(const cvb::CvBlobs &tracked_blobs)
{
      for(auto &blob_p:tracked_blobs)
      {
        stringstream sstm;
        cvb::CvBlob &blob=*blob_p.second;
        sstm << blob.label;
        sstm << "(" << detected_larvae[blob.label].blobs.back().n20 <<")";
        if(detected_larvae[blob.label].round_flag.size()>0 &&
           detected_larvae[blob.label].round_flag.back()==true)
          sstm << "!";
        putText(colorFrame,
            sstm.str(),
            Point2f(blob.centroid.x+12,blob.centroid.y+12),
            FONT_HERSHEY_PLAIN,
            0.8,
            Scalar(255,255,255),
            1,
            CV_AA);

        int PAD=2;
        Mat larvaROI;
        if(!createSimpleROI(colorFrame,
              blob.minx,
              blob.miny,
              blob.maxx,
              blob.maxy,
              ROI_PADDING+PAD,
              larvaROI))
          break;
        Point2f cr(blob.centroid.x-blob.minx,blob.centroid.y-blob.miny);

        if(detected_larvae[blob.label].isCluster==false
           &&
           detected_larvae[blob.label].heads.size()>0)
        {
          circle(larvaROI,
              Point2f(
                detected_larvae[blob.label].heads.back().x+PAD,
                detected_larvae[blob.label].heads.back().y+PAD),
              1,
              Scalar(255,0,0),
              -1);

          circle(larvaROI,
              Point2f(
                detected_larvae[blob.label].tails.back().x+PAD,
                detected_larvae[blob.label].tails.back().y+PAD),
              1,
              Scalar(0,0,255),
              -1);

          createLarvaContour(larvaROI,blob,CV_8UC3,PAD,false,
              Scalar(0,255,0),8);
          if(detected_larvae[blob.label].lrvDistances.size()>0)
            drawSpinePoints(unprocessedFrame,detected_larvae[blob.label]);

          //plotAngle(blob,larvaROI,PAD);
        }
      }
}

/*
 * Function to see if the centre of blob1 matches the centre of blob2.
 * blob1 and blob2 are blobs coming from different frames and the
 * function replies true if the blobs are close enough.
 * Input:
 *  blob1, blob2: the blobs to match
 *  val: the Manhattan distance between centres
 *  factor: defines the threshold for the match (factor*LARVA_OBJECT_LENGTH)
 */
bool centresMatch(
    cvb::CvBlob *blob1,
    cvb::CvBlob *blob2,
    double &val,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1.0)
{
  double objectLength=
      min(max(blob1->maxx-blob1->minx,blob1->maxy-blob1->miny),
               max(blob2->maxx-blob2->minx,blob2->maxy-blob2->miny));

  BOOST_LOG_TRIVIAL(trace) << "CentresMatchS ["<< blob1->label << ", "
                           << blob2->label << "]: Length: " << objectLength
                           << " Difx: "
                           << fabs(blob1->centroid.x - blob2->centroid.x)
                           << " Dify: "
                           << fabs(blob1->centroid.y - blob2->centroid.y)
                           << " Threshold: "
                           << factor*LARVA_OBJECT_LENGTH;

    val=fabs(blob1->centroid.x - blob2->centroid.x) + fabs(blob1->centroid.y - blob2->centroid.y) ;

  if (fabs(blob1->centroid.x - blob2->centroid.x) < factor*LARVA_OBJECT_LENGTH &&
      fabs(blob1->centroid.y - blob2->centroid.y)< factor*LARVA_OBJECT_LENGTH )
    {
      return true;
    }
  else
    {
      return false;
    }
}

/*
 * Function to see if the centres of blobs in "larvae" matches the centre of blob.
 * larvae and blob2 are derived from different frames and the
 * function replies true if the barycentre if blobs in In and blob2
 * are close enough.
 * Input:
 *  In: Complete list of Blobs referenced by the larvae vector
 *  larvae: vector containing indices of larvae to consider for the comparison
 *          with blob. The indices refer to blobs in "In"
 *  blob2: the blob to match with larvae
 *  val: the Manhattan distance between barycentres
 *  factor: defines the threshold for the match (factor*LARVA_OBJECT_LENGTH)
 */
bool centresMatch(
    cvb::CvBlobs &In,
    cvb::CvBlob *blob,
    vector<size_t> &larvae,
    double &val,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1)
{
  double xcomb=0, ycomb=0;
  double objectLength=
      max(blob->maxx-blob->minx,blob->maxy-blob->miny);
  double lrvAreaSum=0;
  if(larvae.size()==1)
  {
    xcomb=In[larvae[0]]->centroid.x;
    ycomb=In[larvae[0]]->centroid.y;
  }
  else
  {
    vector<size_t>::iterator it=larvae.begin();
    while(it!=larvae.end())
    {
      xcomb+=In[*it]->centroid.x*In[*it]->area;
      ycomb+=In[*it]->centroid.y*In[*it]->area;
      lrvAreaSum+=In[*it]->area;
      ++it;
    }
    xcomb=xcomb/lrvAreaSum;
    ycomb=ycomb/lrvAreaSum;
  }
  BOOST_LOG_TRIVIAL(trace) << "CentresMatchP [" << blob->label
                           << ", " << printVector(larvae) << "]: Length: "
                           << objectLength << " Difx: "
                           << fabs(blob->centroid.x - xcomb)
                           << " Dify: " << fabs(blob->centroid.y - ycomb)
                           << " Threshold: " << factor*LARVA_OBJECT_LENGTH;
  val=fabs(blob->centroid.x - xcomb)+fabs(blob->centroid.y - ycomb);
  if (fabs(blob->centroid.x - xcomb)< factor* LARVA_OBJECT_LENGTH &&
      fabs(blob->centroid.y - ycomb)< factor* LARVA_OBJECT_LENGTH )
    {
      return true;
    }
  else
    {
      return false;
    }
}

bool larvaNearRing(cvb::CvBlob &blob)
{
  if(circles.empty())
    return false;
  cv::Point2f cntr(circles[bestCircle][0],circles[bestCircle][1]);
  if(p2fdist(blob.centroid.x,
        blob.centroid.y,
        circles[bestCircle][0],
        circles[bestCircle][1])
      >
      circles[bestCircle][2]*0.95)
  {
    return true;
  }
  return false;
}


// is blob not usable, i.e.  in odor cup or outside petri dish than returns true, otherwise false
// is the larva/blob inside our odor cups or petri dish?
bool larvaToRing(cvb::CvBlob &blob)
{
  // if no cups or dish have been found before, it cannot be inside
  if(circles.empty() && cupContours.empty())
    return false;

  Mat ltrROI; // temporary image

  // create rois inside cupContours using bounding box of blob
  // and store them inside temporary ltrRoi
  // its a rectangular cutout from cupContours in the size blob bounding box
  cupContours(cv::Rect(blob.minx,
        blob.miny,
        blob.maxx-blob.minx+1,
        blob.maxy-blob.miny+1)).copyTo(ltrROI);

  // if this cutout contains only zeroes, there was no overlap
  size_t nz=countNonZero(ltrROI);

  // with all pixels being zero, there is no overlap,
  // the larva/blob cannot be inside cup
  if(nz==0)
    return false;

  // draw filled polygon as blob
  createBlobContour(ltrROI,
                    blob);

  // check weather there have been nonzero pixels outside of blobl/polygon
  // if yes, it "touches" the ring or is completely inside 
  if((size_t) countNonZero(ltrROI) < nz+blob.area)
    return true;
  return false;

}

/*bool larvaToRing(cvb::CvBlob &blob)
{
  if(circles.empty())
    return false;
  cv::Point2f cntr(circles[bestCircle][0],circles[bestCircle][1]);
  if(p2fdist(blob.centroid.x,
             blob.centroid.y,
             circles[bestCircle][0],
             circles[bestCircle][1])
      >
      circles[bestCircle][2]*0.95)
  {
    std::vector<cv::Point2f> p;
    blobToPointVector(blob,p);
    for(size_t i=0;i<p.size();i++)
    {
      Point2f cp=p[i];
      double dst=p2fdist(p[i],cntr);
      if(dst>=circles[bestCircle][2]-2.0)
        return true;
    }
  }
  return false;
}*/

bool larvaToCup(cvb::CvBlob &blob)
{
  Mat larvaROI;
    if(!createSimpleROI(cupContours,
          blob.minx,
          blob.miny,
          blob.maxx,
          blob.maxy,
          2,
          larvaROI))
      return true;

    size_t white;
    white=countNonZero(larvaROI);
    if(white==0)
      return false;
    else
    {
      createLarvaContour(larvaROI,
          blob,
          CV_8UC3,2,true,
          Scalar(255),8);
      if((size_t) countNonZero(larvaROI)<white+blob.area)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

}

/*
 * Quick function to judge if the sizes of BLOB1 and BLOB2 are comparable.
 */
bool blobSizeIsRelevant(
    cvb::CvBlob *BLOB1,
    cvb::CvBlob *BLOB2,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  return (ratio*BLOB1->area > BLOB2->area &&
      ((2-ratio)*BLOB1->area < BLOB2->area));
}

/*
 * Quick function to judge if the sizes of blobs in larvae and BLOB are
 * comparable.
 */
bool blobSizeIsRelevant(
    cvb::CvBlobs &In,
    cvb::CvBlob *BLOB,
    vector<size_t> &larvae,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  vector<size_t>::iterator IT=larvae.begin();
  double areaSUM=0;
  while (IT!=larvae.end())
  {
    areaSUM+=In[*IT]->area;
    ++IT;
  }
  return (ratio*BLOB->area > areaSUM &&
      ((2-ratio)*BLOB->area < areaSUM ));
}

/*
 * Returns all the larvae in the set Blobs within an area around "Blob".
 * The area is determined by the size of the Blob (the longest side
 * of the bounding box multiplied by PADRatio).
 *
 * The vector nearbyLarvae is filled by those found sorted from closest
 * to furthest.
 */
void getNearbyLarvae(cvb::CvBlobs &Blobs, cvb::CvBlob *Blob,
		            vector<size_t> &nearbyLarvae,bool pre=true,
                double PADRatio=1.5)
{
  vector<double> distances;
	double MaxDist = max(Blob->maxx-Blob->minx,Blob->maxy-Blob->miny);
	MaxDist=PADRatio*MaxDist;
  cvb::CvBlobs::iterator It=Blobs.begin();
  while (It!=Blobs.end())
	{
		cvb::CvBlob *cBlob=It->second;
    if (cBlob->centroid.x < (Blob->maxx + MaxDist/2) &&
        cBlob->centroid.x > (Blob->minx - MaxDist/2) &&
        cBlob->centroid.y < (Blob->maxy + MaxDist/2) &&
        cBlob->centroid.y > (Blob->miny - MaxDist/2) &&
        ( ((pre==false) && (assignedNew[It->first].size()<=0)) ||
          ((pre==true)  && (assignedPrevious[It->first].size()<=0)))
        )
    {
      double DIST=fabs(Blob->centroid.x - cBlob->centroid.x) +
        fabs(Blob->centroid.y - cBlob->centroid.y);
      vector<double>::iterator dIt=distances.begin();
      vector<size_t>::iterator nIt=nearbyLarvae.begin();
      if(nearbyLarvae.size()>0)
      {
        while(dIt!=distances.end())
        {
          if(DIST < *dIt)
          {
            distances.insert(dIt,DIST);
            nearbyLarvae.insert(nIt,It->first);
            break;
          }
          ++dIt;
          ++nIt;
        }
        if(dIt==distances.end())
        {
          nearbyLarvae.push_back(It->first);
          distances.push_back(DIST);
        }
      }
      else
      {
        distances.push_back(DIST);
        nearbyLarvae.push_back(It->first);
      }
    }
		++It;
	}
}

/***** Assigning Functions ********

 * These functions are just used as a quick way to set the correct IDs for
 * the assignedNew and assignedPrevious vectors used by the main tracking
 * algorithm. The assume the the IDs have already been correctly identified
 * and they do not change the global larvae structures (e.g. detected_larvae).
*/

/*
 * Function to assign an ID to a larva.
 * preID: is the ID of the larva from the previous frame to which preID
 *         matches
 * postID: is the ID of the blob of the current frame
 */
/*void assign_one(size_t preID,size_t postID)
{
  assignedPrevious[preID].push_back(postID);
  assignedNew[postID].push_back(preID);
  BOOST_LOG_TRIVIAL(trace) << "Assigning: b" << postID << " -> p" << preID;
  assignedPreMap[preID]=1;
}*/

/*
 * Function to match an ID to several IDs
 * preID: is the ID of the blob of the previous frame
 * postID: is the vector with the IDs of the larva in the new frame which
 *         diverged from preID
 */
/*void assign_one(size_t preID,
                vector<size_t> postIDs)
{
  assignedPreMap[preID]=postIDs.size();
  //detected_larvae[preID].isCluster=true;
  for(auto &postID:postIDs)
  {
    size_t NewID=++LARVAE_COUNT;
    assignedNew[postID].push_back(NewID);
    BOOST_LOG_TRIVIAL(trace) << "Assigning: " << postID << " -> " << NewID;
    assignedPrevious[preID].push_back(NewID);
    detected_larvae[preID].diverged_to.push_back(NewID);
    parent_blobs[NewID].push_back(preID);
    children_blobs[preID].push_back(NewID);
  }
  BOOST_LOG_TRIVIAL(debug) << CURRENT_FRAME << ", "
                           <<  preID
                           << " -> "
                           << printVector(assignedPrevious[preID]);
  assignedPreMap[preID]=postIDs.size();
}*/

/*
 * Function to match several IDs to an ID
 * preID: is the vector with the IDs of the larva in the old frame which
 *         collided to form blob with postID
 * postID: is the ID of the blob of the new frame
 */
/*void assign_one(vector<size_t> preIDs,
                size_t postID,size_t newID)
{
  assignedNew[postID].push_back(newID);
  BOOST_LOG_TRIVIAL(debug) << CURRENT_FRAME << ", "
                           << printVector(preIDs) << " -> "
                           << newID;
  for(auto &preID:preIDs)
  {
    assignedNew[postID].push_back(preID);
    assignedPrevious[preID].push_back(postID);
    assignedPreMap[preID]=1;
    children_blobs[preID].push_back(newID);
    parent_blobs[newID].push_back(preID);
  }
}*/

/*
 * Function called whenever clustering is detected. Calls the
 * relevant assign one function and sets the correct details
 * in the detected_clusters vector. Assumes that correct IDs
 * have already been provided.
 * Input:
 *  POST_ID: the ID of the cluster
 *  IDs: the IDs of the blobs before the cluster
 */
/*void assign_clustering(
                      size_t POST_ID,
                      vector<size_t> &IDs
                      )
{
  size_t CLUSTER_ID=++LARVAE_COUNT;
  //newClusters.push_back(CLUSTER_ID);
  //detected_converged[CLUSTER_ID]=IDs;
  parent_blobs[CLUSTER_ID]=IDs;
  for(auto &i:IDs)
    children_blobs[i].push_back(CLUSTER_ID);
  assign_one(IDs,POST_ID,CLUSTER_ID);
}*/

/*
 * Function to calculate the Mahalanobis distance between
 * blob N and larva C.
 * Input:
 *  N: The ID of the Blob in the NEW CvBlobs structure which contains
 *     the latest detected blobs.
 *  C: The ID of the larvae in the detected_larvae structure.
 * Output:
 *  the distance
 */
 double mh_dist(size_t N,size_t C)
{
  Mat candidateCovarMat;
  Mat candidateMeanMat;

  Mat newMeanMat;

  Mat Responses;
  Mat TrainArray;
  /*float size_avg=detected_larvae[C].area_sum/
    detected_larvae[C].area.size();
  float grey_value_avg=detected_larvae[C].grey_value_sum/
    detected_larvae[C].grey_value.size();
  float length_avg=detected_larvae[C].length_sum/
    detected_larvae[C].length.size();
  float perimeter_avg=detected_larvae[C].perimeter_sum/
    detected_larvae[C].perimeter.size();
  float width_avg=detected_larvae[C].width_sum/
    detected_larvae[C].width.size();*/

  float size_avg=avgNVec(detected_larvae[C].area);
  float grey_value_avg=avgNVec(detected_larvae[C].grey_value);
  float length_avg=avgNVec(detected_larvae[C].length);
  float perimeter_avg=avgNVec(detected_larvae[C].perimeter);
  float width_avg=avgNVec(detected_larvae[C].width);

  //float speed_x=avgNVec(detected_larvae[C].midpoint_speed_x);
  //float speed_y=avgNVec(detected_larvae[C].midpoint_speed_y);


  Mat InputArray;
  hconcat(Mat(detected_larvae[C].area),
      Mat(detected_larvae[C].grey_value),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].length),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].perimeter),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].width),
      InputArray);

  //hconcat(InputArray,
  //    Mat(detected_larvae[C].midpoint_speed_x),
  //    InputArray);

  //hconcat(InputArray,
  //    Mat(detected_larvae[C].midpoint_speed_y),
  //    InputArray);

  vector<float> meanVec;
  meanVec.push_back(size_avg);
  meanVec.push_back(grey_value_avg);
  meanVec.push_back(length_avg);
  meanVec.push_back(perimeter_avg);
  meanVec.push_back(width_avg);

  Mat(meanVec).copyTo(candidateMeanMat);
  Mat meanTMat;
  transpose(candidateMeanMat,meanTMat);


  calcCovarMatrix(InputArray,
                  candidateCovarMat,
                  meanTMat,
                  CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
  candidateCovarMat.convertTo(candidateCovarMat,CV_32F);
  invert(candidateCovarMat,candidateCovarMat,DECOMP_SVD);

  Mat newSamplesMat;

  //Setup of new larva
  Mat larvaROI;
  vector<Point2f> newLarvaPoints;
  blobToPointVector(*NEW[N],newLarvaPoints);
  createLarvaContour(larvaROI,(*NEW[N]));
  larvaDistanceMap dstLarva(newLarvaPoints);
  Point2f centroid;
  centroid.x=NEW[N]->centroid.x;
  centroid.y=NEW[N]->centroid.y;
  //computeSpine(*NEW[N],dstLarva,grey_frame);
  fixContour(*NEW[N],
             dstLarva,
             LRVTRACK_CONTOUR_RESOLUTION,
             colorFrame,
             previousFrame);

  float newSize=NEW[N]->area;
  float newGreyval=getGreyValue(larvaROI,*NEW[N],greyFrame);
  float newLength=dstLarva.MaxDist;
  float newPerimeter=getPerimeter(*NEW[N]);
  float newWidth=dstLarva.WidthDist;

  meanVec.clear();
  meanVec.push_back(newSize);
  meanVec.push_back(newGreyval);
  meanVec.push_back(newLength);
  meanVec.push_back(newPerimeter);
  meanVec.push_back(newWidth);

  newMeanMat=Mat(meanVec);

 // cerr << candidateMeanMat << endl;
 // cerr << "==============================" << endl;
 // cerr << newMeanMat << endl;
 // cerr << "==============================" << endl;
 // cerr << candidateCovarMat << endl;

  double ret=Mahalanobis(candidateMeanMat,newMeanMat,candidateCovarMat);
  return ret;
}

/*
 * Function to ask whether assuming that blobP is assigned to blobN
 * is reasonable based on the speed calculated by this assignment.
 */
bool speedMatch(cvb::CvBlob &blobP,
                cvb::CvBlob &blobN,
                //double duration,
                double max_speed,
                double pre_speed_x=0.0,
                double pre_speed_y=0.0)
{
  double frames=(CURRENT_FRAME-detected_larvae[blobP.label].lastFrameWithStats);
  double mduration=frames/VIDEO_FPS;
  double uduration=max_speed*(1.0-(0.25-1/(frames+3)));
  double speedx = (blobP.centroid.x - blobN.centroid.x)/mduration;
  double speedy = (blobP.centroid.y - blobN.centroid.y)/mduration;

  speedx=speedx-pre_speed_x;
  speedy=speedy-pre_speed_y;

  double speed = sqrt(speedx*speedx + speedy*speedy);
  //cerr << "SpeedMatch: P:" << blobP.label << " N:" << blobN.label <<
  //  " M/S: " << uduration << "," << speed << endl;
  if(speed<uduration)
    return true;
  else
    return false;
}

/* Function to check if a Mapping makes sense. Currently uses speed
 * and a quick size check to figure this out. It may use other
 * options in the future.
 */
bool isMappingReasonable(lrvMapping &p
                        //,double duration
                        )
{
      size_t lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
      //dumpDetectedLarvae();
      cvb::CvBlob &blobP=detected_larvae[p.plrv].blobs[lastIdx];
      cvb::CvBlob &blobN=*NEW[p.nlrv];
      //int csx=detected_larvae[p.plrv].centroid_speed_x.size();
      //int csy=detected_larvae[p.plrv].centroid_speed_y.size();
      double pre_speed_x=avgNVec(detected_larvae[p.plrv].centroid_speed_x,3);
      double pre_speed_y=avgNVec(detected_larvae[p.plrv].centroid_speed_y,3);
      if(!speedMatch(
            blobP,
            blobN,
            //duration,
            detected_larvae[p.plrv].max_centroid_speed,
            pre_speed_x,
            pre_speed_y) ||
          !blobSizeIsRelevant(&blobP,&blobN))
       return false;
      else
        return true;
}

/* Function to check if an assignment (several mappings) makes sense.
 * As in the above function we're using speed
 * and a quick size check to figure this out. It may use other
 * options in the future.
 */
bool isMappingReasonable(ltAssignments &m,double duration)
{
  for(size_t i=0;i<m.size();i++)
  {
      lrvMapping &p=m[i];
      size_t lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
      cvb::CvBlob &blobP=detected_larvae[p.plrv].blobs[lastIdx];
      cvb::CvBlob &blobN=*NEW[p.nlrv];
      //int csx=detected_larvae[p.plrv].centroid_speed_x.size();
      //int csy=detected_larvae[p.plrv].centroid_speed_y.size();
      //double pre_speed_x=avgNVec(detected_larvae[p.plrv].centroid_speed_x,3);
      //double pre_speed_y=avgNVec(detected_larvae[p.plrv].centroid_speed_y,3);
      if(!speedMatch(
            blobP,
            blobN,
            duration,
            detected_larvae[p.plrv].max_centroid_speed) ||
          !blobSizeIsRelevant(&blobP,&blobN))
        return false;
  }
  return true;
}

/* Check if assignment contains a mapping of o->n
 */
bool assignmentContainsNewOld(ltAssignments &m,size_t n,size_t o)
{
  for(size_t i=0;i<m.size();i++)
  {
    if(m[i].nlrv==n || m[i].plrv==o)
      return true;
  }
  return false;
}

/* Given a vector of all possible mappings
 * returns the powerset of all possible assignments */
void pair_powersets(vector<lrvMapping> &IN,
    vector<ltAssignments > &OUT)
{
  for (size_t i=1 ; i<=IN.size();i++)
  {
    vector<size_t> pointers;
    for(size_t k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (size_t j=0 ; j<kofn(i,IN.size());j++)
    {
      ltAssignments cvec;
      for(size_t idx=0;idx<i;idx++)
      {
        if(!assignmentContainsNewOld(cvec,
                                  IN[pointers[idx]].nlrv,
                                  IN[pointers[idx]].plrv))
          cvec.push_back(IN[pointers[idx]]);
      }
      if(cvec.size()==i)
        OUT.push_back(cvec);
      for(size_t inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          size_t add=0;
          for(size_t res=inc;res<i;res++)
          {
            add++;
            pointers[res]=pointers[inc-1]+add;
          }
          break;
        }
      }
    }
  }
}

//Adds a set of assignmets following a diverging event
//to the possible assignments map.
void assignMapping(ltAssignments &m,
                   map<size_t,size_t> &newAssignments)
{
  for(size_t i=0;i<m.size();i++)
  {
    newAssignments[m[i].nlrv]=m[i].plrv;
  }
}

//Return the accuracy of a set of mappings
double mappingAccuracy(ltAssignments &a)
{
  double SUM=0;
  for(size_t i=0;i<a.size();i++)
  {
    SUM+=a[i].getDistance();
  }
  return SUM/(1*(1+(a.size()-1)*1.5));
}

//Generate all the available reasonable mappings
void createReasonableMappings(
  vector<size_t> &candidateLarvae,
  vector<size_t> &newLarvae,
  //double duration,
  vector<ltAssignments> &accepted_mappings)
{
    ltAssignments initial_pairs;
    vector<lrvMapping> valid_pairs;
    map<size_t,size_t> mappable_new;
    for(size_t i=0;i<candidateLarvae.size();i++)
    {
      for(size_t j=0;j<newLarvae.size();j++)
      {
        lrvMapping n(newLarvae[j],candidateLarvae[i]);
        //n.candidates=candidateLarvae;
        if(isMappingReasonable(n
              //,duration
              ))
        {
          n.setDistance();
          valid_pairs.push_back(n);
        }
      }
    }
  pair_powersets(valid_pairs,accepted_mappings);
}

/* This is not going to be used in the Off Line version.
 * Main function Handling diverging larvae
 * Input:
 *  candidateLarvae: vector with the IDs of the Larvae that
 *    were known to be involved in the collision
 *  newLarvae: vector with the IDs of the unmapped blobs.
 * Output:
 *  newAssignments: mappints of the type [NEW]=OLD;
 */
void diverge_match(
  vector<size_t> &candidateLarvae,
  vector<size_t> &newLarvae,
  map<size_t, size_t> &newAssignments//,
  //double duration
  )
{
  vector<ltAssignments> valid_mappings;
  createReasonableMappings(
    candidateLarvae,
    newLarvae,
    //duration,
    valid_mappings);

  map<double,size_t> mapping_accuracy;
  for(size_t i=0;i<valid_mappings.size();i++)
  {
    double acc=mappingAccuracy(valid_mappings[i]);
    //cerr << "=========================" << endl;
    //cerr << "Mapping " << i << ": " << endl;
    //for(size_t j=0;j<valid_mappings[i].size();j++)
    //{
    //  valid_mappings[i][j].print();
    //}
    //cerr << "Total: " << acc << endl;
    //cerr << "=========================" << endl;
    mapping_accuracy[acc]=i;
  }

  if(valid_mappings.size()!=0)
    assignMapping(valid_mappings[mapping_accuracy.begin()->second],
                newAssignments);
}


/*
 * For the Offline Version we will be treating every diverging case
 * as if the cluster were always new. The only difference would be
 * the registration of the fact that the new larvae derived from the
 * old larvae.
 */
/*void assign_divergingOL(//cvb::CvBlobs &New,
                      size_t CLUSTER_ID,
                      vector<size_t> &IDs
                      )
{
  // We have the following cases here:
  //  1) CLUSTER_ID Corresponds to no cluster: This means
  //     the cluster is newly found (it started as a cluster)
  //  2) CLUSTER_ID Corresponds to an existing cluster:
  //     we assign all the new blobs as decendants of the previous blob
  map<size_t,vector<size_t> >::iterator dcIT;
  //Not found new cluster NEW IDs to be given to the vector
  //elements
  BOOST_LOG_TRIVIAL(trace) << "Cluster " << CLUSTER_ID
                  << " diverged. Assigning new IDs for diverged larvae";
  assign_one(CLUSTER_ID,IDs);
  //detected_larvae[CLUSTER_ID].isCluster=true;
}
*/
/*int detect_diverging(vector<size_t> &preLarvaeNearby,
                       vector<size_t> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  BOOST_LOG_TRIVIAL(trace) << "Trying to detect diverging clusters";

  BOOST_LOG_TRIVIAL(trace)<< "Size of newLarvaeNearby: "
                          << newLarvaeNearby.size();

  if(newLarvaeNearby.size()<=1)
    return -1; // No diverging is possible
  // Check the complete set first
  vector<size_t>::iterator pIT=preLarvaeNearby.begin();
  while(pIT!=preLarvaeNearby.end())
  {
    //BOOST_LOG_TRIVIAL(trace) << "Checking if nodes "
    //                         << printVector(newLarvaeNearby)
    //                         << " diverged from: " << *pIT;
    double v;
    if(centresMatch(New,Pre[*pIT],newLarvaeNearby,v))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      //BOOST_LOG_TRIVIAL(trace) << "Node " << *pIT << " matches "
      //                         << printVector(newLarvaeNearby);
      assign_divergingOL(//New,
          *pIT,newLarvaeNearby);
      break;
    }
    if(newLarvaeNearby.size()>2)
    {
      //BOOST_LOG_TRIVIAL(trace) << "Checking powersets of "
      //                         << printVector(newLarvaeNearby)
      //                         << " that diverged from: " << *pIT;
      vector<vector<size_t> > pSETS;
      powersets(newLarvaeNearby,pSETS);
      vector<vector<size_t> >::iterator pSIT=pSETS.begin();
      while(pSIT!=pSETS.end())
      {
        if(centresMatch(New,Pre[*pIT],*pSIT,v))
        {
          // Centres of all candidates match with new blob
          // cluster contains all. We can return
          assign_divergingOL(//New,
              *pIT,*pSIT);
          break;
        }
        ++pSIT;
      }
    }
    ++pIT;
  }
  return 0;
}*/

/*int detect_clustering(vector<size_t> &preLarvaeNearby,
                       vector<size_t> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  BOOST_LOG_TRIVIAL(trace) << "Trying to detect clusters";

  if(preLarvaeNearby.size()<=1)
    return -1; // No clustering is possible

  // Check the complete set first
  vector<size_t>::iterator nIT=newLarvaeNearby.begin();
  while(nIT!=newLarvaeNearby.end())
  {
    double v;
    if(centresMatch(Pre,New[*nIT],preLarvaeNearby,v) &&
       blobSizeIsRelevant(Pre,New[*nIT],preLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      //BOOST_LOG_TRIVIAL(trace) << "Centres of new larva: "
      //                         << *nIT << " match centres of "
      //                         << printVector(preLarvaeNearby)
      //                         << ". Assigning clustering";

      assign_clustering(*nIT,preLarvaeNearby);
      [> LarvaFit Test starts here: <]

      [> LarvaFit Test ends here: <]

      //continue;
      break;
    }
    vector<vector<size_t> > pSETS;
    powersets(preLarvaeNearby,pSETS);
    vector<vector<size_t> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      BOOST_LOG_TRIVIAL(trace)
        << "Trying to detect subsets for clustering";
      if(centresMatch(Pre,New[*nIT],*pSIT,v) &&
          blobSizeIsRelevant(Pre,New[*nIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains subset pSIT. We can return
        //BOOST_LOG_TRIVIAL(trace) << "Centres of new larva: "
        //                         << *nIT << " match centres of subset "
        //                         << printVector(*pSIT)
        //                         << ". Assigning clustering"
        //                         << endl;
        assign_clustering(*nIT,*pSIT);
        break;
      }
      ++pSIT;
    }
    ++nIT;
  }
  return 0;
}*/

void findHeadTail(larvaObject &lrv,
                  Point2f &Head,
                  Point2f &Tail//,
                  //bool force_SurroundingValSearch=false
                  )
{
  cvb::CvBlob blob=lrv.blobs.back();
  vector<Point2f> &spine=lrv.lrvDistances.back().Spine;
  Point2f bp(lrv.blobs.back().minx,lrv.blobs.back().miny);
  Point2f sp_front=spine[0]-bp;
  Point2f sp_back=spine.back()-bp;
  {
    double hfdiff=diff(lrv.heads.back(),sp_front);
    double hbdiff=diff(lrv.heads.back(),sp_back);
    double tfdiff=diff(lrv.tails.back(),sp_front);
    double tbdiff=diff(lrv.tails.back(),sp_back);
    if (hfdiff+tbdiff<hbdiff+tfdiff)
    {
      Head=sp_front;
      Tail=sp_back;
    }
    else
    {
      Head=sp_back;
      Tail=sp_front;
    }

  }
}

void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<size_t, larvaObject> &NEW_LARVA)
{
  size_t ID=(*it).first;
  cvb::CvBlob blob=*((*it).second);
  Mat larvaROI,cntPoints;
  createLarvaContour(larvaROI,blob);
  createLarvaContourPoints(cntPoints,blob);

  map<size_t,larvaObject>::iterator curLarva;
  // NEW LARVA OBJECT!
  if ((curLarva=detected_larvae.find(ID))==detected_larvae.end()
      ||
      Prev==In)
  {
    // Create and allocate the new object
    larvaObject newLarva;

    // Set the frame of it's existence
    newLarva.start_frame=CURRENT_FRAME;
    newLarva.end_frame=CURRENT_FRAME;

    // Give the larva the necessary ID
    newLarva.larva_ID=ID;
    newLarva.updated_ID=newLarva.larva_ID;//This may be changed later

    // Add the blob of the larva to its blob history
    newLarva.blobs.push_back(blob);

    // State that the larva is not in a blob
    newLarva.inCluster.push_back(false);

    Point2f centroid=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    Point2f centroidf=Point2f(
        (blob.centroid.x),
        (blob.centroid.y)
        );

    newLarva.capture_times.push_back(CURRENT_FRAME/VIDEO_FPS);

    // Initialize the speed to 0
    newLarva.midpoint_speed_x.push_back(0);
    newLarva.midpoint_speed_y.push_back(0);
    newLarva.max_midpoint_speed=0;

    newLarva.centroid_speed_x.push_back(0);
    newLarva.centroid_speed_y.push_back(0);
    newLarva.max_centroid_speed=0;

    newLarva.centroids.push_back(centroid);
    newLarva.centroids_full.push_back(centroidf);


    ++newLarva.lifetimeWithStats;
    newLarva.lastBlobWithStats=0;

    if(blob.area * LRVTRACK_MPP *LRVTRACK_MPP > 6.0)
      newLarva.isCluster=true;
    else
      newLarva.isCluster=false;

    //Initialize the area values
    newLarva.area.push_back(blob.area);
    newLarva.area_mean=blob.area;
    newLarva.area_sum=blob.area;
    newLarva.area_max=newLarva.area_min=blob.area;
    newLarva.area_min=newLarva.area_min=blob.area;

    if(find(certain_blobs.begin(),certain_blobs.end(),blob.label)!=certain_blobs.end())
    {
      newLarva.lastFrameWithStats=CURRENT_FRAME;

      //detected_larvae[ID]=newLarva;
      //NEW[ID]=newLarva;
      tbb::concurrent_hash_map<size_t,larvaObject>::accessor a;
      NEW_LARVA.insert(a,ID);
      a->second=newLarva;
      return;
    }

    double greyVal=getGreyValue(larvaROI,blob,greyFrame);
    newLarva.grey_value.push_back(greyVal);
    newLarva.grey_value_mean = greyVal;
    newLarva.grey_value_sum= greyVal;
    newLarva.grey_value_max = greyVal;
    newLarva.grey_value_min = greyVal;

    double perimeter=getPerimeter(blob);
    newLarva.perimeter.push_back(perimeter);
    newLarva.perimeter_mean=perimeter;
    newLarva.perimeter_sum=perimeter;
    newLarva.perimeter_max=perimeter;
    newLarva.perimeter_min=perimeter;

    // Initialize the speed to 0
    newLarva.centroid_distance_x.push_back(0);
    newLarva.centroid_distance_y.push_back(0);
    newLarva.centroid_distance_x_sum=0;
    newLarva.centroid_distance_y_sum=0;

    newLarva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));
    newLarva.roundness_mean=newLarva.roundness.back();
    newLarva.roundness_sum=newLarva.roundness.back();
    //cerr << newLarva.larva_ID << ": " << newLarva.roundness.back() << endl;
    // In this block we compute the inner spine for the larva
    vector<Point2f> cntPoints;
    blobToPointVector(blob,cntPoints);
    larvaDistanceMap Distances(cntPoints);
    //computeSpine(blob,Distances,frame);
    //fixContourSimple(blob,Distances,
    if(!newLarva.isCluster)
    {
      fixContour(blob,Distances,
          LRVTRACK_CONTOUR_RESOLUTION,
          colorFrame,previousFrame);

      newLarva.lrvDistances.push_back(Distances);
      newLarva.length.push_back(Distances.MaxDist);
      newLarva.length_mean = Distances.MaxDist;
      newLarva.length_sum= Distances.MaxDist;
      newLarva.length_max = Distances.MaxDist;
      newLarva.length_min = Distances.MaxDist;
      Point2f Head,Tail;
      Point2f bp(newLarva.blobs.back().minx,newLarva.blobs.back().miny);
      Head=Distances.Spine[0]-bp;
      Tail=Distances.Spine.back()-bp;
      newLarva.heads_brightness.push_back(getSurroundingSize(Head,blob,greyFrame));
      newLarva.tails_brightness.push_back(getSurroundingSize(Tail,blob,greyFrame));
      newLarva.angular_speed.push_back(0);

      newLarva.heads.push_back(Head);
      newLarva.tails.push_back(Tail);

      Point2f MP;
      MP.x=Distances.MidPoint.x-newLarva.blobs.back().minx;
      MP.y=Distances.MidPoint.y-newLarva.blobs.back().miny;
      Point2f AxS(MP.x,Tail.y);
      newLarva.headBodyAngle.push_back(angleD(Tail,MP,Head));
      newLarva.orientationAngle.push_back(cvb::cvAngle(&blob));

      newLarva.width.push_back(Distances.WidthDist);
      newLarva.width_mean = Distances.WidthDist;
      newLarva.width_sum= Distances.WidthDist;
      newLarva.width_max = Distances.WidthDist;
      newLarva.width_min = Distances.WidthDist;

      double distMin=min(
          p2fdist(Head,newLarva.heads.back())+p2fdist(Tail,newLarva.tails.back()),
          p2fdist(Tail,newLarva.heads.back())+p2fdist(Head,newLarva.tails.back()));
      newLarva.minHTDist.push_back(distMin);
      newLarva.minHTDist_sum=newLarva.minHTDist.back();
      newLarva.minHTDist_mean=newLarva.minHTDist_sum;

      /*bool curVal = check_roundness(newLarva.larva_ID,
		      newLarva.area.back(),
		      newLarva.length.back(),
		      newLarva.length_mean,
		      newLarva.width.back(),
		      newLarva.width_mean,
		      newLarva.perimeter.back(),
		      newLarva.perimeter_mean,
		      newLarva.roundness.back(),
		      newLarva.roundness.back(),
		      //RndRat,
		      newLarva.minHTDist_mean,
		      distMin,//var,
		      false);*/
      if(Distances.Spine.size()!=2)
	      newLarva.round_flag.push_back(false);
      else
	      newLarva.round_flag.push_back(true);
    }


    newLarva.lastFrameWithStats=CURRENT_FRAME;

    //detected_larvae[ID]=newLarva;
    //NEW[ID]=newLarva;
    tbb::concurrent_hash_map<size_t,larvaObject>::accessor a;
    NEW_LARVA.insert(a,ID);
    a->second=newLarva;
    //check_contour(newLarva);
  }
  // UPDATED LARVA OBJECT!
  else
  {
    //Reference to current larva
    larvaObject &cur_larva=(*curLarva).second;
    //Pointer for the previous blob
    cvb::CvBlob &preBlob = cur_larva.blobs.back();
    cur_larva.end_frame++;

    // Set the ID of the larvaObject to the ID found TODO:Probably unnecessary
    cur_larva.larva_ID=ID;
    cur_larva.updated_ID=ID;

    // Add the current blob to the blobs history of the larva
    cur_larva.blobs.push_back(blob);

    // Create the skeleton of the larva and add it to the skeletons history
    Point2f centroid=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    Point2f centroidf=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    cur_larva.centroids_full.push_back(centroidf);
    cur_larva.centroids.push_back(centroid);

    cur_larva.capture_times.push_back(CURRENT_FRAME/VIDEO_FPS);


    // If not then:
    //  Update area values for larva.

    ++cur_larva.lifetimeWithStats;
    cur_larva.lastBlobWithStats=cur_larva.blobs.size()-1;
    //cur_larva.isCluster=false;

    cur_larva.area.push_back(blob.area);

    //cur_larva.area_mean=(cur_larva.area_mean+blob.area)/2;
    cur_larva.area_sum = cur_larva.area_sum + blob.area;
    cur_larva.area_mean=cur_larva.area_sum/cur_larva.area.size();
    if (cur_larva.area_max < blob.area)
    {
      cur_larva.area_max=blob.area;
    }
    if (cur_larva.area_min > blob.area)
    {
      cur_larva.area_min=blob.area;
    }

    if(find(certain_blobs.begin(),certain_blobs.end(),blob.label)!=certain_blobs.end())
    {
      cur_larva.lastFrameWithStats=CURRENT_FRAME;
      return;
    }
    if(cur_larva.isCluster)
    {
      cur_larva.lastFrameWithStats=CURRENT_FRAME;
      return;
    }
    // Try to find Head and Tail
    //larvaSkel newLarvaSkel(larvaROI,centroid);
    //cur_larva.lrvskels.push_back(newLarvaSkel);

    cur_larva.centroid_distance_x.push_back(
        fabs(blob.centroid.x - preBlob.centroid.x));
    cur_larva.centroid_distance_y.push_back(
        fabs(blob.centroid.y - preBlob.centroid.y));
    cur_larva.centroid_distance_x_sum+=
      fabs(blob.centroid.x - preBlob.centroid.x);
    cur_larva.centroid_distance_x_sum+=
      fabs(blob.centroid.y - preBlob.centroid.y);


    // Point coordinates for head/tail
    Point2f Head,Tail;

    // Map to keep the distances of each point to all the others
    // Pair of points to keep the points with the Maximum distance
    // (i.e. head and tail :) )
    vector<Point2f> cntPoints;
    blobToPointVector(blob,cntPoints);
    larvaDistanceMap Distances(cntPoints);
    // Compute all the inner distances for the larva
    //computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
    //computeSpine(blob,Distances,frame);

    /*fixContourSimple(blob,*/
    fixContour(blob,
        Distances,
        LRVTRACK_CONTOUR_RESOLUTION,
        colorFrame,
        previousFrame,
        &cur_larva.heads,
        &cur_larva.tails,
        &cur_larva.blobs);

    double TIMEFRAME;
    if(cur_larva.inCluster.back()==0)
    {
      TIMEFRAME=1/VIDEO_FPS;
    }
    else
    {
      TIMEFRAME=(CURRENT_FRAME-cur_larva.lastFrameWithStats)/VIDEO_FPS;
    }
    cur_larva.midpoint_speed_x.push_back(
        (Distances.MidPoint.x - cur_larva.lrvDistances.back().MidPoint.x)
        /TIMEFRAME);
    cur_larva.midpoint_speed_y.push_back(
        (Distances.MidPoint.y - cur_larva.lrvDistances.back().MidPoint.y)
        /TIMEFRAME);

    cur_larva.centroid_speed_x.push_back(
        (blob.centroid.x - cur_larva.blobs[cur_larva.blobs.size()-2].centroid.x)
        /TIMEFRAME);
    cur_larva.centroid_speed_y.push_back(
        (blob.centroid.y - cur_larva.blobs[cur_larva.blobs.size()-2].centroid.y)
        /TIMEFRAME);

    double mspeed=sqrt(cur_larva.midpoint_speed_x.back()*
        cur_larva.midpoint_speed_x.back() +
        cur_larva.midpoint_speed_y.back()*
        cur_larva.midpoint_speed_y.back());
    if(fabs(cur_larva.max_midpoint_speed) < fabs(mspeed))
    {
      cur_larva.max_midpoint_speed=mspeed;
    }

    double cspeed=sqrt(cur_larva.centroid_speed_x.back()*
        cur_larva.centroid_speed_x.back() +
        cur_larva.centroid_speed_y.back()*
        cur_larva.centroid_speed_y.back());
    if(fabs(cur_larva.max_centroid_speed) < fabs(cspeed))
    {
      cur_larva.max_centroid_speed=cspeed;
    }

    cur_larva.lrvDistances.push_back(Distances);
    cur_larva.length.push_back(Distances.MaxDist);
    cur_larva.length_mean=(cur_larva.length_mean+Distances.MaxDist)/2;
    cur_larva.length_sum=cur_larva.length_sum+Distances.MaxDist;
    if (cur_larva.length_max < Distances.MaxDist)
    {
      cur_larva.length_max=Distances.MaxDist;
    }
    if (cur_larva.length_min > Distances.MaxDist)
    {
      cur_larva.length_min=Distances.MaxDist;
    }

    cur_larva.width.push_back(Distances.WidthDist);
    cur_larva.width_sum=cur_larva.width_sum+Distances.WidthDist;
    cur_larva.width_mean=cur_larva.width_sum/cur_larva.width.size();
    if (cur_larva.width_max < Distances.WidthDist)
    {
      cur_larva.width_max=Distances.WidthDist;
    }
    if (cur_larva.width_min > Distances.WidthDist)
    {
      cur_larva.width_min=Distances.WidthDist;
    }

    double greyVal=getGreyValue(larvaROI,blob,greyFrame);
    cur_larva.grey_value.push_back(greyVal);
    cur_larva.grey_value_mean=(cur_larva.grey_value_mean+greyVal)/2;
    cur_larva.grey_value_sum=cur_larva.grey_value_sum+greyVal;
    if (cur_larva.grey_value_max < greyVal)
    {
      cur_larva.grey_value_max=greyVal;
    }
    if (cur_larva.grey_value_min > greyVal)
    {
      cur_larva.grey_value_min=greyVal;
    }

    double perimeter=getPerimeter(blob);
    cur_larva.perimeter.push_back(perimeter);
    cur_larva.perimeter_mean=(cur_larva.perimeter_mean+perimeter)/2;
    cur_larva.perimeter_sum=cur_larva.perimeter_sum+perimeter;
    if (cur_larva.perimeter_max < perimeter)
    {
      cur_larva.perimeter_max=perimeter;
    }
    if (cur_larva.perimeter_min > perimeter)
    {
      cur_larva.perimeter_min=perimeter;
    }

    //execute the function to find which is which and assign appropriately
    //to Head/Tail.
    Head=Distances.Spine[0];
    Tail=Distances.Spine.back();
    findHeadTail(cur_larva,Head,Tail);
    double distMin=min(
            p2fdist(Head,cur_larva.heads.back())+p2fdist(Tail,cur_larva.tails.back()),
            p2fdist(Tail,cur_larva.heads.back())+p2fdist(Head,cur_larva.tails.back()));
    cur_larva.minHTDist.push_back(distMin);
    cur_larva.minHTDist_sum+=cur_larva.minHTDist.back();
    cur_larva.minHTDist_mean=cur_larva.minHTDist_sum/cur_larva.minHTDist.size();
    cur_larva.heads_brightness.push_back(getSurroundingSize(Head,blob,greyFrame));
    cur_larva.tails_brightness.push_back(getSurroundingSize(Tail,blob,greyFrame));

    //double preRnd=cur_larva.roundness.back();
    double rndMean=cur_larva.roundness_mean;
    cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));
    cur_larva.roundness_sum+=cur_larva.roundness.back();
    cur_larva.roundness_mean=cur_larva.roundness_sum/cur_larva.roundness.size();
    double curRnd=cur_larva.roundness.back();
    //double RndRat=curRnd/preRnd;
    //double var=Distances.curvatureVariance;

    //if(cur_larva.round_flag.back()==false &&
    /*bool curVal = check_roundness(cur_larva.larva_ID,
	       	       cur_larva.area.back(),
		       cur_larva.length.back(),
		       cur_larva.length_mean,
		       cur_larva.width.back(),
		       cur_larva.width_mean,
		       cur_larva.perimeter.back(),
		       cur_larva.perimeter_mean,
                       rndMean,
                       curRnd,
                       //RndRat,
                       cur_larva.minHTDist_mean,
                       distMin,//var,
                       false);*/
    //{
    if(Distances.Spine.size()!=2)
    	cur_larva.round_flag.push_back(false);
    else
        cur_larva.round_flag.push_back(true);
    //}
    /*else if(cur_larva.round_flag.back()==true &&
            check_roundness(cur_larva.larva_ID,
			    cur_larva.area.back(),
			    cur_larva.length.back(),
			    cur_larva.length_mean,
			    cur_larva.width.back(),
			    cur_larva.width_mean,
			    cur_larva.perimeter.back(),
			    cur_larva.perimeter_mean,
                            rndMean,
                            curRnd,
                            //RndRat,
                            cur_larva.minHTDist_mean,
                            distMin,//var,
                            true))
    {
      cur_larva.round_flag.push_back(false);
    }
    else if(cur_larva.round_flag.back()==true)
    {
      cur_larva.round_flag.push_back(true);
    }
    else if(cur_larva.round_flag.back()==false)
    {
      cur_larva.round_flag.push_back(false);
    }*/

    cur_larva.heads.push_back(Head);
    cur_larva.tails.push_back(Tail);

    Point2f MP;
    MP.x=Distances.MidPoint.x-cur_larva.blobs.back().minx;
    MP.y=Distances.MidPoint.y-cur_larva.blobs.back().miny;
    cur_larva.headBodyAngle.push_back(angleD(Tail,MP,Head));
    Point2f AxS(MP.x,Tail.y);
    cur_larva.orientationAngle.push_back(cvb::cvAngle(&blob));

    //double curAngle=cvb::cvAngle(&blob);
    //double preAngle=cvb::cvAngle(&preBlob);

    //cur_larva.angular_speed.push_back(fast_abs(curAngle-preAngle)*VIDEO_FPS);

    //state that larva is not detected as part of a blob of larvae
    //NOTE: It is important to perform this value setting !AFTER!
    //      searching for the head tail because, the head tail
    //      search behaves differently if the larva was previously
    //      part of a blob.
    cur_larva.inCluster.push_back(0);

    if(LRVTRACK_CSV_OUTPUT)
    {
      csvfile << CURRENT_FRAME/VIDEO_FPS <<
        " , " << cur_larva.larva_ID;
      vector<Point2f> &spine=cur_larva.lrvDistances.back().Spine;
      if(cur_larva.heads.back() == spine.back())
      {
        Point2f &t=spine[0];
        csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM <<
          " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
        vector<Point2f>::iterator it=spine.begin();
        it+=2;
        for(;it!=spine.end();it++)
        {
          //Point2f cp=*it-t;
          Point2f cp=*it;
          csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , "
            << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
        }

      }
      else
      {
        Point2f &t=spine.back();
        csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM <<
          " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
        vector<Point2f>::reverse_iterator it=spine.rbegin();
        it+=2;
        for(;it!=spine.rend();it++)
        {
          //Point2f cp=*it-t;
          Point2f cp=*it;
          csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , "
            << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
        }

      }
      //csvfile << " , " << cur_larva.centroids.back().x + blob.minx <<
      //" , " << cur_larva.centroids.back().y + blob.miny <<
      //" , " << cur_larva.inCluster.back() <<
      csvfile << endl;
    }

    cur_larva.lastFrameWithStats=CURRENT_FRAME;
    //check_contour(cur_larva);
  }
}

void bg_without_larvae(Mat &fgImg)
{
  IplImage *fflabelImg;

  cvb::CvBlobs blobs;
  IplImage ipl_thresholded = fgImg;

  fflabelImg=cvCreateImage(
      cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);

  cvLabel(&ipl_thresholded, fflabelImg, blobs);
  cvb::cvFilterByArea(blobs, LRVTRACK_MIN_OBJ_SIZE, LRVTRACK_MAX_OBJ_SIZE);
  cvb::CvBlobs::iterator it=blobs.begin();
  //double c,max=0.0;
  while (it!=blobs.end())
  {
    Mat ROI;
    cvb::CvBlob &blob=*(it->second);
    createLarvaROI(fgImg, ROI, blob);
    createLarvaContour(ROI,
                        blob,
                        CV_8UC1,
                        0,
                        true,
                        Scalar(0),
                        8);

    ++it;
  }
  cvReleaseImage(&labelImg);
}


/* The function removes the background from the new frame
   and applies thresholding to derive the image in which we look
   for blobs.
    Input:
       * newFrame: the new frame from the get_next_frame function
       * greyBgFrame: a grey image of the background
       * showFrame: the BGR image returned to the user
    Output:
       * processedFrame: the frame containing only the greyscale shapes
                         of the larvae
*/
void process_frame(Mat &newFrame,
                  Mat &greyBgFrame,
                  Mat &showFrame,
                  Mat &processedFrame)
{
    //int THRESH=THRESH_BINARY;
    //size_t CMAX=255;
    //size_t NB=137;
    //int THRESHOLD=-30;
    //int METHOD=ADAPTIVE_THRESH_MEAN_C;
    //int METHOD=ADAPTIVE_THRESH_GAUSSIAN_C;
    Mat fgFrame; //Foreground frame
    Mat fgROI; // Foreground region of interest (only the part we want)
    Mat fgImage; //Final image where the foreground is grey and (
                  // hopefully everything else is black

    // fgFrame is created by the absolute difference of the
    // newFrame and the background. This gives us the best
    // flexibility also when the larva are over areas that we
    // consider background.
    // TODO: Check if this is true:
    // The addWeighted would provide us with a fuzzier
    // removal
    // .
    //absdiff(newFrame,greyBgFrame,fgFrame);
    if(LRVTRACK_EXTRACT_OFFLINEBG)
      //absdiff(newFrame,greyBgFrame,fgFrame);
      addWeighted(newFrame, 1.0, greyBgFrame, -1.0, 0.0, fgFrame);
    else
      addWeighted(newFrame, 1.0, greyBgFrame, -4.0, 0.0, fgFrame);
    //equalizeHist(fgFrame,fgFrame);
    normalize(fgFrame,fgFrame,0,255,CV_MINMAX);


    //imshow("foreground",fgFrame);
    //waitKey(-1);
    fgROI=Mat(fgFrame.rows , fgFrame.cols,fgFrame.depth(),Scalar(0));

    //Use the registered bestCircle to get construct the ROI as a circle
    //where all things outside the circle (petri-dish) are black.
    if(circles.size()>0)
    {
      circle(fgROI,
          Point2f(circles[bestCircle][0],circles[bestCircle][1]),
          int(circles[bestCircle][2]),Scalar(255),-1);

      if(!showFrame.empty())
      {
        // Add the circle as a marker to the showFrame
        circle(showFrame,
            Point2f(circles[bestCircle][0],circles[bestCircle][1]),
            int(circles[bestCircle][2]),
            Scalar(0,255,0),1);
        //if(cups.em
      }
    }
    else
    {
      // This is a case where there is no circle found
      fgROI=Mat(greyBgFrame.rows ,
          greyBgFrame.cols,
          greyBgFrame.depth(),
          Scalar(255));
    }

    for(auto &b:cupBlobs)
    {
      Mat cupROI;
      createBlobContour(showFrame,
          *b.second,
          CV_8UC3,
          2,
          false,
          Scalar(0,0,255),
          8);
    }

    // Same process for odor cups
    /*if(cups.size()>0)
    {
      for(size_t i=0; i<cups.size(); ++i)
      {
        circle(fgROI,
            Point2f(cups[i][0],cups[i][1]),
            int(cups[i][2]),
            Scalar(0),
            -1);
        circle(showFrame,
            Point2f(cups[i][0],cups[i][1]),
            int(cups[i][2]),
            Scalar(0,0,255),
            1);
      }
    }*/

    // Construct a complete image of the BW fgROI and the whole frame
    fgImage=fgFrame&fgROI;
    Mat test(fgImage.rows,fgImage.cols,fgImage.type());
    // Apply a thresholding to extract those elements that escape.
    //equalizeHist(fgImage,fgImage);
    //normalize(fgImage,fgImage,0,255,CV_MINMAX);
    //lrvTrackBilateral(fgImage,test,21,10,10);
    //adaptiveBilateralFilter(fgImage,test,Size(51,51),90,90);
    //GaussianBlur(fgImage,fgImage,Size(5,5),0);
    threshold(fgImage,processedFrame,0,255,THRESH_OTSU+THRESH_BINARY);
    //imshow("processFrame", fgImage);
    /*adaptiveThreshold(fgImage,
        processedFrame,
        CMAX,
        METHOD,
        THRESH,
        NB,
        THRESHOLD);*/
    Mat element = getStructuringElement(MORPH_CROSS, Size(3, 3));
    // Here we dilate since our thresholding effort will almost always
    // removes outer parts of the contour
    dilate(processedFrame,processedFrame,element);
    //dilate(processedFrame,processedFrame,element);
    // The processedFrame is the outcome of the good image and filtered for
    // noise by the processedFrame
    processedFrame=fgImage&processedFrame;
    //imshow("PF",processedFrame);
    //processedFrame=processedFrame*2;
    //normalize(processedFrame,processedFrame,0,255,CV_MINMAX);
    //equalizeHist(processedFrame, processedFrame);
    //
	
}

//Match lost larvae
/*size_t matchLostLarvae(size_t newLabel)
{
  for(vector<size_t>::iterator lIT=lost_blobs.begin();
      lIT!=lost_blobs.end();lIT++)
  {
      cvb::CvBlob &blobN=*NEW[newLabel];
      if(detected_larvae[*lIT].blobs.size()==0)
        return 0;
      cvb::CvBlob &blobP=detected_larvae[*lIT].blobs.back();
      double duration=(CURRENT_FRAME
             - detected_larvae[*lIT].lastFrameWithStats)/VIDEO_FPS;
      if(speedMatch(
            blobP,
            blobN,
            duration,
            detected_larvae[*lIT].max_centroid_speed) &&
          blobSizeIsRelevant(&blobP,&blobN))
      {
        size_t ret=*lIT;
        lost_blobs.erase(lIT);
        return ret;
      }
  }
  return 0;
}*/

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  //cvb::CvBlobs::iterator it=In.begin();
  tbb::concurrent_hash_map<size_t, larvaObject> NEW_LARVA;
  larvaeUpdateBody body(In,Prev,NEW_LARVA);
  if(LRVTRACK_PARALLEL)
    parallel_for_(Range(0, In.size()), body,LRVTRACK_THREADS);
  else
  {
    for(cvb::CvBlobs::iterator it=In.begin();it!=In.end();it++)
    {
      updateOneLarva(In,Prev,it,NEW_LARVA);
    }
    //parallel_for_(Range(0, In.size()), body,1.0);
  }

  tbb::concurrent_hash_map<size_t, larvaObject>::iterator nit=
    NEW_LARVA.begin();
  while (nit!=NEW_LARVA.end())
  {
    detected_larvae[nit->first]=nit->second;
    is_larva(nit->second);
    ++nit;
    }
}


/* larvae tracking: does pairwise mapping of all larvae from previous to current frame and vice-versa
 *
 * next considers all mapping cases:
 *  i)   1 to 1:    unique match, very good
 *  ii)  many to 1: multiple larvae/blobs have fused
 *  iii) 1 to many: one blob has split into multiple laevae
 *  iv)  0 to 1:    larva popped up, e.g. first frame or entered from cup
 *  v)   1 to 0:    larva has vanished, e.g. in cups or out of petri dish
 *
  */

void contourLarvaeTrack(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs &out)
{

  //new_matches_pre: contains for each cvblobID returned by the blob extraction
  //                 the matching larvaIDs from the previous frame
  std::map<size_t, vector<size_t> > new_matches_pre;

  //pre_matched_by_new: contains for each larvaID from the previous frame all
  //                    those cvblobIDs that matched it
  std::map<size_t, vector<size_t> > pre_matched_by_new;

  map<size_t, size_t> new_assigned;

  certain_blobs.empty();

  //First step:
  // Check larvae from the new frame matching contours of larvae in the
  // previous frame.
  // If the match we:
  //  1) for the cvblobID add the matching previous larvaIDs in 'new_matches_pre'
  //  2) for each larvaID matched push cvblobID in 'pre_matched_by_new'
  //  3) initialize new_assigned map keeping the assignmend of each cvBlobID to a
  //      given larvaID
  for(auto &new_blob_p:In)
  {
    cvb::CvBlob &new_blob=*new_blob_p.second;
    new_assigned[new_blob.label]=0;
    for(auto &pre_blob_p:Prev)
    {
      cvb::CvBlob &pre_blob=*pre_blob_p.second;
      //Quick screening if the bounding boxes don't even touch we skip
      if( new_blob.minx > pre_blob.maxx ||
          new_blob.miny > pre_blob.maxy ||
          new_blob.maxx < pre_blob.minx ||
          new_blob.maxy < pre_blob.miny )
      {
        continue;
      }
      //Ratio of pixels belonging to new_blob that also belong in the old blob
      double n1;
      //Ratio of pixels belonging to old_blob that also belong in the new blob
      double o1;
      blobToBlobOverlap(new_blob,
          pre_blob,
          n1,
          o1);
      //If the overlap is more than 20% new_blob is derived from pre_blob
      if(n1>0.1) // CAUTION: frame rate dependent magic constant
      {
        new_matches_pre[new_blob.label].push_back(pre_blob.label);
        pre_matched_by_new[pre_blob.label].push_back(new_blob.label);
      }
    }
  }

  /*//Second step:
  // Check larvae from the previous frame matching contours of larvae in the
  // new frame.
  // If the match we:
  //  1) for the larvaID add the matching previous larvaIDs in 'pre_matches_new'
  //  2) for each cvblobID matched push cvblobID in 'new_matched_by_pre'
  for(auto &pre_blob_p:Prev)
  {
    cvb::CvBlob &pre_blob=*pre_blob_p.second;
    for(auto &new_blob_p:Prev)
    {
      cvb::CvBlob &new_blob=*new_blob_p.second;
      //Quick screening if the bounding boxes don't even touch we skip
      if( pre_blob.minx > new_blob.maxx ||
          pre_blob.miny > new_blob.maxy ||
          pre_blob.maxx < new_blob.minx ||
          pre_blob.maxy < new_blob.miny )
      {
        continue;
      }
      //Ratio of pixels belonging to pre_blob that also belong in the old blob
      double n1;
      //Ratio of pixels belonging to old_blob that also belong in the new blob
      double o1;
      blobToBlobOverlap(pre_blob,
          new_blob,
          n1,
          o1);
      //If the overlap is more than 20% pre_blob is derived from new_blob
      if(n1>0.2)
      {
        pre_matches_new[pre_blob.label].push_back(new_blob.label);
        new_matched_by_pre[new_blob.label].push_back(pre_blob.label);
      }
    }
  }*/
  //Here we assign:
  for(auto &new_blob_p:In)
  {
    size_t cid=new_blob_p.second->label;
    cvb::CvBlob &new_blob=*new_blob_p.second;
    //If the current cvblobID is matched by only one larva in the previous frame
    if(new_matches_pre[cid].size()==1)
    {
      //TODO: Check whether there are other new that match the old
      //if not check whether there's good correcpondence between objects
      //if not give new number.
      size_t pid=new_matches_pre[cid][0];
      if(pre_matched_by_new[pid].size()==1 && pre_matched_by_new[pid][0]==cid)
      {
        double n1;
        double o1;
        blobToBlobOverlap(*In[cid],*Prev[pid],n1,o1);
        if(fabs(n1-o1)>0.3 && larvaNearRing(new_blob))
        {
          //Check for split and one goes to border
          size_t BLOB_ID=++LARVAE_COUNT;
          new_blob.n20=new_blob.label;
          new_blob.label=BLOB_ID;
          new_assigned[new_blob.n20]=new_blob.label;
          if(n1>o1)
          {
            //Split and border
            BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: "
              << CURRENT_FRAME
              << " WARN,Check for split to border on frame: " << CURRENT_FRAME
              << ", pre-new, " << pid << "-" << BLOB_ID;
            partial_lost_blobs[CURRENT_FRAME].push_back(pid);
          }
          else
          {
            //Merge and border
            BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: "
              << CURRENT_FRAME
              << " WARN,Check for merge from border on frame: " << CURRENT_FRAME
              << ", pre-new, " << pid << "-" << BLOB_ID;
          }
        }
        else
        {
          new_blob.n20=new_blob.label;
          new_blob.label=pid;
          new_assigned[new_blob.n20]=new_blob.label;
          /*BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: 1->1," << pid
                                   << "->" << new_blob.n20<<"("
                                   << new_blob.label <<")";*/
        }
      }
    }
    //If nothing was really matching the previous frame we have a new object
    //we give it a new ID
    else if(new_matches_pre[cid].size()==0)
    {
      size_t CLUSTER_ID=++LARVAE_COUNT;
      int OUTVAL;
      new_blob.n20=new_blob.label;
      new_blob.label=CLUSTER_ID;
      new_assigned[new_blob.n20]=new_blob.label;
      Point2f pcentroid;
      pcentroid.x=new_blob.centroid.x;
      pcentroid.y=new_blob.centroid.y;
      Point2f center(circles[bestCircle][0],circles[bestCircle][1]);
      Point2f cupleft;
      Point2f cupright;
      vector<size_t> cIds;
      for (auto &i : cupBlobs)
      {
	      cIds.push_back(i.first);
      }
      if(cupBlobs.size() > 0)
      {
	      std::cout << "DBG: pcentroid,center:" << pcentroid << " | " << center <<  endl;
	      if(cupBlobs[cIds[0]]->centroid.x<cupBlobs[cIds[1]]->centroid.x)
	      {
		      std::cout << "DBG: 0 is left"<<  endl;
		      cupleft.x=cupBlobs[cIds[0]]->centroid.x;
		      cupleft.y=cupBlobs[cIds[0]]->centroid.y;
		      cupright.x=cupBlobs[cIds[1]]->centroid.x;
		      cupright.y=cupBlobs[cIds[1]]->centroid.y;
	      }
	      else
	      {
		      std::cout << "DBG: 1 is left"<<  endl;
		      cupleft.x=cupBlobs[cIds[1]]->centroid.x;
		      cupleft.y=cupBlobs[cIds[1]]->centroid.y;
		      cupright.x=cupBlobs[cIds[0]]->centroid.x;
		      cupright.y=cupBlobs[cIds[0]]->centroid.y;
	      }
	      std::cout << "DBG: cupleft, cupright: " << cupleft << " | " << cupright << endl;

	      if(p2fdist(pcentroid,cupleft) < p2fdist(pcentroid,cupright))
	      {
		      //left side
		      std::cout << "DBG: left side" << endl;
		      if(p2fdist(pcentroid,cupleft)<
				      (circles[bestCircle][2]-p2fdist(pcentroid,center)))
		      {
			      std::cout << "DBG: left cup" << endl;
			      //left cup
			      OUTVAL=-1;
		      }
		      else
		      {
			      std::cout << "DBG: ring" << endl;
			      //ring
			      OUTVAL=0;
		      }
	      }
	      else
	      {
		      std::cout << "DBG: right side" << endl;
		      //right side
		      if(p2fdist(pcentroid,cupright)<
				      (circles[bestCircle][2]-p2fdist(pcentroid,center)))
		      {
			      //right cup
			      OUTVAL=-2;
		      }
		      else
		      {
			      //ring
			      OUTVAL=0;
		      }
	      }
      }
      else
      {
	      OUTVAL=0;
      }
      BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: " << CURRENT_FRAME << ",NEW ," << OUTVAL
                               << "->" << cid<<"(" << new_blob.label <<")";
    }
    //If the current cvblobID matched more than one larva in the previous frame
    //we give the cvblobID a new larvaID and update the parent_blobs of the newID
    //and the children_nodes of each of the previous IDs
    else
    {
      size_t CLUSTER_ID=++LARVAE_COUNT;
      new_blob.n20=new_blob.label;
      new_blob.label=CLUSTER_ID;
      new_assigned[new_blob.n20]=new_blob.label;
      for(auto &pid:new_matches_pre[cid])
      {
        if(find(children_blobs[pid].begin(),
              children_blobs[pid].end(),CLUSTER_ID)
            ==children_blobs[pid].end())
        {
          children_blobs[pid].push_back(CLUSTER_ID);
        }
        if(find(parent_blobs[CLUSTER_ID].begin(),
              parent_blobs[CLUSTER_ID].end(),pid)
            ==parent_blobs[CLUSTER_ID].end())
        {
          parent_blobs[CLUSTER_ID].push_back(pid);
        }
      }
      BOOST_LOG_TRIVIAL(debug) <<  "CNTTRACK: " << CURRENT_FRAME << ",N->1,"
        << printVector(new_matches_pre[cid]) << "->" << cid<<"("
                                   << new_blob.label <<")";
      certain_blobs.push_back(new_blob.label);
      for(auto &num:new_matches_pre[cid])
        BOOST_LOG_TRIVIAL(debug) <<  "DOT: " << num << "->" << new_blob.label << ";";
    }
  }
  //Now from the point of view of the previous frame
  for(auto &pre_blob_p:Prev)
  {
    size_t pid=pre_blob_p.second->label;
    int OUTVAL;
    if(pre_matched_by_new[pid].size()==0)
    {
      std::cout << "DBG: Larva Lost." << endl;
      //Larva is lost so we add its id to the lost_blobs vector
      Point2f pcentroid;
      pcentroid.x=pre_blob_p.second->centroid.x;
      pcentroid.y=pre_blob_p.second->centroid.y;
      Point2f center(circles[bestCircle][0],circles[bestCircle][1]);
      Point2f cupleft;
      Point2f cupright;
      vector<size_t> cIds;
      for (auto &i : cupBlobs)
      {
	      cIds.push_back(i.first);
      }
      if(cupBlobs.size() > 0)
      {
	      std::cout << "DBG: pcentroid,center:" << pcentroid << " | " << center <<  endl;
	      if(cupBlobs[cIds[0]]->centroid.x<cupBlobs[cIds[1]]->centroid.x)
	      {
		      std::cout << "DBG: 0 is left"<<  endl;
		      cupleft.x=cupBlobs[cIds[0]]->centroid.x;
		      cupleft.y=cupBlobs[cIds[0]]->centroid.y;
		      cupright.x=cupBlobs[cIds[1]]->centroid.x;
		      cupright.y=cupBlobs[cIds[1]]->centroid.y;
	      }
	      else
	      {
		      std::cout << "DBG: 1 is left"<<  endl;
		      cupleft.x=cupBlobs[cIds[1]]->centroid.x;
		      cupleft.y=cupBlobs[cIds[1]]->centroid.y;
		      cupright.x=cupBlobs[cIds[0]]->centroid.x;
		      cupright.y=cupBlobs[cIds[0]]->centroid.y;
	      }
	      std::cout << "DBG: cupleft, cupright: " << cupleft << " | " << cupright << endl;

	      if(p2fdist(pcentroid,cupleft) < p2fdist(pcentroid,cupright))
	      {
		      //left side
		      std::cout << "DBG: left side" << endl;
		      if(p2fdist(pcentroid,cupleft)<
				      (circles[bestCircle][2]-p2fdist(pcentroid,center)))
		      {
			      std::cout << "DBG: left cup" << endl;
			      //left cup
			      OUTVAL=-1;
		      }
		      else
		      {
			      std::cout << "DBG: ring" << endl;
			      //ring
			      OUTVAL=0;
		      }
	      }
	      else
	      {
		      std::cout << "DBG: right side" << endl;
		      //right side
		      if(p2fdist(pcentroid,cupright)<
				      (circles[bestCircle][2]-p2fdist(pcentroid,center)))
		      {
			      //right cup
			      OUTVAL=-2;
		      }
		      else
		      {
			      //ring
			      OUTVAL=0;
		      }
	      }
      }
      else
      {
	      OUTVAL=0;
      }
      lost_blobs[CURRENT_FRAME].push_back(pid);
      BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: " << CURRENT_FRAME
                               << ",LOST," << pid << "->" << OUTVAL ;
    }
    else if(pre_matched_by_new[pid].size()==1)
    {
      size_t cid=pre_matched_by_new[pid][0];
      if(new_assigned[cid]==0)
      {
        BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: " << CURRENT_FRAME
          << ",WARN," << cid
          << " unassigned and " << pid
          << " is pointing to it";
      }
      if(new_matches_pre[cid].size()==0)
      {
        BOOST_LOG_TRIVIAL(debug) <<  "CNTTRACK: " << CURRENT_FRAME
          << ",WARN," << pid << ">>" << cid
          << " but " << cid << " >> 0";
      }
    }
    //If the current preblobID was matched by more than one larva in the new
    //frame, then
    else if(pre_matched_by_new[pid].size()>1)
    {
      stringstream sstm;
      sstm <<  "CNTTRACK: "  << CURRENT_FRAME
        << ",1->N,"
        << pid << "->"
        << printVector(pre_matched_by_new[pid]) << "(";
      for(auto &nid:pre_matched_by_new[pid])
      {
        if(new_assigned[nid]==0)
        {
          cvb::CvBlob &new_blob=*In[nid];
          size_t CLUSTER_ID=++LARVAE_COUNT;
          new_blob.n20=new_blob.label;
          new_blob.label=CLUSTER_ID;
          if(find(parent_blobs[CLUSTER_ID].begin(),
                  parent_blobs[CLUSTER_ID].end(),pid)
              ==parent_blobs[CLUSTER_ID].end())
          {
            parent_blobs[CLUSTER_ID].push_back(pid);
          }
          if(find(children_blobs[pid].begin(),
                  children_blobs[pid].end(),CLUSTER_ID)
              ==children_blobs[pid].end())
          {
            children_blobs[pid].push_back(CLUSTER_ID);
          }
          new_assigned[new_blob.n20]=new_blob.label;
          sstm << new_blob.label;
          BOOST_LOG_TRIVIAL(debug) <<  "DOT: " << pid << "->" << new_blob.label << ";";
        }
        else
        {
          if(find(parent_blobs[new_assigned[nid]].begin(),
                  parent_blobs[new_assigned[nid]].end(),pid)
              ==parent_blobs[new_assigned[nid]].end())
          {
            parent_blobs[new_assigned[nid]].push_back(pid);
          }
          if(find(children_blobs[pid].begin(),
                  children_blobs[pid].end(),new_assigned[nid])
              ==children_blobs[pid].end())
          {
          children_blobs[pid].push_back(new_assigned[nid]);
          }
          sstm << new_assigned[nid];
          BOOST_LOG_TRIVIAL(debug) <<  "DOT: " << pid << "->" << new_assigned[nid] << ";";
        }
        sstm << " ";
      }
      BOOST_LOG_TRIVIAL(debug) << sstm.str() << ")";
    }
  }

  for(auto &cvBlobEntry: new_assigned)
  {
    size_t cvBlobID=cvBlobEntry.first;
    size_t larvaID=cvBlobEntry.second;
    if(larvaID==0)
    {
      BOOST_LOG_TRIVIAL(debug) << "CNTTRACK: " << CURRENT_FRAME
                               << ",ERR ,No assignment for cvblobID "
                               << cvBlobID;
    }
    out[larvaID]=In[cvBlobID];
  }
  updateLarvae(out,Prev);
}

/* Function to extract and process each frame so that the background is dark
 * and the forground white. Working with UMat
 *   Input:
 *       * capture: the capture device
 *   Output:
 *       * output: the processed frame
 *       * origFrame: the RGB version of the image
 */
bool u_get_next_frame(VideoCapture &capture, UMat &output, UMat &colorOutput,size_t step=1)
{
  if(TOTAL_FRAMES<CURRENT_FRAME+step && TOTAL_FRAMES!=0)
  {
    previousOrigFrame=Mat();
    return false;
  }
  if(step!=1)
    for(int i=0;i<(int) step-1;i++)
    {
      bool res=capture.grab();
      if(!res)
      {
        previousOrigFrame=Mat();
        return res;
      }
    }
  CURRENT_FRAME+=step;

  bool retval=capture.read(unprocessedFrame);
  unprocessedFrame.copyTo(output);
  if(retval==false)
  {
    previousOrigFrame=Mat();
    return retval;
  }

  if(step==1)
  {
    if( !previousOrigFrame.empty())
    {
      UMat dif;
      absdiff(output,previousOrigFrame,dif);
      //cout << "NORM: " << norm(dif) << endl;
      if(norm(dif)>(output.cols*output.rows*0.095))
      {
        //Frame corrupted, we keep the previous
        previousOrigFrame.copyTo(output);
      }
      else
        output.copyTo(previousOrigFrame);
    }
    else{
      output.copyTo(previousOrigFrame);
    }
  }


  UMat xchange;
  if(output.channels()==3)
  {
    cvtColor(output,output,CV_BGR2GRAY);
  }

  if(LRVTRACK_INVERT==true)
  {
    UMat ctout;
    //double MF=1.3;
    //size_t PF=0;
    int THRESHOLD=255;
    //double cpow=1.25;
    //output=output*MF+PF;
    addWeighted(output,1.0,output,-2.0,THRESHOLD,output);
    output.convertTo(output,CV_32F);
    //pow(output,cpow,output);
    convertScaleAbs(output,output,1,0);
  }
  /*brightnessContrastGamma(//output,
                          output,
                          LRVTRACK_BRIGHTNESS,
                          LRVTRACK_CONTRAST,
                          LRVTRACK_GAMMA);*/

  cvtColor(output,colorOutput,CV_GRAY2BGR);

  //imshow("OUTPUT",output);
  //waitKey(1);
  return retval;
//#endif
}


/* Function to extract and process each frame so that the background is dark
 * and the forground white.
 *   Input:
 *       * capture: the capture device
 *   Output:
 *       * output: the processed frame
 *       * origFrame: the RGB version of the image
 */
bool get_next_frame(VideoCapture &capture, Mat &output, Mat &colorOutput,size_t step=1)
{
  if(TOTAL_FRAMES<CURRENT_FRAME+step && TOTAL_FRAMES!=0)
  {
    previousOrigFrame=Mat();
    return false;
  }
  if(step!=1)
    for(int i=0;i<(int) step-1;i++)
    {
      bool res=capture.grab();
      if(!res)
      {
        previousOrigFrame=Mat();
        return res;
      }
    }
  CURRENT_FRAME+=step;

  bool retval=capture.read(unprocessedFrame);
  unprocessedFrame.copyTo(output);
  if(retval==false)
  {
    previousOrigFrame=Mat();
    return retval;
  }

  if(step==1)
  {
    if( !previousOrigFrame.empty())
    {
      Mat dif;
      absdiff(output,previousOrigFrame,dif);
      //cout << "NORM: " << norm(dif) << endl;
      if(norm(dif)>(output.cols*output.rows*0.095))
      {
        //Frame corrupted, we keep the previous
        previousOrigFrame.copyTo(output);
      }
      else
        output.copyTo(previousOrigFrame);
    }
    else{
      output.copyTo(previousOrigFrame);
    }
  }


  Mat xchange;
  if(output.channels()==3)
  {
    cvtColor(output,output,CV_BGR2GRAY);
  }

  if(LRVTRACK_INVERT==true)
  {
    Mat ctout;
    //double MF=1.3;
    //size_t PF=0;
    int THRESHOLD=255;
    //double cpow=1.25;
    //output=output*MF+PF;
    addWeighted(output,1.0,output,-2.0,THRESHOLD,output);
    output.convertTo(output,CV_32F);
    //pow(output,cpow,output);
    convertScaleAbs(output,output,1,0);
  }
  brightnessContrastGamma(//output,
                          output,
                          LRVTRACK_BRIGHTNESS,
                          LRVTRACK_CONTRAST,
                          LRVTRACK_GAMMA);
  Scalar meanS = mean(output);
  double meanVal = meanS.val[0];
  int nRows = output.rows;
  int nCols = output.cols;
  uchar *p;
  for ( int i=0; i< nRows; ++i)
  {
	  p = output.ptr<uchar>(i);
	  for(int j=0; j<nCols; ++j)
	  {
		  if(p[j]<meanVal)
			  p[j]=meanVal;
	  }
  }
  normalize(output,output,0,255,CV_MINMAX);
  cvtColor(output,colorOutput,CV_GRAY2BGR);

  //imshow("OUTPUT",output);
  //waitKey(0);
  return retval;
//#endif
}

/*
 * Function to take care of the various input possibilities.
 * Set up the parameters in the capture device for direct camera
 * input/file input and correct FPS.
 */
int setup_capture_input(VideoCapture &capture)
{
  // Check whether we have camera input or file
  if(LRVTRACK_CAMERA_INPUT != -2)
    {
      //capture.open(CV_CAP_DC1394);
      // These settings are for our setup.
      // TODO: Autodetection is easy for windows but still need to look into this.
      capture.open(300);
      //capture.set(CV_CAP_PROP_FRAME_WIDTH,1280);
      //capture.set(CV_CAP_PROP_FRAME_HEIGHT,1024);
      capture.set(CV_CAP_PROP_FPS,24);
    }
  else if (LRVTRACK_FILE_INPUT != "")
    {
      capture.open(LRVTRACK_FILE_INPUT);
    }
  capture.set(CV_CAP_PROP_FORMAT,CV_8U);
  TOTAL_FRAMES=capture.get(CV_CAP_PROP_FRAME_COUNT);
  if (capture.get(CV_CAP_PROP_FPS)==0)
    {
      VIDEO_FPS=24.1;
    }
  else
    {
      // Funny case where the FPS are mentioned to be double what they actually are.
      if ( capture.get(CV_CAP_PROP_FPS) > 30 )
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS)/2;
          LRVTRACK_FRAME_WIDTH=capture.get(CV_CAP_PROP_FRAME_WIDTH);
          LRVTRACK_FRAME_HEIGHT=capture.get(CV_CAP_PROP_FRAME_HEIGHT);
          //cout << VIDEO_FPS << endl;
        }
      else
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
          //cout << VIDEO_FPS << endl;
        }
    }

  if (!capture.isOpened())
    return -1;

  return 0;
}

void extract_background_offline(VideoCapture &capture,
                       Mat &greyBgFrame)
{
  Mat tmpFrame;
  Mat origFrame;
  Mat tmpFrame32f;
  size_t width=capture.get(CV_CAP_PROP_FRAME_WIDTH);
  size_t height=capture.get(CV_CAP_PROP_FRAME_HEIGHT);
  LRVTRACK_EXTRACT_OFFLINEBG_MIN=true;
  Mat resultframe;
  if(LRVTRACK_EXTRACT_OFFLINEBG_MIN)
    resultframe= Mat(height,width,CV_8UC1,Scalar(255));
  else
    resultframe = Mat::zeros(height,width,CV_32FC1);

  size_t count=0;//capture.get(CV_CAP_PROP_FRAME_COUNT);
  size_t total=capture.get(CV_CAP_PROP_FRAME_COUNT);
  cerr << "Background computation" << endl;
  while(get_next_frame(capture,tmpFrame,origFrame,10) && count < 100)
  {
    if(LRVTRACK_EXTRACT_OFFLINEBG_MIN)
      min(resultframe,tmpFrame,resultframe);
    else
    {
      tmpFrame.convertTo(tmpFrame32f,CV_32FC1);
      add(resultframe,tmpFrame32f,resultframe);
    }
    std::cerr << "B: " << CURRENT_FRAME << "/" << total << "\r";
    count++;
  }
  cerr << endl;
  CURRENT_FRAME=0;
  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error)
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  if(LRVTRACK_EXTRACT_OFFLINEBG_MIN)
    resultframe.convertTo(greyBgFrame,CV_8UC1);
  else
  {

    //resultframe*=(1.0/count);
    //resultframe.mul(resultframe,1.0/count);
    resultframe.convertTo(resultframe,resultframe.type(),1.0/count);
    normalize(resultframe,greyBgFrame,0,255,CV_MINMAX);
  }
  greyBgFrame.convertTo(greyBgFrame,CV_8UC1);

  // Find the petri-dish if not already defined from command line
  if(circles.size()==0)
  {
    //threshold(greyBgFrame,t,240,255,THRESH_BINARY|THRESH_OTSU);
    /*adaptiveThreshold(ctout,
      ctout,
      CMAX,
      METHOD,
      THRESH,
      NB,
      THRESHOLD);*/

    //addWeighted(greyBgFrame,0,greyBgFrame,2.0,0,ctout);
    Mat ctout;

    //int votes=270;
    double bgNorm=norm(greyBgFrame,NORM_L1)/(greyBgFrame.rows * greyBgFrame.cols);
    cout << "NORM: " << bgNorm << endl;
    //GaussianBlur(greyBgFrame, ctout, Size(0, 0), (2*(int)bgNorm), (2*(int)bgNorm));
    //bilateralFilter(greyBgFrame,ctout,128,128,128);
    //bilateralFilter(greyBgFrame,ctout,164,350,350);
    /*Mat element = getStructuringElement(MORPH_CROSS, Size(11, 11));
    dilate(ctout,ctout,element);
    dilate(ctout,ctout,element);
    dilate(ctout,ctout,element);*/
    //GaussianBlur(greyBgFrame, ctout, Size(0, 0), 32, 32 );
    //greyBgFrame.copyTo(ctout);
    //equalizeHist(ctout,ctout);
    greyBgFrame.copyTo(ctout);
    normalize(ctout,ctout,0,255,CV_MINMAX);
    //threshold(ctout,ctout,bgNorm-20,255,THRESH_BINARY_INV);
    adaptiveThreshold(ctout,ctout,255,ADAPTIVE_THRESH_GAUSSIAN_C,THRESH_BINARY_INV,355,0);
     // imshow("DISH",ctout);
     // waitKey(-1);
    //threshold(ctout,ctout,70,255,THRESH_OTSU);
    //bitwise_not(ctout, ctout);
      IplImage ipl= ctout;
      labelImg=cvCreateImage(
          cvGetSize(&ipl), IPL_DEPTH_LABEL, 1);
      cvLabel(&ipl, labelImg, dishBlob);
      Mat ctout_col;
      cvb::cvFilterByArea(dishBlob, 900400 ,6292529);
      int bestDishIdx=0;
      double bestRoundness=DBL_MAX;
      cout << "DishBlob Size: " << dishBlob.size() << endl;
      if(dishBlob.size() !=1 )
      {
        cout << dishBlob.size() << " dish blobs found!!!" << endl;
      }
      for(auto &dish:dishBlob)
      {
	double per = getPerimeter(*dish.second);
	double area = dish.second->area;
	double roundness = (per*per)/(2*CV_PI*area);
	if(roundness<bestRoundness)
          bestDishIdx=dish.first;
      }
      cerr << "We choose the best one : " << bestDishIdx << endl;
      cerr << "With Area: " << dishBlob[bestDishIdx]->area << endl;
      CvBlob *pdish=dishBlob[bestDishIdx];
      float cx=pdish->minx+((pdish->maxx-pdish->minx)/2);
      cerr << "cx: " << cx << endl;
      float cy=pdish->miny+((pdish->maxy-pdish->miny)/2);
      cerr << "cy: " << cy << endl;
      float radius=min((pdish->maxx-pdish->minx), (pdish->maxx-pdish->minx));
      radius = 0.48 * radius;
      if ( (int) radius % 2 == 1 )
        radius++;
      Vec3f petridish(cx,cy,radius);
      circles.push_back(petridish);
    //Loop until we find a set of circles we can work with :)
    //while (circles.size()==0 && votes >= 10)
    //{
      //HoughCircles(ctout, circles, CV_HOUGH_GRADIENT,
          //4,   // accumulator resolution (size of the image / 2)
          //10,  // minimum distance between two circles
          //100, // Canny high threshold
          //votes, // minimum number of votes
          //greyBgFrame.rows/3.0, greyBgFrame.rows/2.0); // min and max radiusV
    //  votes-=20;
    //}
    cout << "Petri dish circles size: "  << circles.size() << endl;
    cout << "X,Y,R: " << circles[0][0] << ", " << circles[0][1] << ", " << circles[0][2] << endl;

    cvtColor(ctout,ctout,CV_GRAY2BGR);
    for(size_t i=0; i<circles.size();i++)
    {
      circle(ctout,
          Point2f(circles[i][0],circles[i][1]),
          circles[i][2],
          Scalar(0,255,0),1);
    }

    // Once we have a set of circles we try to get the one that will give
    // us the best ratio brigthness/size. Assign its ID to bestCircle
    double cutlim;
    cutlim=DBL_MAX;
    //size_t mysz=circles.size();
    /*for(size_t i=0; i<circles.size();i++)
    {
      Mat cutout=Mat::zeros(ctout.rows,ctout.cols, ctout.type());
      circle(cutout,
          Point2f(circles[i][0],circles[i][1]),
          circles[i][2],
          Scalar(255),-1);
      cutout=cutout&greyBgFrame;
      double val=((double) sum(cutout)[0])/(double)(CV_PI*circles[i][2]*circles[i][2]);
      //val=(val*val)/(CV_PI*circles[i][2]*circles[i][2]);
      //(CV_PI*circles[i][2]*circles[i][2]);
      if(val<cutlim)
      {
        cutlim=val;
        bestCircle=i;
      }
    }*/
    bestCircle=0;
    /*Mat copyMat;
    greyBgFrame.copyTo(copyMat);
    cvtColor(copyMat,copyMat,CV_GRAY2BGR);
    circle(copyMat,
        Point2f(circles[0][0],circles[0][1]),
        circles[0][2],
        Scalar(0,255,0),1);

    imshow("BG with circle", copyMat);
    waitKey(1000000000);*/
  }
  else{
    bestCircle=0;
  }

  // Find the odour cups:
  if (LRVTRACK_ODOUR_CUPS>0 || LRVTRACK_ODOUR_CUPS==-1)
    {
      Mat fgROI;
      fgROI=Mat::zeros(greyBgFrame.rows ,
                           greyBgFrame.cols,
                           greyBgFrame.depth());

      if(circles.size()>0)
        {
          circle(fgROI,
                     Point2f(circles[bestCircle][0],circles[bestCircle][1]),
                     int(circles[bestCircle][2]),
                     Scalar(255),
                     -1);
        }
      else
        fgROI=Mat::ones(greyBgFrame.rows,
			greyBgFrame.cols,
			greyBgFrame.depth());

      Mat thr;
      //thr=greyBgFrame&fgROI;
      greyBgFrame.copyTo(thr);
      thr=fgROI&thr;
      morphologyEx(thr,thr,MORPH_OPEN,Mat(),Point2f(-1,-1),9);
      threshold(thr,thr,0,255,THRESH_OTSU+THRESH_BINARY);
      IplImage ipl= thr;
      labelImg=cvCreateImage(
          cvGetSize(&ipl), IPL_DEPTH_LABEL, 1);
      cvLabel(&ipl, labelImg, cupBlobs);
      double minArea=1000*(613/LRVTRACK_PETRIDISH);
      double maxArea=1200*(3090/LRVTRACK_PETRIDISH);
      cout << "CupBlobs Size before size filter: " << cupBlobs.size() << endl;
      cout << "minArea: " << minArea << endl;
      cout << "maxArea: " << maxArea << endl;
      for (auto &c: cupBlobs)
      {
	      cout << "C: " << c.second->area << endl;
      }
      cvb::cvFilterByArea(cupBlobs, minArea, maxArea);
      cout << "CupBlobs Size: " << cupBlobs.size() << endl;
      //cout << "Blob Cup1" << endl;
  //imshow("background",greyBgFrame);
  //waitKey(10000);
  }

  if(circles.empty() && cupBlobs.empty())
    return;

  cupContours=Mat(greyBgFrame.rows,greyBgFrame.cols, greyBgFrame.type(),Scalar(255));
  if(!circles.empty())
      circle(cupContours,
          Point2f(circles[bestCircle][0],circles[bestCircle][1]),
          int(circles[bestCircle][2]),
          Scalar(0),
          -1);

  cv::Point2f cc=Point2f(circles[bestCircle][0],
      circles[bestCircle][1]);

  double ppm;
  if(LRVTRACK_MPP==0.0)
  {
    ppm=LRVTRACK_PETRIDISH/(2*circles[bestCircle][2]);
    LRVTRACK_MPP=ppm;
  }
  else
    ppm=LRVTRACK_MPP;

  for(auto &b:cupBlobs)
  {
    if(LRVTRACK_ODOR_LR == "left")
    {
      // BUG: This was hard to find but horrible!
      if(b.second->centroid.x < greyBgFrame.cols/2)
	writeIni(ppm*(b.second->centroid.x-cc.x),ppm*(cc.y-b.second->centroid.y));
	//writeIni(ppm*(b.second->centroid.x-cc.x),ppm*(b.second->centroid.y-cc.y));
    }
    if(LRVTRACK_ODOR_LR == "right")
    {
      // BUG: This was hard to find but horrible!
      if(b.second->centroid.x > greyBgFrame.cols/2)
        writeIni(ppm*(b.second->centroid.x-cc.x),ppm*(cc.y-b.second->centroid.y));
        //writeIni(ppm*(b.second->centroid.x-cc.x),ppm*(b.second->centroid.y-cc.y));
    }
    createBlobContour(cupContours,
        *b.second,
        CV_8UC3,0,true,
        Scalar(255),8);
  }
  //cupContoursWhitePix=countNonZero(cupContours);
}

/*
 * Function to extract the background
 * The function retrieves a frame from the capture device and
 *  returns the background into the bgFrame image.
 *    Input:
 *     * capture: the video capture device
 *    Output:
 *     * greyBgFrame: the Suggested background
 *
 */
void extract_background(VideoCapture &capture,
                       Mat &greyBgFrame)
{
//#ifdef LRVTRACK_WITH_OPENCL
//#else
  Mat origFrame;
  // Grab a frame from the stream
  if(!get_next_frame(capture,greyBgFrame,origFrame))
  {
    //TODO: Error handling
    exit(1);
  }
  greyBgFrame.copyTo(origFrame);
  int votes=140; //For the HoughCircles call
  //int THRESH=THRESH_BINARY;

  //Mat ctout=greyBgFrame;
  Mat ctout;
  Mat t;
  //size_t CMAX=255;
  //size_t NB=275;
  //int THRESHOLD=0;
  //int METHOD=ADAPTIVE_THRESH_MEAN_C;
  // Thresholding here helps circle detection (trial and error)
  // TODO: Check if it's true.
  if(circles.size()==0)
  {
    //threshold(greyBgFrame,t,240,255,THRESH_BINARY|THRESH_OTSU);
    /*adaptiveThreshold(ctout,
      ctout,
      CMAX,
      METHOD,
      THRESH,
      NB,
      THRESHOLD);*/

    //addWeighted(greyBgFrame,0,greyBgFrame,2.0,0,ctout);
    //GaussianBlur(ctout, ctout, Size(15, 15), 4, 4 );
    GaussianBlur(greyBgFrame, ctout, Size(27, 27), 12, 12 );
    //Loop until we find a set of circles we can work with :)
    while (circles.size()==0 && votes >= 10)
    {
      HoughCircles(ctout, circles, CV_HOUGH_GRADIENT,
          8,   // accumulator resolution (size of the image / 2)
          300,  // minimum distance between two circles
          60, // Canny high threshold
          votes, // minimum number of votes
          greyBgFrame.rows/2.5, greyBgFrame.rows/1.95); // min and max radiusV
      votes-=20;
    }

    // Once we have a set of circles we try to get the one that will give
    // us the best ratio brigthness/size. Assign its ID to bestCircle
    double cutlim;
    cutlim=DBL_MAX;
    //size_t mysz=circles.size();
    /*for(size_t i=0; i<circles.size();i++)
    {
      Mat cutout=Mat::zeros(ctout.rows,ctout.cols, ctout.type());
      circle(cutout,
          Point2f(circles[i][0],circles[i][1]),
          circles[i][2],
          Scalar(255),-1);
      cutout=cutout&greyBgFrame;
      double val=((double) sum(cutout)[0])/(double)(CV_PI*circles[i][2]*circles[i][2]);
      //val=(val*val)/(CV_PI*circles[i][2]*circles[i][2]);
      //(CV_PI*circles[i][2]*circles[i][2]);
      if(val<cutlim)
      {
        cutlim=val;
        bestCircle=i;
      }
    }*/
    bestCircle=0;
  }
  else{
    bestCircle=0;
  }
  Mat circ;
  origFrame.copyTo(circ);
  cvtColor(circ,circ,CV_GRAY2BGR);
  circle(circ,Point2f(circles[bestCircle][0],circles[bestCircle][1]),
      int(circles[bestCircle][2]),
      Scalar(0,0,255),
      1);
  //imshow("Circles",circ);
  //waitKey(1);
  //int z=0;
  // Find the odour cups:
  if (LRVTRACK_ODOUR_CUPS>0 || LRVTRACK_ODOUR_CUPS==-1)
    {
      Mat fgROI;
      fgROI=Mat::zeros(greyBgFrame.rows ,
                           greyBgFrame.cols,
                           greyBgFrame.depth());

      if(circles.size()>0)
        {
          circle(fgROI,
                     Point2f(circles[bestCircle][0],circles[bestCircle][1]),
                     int(circles[bestCircle][2]),
                     Scalar(255),
                     -1);
        }

      Mat thr;
      //thr=greyBgFrame&fgROI;
      greyBgFrame.copyTo(thr);
      double thresholdlow=0;
      double thresholdhigh=0;
      morphologyEx(thr,thr,MORPH_OPEN,Mat(),Point2f(-1,-1),9);
      /*
      threshold(thr,
                    thr,
                    thresholdlow,
                    thresholdhigh,
                    THRESH_BINARY|THRESH_OTSU);
                    */
      HoughCircles(greyBgFrame, cups, CV_HOUGH_GRADIENT,
          2,   // accumulator resolution (size of the image / 2)
          500,  // minimum distance between two circles
          120, // Canny high threshold
          70, // minimum number of votes
          60, 90); // min and max radiusV
      /*for(auto &cup: cups)
      {

        cupROI=Mat::zeros(greyBgFrame.rows ,
                           greyBgFrame.cols,
                           greyBgFrame.depth());

      }*/
    }

  //Initialize a first background frame. Everything that is in the circle
  //is black (i.e. will be excluded from the background).
  //Outside the circle is white.
  if(circles.size()>0)
  {
    circle(greyBgFrame,
        Point2f(circles[bestCircle][0], circles[bestCircle][1]),
        circles[bestCircle][2],
        Scalar(),
        -1);
  }
  else
  {
    greyBgFrame=Mat(greyBgFrame.rows,
        greyBgFrame.cols,
        greyBgFrame.depth(),
        Scalar(0));
  }
  Mat tempFrame;
  Mat showFrame;
  //origFrame.copyTo(showFrame);
  process_frame(origFrame,
               greyBgFrame,
               showFrame,
               tempFrame);
  // We remove the larvae from the background
  bg_without_larvae(tempFrame);
  //imshow("Background",tempFrame);

  // We add the greyBgFrame with the tempFrame to get the full
  // background image
  add(greyBgFrame,tempFrame,greyBgFrame);
//#endif

 return;

}

void setROI(string &input)
{
  boost::char_separator<char> sep(", ");
  boost::tokenizer<boost::char_separator <char> > tokens(input,sep);
  size_t i=0;
  Vec3f b;
  for(const auto& t: tokens)
  {
    b[i]=strtof(t.c_str(),NULL);
    i++;
  }
  circles.push_back(b);
}

// Function to handle command-line arguments
int handle_args(int argc, char* argv[])
{
  char cDate[16];
  time_t curtime;
  struct tm *loctime;

  /* Get the current time. */
  curtime = time (NULL);

  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  strftime(cDate,16,"%Y%m%d_%H%M%S",loctime);

  string DATE=string(cDate);
  LRVTRACK_DATE=DATE;
  string VIDEO_NAME=DATE + VIDEO_TYPE;
  string PROCESSED_VIDEO_NAME=DATE + "_PROC"+VIDEO_TYPE;

  try
    {
      // allowed only on command line
      po::options_description generic("Generic options");
      generic.add_options()
      ("version", "Print version string")
      ("help,h", "Produce help message")
      ("list-cameras,l","List available cameras.")
      ("verbose,v",
       po::value<int>(&LRVTRACK_VERBOSE_LEVEL),
       "Verbosity (0-3).")
      ("dstep,d",
       po::value<int>(&LRVTRACK_DSTEP),
       "Dstep.")
      ("wstep,w",
       po::value<double>(&LRVTRACK_WSTEP),
       "Wstep.")
      ;

      po::options_description IOopts("Input Output Options");

      IOopts.add_options()

      ("results-folder,r",
       po::value<string>(&LRVTRACK_RESULTS_FOLDER)->implicit_value("."),
       "Main folder where the results will be recorded. The experiment folder for each experiment (based on the date and time) will be created under this folder.(Default: Current folder)"
      )

      ("file-suffix,f",
       po::value<string>(&LRVTRACK_NAME)->implicit_value(""),
       "Suffix to use for the experiment output files under the experiment folder. The experiment folder name will be DATE-<file-suffix>.(Default: empty)"
      )

      ("save-video,s",
       po::value<string>
       (&LRVTRACK_SAVE_VIDEO)->implicit_value(DATE+VIDEO_TYPE),
       "Save original input as video. \
                         The option is disabled when the input is from a video \
                         source.(Default: <date-filesuffix>)"
      )

      ("save-tracked-video,t",
       po::value<string>
       (&LRVTRACK_SAVE_PROCESSED_VIDEO)->implicit_value(
         DATE+"_PROC"+VIDEO_TYPE),
       "Save resulting output as video with the filename \
                         specified.(Default: <date-filesuffix>_PROC)"
      )
      ("file-input,i",
       po::value<string>(&LRVTRACK_FILE_INPUT)->implicit_value(""),
       "Filename of video to be used as input for tracker."
      )

      ("metadata-file",
       po::value<string>(&LRVTRACK_INPUT_METADATA)->implicit_value(""),
       "Filename of metadata for video to be used as input for tracker."
      )

      ("camera-input,c",
       po::value<int>(&LRVTRACK_CAMERA_INPUT)->implicit_value(0),
       "Camera feed to be used as input for tracker. For now unfortunately \
                         one can only use numbers and no listing is possible."
      )

      ("csv-output,u",
       "Output will be csv file named as \"date_time@<suffix>.csv\". The file will be \
       saved under the experiment folder."
      )

      ("choreography-output,j",
       "Output summary and blob files to be processed by choreography."
      );

      po::options_description setupOpts("Setup Options");
      setupOpts.add_options()

      ("parallel,p",
       "When the flag is specified lrvTrack will try to make use of multiple cores."
      )
      ("thread-count",
       po::value<size_t> (&LRVTRACK_THREADS)->implicit_value(4),
       "Number of threads to use for parallel computation."
      )

      ("invert,e",
       "When the flag is specified lrvTrack will assume that the larvae are the dark objects in the image."
      )

      ("roi",
       po::value<string>(&LRVTRACK_ROI_INPUT)->implicit_value(""),
       "Set a circle as Region of Interest. Input: x,y,radius"
      )

      ("mpp",
       po::value<double>(&LRVTRACK_MPP)->implicit_value(0.0),
       "Milimeters per pixel. If not set will be calculated from region of interest and size of dish (see -z)."
      )


      ("gamma,g",
       po::value<double> (&LRVTRACK_GAMMA)->implicit_value(1.0),
       "Gamma Correction value."
      )
      ("brightness,b",
       po::value<double> (&LRVTRACK_BRIGHTNESS)->implicit_value(0.0),
       "Brightness value."
      )
      ("contrast",
       po::value<double> (&LRVTRACK_CONTRAST)->implicit_value(1.0),
       "Contrast value."
      )
      ("min-obj-size",
       po::value<size_t> (&LRVTRACK_MIN_OBJ_SIZE)->implicit_value(20),
       "Minimum size of object to consider (in pixels)."
      )
      ("max-obj-size",
       po::value<size_t> (&LRVTRACK_MAX_OBJ_SIZE)->implicit_value(1200),
       "Maximum size of object to consider (in pixels)."
      )

      ("petri-dish-size,z",
       po::value<size_t> (&LRVTRACK_PETRIDISH)->implicit_value(90),
       "Size of petri-dish (actually the circle enclosed by the ROI) in mm."
      )

      ("odour-cups,o",
       po::value<int> (&LRVTRACK_ODOUR_CUPS)->implicit_value(10),
       "When the flag is specified lrvTrack will look for odour cups. If specified\
                         with no arguments (e.g. -o) lrvTrack will look for as many as it can find.\
                         The number of odour cups can be specified as a value for this option."
      );

      po::options_description procOpts("Processing Options");
      procOpts.add_options()

      ("smoothing,m",
       po::value<double> (&LRVTRACK_SMOOTHING)->implicit_value(2.0),
       "Smoothing value for the contour reconstruction."
      )

      ("use-model,k",
       "Use model to resolve short collisions."
      )
      ("model-duration,",
       po::value<double> (&LRVTRACK_MODEL_DURATION)->implicit_value(3.5),
       "Duration of collisions on which to use the model in seconds."
      )

      ("normalize,n",
       "Setting this flag will disable normalization of the brightness values of \
                         each frame. This degrades the accuracy of results but reduces \
                         cpu utilization. (Default: Not set)"
      );

      po::options_description displayOpts("Display Options");
      displayOpts.add_options()

      ("show-liveview,V",
       "Show analysis live (Default: False)"
      )
      ("extract_offline,x",
       "Extract background from the whole video (Default: False)"
      )

      ("extract_offline_min",
       "Extract minimum background from the whole video (Default: False)"
      )

      ("show-skeleton,S",
       po::value<bool> (&LRVTRACK_SHOW_SKELETON),
       "Show skeleton for detected larvae (Default: False)"
      )

      ("show-contour,C",
       po::value<bool> (&LRVTRACK_SHOW_CONTOUR),
       "Show contour for detected larvae (Default: False)"
      )

      ("show-orientation,O",
       po::value<bool> (&LRVTRACK_SHOW_ORIENTATION),
       "Show orientation for detected larvae (Default: False)"
      )

      ("show-centroid,Z",
       po::value<bool> (&LRVTRACK_SHOW_CENTROID),
       "Show centroid for detected larvae (Default: True)"
      )

      ("show-head-tail,H",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show head and tail for detected larvae (Default: True)"
      )

      ("show-tags,T",
       po::value<bool> (&LRVTRACK_SHOW_TAGS),
       "Show larvae tags (Default: True)"
      )

      ("batch-testing,B",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show no output except from the matchings. Used for testing internal parameters."
      );

      po::options_description cmdline_options;
      cmdline_options.add(generic).add(IOopts).add(setupOpts).add(procOpts).add(displayOpts);

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
      po::notify(vm);

      if (vm.count("help"))
        {
          cout << cmdline_options << "\n";
          exit(1);
        }
      if (vm.count("version"))
        {
          cout << LRVTRACK_VERSION_MAJOR << "." << LRVTRACK_VERSION_MINOR << "\n";
          exit(1);
        }
      if (vm.count("results-folder"))
        {
          LRVTRACK_RESULTS_FOLDER=vm["results-folder"].as<string>();
        }
      else
        {
          LRVTRACK_RESULTS_FOLDER=".";
        }

      if(vm.count("csv-output"))
        {
          LRVTRACK_CSV_OUTPUT=true;
        }
      if(vm.count("choreography-output"))
        {
          LRVTRACK_CHOREOGRAPHY_OUTPUT=true;
        }
      if(vm.count("use-model"))
        {
          LRVTRACK_USE_MODEL=true;
        }
      if(vm.count("normalize"))
        {
          LRVTRACK_NORMALIZE=true;
        }
      if(vm.count("odour-cups"))
        {
          LRVTRACK_ODOUR_CUPS=vm["odour-cups"].as<int>();
        }

      if (vm.count("file-suffix"))
        {
          LRVTRACK_NAME=LRVTRACK_NAME;
        }
      else
        {
          LRVTRACK_NAME=DATE;
        }

      if (vm.count("save-video"))
        {
          LRVTRACK_SAVE_VIDEO=LRVTRACK_NAME+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_VIDEO="";
        }

      if (vm.count("save-tracked-video"))
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO=LRVTRACK_NAME+"_PROC"+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO="";
        }

      if (vm.count("list-cameras"))
        {
          BOOST_LOG_TRIVIAL(error) << "Camera listing not yet implemented." ;
          exit(1);
        }

      if (vm.count("file-input")<=0 && vm.count("camera-input")<=0)
        {
          BOOST_LOG_TRIVIAL(error) << "Error: No Input Specified. Exiting..." ;
          exit(2);
        }
      else if (vm.count("file-input")>0 && vm.count("camera-input")>0)
        {
          BOOST_LOG_TRIVIAL(error)
            << "Error: Ambiguous Input Specified. Exiting...";
          exit(2);
        }

      if (vm.count("file-input")>0)
        {
          LRVTRACK_FILE_INPUT=vm["file-input"].as<string>();
          if (LRVTRACK_FILE_INPUT=="")
            {
              BOOST_LOG_TRIVIAL(error)
              << "Error: Input file flag given but no file specified. Exiting...";
              exit(3);
            }
        }
      if (vm.count("metadata-file")>0)
        {
          LRVTRACK_INPUT_METADATA=vm["metadata-file"].as<string>();
          if (LRVTRACK_INPUT_METADATA=="")
            {
              cout << "METADATA file from: " << LRVTRACK_FILE_INPUT << endl;
              boost::filesystem::path p(LRVTRACK_FILE_INPUT);
              boost::filesystem::path dir = p.parent_path();
              if(dir.string() == "")
                LRVTRACK_INPUT_METADATA="./metadata.txt" ;
              else
                LRVTRACK_INPUT_METADATA=dir.string() + "/metadata.txt" ;
              if(boost::filesystem::exists(LRVTRACK_INPUT_METADATA) )
              {
                cout << "METADATA file found: " << LRVTRACK_INPUT_METADATA << endl;
                readIni();
              }
              else
              {
                LRVTRACK_INPUT_METADATA="";
                cout << "Could not fine metadata file" << endl;
              }
            }
        }
      if (vm.count("roi")>0)
        {
          LRVTRACK_ROI_INPUT=vm["roi"].as<string>();
          setROI(LRVTRACK_ROI_INPUT);
        }
      if (vm.count("mpp")==0)
        {
          LRVTRACK_MPP=0.0;
        }

      if (vm.count("camera-input")>0)
        {
          LRVTRACK_CAMERA_INPUT=vm["camera-input"].as<int>();
        }
      else
        {
          LRVTRACK_CAMERA_INPUT=-2;
        }
      if(vm.count("parallel")>0)
      {
        LRVTRACK_PARALLEL=true;
      }
      else
        LRVTRACK_PARALLEL=false;
      if(vm.count("show-liveview")>0)
      {
        LRVTRACK_SHOW_LIVE=true;
      }
      if(vm.count("extract_offline")>0)
      {
        LRVTRACK_EXTRACT_OFFLINEBG=true;
      }
      if(vm.count("extract_offline_min")>0)
      {
        LRVTRACK_EXTRACT_OFFLINEBG_MIN=true;
      }
      if(vm.count("invert")>0)
      {
        LRVTRACK_INVERT=true;
      }
    }
  catch(po::error &e)
    {
      BOOST_LOG_TRIVIAL(error) << "Problem parsing options: " << e.what();
      exit(1);
    }

  return 0;
}

/* Second Pass
 *
 * The functions below belong to the second pass of data
 *  - Verify larvae and non-larvae objects
 *  - Detect heads/tails
 *
*/

size_t evalPartitionSize(std::vector<std::vector<size_t> > &v,
                          std::vector<size_t> &blobs)
{
  //size_t vsz=v.size();
  //size_t bsz=blobs.size();
  size_t totalError=0;
  size_t b=0;
  for(auto &i:v)
  {
    size_t expected_size=0;
    for(auto &l:i)
      expected_size+=detected_larvae[l].area.back();
    totalError+=abs(detected_larvae[blobs[b++]].area.back()-
                    expected_size);
  }
  return totalError;
}

void findEndObjectsFwHelper(larvaObject &t,
                            size_t c,
                            size_t objnum,
                            std::map<size_t,size_t> &blobSizeMap,
                            std::vector<size_t> &e)
{
  if(e.size()==objnum)
    return;
  BOOST_LOG_TRIVIAL(debug) << "Verify: FindFW: Main Obj: " <<
      t.larva_ID << " CurObj: " << c << " Objnum: " << objnum;
  if(blobSizeMap[c]==1)
  {
    BOOST_LOG_TRIVIAL(debug) <<
      "Verify: FindFW: blobSize 1 and not Cluster so we push: " << c;
    e.push_back(c);
  }
  else
  {
    if(!children_blobs[c].empty())
    {
      if(e.size()==objnum)
        return;
      BOOST_LOG_TRIVIAL(debug) <<
        "Verify: FindFW: children not empty: " << c;
      BOOST_LOG_TRIVIAL(debug) << "Verify: FindFW: "
        << printVector(children_blobs[c]);
      for(auto &n: children_blobs[c])
      {
        if(e.size()==c)
          break;
        findEndObjectsFwHelper(t,n,objnum,blobSizeMap,e);
      }
    }
  }
}

/*void findEndObjectsBwHelper(larvaObject &t,
                            size_t c,
                            size_t objnum,
                            std::map<size_t,size_t> &blobSizeMap,
                            std::vector<size_t> &e)
{
  if(e.size()==objnum)
    return;
  BOOST_LOG_TRIVIAL(debug) << "Verify: FindBW: Main Obj: " <<
      t.larva_ID << " CurObj: " << c << " Objnum: " << objnum;
  if(blobSizeMap[c]==1)
  {
    BOOST_LOG_TRIVIAL(debug) <<
      "Verify: FindBW: blobSize 1 and not Cluster so we push: " << c;
    e.push_back(c);
  }
  else
  {
    if(!parent_blobs[c].empty())
    {
      BOOST_LOG_TRIVIAL(debug) <<
        "Verify: FindBW: parents not empty: " << c;
      BOOST_LOG_TRIVIAL(debug) << "Verify: FindBW: "
        << printVector(parent_blobs[c]);
      for(auto &n: parent_blobs[c])
      {
        findEndObjectsBwHelper(t,n,objnum,blobSizeMap,e);
      }
    }
  }
}*/

/*void findEndObjects(larvaObject &t,
                    std::vector<size_t> &e,
                    std::map<size_t,size_t> &blobSizeMap,
                    size_t obj_num,
                    bool forward)
{
  if(forward)
    findEndObjectsFwHelper(t,t.larva_ID,obj_num,blobSizeMap,e);
  else
    findEndObjectsBwHelper(t,t.larva_ID,obj_num,blobSizeMap,e);

  if(e.size()!=(size_t) t.blobSize)
  {
    BOOST_LOG_TRIVIAL(debug) << "WARNING: Cannot find end Objects "
      << "Larva: " << t.larva_ID << " with size " << t.blobSize;
  }
}*/

void clarifyLarvae(larvaObject &l,
                   std::vector<size_t> &neighbors,
                   size_t neighborsize,
                   std::vector<size_t> &allfound,
                   std::map<size_t,size_t> &blobSizeMap,
                   bool &change)
{
    if(neighborsize>allfound.size())
    {
      BOOST_LOG_TRIVIAL(debug) << "Clarify: For Larva: " <<
        l.larva_ID << " WARN Nsize=" <<
        neighbors.size() << ", ALL=" << allfound.size();
      return;
    }
    partition_generator<size_t> pgs(allfound,neighbors.size());
    pgs.m_partitions_of_n();
    vector<vector <size_t> > &BestPartition=pgs.RES[0];
    size_t minError=SIZE_MAX;
    for(auto &pg: pgs.RES)
    {
      size_t ne;
      ne=evalPartitionSize(pg,neighbors);
      if(ne<minError)
      {
        minError=ne;
        BestPartition=pg;
        change=1;
      }
    }
    size_t newsz=0;
    for(auto i=0;i<(int)neighbors.size();i++)
    {
      size_t sz=0;
      for(auto &s:BestPartition[i])
        sz+=detected_larvae[s].blobSize;

      newsz+=sz;
      detected_larvae[neighbors[i]].blobSize=sz;
      blobSizeMap[neighbors[i]]=sz;
      detected_larvae[neighbors[i]].isCluster=(sz-1);
      BOOST_LOG_TRIVIAL(debug) << "Verify: FIX: OBJECT"
        << neighbors[i]
        << " updated with blobSize: " << sz;
      change=true;
    }
    //l.blobSize=newsz;
    //if(newsz>1)
    //  l.isCluster=true;
    //BOOST_LOG_TRIVIAL(debug) << "Verify: FIX: OBJECT"
    //  << l.larva_ID
    //  << " updated with blobSize: " << newsz;
}

/*void verifyLarva(larvaObject &f,bool &change)
{
  if(!f.diverged_to.empty())
  {
    size_t sum=0;
    for(auto &v: f.diverged_to)
      sum+=detected_larvae[v].blobSize;
    BOOST_LOG_TRIVIAL(debug) << "Verify: Checking: CLUSTER"
      << f.larva_ID
      << " has blobSize: " << f.blobSize
      << " and the sum of those diverging from it is: " << sum;
    if(f.blobSize!=sum)
    {
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: CLUSTER"
        << f.larva_ID
        << " has blobSize: " << f.blobSize
        << " and the sum of those diverging from it is: " << sum;
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: From: "
        << printVector(f.diverged_to);
      std::vector<size_t> endObjs;
      f.blobSize=std::max(sum,(size_t) f.blobSize);

      if(f.blobSize>1)
        f.isCluster=true;
      findEndObjects(f,endObjs,1);
      if(endObjs.size()!=f.blobSize)
      {
        endObjs.clear();
        findEndObjects(f,endObjs,0);
      }
      if(endObjs.size()!=f.blobSize)
      {
        BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: for " << f.larva_ID
          << " we cannot find enough end objects." ;
        return;
      }

      BOOST_LOG_TRIVIAL(debug) << "Verify: Clarifying: " << f.larva_ID
        << " with end objects: "<< printVector(endObjs)
        << " for ObjectSize: " << f.blobSize;
      clarifyLarvae(f,f.diverged_to,sum,endObjs,change);
    }
  }
  if(parent_blobs[f.larva_ID].size()>1)
  {
    size_t sum=0;
    for(auto &v: parent_blobs[f.larva_ID])
      sum+=detected_larvae[v].blobSize;
    BOOST_LOG_TRIVIAL(debug) << "Verify: Checking: CLUSTER"
      << f.larva_ID
      << " has blobSize: " << f.blobSize
      << " and the sum of those diverging from it is: " << sum;
    if(f.blobSize!=sum)
    {
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: CLUSTER"
        << f.larva_ID
        << " has blobSize: " << f.blobSize
        << " and the sum of those diverging from it is: " << sum;
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: From: "
        << printVector(parent_blobs[f.larva_ID]);
      std::vector<size_t> endObjs;
      f.blobSize=std::max(sum,(size_t) f.blobSize);

      if(f.blobSize>1)
        f.isCluster=true;
      findEndObjects(f,endObjs,0);
      if(endObjs.size()!=f.blobSize)
      {
        endObjs.clear();
        findEndObjects(f,endObjs,0);
      }
      if(endObjs.size()!=f.blobSize)
      {
        BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: for " << f.larva_ID
          << " we cannot find enough end objects." ;
        return;
      }

      BOOST_LOG_TRIVIAL(debug) << "Verify: Clarifying: " << f.larva_ID
        << " with end objects: "<< printVector(endObjs)
        << " for ObjectSize: " << f.blobSize;
      clarifyLarvae(f,parent_blobs[f.larva_ID],sum,endObjs,change);
    }
  }
}*/
void standardize1d(std::vector<double> &v,std::vector<double> &out,double scale=1.0)
{
    std::vector<double> m,s,tmp;
    meanStdDev(v,m,s);

    out.resize(v.size());
    for(int i=0;i<out.size();i++)
        out[i]=scale*((v[i]-m[0])/s[0]);

}

void determineHeadTailBreaks(larvaObject &f)
{
    //stringstream d;
    //Construct function scaled(RND) + scaled(length) - scaled(width) of larva f
    if(f.length.size() < 5)
        return;

    // TS: transform to standard normal distribution
    std::vector<double> l_scaled, w_scaled, r_scaled, p_scaled;
    standardize1d(f.length, l_scaled);
    standardize1d(f.width, w_scaled, -1.0);
    standardize1d(f.roundness, r_scaled);
    standardize1d(f.perimeter, p_scaled);

    std::vector<double> c, c_std;
    c.resize(f.perimeter.size());
    for(int i=0;i<c.size();i++)
    {
      // FRAGE: wie kommt diese formel zustande? lookup online?
        c[i]=l_scaled[i]+w_scaled[i]+r_scaled[i]+p_scaled[i];
    }
    //standardize the result
    standardize1d(c,c_std);
    //Smoothing with filter of 3
    c[0]=(c_std[0]+c_std[1])/2;
    c[c.size()-1] = (c_std[c_std.size()-1]+c_std[c_std.size()-2])/2.0;
    for (int i=1; i<c.size()-1;i++)
    {
        c[i]=(c_std[i-1]+c_std[i]+c_std[i+1])/3.0;
        //BOOST_LOG_TRIVIAL(debug) << "BREAKCHECK: " << f.larva_ID << "," << c[i];

    }
    //update round_flag vector.
    for (int i=0; i<c.size();i++)
    {
      if(c[i]<-3.0)// FRAGE: warum -3? empirically
            f.round_flag[i]=true;
    }

}

void determineHeadTail(larvaObject &f)
{
  stringstream d;
  size_t F=24;
  size_t sframe=f.start_frame;
  if (f.tails.size() == 0 )
  {
      	d << "|=============== Blob " << f.larva_ID << " has 0 info for tail: turning to cluster ========";
	  f.isCluster=true;
	  return;
  }
  if (f.heads.size() == 0 )
  {
      	d << "|=============== Blob " << f.larva_ID << " has 0 info for head: turning to cluster ========";
	  f.isCluster=true;
	  return;
  }
  try
  {
    determineHeadTailBreaks(f);
    d << endl;
    d << "|=============== Determine head/tail of: " << f.larva_ID
      << " ============= " << endl;
    if(f.isCluster==true)
    {
      d << "|=============== Blob " << f.larva_ID << " is a cluster ========"
        << endl;
      d << "|========================================= " << endl;
      BOOST_LOG_TRIVIAL(debug) << d.str();
      return;
    }
    size_t range=f.blobs.size();
    std::vector<size_t> breaks;
    bool INBREAK=false;
    // FRAGE: TS: abschnitte von pausen detektieren, flanken finden?
    for(size_t i=0;i<range;i++)
    {
      if(f.round_flag[i]==true && INBREAK==false)
      {
        breaks.push_back(i);
        INBREAK=true;
      }
      else if(f.round_flag[i]==false && INBREAK==true)
      {
        if(f.round_flag.size()>=i && f.round_flag[i+1]==true)
          continue;
        breaks.push_back(i);
        INBREAK=false;
      }
    }

    breaks.push_back(range-1);
    d << "|  Range: " << range << endl;
    d << "|  Breaks: " << breaks.size() << " " << printVector(breaks) << endl;

    for(size_t i=0;i<breaks.size();i++)
    {
      d << "|    Break " << i << ":" << endl;
      size_t endBreak=breaks[i];
      size_t preBreak;
      if(i>0)
        preBreak=breaks[i-1]+1;
      else
        preBreak=0;

      Point2f startPos(f.blobs[preBreak].centroid.x,
          f.blobs[preBreak].centroid.y);
      Point2f endPos(f.blobs[endBreak].centroid.x,
          f.blobs[endBreak].centroid.y);

      Point2f startbp(f.blobs[preBreak].minx,f.blobs[preBreak].miny);
      Point2f endbp(f.blobs[endBreak].minx,f.blobs[endBreak].miny);

      std::cerr << "Prebreak: " << preBreak << " Tail size: " << f.tails.size() << std::endl;
      d << "|      Start[ " << preBreak+sframe <<"]: Centroid" << startPos << " "
        << "IHead" << f.heads[preBreak]+startbp << " "
        << "ITail" << f.tails[preBreak]+startbp
        << endl;
      d << "|      End[ " << endBreak+sframe <<"]: Centroid" << endPos << " "
        << "IHead" << f.heads[endBreak]+endbp << " "
        << "ITail" << f.tails[endBreak]+endbp
        << endl;

      //double sumHeadDist=0.0; //Distance to endpoint
      //double sumTailDist=0.0;
      double sumHeadGrey=0.0; //Grey value sum
      double sumTailGrey=0.0;
      double sumHeadDiff=0.0; //Difference of current centroid
      //from previous endpoints
      double sumTailDiff=0.0;
      double sumHeadCnt=0.0;  //Distance covered by endpoint
      double sumTailCnt=0.0;

      for(size_t j=preBreak;j<=endBreak;j++)
      {
        Point2f bpPoint(f.blobs[j].minx,f.blobs[j].miny);
        Point2f hPoint=filterPoint(f.heads,j,F,f.blobs);
        Point2f tPoint=filterPoint(f.tails,j,F,f.blobs);
	//if(isfinite(f.heads_brightness[j]) &&
	//   isfinite(f.tails_brightness[j]))
	if(f.heads_brightness[j]>0 && f.heads_brightness[j] < 9999 &&
	   f.tails_brightness[j]>0 && f.tails_brightness[j] < 9999)
	{
		sumHeadGrey+=f.heads_brightness[j];
		sumTailGrey+=f.tails_brightness[j];
	}
        if(j>preBreak+1)
        {
          Point2f prebpPoint(f.blobs[j-2].minx,f.blobs[j-2].miny);
          Point2f prehPoint=filterPoint(f.heads,j-2,F,f.blobs);
          Point2f pretPoint=filterPoint(f.tails,j-2,F,f.blobs);
          //Point2f centroid=filterCentroid(f.blobs,j,F);
          Point2f midpoint=filterMidpoint(f.lrvDistances,j,F);
	  double tmpd1,tmpd2,tmpd3,tmpd4;
	  tmpd1 = p2fdist(midpoint,prehPoint);
	  if(tmpd1>0 & tmpd1<100)
          	sumHeadDiff+=tmpd1;
	  tmpd2 = p2fdist(midpoint,pretPoint);
	  if(tmpd2>0 & tmpd2<100)
          	sumTailDiff+=tmpd2;
	  tmpd3 = p2fdist(hPoint,prehPoint);
	  if(tmpd3>0 & tmpd3<100)
          	sumHeadCnt+=tmpd3;
	  tmpd4 = p2fdist(tPoint,pretPoint);
	  if(tmpd4>0 & tmpd4<100)
          	sumTailCnt+=tmpd4;

// FRAGE: warum 100? Cnt=Count? aber eigentlich summierte distanz die kopf zuruecklegt
// lookup unit of p2fdist, maybe pixels
	  /*BOOST_LOG_TRIVIAL(debug) << "HTCHECK: " << f.larva_ID << "," << j <<
		  "," << f.heads_brightness[j] <<
		  "," << f.tails_brightness[j] <<
		  "," << sumHeadGrey <<
		  "," << sumTailGrey <<
		  "," << tmpd1 <<
		  "," << tmpd2 <<
		  "," << tmpd3 <<
		  "," << tmpd4;*/
        }
	else
	{
	  /*BOOST_LOG_TRIVIAL(debug) << "HTCHECK: " << f.larva_ID << "," << j <<
	  	"," << f.heads_brightness[j] <<
	  	"," << f.tails_brightness[j] <<
		"," << sumHeadGrey <<
		"," << sumTailGrey <<
		"," << 0 <<
		"," << 0 <<
		"," << 0 <<
		"," << 0 ;*/
	}
	if(!isfinite(sumHeadGrey) || !isfinite(sumTailGrey))
	{
		BOOST_LOG_TRIVIAL(debug) << "VALUE GETS TO NaN here!!!" << endl;
	}
      }
      double vote_for_nochange=0;
      size_t duration=endBreak-preBreak+1;

      if(sumHeadDiff/duration<sumTailDiff/duration) //If head was infront most of the time vote for nochange
	vote_for_nochange++;
        //vote_for_nochange+=fabs(sumHeadDiff-sumTailDiff)/(sumHeadDiff+sumTailDiff);
      else
	vote_for_nochange--;
        //vote_for_nochange-=fabs(sumHeadDiff-sumTailDiff)/(sumHeadDiff+sumTailDiff);

      if(sumHeadGrey!=0 && sumTailGrey!=0 && isfinite(sumHeadGrey) && isfinite(sumTailGrey)
		      && duration>1 )
      {
        if(sumHeadGrey/duration<sumTailGrey/duration) //If head was darker vote for no change
	  vote_for_nochange++;
          //vote_for_nochange+=fabs(sumHeadGrey-sumTailGrey)/(sumHeadGrey+sumTailGrey);
        else
	  vote_for_nochange--;
          //vote_for_nochange-=fabs(sumHeadGrey-sumTailGrey)/(sumHeadGrey+sumTailGrey);
      }

      if(sumHeadCnt!=0 && sumTailCnt!=0 && isfinite(sumHeadCnt) && isfinite(sumTailCnt)
		      && duration>1 )
      {
        if(sumHeadCnt/duration>sumTailCnt/duration) //If head made more distance vote for no change
	  vote_for_nochange++;
          //vote_for_nochange+=fabs(sumHeadCnt-sumTailCnt)/(sumHeadCnt+sumTailCnt);
        else
	  vote_for_nochange--;
          //vote_for_nochange-=fabs(sumHeadCnt-sumTailCnt)/(sumHeadCnt+sumTailCnt);
      }

      d << "|" << endl;
      d << "|      LarvaID, GHead, HeadForward, HeadDistance, Duration, Votes" << endl;
      d << "|      " <<f.larva_ID*1000 + i << "0, ";
      d << sumHeadGrey/duration<< ", ";
      d << sumHeadDiff/duration<< ", ";
      d << sumHeadCnt/duration<< ", " << endBreak-preBreak <<", ";
      d << (vote_for_nochange>=0) << ",0" << endl;
      d << "|      LarvaID, GTail, TailForward, TailDistance, Duration, Votes" << endl;
      d << "|      " <<f.larva_ID*1000 + i << "1, ";
      d << sumTailGrey/duration<< ", ";
      d << sumTailDiff/duration<< ", ";
      d << sumTailCnt/duration<< ", " << endBreak-preBreak <<", ";
      d << (vote_for_nochange<0) << ",1" << endl;
      d << "|" << endl;

      if(vote_for_nochange<0 && breaks.size()==1)
      {
        f.heads.swap(f.tails);
      }
      else if(vote_for_nochange<0)
      {
        for(size_t k=preBreak;k<=endBreak;k++)
        {
          cv::Point2f b=f.heads[k];
          f.heads[k]=f.tails[k];
          f.tails[k]=b;
        }
      }
      d << "|      Start[ " << preBreak <<"]: Centroid" << startPos << " "
        << "RHead" << f.heads[preBreak]+startbp << " "
        << "RTail" << f.tails[preBreak]+startbp
        << endl;
      d << "|      End[ " << endBreak <<"]: Centroid" << endPos << " "
        << "RHead" << f.heads[endBreak]+endbp << " "
        << "RTail" << f.tails[endBreak]+endbp
        << endl;
      d << "|========================================= " << endl;
    }
    BOOST_LOG_TRIVIAL(debug) << d.str();
  }
  catch(...)
  {
    BOOST_LOG_TRIVIAL(debug) << d.str();
    exit(0);
  }
}

#define blobSum(l,g) accumulate(l.begin(), \
                              l.end(),0,\
                              [&g](size_t total, size_t cur)\
                              {return g[cur]+total;} );

/*void updateBlobSize(larvaObject &l,
                 map<size_t,size_t> &blobSizeMap,
                 bool &change)
{
  size_t lID=l.larva_ID;
  size_t parents=blobSum(parent_blobs[lID],blobSizeMap);
  size_t children=blobSum(children_blobs[lID],blobSizeMap);
  int size_as_child=0;
  int size_as_parent=0;

  if(parents==0 && children==0)
  {
    if(blobSizeMap[lID]==0)
    {
      change=true;
      blobSizeMap[lID]=1;
    }
    return;
  }

  //Node as a child
  if(parents!=0)
  {
    //One parent and being the onlychild
    // Weird case that perhaps doesn't occur.
    if(parent_blobs[lID].size()==1 &&
       children_blobs[parent_blobs[lID].front()].size()==1)
    {
      size_as_child=blobSizeMap[parent_blobs[lID].front()];
    }
    else if(parent_blobs[lID].size()==1 &&
       children_blobs[parent_blobs[lID].front()].size()>1)
    {
      [>size_as_child=blobSizeMap[parent_blobs[lID].front()];
      for(auto &pc: children_blobs[parent_blobs[lID].front()])
      {
        if(pc!=lID)
          size_as_child-=blobSizeMap[pc];
      }<]
      size_as_child=0;
    }
    else if(parent_blobs[lID].size()>1)
    {
      bool onlychild=true;
      for(auto &pp: parent_blobs[lID])
      {
        if(children_blobs[pp].size()>1)
        {
          onlychild=false;
          break;
        }
        if(children_blobs[pp].front()!=lID)
        {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning one parent but different child";
          onlychild=false;
          break;
        }
      }
      if(onlychild)
      {
        size_as_child=blobSum(parent_blobs[lID],blobSizeMap);
      }
      else
      {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning two parents with foster siblings. Is unhandled";
      }
    }
  }

  //Node as parent
  if(children!=0)
  {
    //One child and being the onlyparent
    // Weird case that perhaps doesn't occur.
    if(children_blobs[lID].size()==1 &&
       parent_blobs[children_blobs[lID].front()].size()==1)
    {
      size_as_parent=blobSizeMap[children_blobs[lID].front()];
    }
    else if(children_blobs[lID].size()==1 &&
       parent_blobs[children_blobs[lID].front()].size()>1)
    {
      [>size_as_parent=blobSizeMap[children_blobs[lID].front()];
      for(auto &pc: parent_blobs[children_blobs[lID].front()])
      {
        if(children_blobs[pc].size()>1)
        {
          size_as_parent=0;//A bit of a mess :)
          break;
        }
        //if(pc!=lID) //Also unhandled
         // size_as_parent-=blobSizeMap[pc];
      }<]
      size_as_parent=0;
    }
    else if(children_blobs[lID].size()>1)
    {
      bool singleparent=true;
      for(auto &pp: children_blobs[lID])
      {
        if(parent_blobs[pp].size()>1)
        {
          singleparent=false;
          break;
        }
        if(parent_blobs[pp].front()!=lID)
        {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning one child but different parent";
          singleparent=false;
          break;
        }
      }
      if(singleparent)
      {
        size_as_parent=blobSum(children_blobs[lID],blobSizeMap);
      }
      else
      {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning two children with other parents. Is unhandled";
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug) << "Verify: L"
    << lID
    << " with blobSize: " << blobSizeMap[lID]
    << " parents: " << printVector(parent_blobs[lID])
    << " as parent: " << size_as_parent
    << " children: " << printVector(children_blobs[lID])
    << " as child: " << size_as_child;

  if(size_as_parent>size_as_child && size_as_child>0)
  {
    BOOST_LOG_TRIVIAL(debug) << "Verify: L"
      << lID
      << " parents>children clarifying";
    vector<size_t> allFwBlobs;
    findEndObjects(l,allFwBlobs,blobSizeMap,size_as_parent,1);
    BOOST_LOG_TRIVIAL(debug) << "Verify: L"
      << lID
      << " found end Objects:" << printVector(allFwBlobs);
    clarifyLarvae(l,
        parent_blobs[lID],
        parents,
        allFwBlobs,
        blobSizeMap,
        change);
    if(blobSizeMap[lID]<(size_t) size_as_parent)
    {
      blobSizeMap[lID]=size_as_parent;
      change=true;
    }
  }

  else if(size_as_parent<size_as_child && size_as_parent>0)
  {
    BOOST_LOG_TRIVIAL(debug) << "Verify: L"
      << lID
      << " parents<children clarifying";
    vector<size_t> allBwBlobs;
    findEndObjects(l,allBwBlobs,blobSizeMap,size_as_child,0);
    BOOST_LOG_TRIVIAL(debug) << "Verify: L"
      << lID
      << " found end Objects:" << printVector(allBwBlobs);
    clarifyLarvae(l,
        children_blobs[lID],
        children,
        allBwBlobs,
        blobSizeMap,
        change);
    if(blobSizeMap[lID]<(size_t) size_as_child)
    {
      blobSizeMap[lID]=size_as_child;
      change=true;
    }
  }
  else
  {
    if(max(size_as_parent,size_as_child)>(int) blobSizeMap[lID])
    {
      blobSizeMap[lID]=max(size_as_parent,size_as_child);
      change=true;
    }
  }
}*/

bool modelHeadTailRepair(larvaObject &l, lrvFit &f)
{
  Point2f bp(l.blobs[0].minx,l.blobs[0].miny);
  double hh=p2fdist(f.spine[0],l.heads[0]+bp);
  double ht=p2fdist(f.spine[0],l.tails[0]+bp);
  double tt=p2fdist(f.spine.back(),l.tails[0]+bp);
  double th=p2fdist(f.spine.back(),l.heads[0]+bp);
  bool switched=false;
  if(hh+tt>ht+th)
    switched=true;

  if(switched)
  {
    //Find longest consistent track and
    //change the shortest to match that

    //Length of track before model
    larvaObject &pl=(reincarnations[f.ID].size()!=0) ?
      detected_larvae[reincarnations[f.ID].back()] :
      detected_larvae[f.ID];

    //TODO check the case where the Head/Tail are switched
    //within model
    size_t prelength=0;
    auto preItHead=pl.heads.end();
    auto preItTail=pl.tails.end();
    for(auto b=pl.round_flag.rbegin();
        b!=pl.round_flag.rend();b++)
    {
      preItHead--;
      preItTail--;
      if(!*b)
        prelength++;
      else
        break;
    }

    size_t postlength=0;
    auto postItHead=l.heads.begin();
    for(auto b=l.round_flag.begin();b!=l.round_flag.end();b++)
    {
      postItHead++;
      if(!*b)
        postlength++;
      else
        break;
    }

    if(postlength>prelength && prelength < 2*VIDEO_FPS)
    {
      BOOST_LOG_TRIVIAL(debug)
        << "Head Tail doesn't match. Switching based on Post! Post Larva: " << l.larva_ID
        << " model of: " << f.ID << " EndModel: " << l.start_frame
        << endl;
      std::reverse(f.spine.begin(),f.spine.end());
      swap_ranges(preItHead,pl.heads.end(),preItTail);
    }
    if(prelength>postlength  && postlength < 2*VIDEO_FPS)
    {
      BOOST_LOG_TRIVIAL(debug)
        << "Head Tail doesn't match. Switching based on Pre! Post Larva: " << l.larva_ID
        << " model of: " << f.ID << " EndModel: " << l.start_frame
        << endl;
      swap_ranges(l.heads.begin(),postItHead,l.tails.begin());
    }

  }
  return true;
}

void reAssignID()
{
  //for(auto &p:detected_larvae)
  //{
    //if cluster
      //get output objects of cluster
        //for all larvae output of cluster
          //if larvae is already assigned (has an updated id != id ) then
            //check if the original larva with larva_id = updated_id has an updated_id
            //and replace value
          //else
            //See if paths are long enough
            //get contents of the cluster
  //}
}
// FRAGE: wie funktioniert das alles?
void collisionSearch()
{
  int cid=0;
  size_t COLLISIONS_RAW=0;
  size_t COLLIDED_SINGLE_LARVA_IN_OBJECTS=0;
  size_t COLLIDED_SINGLE_LARVA_OUT_OBJECTS=-1;
  //size_t COLLIDED_SINGLE_LARVAE_OUT=0;
  size_t ASSIGNABLE_SINGLE_LARVAE=0;
  size_t RESOLVED_BY_MODEL=0;
  size_t PARTIAL_RESOLVED_BY_MODEL=0;
  size_t ASSIGNMENTS_BY_MODEL=0;
int blobid = 0;

  for(auto &p:detected_larvae)
  {
    blobid++;

    if(!p.second.isCluster)
        ASSIGNABLE_SINGLE_LARVAE++;

    if(p.second.isCluster && parent_blobs[p.first].size()>0)
    {
      if(parent_blobs[p.first].size()>1)
        COLLISIONS_RAW++;
      if(any_of(parent_blobs[p.first].begin(),
                parent_blobs[p.first].end(),
                [](size_t i){return detected_larvae[i].isCluster==false;}))
        COLLIDED_SINGLE_LARVA_IN_OBJECTS++;

      // This part below is a mistake. It counts the
      // 	blob->blob->single larva first appearance
      // as well
      // The best way to calculate the possible assignments is to calculate
      // all single detected larvae minus the total number of larvae on the plate
      /*size_t newlarvaeout=count_if(children_blobs[p.first].begin(),
                children_blobs[p.first].end(),
                [](size_t i){return detected_larvae[i].isCluster==false;});
        COLLIDED_SINGLE_LARVAE_OUT+=newlarvaeout;

        if(newlarvaeout>0)
          COLLIDED_SINGLE_LARVA_OUT_OBJECTS++;*/

      //exclude cases with more than 2 in/out larvae
      if(parent_blobs[p.first].size()==1)
        std::cout << "Parentblobs of: " << p.first << " has size 1!!!" << endl;
      if(children_blobs[p.first].size()==1)
        std::cout << "Childrenblobs of: " << p.first << " has size 1!!!" << endl;

      if(parent_blobs[p.first].size()!=2 || children_blobs[p.first].size()!=2)
        continue;
      //exclude those collisions between clusters, we don't do this
      if(all_of(parent_blobs[p.first].begin(),
      //if(any_of(parent_blobs[p.first].begin(),
                parent_blobs[p.first].end(),
                [](size_t i){return detected_larvae[i].isCluster==true;})
        )
      {
        continue;
      }
      larvaObject &l=p.second;
      //exclude long collisions
      if(l.lastFrameWithStats-l.start_frame > LRVTRACK_MODEL_DURATION*VIDEO_FPS)
        continue;
      //Eliminate those cases with wrong initial shape
      if(any_of(parent_blobs[p.first].begin(),
                parent_blobs[p.first].end(),
                [](size_t i){return (detected_larvae[i].round_flag.size()>0 &&
                                    detected_larvae[i].round_flag.back())==true;})
        )
          continue;
      //Try to model the rest
      stringstream n;
      //n << "modelOutput" << cid++ << ".avi";
      vector<size_t> parent_larvae(parent_blobs[p.first].size());
      auto it = std::copy_if (parent_blobs[p.first].begin(),
                              parent_blobs[p.first].end(),
                              parent_larvae.begin(),
                              [](size_t i){
                              return !detected_larvae[i].isCluster;
                              } );
      parent_larvae.resize(std::distance(parent_larvae.begin(),it));
      larvaeModels.emplace_back(parent_larvae,
      //larvaeModels.emplace_back(parent_blobs[p.first],
                       l.blobs[0],
                       l.start_frame,
                       l.lastFrameWithStats);
      collisionModel &C=larvaeModels.back();
      size_t modelFrame;
      VideoWriter modelOut;
      char fname[1024];


      printf("colsearch blobid: %d  is cluser of two!!\n", blobid);

      for(modelFrame=l.start_frame+1;modelFrame<=l.end_frame;modelFrame++)
      {
        cerr << "Model: Frame: " << cid << "(" << l.start_frame << ")"
          << " F:" << modelFrame-l.start_frame << "/"
          << l.end_frame-l.start_frame ;
        if(modelFrame==l.end_frame)
          cerr << " Error: " <<  accumulate(C.frameError.begin(),C.frameError.end(),0.0);
        cerr  << "\r";
        Mat ret;
        // TS: fit vector version of larva contours onto pixel contours
        C.updateModel(l.blobs[modelFrame-l.start_frame],modelFrame,ret);

        // parameter space:
        //    whole larvae moves within 3x3 pixel neighborbood translation
        //    bending of head by angles (discretized)
        //imshow("Model",ret);
        //waitKey(1);
/*        copyMakeBorder(ret, ret, 15, 15 , 15, 15,BORDER_CONSTANT , Scalar(255,255,255));
        cvtColor(ret,ret,CV_GRAY2BGR);
        if(!modelOut.isOpened())
        {
          modelOut.open(n.str(),
              VIDEO_CODEC,
              VIDEO_FPS,
              ret.size());
          if (!modelOut.isOpened())
          {
            std::cerr << "Error opening video output file. " <<
              "Problems with video codec. Exiting..." << std::endl;
            exit(5);
          }
        }
        modelOut << ret;*/
      }
//      modelOut.release();
      //Temporary assignment for testing
      std::cerr << endl;
      if(C.larvae_models.size()==2)
      {
        lrvFit &fit1 = C.larvae_models[0].back();
        lrvFit &fit2 = C.larvae_models[1].back();
      	if( accumulate(C.frameError.begin(),C.frameError.end(),0.0)
			>
			C.frameError.size() * 8.0  * (detected_larvae[fit1.ID].area_sum
			 /
			 detected_larvae[fit1.ID].area.size()))
	{
		cerr << "Skipping this assignment. Error too big" << endl;
	  	continue;
	}

        CvPoint2D64f cl1=detected_larvae[children_blobs[p.first][0]].blobs[0].centroid;
        CvPoint2D64f cl2=detected_larvae[children_blobs[p.first][1]].blobs[0].centroid;
        double da1=p2fdist(fit1.spine[5],cl1); // 5 == centriod
        double da2=p2fdist(fit2.spine[5],cl2);
        double db1=p2fdist(fit1.spine[5],cl2);
        double db2=p2fdist(fit2.spine[5],cl1);

        size_t replacingID1=detected_larvae[fit1.ID].updated_ID;
        size_t replacingID2=detected_larvae[fit2.ID].updated_ID;
        size_t replacedID1=children_blobs[p.first][0];
        size_t replacedID2=children_blobs[p.first][1];

        if(da1+da2<db1+db2)
        {

          reincarnations[replacingID1].push_back(replacedID1);
          reincarnations[replacingID2].push_back(replacedID1);

          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][0]],fit1);
          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][1]],fit2);

          detected_larvae[replacedID1].updated_ID=replacingID1;
          detected_larvae[replacedID2].updated_ID=replacingID2;
	  cout << "MODEL: ID: " << replacedID1 << ">" << replacingID1 << endl;
	  cout << "MODEL: ID: " << replacedID2 << ">" << replacingID2 << endl;
        }
        else
        {
          reincarnations[replacingID2].push_back(replacedID1);
          reincarnations[replacingID1].push_back(replacedID2);
          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][1]],fit1);
          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][0]],fit2);
          detected_larvae[replacedID1].updated_ID=replacingID2;
          detected_larvae[replacedID2].updated_ID=replacingID1;
	  cout << "MODEL: ID: " << replacedID1 << ">" << replacingID2 << endl;
	  cout << "MODEL: ID: " << replacedID2 << ">" << replacingID1 << endl;
        }

        RESOLVED_BY_MODEL++;
        ASSIGNMENTS_BY_MODEL+=2;
        C.SUCCESS=true;
      }
      if(C.larvae_models.size()==1) // TS: identify shortly lost (touching cluster, or cups? ) single larva
      {

        lrvFit &fit1 = C.larvae_models[0].back();
      	if( accumulate(C.frameError.begin(),C.frameError.end(),0.0)
			>
			C.frameError.size() * 4.0  * (detected_larvae[fit1.ID].area_sum
			 /
			 detected_larvae[fit1.ID].area.size()))
	{
		cerr << "Skipping this assignment. Error too big" << endl;
	  	continue;
	}
        CvPoint2D64f cl1=detected_larvae[children_blobs[p.first][0]].blobs[0].centroid;
        CvPoint2D64f cl2=detected_larvae[children_blobs[p.first][1]].blobs[0].centroid;
        double da1=p2fdist(fit1.spine[5],cl1);
        double db1=p2fdist(fit1.spine[5],cl2);

        size_t replacingID1=detected_larvae[fit1.ID].updated_ID;
        size_t replacedID1=children_blobs[p.first][0];
        size_t replacedID2=children_blobs[p.first][1];

        if(da1<db1 && !detected_larvae[children_blobs[p.first][0]].isCluster)
        {
          reincarnations[replacingID1].push_back(replacedID1);
          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][0]],fit1);
          detected_larvae[replacedID1].updated_ID=replacingID1;
	  cout << "MODEL: ID: " << replacedID1 << ">" << replacingID1 << endl;
          PARTIAL_RESOLVED_BY_MODEL++;
          ASSIGNMENTS_BY_MODEL++;
          C.SUCCESS=true;
        }
        else if ( da1>db1 && !detected_larvae[children_blobs[p.first][1]].isCluster)
        {
          reincarnations[replacingID1].push_back(replacedID2);
          //modelHeadTailRepair(detected_larvae[children_blobs[p.first][1]],fit1);
          detected_larvae[replacedID2].updated_ID=replacingID1;
	  cout << "MODEL: ID: " << replacedID2 << ">" << replacingID1 << endl;
          PARTIAL_RESOLVED_BY_MODEL++;
          ASSIGNMENTS_BY_MODEL++;
          C.SUCCESS=true;
        }
        else
        {
          C.SUCCESS=false;
        }
      }
      cid++;
      std::cerr << endl;
    }
  }
  ASSIGNABLE_SINGLE_LARVAE -= MAX_LARVAE_DETECTED;
  BOOST_LOG_TRIVIAL(debug) << "COLLISIONS RAW: " <<  COLLISIONS_RAW << endl
    << "COLLIDED_SINGLE_LARVA_IN_OBJECTS:" << COLLIDED_SINGLE_LARVA_IN_OBJECTS << endl
    << "COLLIDED_SINGLE_LARVA_OUT_OBJECTS:" << COLLIDED_SINGLE_LARVA_OUT_OBJECTS << endl
    << "ASSIGNABLE_SINGLE_LARVAE:" << ASSIGNABLE_SINGLE_LARVAE << endl
    << "RESOLVED_BY_MODEL:" << RESOLVED_BY_MODEL << endl
    << "PARTIAL_RESOLVED_BY_MODEL:" << PARTIAL_RESOLVED_BY_MODEL << endl
    << "ASSIGNMENTS_BY_MODEL:" << ASSIGNMENTS_BY_MODEL << endl;
}

/*bool checkForInconsistentNodes(map<size_t,size_t> &blobSizeMap,
    vector<size_t> &inconsistent_nodes)
{
}*/

void constructCandidates()
{
  map<size_t,vector<size_t> > candidates;
  for (auto &p: detected_larvae)
  {
    larvaObject &l=p.second;
    size_t lID=l.larva_ID;
    if(l.updated_ID!=l.larva_ID)
    {
      if(candidates[l.updated_ID].size()>0)
      {
        for(auto &c:candidates[l.updated_ID])
          candidates[lID].push_back(c);
      }
      else
      {
        candidates[lID].push_back(l.updated_ID);
      }
      continue;
    }

    //If we have parents
    //TODO: Two cases:
    //        - Parents were all in the ROI
    //        - Some of the parents were in the ROI
    if(parent_blobs[lID].size()!=0)
    {
    }
  }
}

void secondPass()
{
  /* We loop over each larva object and try to
   * identify head/tail and whether it's a larva or not */

  BOOST_LOG_TRIVIAL(debug) <<
    "Starting Second Pass: Larvae on Dish guess: " << MAX_LARVAE_DETECTED;

  if(detected_larvae.size()>10000)
    cerr << "More than 10000 objects detected... Quitting!" << endl;
  std::map<size_t,size_t> blobSizeMap;
  std::map<size_t, std::vector<size_t> > detailed_clusters;

  //Solve this using Simplex with guarantee that the minimum solution
  //is integral... TODO check this!!!
  //Mat mConstr;
  //Mat mf;
  //Mat res;
  vector<size_t> indexToID;
  map<size_t,size_t> IDtoIndex;
  size_t nidx=0;
  cout << "MPP: " << LRVTRACK_MPP << endl;
  for(auto &p: detected_larvae)
  {
    IDtoIndex[p.second.larva_ID]=nidx++;

    vector<double> clLenStDevArr;
    vector<double> clLenMeanArr;
    vector<double> clWidStDevArr;
    vector<double> clWidMeanArr;
    vector<double> clPerStDevArr;
    vector<double> clPerMeanArr;
    meanStdDev(p.second.length,clLenMeanArr,clLenStDevArr);
    meanStdDev(p.second.width,clWidMeanArr,clWidStDevArr);
    meanStdDev(p.second.perimeter,clPerMeanArr,clPerStDevArr);
    cout << "AVG LARVA SIZE: " << p.second.larva_ID << ", "
         << p.second.area_sum/p.second.area.size() << ", "
         << clLenStDevArr[0] << ", "
         << clWidStDevArr[0] << ", "
         << clPerStDevArr[0]
         << endl;
    //is_larva(p.second);
  }
  //LPconstr(mConstr,indexToID,IDtoIndex);
  //LPmaxFunction(mf);
  //cout << "MF: " << endl << mf << endl;
  //cout << "CONSTR: " << endl << mConstr << endl;
  //int LPout=solveLP(mf,mConstr,res);

  //res.convertTo(res,CV_8UC1);
  vector<double> solution;
  int LPout=LPconstr(solution,indexToID,IDtoIndex);
  if(LPout>=0)
    cout << "RES: OK ";// << res.t() << endl;
  else
  {
    cerr << "LPout: " << LPout << endl;
    exit(100);
  }
  cout << "Sizes: " ;
  size_t LAST_FRAME=0;
  for (auto &p: detected_larvae)
  {
    char sz=(size_t) (solution[IDtoIndex[p.second.larva_ID]]+0.5);
    size_t size_of_node=(size_t) sz;
    cout << size_of_node << ", ";
    if(size_of_node>1)
    {
      p.second.isCluster=1;
      p.second.blobSize=size_of_node;
    }
    if(p.second.end_frame>LAST_FRAME)
      LAST_FRAME=p.second.lastFrameWithStats;
  }
  cout << endl;
  size_t single_larvae_objects=0;
  double single_larvae_total=0.0;
  for(size_t f=1;f<=LAST_FRAME;f++)
  {
    size_t FSZ=0;
    for (auto &p: detected_larvae)
    {
      if(p.second.blobSize==1)
      {
        single_larvae_total+=p.second.area_mean;
        single_larvae_objects++;
      }
      if(p.second.start_frame < f &&
        p.second.end_frame > f)
        FSZ+=p.second.blobSize;
    }
    if(FSZ>MAX_LARVAE_DETECTED)
      MAX_LARVAE_DETECTED=FSZ;
  }
  double AVG_LARVA_SIZE=((double)single_larvae_total/(single_larvae_objects*1.02));
  BOOST_LOG_TRIVIAL(debug) <<
    "Second Pass: Updated Larvae on Dish guess: " << MAX_LARVAE_DETECTED;

  solution.clear();
  //LPout=LPconstr(solution,indexToID,IDtoIndex,AVG_LARVA_SIZE);
  LPout=LPconstrAvg(solution,indexToID,IDtoIndex,AVG_LARVA_SIZE);

  if(LPout>=0)
    cout << "RES: OK ";// << res.t() << endl;
  else
  {
    cerr << "LPout: " << LPout << endl;
    exit(100);
  }
  cout << "Sizes: " ;
  for (auto &p: detected_larvae)
  {
    char sz=(size_t) (solution[IDtoIndex[p.second.larva_ID]]+0.5);
    size_t size_of_node=(size_t) sz;
    cout << size_of_node << ", ";
    if(size_of_node>1)
    {
      p.second.isCluster=true;
      p.second.blobSize=size_of_node;
    }
    if(size_of_node<=1)
    {
      p.second.isCluster=false;
      p.second.blobSize=1;
    }
  }
  cout << endl;

  std::map<size_t,larvaObject>::iterator dlit;
  dlit=detected_larvae.begin();
  for(;dlit!=detected_larvae.end();++dlit)
  {
    determineHeadTail(dlit->second);
  }
  if(LRVTRACK_USE_MODEL)
    collisionSearch();
}

int main(int argc, char* argv[])
{
  bool SHOWTAGS=true;
  //bool TRACK=false;
  //bool STEP=true;
  bool STOP=false;
  VideoWriter vidOut,vidPOut;

  /*if (!combOut.isOpened())
  {
    combOut.open("combMov.avi",
        VIDEO_CODEC,
        16.0,
        Size(168,168));
    if (!combOut.isOpened())
    {
      std::cerr << "Error opening video output file. " <<
        "Problems with video codec. Exiting..." << std::endl;
      exit(5);
    }
  }*/

  log_init();

#ifdef LRVTRACK_WITH_CUDA
  gpu::printShortCudaDeviceInfo(gpu::getDevice());
//#elif defined(LRVTRACK_WITH_OPENCL)
#endif
  std::cout << "Using OPENCL: " <<  ocl::useOpenCL() << endl;
  if( ocl::useOpenCL())
  {
    ocl::Device oclD;
    oclD=oclD.getDefault();
    cout << "OpenCL device: " << oclD.name() << endl;
  }
//#endif
  // Initialize video capture device
  VideoCapture capture;
  // Handle command line arguments
  handle_args(argc,argv);
  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error)
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  if(LRVTRACK_EXTRACT_OFFLINEBG)
    extract_background_offline(capture,bgFrame);
  else
    extract_background(capture,bgFrame);


  //Initialize the frame size
  FRAME_COLS=bgFrame.cols;
  FRAME_ROWS=bgFrame.rows;

  cvb::CvBlobs preBlobs;
  //larvaFit lf;
  double totalFitError=0.0;
  //int TOGGLE=0;
  while (!STOP)
  {
	  frameCount++;
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    process_frame(greyFrame,bgFrame,colorFrame,processedFrame);
    //imshow("Orig",greyFrame);
    //imshow("Processed" ,processedFrame);
    //waitKey(10000000);
    //imshow("PF",processedFrame);
    /*if(CURRENT_FRAME==3)
      TOGGLE=1;
    char k=waitKey(1);
    if(TOGGLE==1)
    {
      while(k!=' ' && k!='.')
      {
        k=waitKey(1);
      }
    }

    if(k==' ' && TOGGLE==0)
      TOGGLE=1;
    if(k==' ' && TOGGLE==1)
      TOGGLE=0;
    if(k=='.')
      TOGGLE=1;*/

    IplImage ipl_thresholded = processedFrame;
    labelImg=cvCreateImage(
        cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);
    NEW.clear();	// blob structre
    cvLabel(&ipl_thresholded, labelImg, NEW);

//    IplImage cf = colorFrame;
//
    IplImage *tmp_P = cvCreateImage(cvGetSize(&ipl_thresholded),IPL_DEPTH_8U,3);
    cvZero(tmp_P);
 //   cvtColor(labelImg,tmp_P,CV_GRAY2BGR,1);
//    cvMerge(labelImg, NULL, NULL, tmp_P);
    cvRenderBlobs(labelImg, NEW, tmp_P, tmp_P);

    cvb::cvFilterByArea(NEW, LRVTRACK_MIN_OBJ_SIZE, LRVTRACK_MAX_OBJ_SIZE);

    cvZero(tmp_P);
    cvRenderBlobs(labelImg, NEW, tmp_P, tmp_P);

    cvb::CvBlobs tracked_blobs;
    cvb::CvBlobs blob1;
    cvb::CvBlobs::iterator blIT=NEW.begin();

    // dump all blobs in cups or outside dish
    for(;blIT!=NEW.end();)
    {
      if(larvaToRing(*blIT->second))
      {
        NEW.erase(blIT++);
      }
      else
        blIT++;
    }

    cvZero(tmp_P);
    cvRenderBlobs(labelImg, NEW, tmp_P, tmp_P);

    //cvb::CvBlobs::iterator it=NEW.begin();
    if(preBlobs.size()>0)
    {
      if(NEW.size()>MAX_LARVAE_DETECTED)
        MAX_LARVAE_DETECTED=NEW.size();
      //newLarvaeTrack(NEW,preBlobs,tracked_blobs);

// "connect" this previous frames blobs with current frames blobs for tracking 
      contourLarvaeTrack(NEW,preBlobs,tracked_blobs);

      if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
      {
        //printSummary(preBlobs,NEW,false);
      }
      preBlobs=tracked_blobs;
      //stringstream sstr;
      //sstr<<lf;
      //BOOST_LOG_TRIVIAL(debug) << sstr;
      //unsigned long long a=getmsofday();
      //totalFitError+=lf.optimize(detected_larvae[LRVTRACK_VERBOSE_LEVEL].blobs.back());
      //BOOST_LOG_TRIVIAL(debug) << "FitTime: " << getmsofday()-a << "micros";
      //lf.showContour();
    }
    else
    {
      preBlobs.clear();
      MAX_LARVAE_DETECTED=NEW.size();
      //normalize blobs for next use
      cvb::CvBlobs::iterator it=NEW.begin();
      int i=1;
      while (it!=NEW.end())
      {

        preBlobs[i]=(*it).second;
        preBlobs[i]->label=i;
        it++;
        i++;
      }
      updateLarvae(preBlobs,preBlobs);
      //initialize(detected_larvae[3]);
      LARVAE_COUNT=preBlobs.size();
      //lf.setup(detected_larvae[LRVTRACK_VERBOSE_LEVEL]);
      //lf.setloops(LRVTRACK_DSTEP,LRVTRACK_WSTEP);
      //lf.optimize(detected_larvae[LRVTRACK_VERBOSE_LEVEL].blobs.back());
    }

    //dumpDetectedLarvae();
    if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
    {
      //printSummary(preBlobs,NEW,true);
    }


   if (LRVTRACK_SHOW_LIVE)
   {
     showTags(tracked_blobs);
   }

    cvZero(tmp_P);
    cvRenderBlobs(labelImg, tracked_blobs, tmp_P, tmp_P);

    cvReleaseImage(&tmp_P);
    cvReleaseImage(&labelImg);

    if(LRVTRACK_SHOW_LIVE)
    {
      //make_Video(1,15,17);
      //make_Video(1,7,9);
      imshow("Extracted Frame",colorFrame);
      char k=waitKey(1);
      if(k=='x')
        exit(0);
    }
    //vidPOut << colorFrame;
    //handleKeys(STEP,STOP,TRACK,SHOWTAGS);
    std::cerr << "F: " << CURRENT_FRAME << "/" << TOTAL_FRAMES << "\r";
  }
  std::cerr << endl;
  BOOST_LOG_TRIVIAL(debug) << "Total Fit Error:" << totalFitError ;

  //Second Pass
  secondPass();

  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error)
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  CURRENT_FRAME=0;
  /*while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    //apply_model();
  }*/

  // Third pass
  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error)
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  CURRENT_FRAME=0;
  if (LRVTRACK_SAVE_PROCESSED_VIDEO!="" && !vidPOut.isOpened())
  {
    vidPOut.open(LRVTRACK_SAVE_PROCESSED_VIDEO,
        VIDEO_CODEC,
        VIDEO_FPS,
        bgFrame.size());
    if (!vidPOut.isOpened())
    {
      std::cerr << "Error opening video output file. " <<
        "Problems with video codec. Exiting..." << std::endl;
      exit(5);
    }
  }
  while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    /*cv::circle(colorFrame,
               Point(circles[bestCircle][0],circles[bestCircle][1]),
               circles[bestCircle][2],
               Scalar(255,100,100));*/
    if (SHOWTAGS!=0)
    {
      showTags2();
    }

    stringstream d;
    d<<CURRENT_FRAME;
    Scalar frame_color;
    if(LRVTRACK_INVERT==true)
    {
      frame_color=Scalar(0,0,0);
    }
    else{
      frame_color=Scalar(200,200,200);
    }
    putText(unprocessedFrame,
        d.str(),
        Point2f(20,50),
        FONT_HERSHEY_PLAIN,
        2.0,
        frame_color,
        2,
        CV_AA);

    stringstream n;
    if(LRVTRACK_SAVE_PROCESSED_VIDEO=="")
    {
      vector<int> compression_params;
#if CV_VERSION_MAJOR<3
      compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
#else
      compression_params.push_back(IMWRITE_PNG_COMPRESSION);
#endif
      compression_params.push_back(0);

      n << "img/ff" << CURRENT_FRAME << ".png";
      try {
        imwrite(n.str(), unprocessedFrame, compression_params);
      }
      catch (runtime_error& ex) {
        fprintf(stderr,
            "Exception converting image to PNG format: %s\n",
            ex.what());
        return 1;
      }
      std::cerr << "R: " << CURRENT_FRAME << "/" << TOTAL_FRAMES << "\r";
    }
    else
    {
      //vidPOut << colorFrame;
      vidPOut << unprocessedFrame;
      std::cerr << "R: " << CURRENT_FRAME << "/" << TOTAL_FRAMES << "\r";
    }
  }
  cerr << "HELLO World" << endl;
  cerr << "FINISHED!!!" << endl;
  std::cerr << endl;
  cupContours.release();
  processedFrame.release();
  greyFrame.release();
  unprocessedFrame.release();
  colorFrame.release();
  bgFrame.release();
  previousFrame.release();
  previousOrigFrame.release();
  toRingTemplate.release();
  return 0;
}
