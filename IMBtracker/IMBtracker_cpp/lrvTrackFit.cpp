#include "lrvTrackDebug.hpp"
#include "lrvTrackFit.hpp"
#include "blobUtils.hpp"
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

using namespace lrvTrack;

double modelError(Mat &real,Mat &model)
{
  Mat Diff;
  if(real.cols != model.cols || real.rows != model.rows)
	  return INT_MAX;
  bitwise_xor(real,model*2.5,Diff);
  imwrite("err.png", Diff);
  size_t nz=countNonZero(Diff);
  return nz;
}

class optimizeBody
{
  private:
    std::vector<fitData> &fitSpace;
    const lrvFit &lf;
    Mat &ref;
    Mat &bg;
  public:
    double minerr;
    int mini;
    optimizeBody(std::vector<fitData> &fs,
        lrvFit &l,
        Mat &aref,
        Mat &abg
        ):
      fitSpace(fs),
      lf(l),
      ref(aref),
      bg(abg),
      minerr(DBL_MAX),
      mini(-1)
  {
  }

    optimizeBody(optimizeBody &x, tbb::split) :
    fitSpace(x.fitSpace),
    lf(x.lf),
    ref(x.ref),
    bg(x.bg),
    minerr(DBL_MAX),
    mini(-1)
  {}
    void join(const optimizeBody &y)
    {
      if(y.minerr<=minerr)
      {
        minerr=y.minerr;
        mini=y.mini;
      }
    }

    void operator ()(const tbb::blocked_range<size_t>& range)
    {
      for (int i = range.begin(); i < range.end(); ++i)
      {
        lrvFit nl=lf;
        nl.larvaFitData=fitSpace[i];
        double r;
        if(bg.empty())
          r=nl.errorFunc(ref);
        else
          r=nl.errorFunc(ref,bg);
        if(r<=minerr)
        {
          minerr=r;
          mini=i;
        }
      }
    }
};

/*void collisionModel::csvLine(size_t CURRENT_FRAME, 
                          size_t VIDEO_FPS, 
                          cv::Point2f &cc, 
                          double ppm, 
                          std::map<size_t , std::string> &csvlines)
{
  if(SFRAME<CURRENT_FRAME || EFRAME>CURRENT_FRAME)
    return;
  for(auto &l: larvae_models)
  {
    size_t c_index=CURRENT_FRAME-SFRAME;
    stringstream data;

    data << (float) CURRENT_FRAME/VIDEO_FPS << ",";
    int i=0;
    for (auto &s : l[c_index].spine)
    {
      data << (s.x-cc.x)*ppm << "," << (-s.y+cc.y)*ppm <<",";
    }
    for (auto &s : l[c_index].cpoints)
    {
        data << (s.x-cc.x)*ppm << "," << (-s.y+cc.y)*ppm <<",";
    }
    data << (l[c_index].spine[5].x-cc.x)*ppm << ","
      << (-l[c_index].spine[5].y+cc.x)*ppm << ",";
    data << 2;
    data << endl;
    csvlines[l[c_index].ID] = data.str();
  }
}*/

collisionModel::collisionModel(std::vector<size_t> clarvae,
                                    cvb::CvBlob &blob,
                                    size_t START_FRAME,
                                    size_t END_FRAME)
{
  SFRAME=START_FRAME;
  EFRAME=END_FRAME;
  larvae_models.resize(clarvae.size());
  frameError.resize(abs((long) (START_FRAME-END_FRAME))+1);
  size_t lrv=0;
  larvae_models[0].resize(clarvae.size());
  Point2f bp(blob.minx-PADDING-1,blob.miny-PADDING-1);
  Mat CollisionROI; 
  Mat CollisionObj; 

  createLarvaContour_custom(CollisionObj,
      blob,
      CV_8UC1,
      blob.minx,
      blob.maxx,
      blob.miny,
      blob.maxy,
      PADDING+1,
      true,
      cv::Scalar(255),
      8,
      cv::Scalar(0),
      false);


  char name[1024];
  sprintf(name, "56_collisionModel_Obj_%d_%d.png", blob.label, START_FRAME);
  	imwrite(name, CollisionObj);

  for(auto &cl : clarvae)
  {
    larvae_models[lrv].resize(abs((long) (START_FRAME-END_FRAME))+1);
    lrvFit &c=larvae_models[lrv][0];
    larvaObject &o=detected_larvae[cl];
    size_t idx=START_FRAME-o.start_frame-1;
    c.setup(o,blob.minx,blob.maxx,blob.miny,blob.maxy,START_FRAME);
    Mat f;//(blob.maxx-blob.minx+PADDING,blob.maxy-blob.miny+PADDING,CV_8UC1);
    Mat sp;//(blob.maxx-blob.minx+PADDING,blob.maxy-blob.miny+PADDING,CV_8UC1);

    createLarvaContour_custom(f,
        o.blobs[idx],
        CV_8UC1,
        blob.minx,
        blob.maxx,
        blob.miny,
        blob.maxy,
        PADDING+1,
        true,
	cv::Scalar(255),
	8,
	cv::Scalar(0),
	false);
    f.copyTo(sp);
    //Create copy of f with spine marks
    for(auto &spp : o.lrvDistances[idx].Spine)
    {
      cv::circle(sp,spp-bp,0,cv::Scalar(200),-1);
    }
    Mat r;
    c.optimizeAndReturn(f,r);

    cerr << "CollisionROI: " << CollisionROI.rows << " x " << CollisionROI.cols << endl;
    cerr << "CollisionROI: " << CollisionROI.type() << endl;
    cerr << "r: " << r.rows << " x " << r.cols << endl;
    cerr << "r: " << r.type() << endl;
    
    if(CollisionROI.empty())
      sp.copyTo(CollisionROI);
    else
      CollisionROI|=r; //collisionROI = CoillisionROI | R


    lrv++;
  }
    sprintf(name, "57_collisionModel_R_%d_%d.png", blob.label, START_FRAME);
  	imwrite(name, CollisionROI);

    frameError[0]=modelError(CollisionROI,CollisionObj);
  //Mat exp;
  //tile2same(CollisionROI,CollisionObj,exp);
  /*imshow("Model Outcome",exp);
  waitKey(1);*/
}

void collisionModel::createBg(int exclude,size_t FRAME,Mat &bg)
{
  int idx=FRAME-SFRAME;
  for(auto &cm: larvae_models)
  {
    if(cm[idx].ID==exclude)
      continue;
    if(bg.empty())
    {
      cm[idx].createMatfromFit(bg);
    }
    else
    {
      Mat n;
      cm[idx].createMatfromFit(n);
      bitwise_or(n,bg,bg);
    }
  }
}

void collisionModel::updateModel(cvb::CvBlob &blob,size_t FRAME,Mat &ret)
{
  Mat result;
  Mat f;
  createLarvaContour_custom(f,
      blob,
      CV_8UC1,
      blob.minx,
      blob.maxx,
      blob.miny,
      blob.maxy,
      PADDING+1);

  size_t idx=0;
  for(auto &cm: larvae_models)
  {
    lrvFit &c=cm[FRAME-SFRAME];
    lrvFit &pc=cm[FRAME-SFRAME-1];
    c=pc;
    c.changeBaseCoords(blob.minx,blob.maxx,blob.miny,blob.maxy);
  }
  for(auto &cm: larvae_models)
  {
    Mat pl;
    Mat ret;
    Mat bg;
    lrvFit &c=cm[FRAME-SFRAME];
    createBg(c.ID,FRAME,bg);

    if(bg.empty())
    {
      Mat fg;
      f.copyTo(bg);
      c.createMatfromFit(fg);
      dilate(fg,fg,Mat(),Point(-1,-1),2);
      bitwise_not (fg,fg); // invert
      bitwise_and(bg,fg,bg); // masking
    }
    //c.createMatfromFit(pl,bg);
    c.optimize(f,bg,ret);
    c.createMatfromFit(pl);
    c.createMatfromFit(pl,bg);
    result=pl;
    idx++;
  }
  //Mat exp;
  frameError[FRAME-SFRAME]=modelError(result,f);
  tile2same(result,f,ret);
  //imshow("Model Outcome",ret);
  //waitKey(10000);
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

/*double optimize(cvb::CvBlob &blob, std::vector<larvaObject> &l)
{
  for(auto i=0;i<l.size();i++)
  {
    larvaFit f(l[i]);
    f.optimize(blob);
  }
}*/

ostream& operator<<(ostream& os,lrvFit &l)
{
  os << "[l:" << l.larvaFitData.length << ", w:" <<
        l.larvaFitData.width << ", mt:" <<
        l.larvaFitData.midtail << ", ga:" <<
        l.larvaFitData.global_angle << ", as:" <<
        printVector(l.larvaFitData.angles) << "]";
  return os;
}

void lrvFit::csvLine(size_t CURRENT_FRAME, 
                          size_t SFRAME,
                          size_t EFRAME,
                          size_t VIDEO_FPS, 
                          cv::Point2f &cc, 
                          double ppm, 
                          std::string &csvline)
{
  calculateContour2f();
  if(SFRAME>CURRENT_FRAME || EFRAME<CURRENT_FRAME)
    return;
  //size_t c_index=CURRENT_FRAME-SFRAME;
  stringstream data;

  //data << (float) CURRENT_FRAME/VIDEO_FPS << ",";
  data << (float) CURRENT_FRAME << ",";
  //int i=0;
  for (auto &s : spine)
  {
    data << (s.x-cc.x)*ppm << "," << (-s.y+cc.y)*ppm <<",";
  }
  for (auto &s : contour)
  {
    data << (s.x-cc.x)*ppm << "," << (-s.y+cc.y)*ppm <<",";
  }
  data << (spine[5].x-cc.x)*ppm << ","
    << (-spine[5].y+cc.x)*ppm << ",";
  data <<  ",";
  data <<  ",";
  data <<  ",";
  data <<  ",";
  data <<  ",";
  data <<  ",";
  data << 2;
  data << endl;
  csvline = data.str();
}

void lrvFit::paintPoly(Mat &ROI, vector<Point> f,size_t fsize)
{
  fillConvexPoly(ROI,&f[0],fsize,Scalar(255));
}

double lrvFit::errorFunc(Mat &b1C)
{
  Mat l1C,Diff;

  createMatfromFit(l1C);
  //std::cout << "l1C size: " << l1C.cols << "x" << l1C.rows << std::endl;
  //std::cout << "b1C size: " << b1C.cols << "x" << b1C.rows << std::endl;
  if(l1C.cols != b1C.cols || l1C.rows != b1C.rows)
  {
	  //std::cout << "DIFFERENCE BETWEEN SIZES OF MODEL AND REFERENCE!!!" << std::endl;
	  return INT_MAX;
  }
  bitwise_xor(b1C,l1C*2.5,Diff);
  size_t nz=countNonZero(Diff);
  return nz;
}

double lrvFit::errorFunc(Mat &real,
                           Mat &bg
                           )
{
  Mat l1C,Diff;

  createMatfromFit(l1C,bg,true);
  if(l1C.cols != real.cols || l1C.rows != real.rows)
  {
	  std::cout << "DIFFERENCE BETWEEN SIZES OF MODEL AND REFERENCE!!!" << std::endl;
	  return INT_MAX;
  }
  bitwise_xor(real,l1C*2.5,Diff);
  size_t nz=countNonZero(Diff);
  return nz;
}

/*
 * From point1 with global angle <angle> and distance <d> we return point p2
 */
void lrvFit::pointToPoint(cv::Point2f &p1, double angle, double d, cv::Point2f &p2)
{
  double qpi=0.5*CV_PI;
  double dpi=1.5*CV_PI;

  if(angle>=2*CV_PI)
    angle-=2*CV_PI;

  if(angle<qpi && angle>=0)
  {
    double k=tan(angle);
    double xt=sqrt((d*d)/(1+k*k));
    p2.x=p1.x+xt;
    double yt=k*xt;
    p2.y=p1.y+yt;
  }
  else if(angle>qpi && angle<dpi)
  {
    double k=tan(angle);
    double xt=-sqrt((d*d)/(1+k*k));
    p2.x=p1.x+xt;
    double yt=k*xt;
    p2.y=p1.y+yt;
  }
  else if(angle>dpi && angle<2*CV_PI)
  {
    double k=tan(angle);
    double xt=sqrt((d*d)/(1+k*k));
    p2.x=p1.x+xt;
    double yt=k*xt;
    p2.y=p1.y+yt;
  }
  else
  {
    p2.x=p1.x;
    p2.y=p1.y+(((angle-qpi)/qpi)-1)*d;
  }
}

void lrvFit::setloops(int ds,double ws)
{
  degree_step=ds*0.0174532925199432957692369076848;
  wstep=ws;
  setloops();
}

void lrvFit::setloops()
{
    ga.clear();
    //ga.push_back(-2*degree_step);
    ga.push_back(-degree_step);
    ga.push_back(0);
    ga.push_back(degree_step);
    //ga.push_back(2*degree_step);

    a1.clear();
    a1.push_back(-2*degree_step);
    a1.push_back(-degree_step);
    a1.push_back(0);
    a1.push_back(degree_step);
    a1.push_back(2*degree_step);

    a2.clear();
    //a2.push_back(-2*degree_step);
    //a2.push_back(-degree_step);
    a2.push_back(0);
    //a2.push_back(degree_step);
    //a2.push_back(2*degree_step);

    //a1=ga;
    //a2=ga;
    a3=ga;
    a4=a2;


    wl.clear();
    //wl.push_back(-2*wstep);
    //wl.push_back(-wstep);
    wl.push_back(0);
    //wl.push_back(wstep);
    //wl.push_back(2*wstep);

}

void lrvFit::filterAngle(std::vector<double> &a,
                           double &val,
                           double lim,
                           double add)
{
  if(val+a[0]<=lim)
    transform(a.begin(), a.end(), a.begin(),
        bind2nd(std::plus<double>(), add));
  else if(val+a.back()>=2*CV_PI-lim)
    transform(a.begin(), a.end(), a.begin(),
        bind2nd(std::minus<double>(), add));
}

void lrvFit::generate(std::vector<fitData> &fitSpace)
{

  filterAngle(a1,larvaFitData.angles[0],2*degree_step,2*degree_step);
  filterAngle(a2,larvaFitData.angles[1],2*degree_step,2*degree_step);
  filterAngle(a3,larvaFitData.angles[2],2*degree_step,degree_step);
  filterAngle(a4,larvaFitData.angles[3],2*degree_step,degree_step);
  //filterAngle(ga,larvaFitData.global_angle,degree_step);

  fitSpace=std::vector<fitData>((a1.size()*
                                a2.size()*
                                a3.size()*
                                a4.size()*
                                ga.size()*wl.size()*81)+1);
  //fitSpace=std::vector<fitData>(6076);

  fitData &l=larvaFitData;
  //BOOST_LOG_TRIVIAL(debug) << printVector(l.angles);

 //Find midpoint tail equation to minimize sideways motion
  fitSpace[0]=larvaFitData;
  double mtx = spine[8].x; //midtail
  double mty = spine[8].y;
  double tx = spine[11].x; //tail
  double ty = spine[11].y;
  /*double a,b,br;
  double bw=2.0; //Width of zone 
  bool vert=false;
  bool hor=false;
  bool fail=false;
  if(mtx-tx!=0 && mty-ty!=0)
  {
    a = (mty-ty)/(mtx-tx) ; //slope
    br = abs(a)*bw;
  }
  else if( mty-ty == 0 && mtx-tx==0)
  {
    fail=true;
    BOOST_LOG_TRIVIAL(debug) << "FAILED FINDING OUT SLOPE WHEN RESTRICTED MODEL'S SIDEWAYS MOTION" << endl;
  }
  else if (mty-ty==0)
  {
    hor=true;
  }
  else if (mtx-tx==0)
  {
    vert=true;
  }*/
  //We want to have about 6 search steps
  //	larvae can MAXIMALLY! move about 8mm per second in one direction
  //    search space: +- 2/FPS mm
  //    search space in px: +-(2/FPS)/MilimetersPerPixel
  double max_pix_per_frame_f = (12/VIDEO_FPS)/LRVTRACK_MPP;
  int coord_step = (int) max_pix_per_frame_f/6;
  if (coord_step < 1)
    coord_step = 1;
  int max_pix_per_frame = max_pix_per_frame_f;
  // BOOST_LOG_TRIVIAL(debug) << "MODEL STATS: coord_step=" << coord_step << endl;
  // BOOST_LOG_TRIVIAL(debug) << "MODEL STATS: LRVTRACK_MPP=" << LRVTRACK_MPP << endl;
  
  size_t i=1;
  for(auto iag=0;iag<ga.size();++iag)
  {
    for(auto ia1=0; ia1<a1.size() ; ++ia1)
    {
      for(auto ia2=0; ia2<a2.size() ; ++ia2)
      {
        for(auto ia3=0; ia3<a3.size() ; ++ia3)
        {
          for(auto ia4=0; ia4<a4.size() ; ++ia4)
          {
            for(auto wp=0;wp<wl.size();wp++)
            {
	      for(auto mpx=-4*coord_step;mpx<=4*coord_step;mpx+=coord_step)
              {
              	for(auto mpy=-4*coord_step;mpy<=4*coord_step;mpy+=coord_step)
                {
                  Point2f mp=Point2f(larvaFitData.midtail.x+mpx,
                              larvaFitData.midtail.y+mpy);
		  /*if(!fail && !hor && !vert)
		  {
	            double yminusax=mpy-a*mpx;
		    if(yminusax>br || yminusax<-br)
		    {
                      mp=Point2f(larvaFitData.midtail.x,
                              larvaFitData.midtail.y);
		    }
		  }
		  if(!fail && !hor && vert)
		  {
		    if(mpx>bw || mpx<-bw)
		    {
                      mp=Point2f(larvaFitData.midtail.x,
                              larvaFitData.midtail.y);
		    }
		  }
		  if(!fail && !vert && hor)
		  {
		    if(mpy>bw || mpy<-bw)
		    {
                      mp=Point2f(larvaFitData.midtail.x,
                              larvaFitData.midtail.y);
		    }
		  }*/
                  fitSpace[i++]=fitData(mp,
                                        orig_length*(1.0+wl[wp]),
                                        orig_width*(1.0-wl[wp]),
                                        l.angles[0]+a1[ia1],
                                        l.angles[1]+a2[ia2],
                                        l.angles[2]+a3[ia3],
                                        l.angles[3]+a4[ia4],
                                        l.global_angle+ga[iag]
                                        );
                }
              }
            }
          }
        }
      }
    }
  }
}

double lrvFit::optimize(Mat &ref)
{
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  double minerr=DBL_MAX;

  if(LRVTRACK_PARALLEL)
  {
    Mat empty;
    optimizeBody body(fitSpace,*this,ref,empty);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, fitSpace.size()), body);
    larvaFitData=fitSpace[body.mini];
    minerr=body.minerr;
  }
  else
  {
    std::vector<fitData>::iterator min=fitSpace.begin()+1;

    for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
    {
      larvaFitData=*i;
      //cout << *this << endl;
      double r=errorFunc(ref);

      if(r<minerr)
      {
        minerr=r;
        min=i;
      }
    }
    larvaFitData=*min;
  }

  //Mat l1C;
  //createMatfromFit(l1C);
  //cv::imshow("Fit Contour", l1C);
  //cv::imshow("Old Frame", greyFrame);
  //cv::imshow("Old Contour", ref);
  //cv::moveWindow("Old Contour",200,200);
  //cv::moveWindow("Fit Contour",200,130);
  //cv::waitKey(1);
  return minerr;

}

double lrvFit::optimizeAndReturn(Mat &ref,Mat &ret)
{
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  double minerr=DBL_MAX;
  if(LRVTRACK_PARALLEL)
  {
    Mat empty;
    optimizeBody body(fitSpace,*this,ref,empty);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, fitSpace.size()), body);
    larvaFitData=fitSpace[body.mini];
    minerr=body.minerr;
  }
  else
  {
    std::vector<fitData>::iterator min=fitSpace.begin()+1;

    for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
    {
      larvaFitData=*i;
      //cout << *this << endl;
      double r=errorFunc(ref);

      if(r<minerr)
      {
        minerr=r;
        min=i;
      }
    }
    larvaFitData=*min;
  }
  createMatfromFit(ret);
  return minerr;
}

double lrvFit::optimize(Mat &ref,Mat &bg, Mat &ret)
{
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  double minerr=DBL_MAX;
  if(LRVTRACK_PARALLEL)
  {
    optimizeBody body(fitSpace,*this,ref,bg);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, fitSpace.size()), body);
    larvaFitData=fitSpace[body.mini];
    minerr=body.minerr;
  }
  else
  {
    std::vector<fitData>::iterator min=fitSpace.begin()+1;

    for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
    {
      larvaFitData=*i;
      //cout << *this << endl;
      double r=errorFunc(ref,bg);

      if(r<minerr)
      {
        minerr=r;
        min=i;
      }
    }
    larvaFitData=*min;
  }

  //Mat l1C;
  createMatfromFit(ret);
  //cv::imshow("Old Frame", greyFrame);
  /*cv::imshow("Old Contour", ref);
  cv::imshow("Model Contour", ret);
  cv::waitKey(1);*/
  /*cv::imshow("Fit Contour", l1C);
  cv::moveWindow("Old Contour",200,200);
  cv::moveWindow("Fit Contour",200,130);*/
  return minerr;
}

double lrvFit::optimize(Mat &ref,Mat &bg)
{
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  double minerr=DBL_MAX;
  if(LRVTRACK_PARALLEL)
  {
    optimizeBody body(fitSpace,*this,ref,bg);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, fitSpace.size()), body);
    larvaFitData=fitSpace[body.mini];
    minerr=body.minerr;
  }
  else
  {
    std::vector<fitData>::iterator min=fitSpace.begin()+1;

    for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
    {
      larvaFitData=*i;
      //cout << *this << endl;
      double r=errorFunc(ref,bg);

      if(r<minerr)
      {
        minerr=r;
        min=i;
      }
    }
    larvaFitData=*min;
  }

  /*Mat l1C;
  createMatfromFit(l1C);
  cv::imshow("Fit Contour", l1C);
  cv::imshow("Old Contour", ref);
  cv::waitKey(1);*/
  //cv::imshow("Old Frame", greyFrame);
  //cv::moveWindow("Old Contour",200,200);
  //cv::moveWindow("Fit Contour",200,130);
  return minerr;
}

void lrvFit::showContour()
{
  Mat l1C;
  createMatfromFit(l1C);
  cv::imshow("Contour", l1C);
  cv::imshow("Picture", greyFrame);
  cv::waitKey(1);
}

lrvFit::lrvFit(larvaObject &l,int dstep, double wstep,
               size_t minx,size_t maxx,size_t miny,size_t maxy, size_t FRAME)
{
  lrvFit(l,minx,maxx,miny,maxy,FRAME);
  setloops(dstep,wstep);
}

lrvFit::lrvFit(larvaObject &l,
               size_t a_minx,
               size_t a_maxx,
               size_t a_miny,
               size_t a_maxy, 
               size_t FRAME)
{
  setup(l,a_minx,a_maxx,a_miny,a_maxy,FRAME);
}

void lrvFit::setup(larvaObject &l,
               size_t a_minx,
               size_t a_maxx,
               size_t a_miny,
               size_t a_maxy, 
               size_t FRAME)
{
  //ID=l.larva_ID;
  ID=l.updated_ID;
  spine=l.lrvDistances.back().Spine;
  Point2f bp(l.blobs.back().minx,l.blobs.back().miny);
  if(l.heads.back()+bp==spine.back())
    reverse(spine.begin(),spine.end());
  intSpine=vector<Point2f>(IPOL*(spine.size()-1)+1,Point2f(-1,-1));
  cpoints=vector<Point>(2*(intSpine.size()-1));
  double width=l.width.back();
  orig_width=width/2;
  double length=l.length.back();
  orig_length=length;
  minx=a_minx;
  maxx=a_maxx;
  miny=a_miny;
  maxy=a_maxy;

  Point2f V1=spine[6]-spine[8];
  double angleg=angleD(V1,Point2f(0,0),Point2f(1,0));
  if(angleg<0)
    angleg=2*CV_PI+angleg;

  vector<double> angles(4,0.0);
  /*angles[0]=((l.lrvDistances.back().Angles[1]));
  angles[1]=((l.lrvDistances.back().Angles[3]));
  angles[2]=((l.lrvDistances.back().Angles[5]));
  angles[3]=((l.lrvDistances.back().Angles[7]));*/
  angles[0]=angleD(spine[4],
                   spine[2],
                   spine[0]);
  angles[1]=angleD(spine[6],
                   spine[4],
                   spine[2]);
  angles[2]=angleD(spine[8],
                   spine[6],
                   spine[4]);
  angles[3]=angleD(spine[6],
                   spine[8],
                   spine[11]);
  larvaFitData=fitData(spine[8],
                        length,
                        width/2,
                        angles,
                        angleg);
  Mat test;
  createMatfromFit(test);
  setloops();
}

/*
 * Sets up spine from larvaFitData
 */
void lrvFit::setupSpine()
{
  double seglen=larvaFitData.length/11;
  spine[8].x=larvaFitData.midtail.x;
  spine[8].y=larvaFitData.midtail.y;

  //For point spine[6]
  pointToPoint(spine[8],larvaFitData.global_angle,seglen*2,spine[6]);
  //For point spine[11]
  double phi=larvaFitData.global_angle+(2*CV_PI-larvaFitData.angles[3]);
  //double phi=larvaFitData.angles[3];
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[8],phi,seglen*3,spine[11]);
  //For point spine[10] we interpolate
  spine[9]=((spine[11]-spine[8])*(1.0/3.0))+spine[8];
  spine[10]=(spine[9]+spine[11])*0.5;

  //For point spine[4]
  phi=(2*CV_PI-larvaFitData.angles[2])-(CV_PI-larvaFitData.global_angle);
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[6],phi,seglen*2,spine[4]);

  //For point spine[2]
  phi=(2*CV_PI-larvaFitData.angles[1])-(CV_PI-phi);
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[4],phi,seglen*2,spine[2]);

  //For point spine[0]
  phi=(2*CV_PI-larvaFitData.angles[0])-(CV_PI-phi);
  if(phi>2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[2],phi,seglen*2,spine[0]);

  spine[7]=(spine[6]+spine[8])*0.5;
  spine[5]=(spine[4]+spine[6])*0.5;
  spine[3]=(spine[2]+spine[4])*0.5;
  spine[1]=(spine[0]+spine[2])*0.5;

  //BOOST_LOG_TRIVIAL(debug) << "Setup Spine: Computed Spine: " << printVector(spine) << endl;
}
void lrvFit::changeBaseCoords(size_t new_minx,
    size_t new_maxx,
    size_t new_miny,
    size_t new_maxy)
{
  int dminx=new_minx-minx;
  int dminy=new_miny-miny;
  
  Point bp(dminx,dminy);
  Point2f fbp(dminx,dminy);
  for(auto &p:cpoints)
  {
    p=p-bp;
  }
  for(auto &p:spine)
  {
    p=p-fbp;
  }
  for(auto &p:intSpine)
  {
    p=p-fbp;
  }
  larvaFitData.midtail-=fbp;
  minx=new_minx;
  maxx=new_maxx;
  miny=new_miny;
  maxy=new_maxy;
}

void lrvFit::createMatfromFit(Mat &larvaFitContour,
                                Mat &fitBase,
                                bool verbose)
{
  setupSpine();
  Mat tmp=cv::Mat(maxy-miny+1+(2*PADDING+2),
        maxx-minx+1+(2*PADDING+2),
        CV_8UC1,Scalar(0));

  Point2f bp(minx-PADDING-1,miny-PADDING-1);

  if(fitBase.rows!=tmp.rows || fitBase.cols!=tmp.cols)
    return;

  tmp=Scalar(0);

  for(size_t i=0; i<spine.size()-1;i++)
  {
    size_t v=i*IPOL;
    for(auto j=0;j<IPOL;j++)
    {
      intSpine[v+j]=(spine[i+1]-spine[i])*((double)j/IPOL)+spine[i];
    }
  }
  intSpine.back()=spine.back();

  cpoints[0]=intSpine[0]-bp;
  cpoints[cpoints.size()/2]=intSpine.back()-bp;
  for(int i=1;i<(int)spine.size()-1;i++)
  {
    circle(tmp,
           0.5*(spine[i-1]+spine[i])-bp,
           w((double)(i)/(double)spine.size())*(larvaFitData.width),
           Scalar(255),
           -1);
    circle(tmp,
           spine[i]-bp,
           w((double)i/(double)spine.size())*(larvaFitData.width),
           Scalar(255),-1);
  }
  /*for(auto i=1;i<intSpine.size()-1;i++)
  {
    calculateContourPoints(intSpine[i-1],
                           intSpine[i],
                           intSpine[i+1],
                           bp,
                           //2.0*(i+1.0)/intSpine.size(),
                           (double) i/intSpine.size(),
                           larvaFitData.width,
                           cpoints[i],
                           cpoints[cpoints.size()-i]);
  }*/
  //size_t csz=cpoints.size();
  //vector<vector<Point> > C;
  //C.push_back(cpoints);
  //fillPoly(tmp,C,Scalar(255));
  /*for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(100));
  }*/
  larvaFitContour = tmp | fitBase;
}

void lrvFit::createMatfromFit(Mat &larvaFitContour,
                                bool verbose)
{
  setupSpine();
  Mat tmp=cv::Mat(maxy-miny+1+(2*PADDING+2),
        maxx-minx+1+(2*PADDING+2),
        CV_8UC1,Scalar(0));

  Point2f bp(minx-PADDING-1,miny-PADDING-1);

  for(size_t i=0; i<spine.size()-1;i++)
  {
    size_t v=i*IPOL;
    for(auto j=0;j<IPOL;j++)
    {
      intSpine[v+j]=(spine[i+1]-spine[i])*((double)j/IPOL)+spine[i];
    }
  }
  intSpine.back()=spine.back();

  cpoints[0]=intSpine[0]-bp;
  cpoints[cpoints.size()/2]=intSpine.back()-bp;
  for(int i=1;i<(int)spine.size()-1;i++)
  {
    circle(tmp,
           0.5*(spine[i-1]+spine[i])-bp,
           w((double)(i)/(double)spine.size())*(larvaFitData.width),
           Scalar(255),
           -1);
    circle(tmp,
           spine[i]-bp,
           w((double)i/(double)spine.size())*(larvaFitData.width),
           Scalar(255),-1);
  }
  /*for(auto i=1;i<intSpine.size()-1;i++)
  {
    calculateContourPoints(intSpine[i-1],
                           intSpine[i],
                           intSpine[i+1],
                           bp,
                           //2.0*(i+1.0)/intSpine.size(),
                           (double) i/intSpine.size(),
                           larvaFitData.width,
                           cpoints[i],
                           cpoints[cpoints.size()-i]);
  }*/

  //size_t csz=cpoints.size();
  //const Point *p=&cpoints[0];
  //fillPoly(tmp,(const Point **) &p,(int *) &csz,1,Scalar(255));
  /*for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(100));
  }*/
  tmp.copyTo(larvaFitContour);
}

void lrvFit::calculateContourPoints(Point2f &a,
                            Point2f &b,
                            Point2f &c,
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
  width=w(b_index)*width;
  if(fabs(sp.x)<0.0001 && fabs(sp.y)<0.0001) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=-r.x/r.y;
      double xc1=sqrt((width*width)/(1+pslope*pslope));
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
      lPoint.x=width;
      lPoint.y=0;
      rPoint.x=-width;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=width;
      rPoint.x=0;
      rPoint.y=-width;
    }

  }
  else
  {
    double ratio=width/sqrt(sp.x*sp.x+sp.y*sp.y);
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
    cl=lPoint+b;
    cr=rPoint+b;
  }
  else if((avals[2]<=avals[0] || avals[2]>=avals[1]) &&
          (avals[3]>=avals[0] || avals[3]<=avals[1]))
  {
    //points are reverse lPoint is to the left, rPoint is to the right
    cl=rPoint+b;
    cr=lPoint+b;
  }
}

void lrvFit::calculateContourPoints(Point2f &a,
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
  width=w(b_index)*width;
  if(fabs(sp.x)<0.0001 && fabs(sp.y)<0.0001) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=-r.x/r.y;
      double xc1=sqrt((width*width)/(1+pslope*pslope));
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
      lPoint.x=width;
      lPoint.y=0;
      rPoint.x=-width;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=width;
      rPoint.x=0;
      rPoint.y=-width;
    }

  }
  else
  {
    double ratio=width/sqrt(sp.x*sp.x+sp.y*sp.y);
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
}

void lrvFit::calculateContour2f()
{
  contour.resize((spine.size()*2)-2);
  contour[0]=spine[0];
  contour[spine.size()-1]=spine.back();
  Point2f bp(minx-PADDING-1,miny-PADDING-1);
  for(auto i=1;i<spine.size()-1;i++)
  {
    calculateContourPoint2f(spine[i-1],
                           spine[i],
                           spine[i+1],
                           (double) i/spine.size(),
                           larvaFitData.width,
                           contour[i],
                           contour[contour.size()-i]);
  }

}

void lrvFit::calculateContourPoint2f(Point2f &a,
                            Point2f &b,
                            Point2f &c,
                            double b_index,
                            double width,
                            Point2f &cl,
                            Point2f &cr)
{
  //Make b our reference point
  Point2f r=a-b;
  Point2f f=c-b;
  Point2f sp=r+f;
  Point2f lPoint,rPoint;
  width=w(b_index)*width;
  if(fabs(sp.x)<0.0001 && fabs(sp.y)<0.0001) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=-r.x/r.y;
      double xc1=sqrt((width*width)/(1+pslope*pslope));
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
      lPoint.x=width;
      lPoint.y=0;
      rPoint.x=-width;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=width;
      rPoint.x=0;
      rPoint.y=-width;
    }

  }
  else
  {
    double ratio=width/sqrt(sp.x*sp.x+sp.y*sp.y);
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
    cl=lPoint+b;
    cr=rPoint+b;
  }
  else if((avals[2]<=avals[0] || avals[2]>=avals[1]) &&
          (avals[3]>=avals[0] || avals[3]<=avals[1]))
  {
    //points are reverse lPoint is to the left, rPoint is to the right
    cl=rPoint+b;
    cr=lPoint+b;
  }
}
