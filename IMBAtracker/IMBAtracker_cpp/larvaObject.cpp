/*
 * Class containing information about each larva and its history on the plate.
 */
#include "larvaObject.hpp"
#include "lrvTrackBase.hpp"
#include "lrvTrackDebug.hpp"
#include <iostream>
#include <sstream>

using namespace std;
using namespace lrvTrack;

void larvaObject::dump() const
{
  cerr << endl;
  cerr << "============== Larva ID:" << larva_ID << "===============" <<  endl;
  cerr << "start_frame: " << start_frame << endl;
  cerr << "lifetimeWithStats: " << lifetimeWithStats << endl;
  cerr << "lastBlobWithStats: " << lastBlobWithStats << endl;
  cerr << "lastFrameWithStats: " << lastFrameWithStats << endl;
  cerr << "larva_ID: " << larva_ID << endl;
  cerr << "updated_ID: " << updated_ID << endl;
  //cerr << "parentBlobID: " << parentBlobID << endl;
  cerr << "isCluster: " << isCluster << endl;
  cerr << "vector capture_times: " << endl;
  cerr << "  ";
  cerr << "size: " << capture_times.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(capture_times) << endl;
  cerr << "vector blobs: " << endl;
  cerr << "  ";
  cerr << "size: " << capture_times.size() << endl;

  cerr << "vector area: " << endl;
  cerr << "  ";
  cerr << "size: " << area.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(area) << endl;
  cerr << "area_mean: " << area_mean << endl;
  cerr << "area_sum: " << area_sum << endl;
  cerr << "area_max: " << area_max << endl;
  cerr << "area_min: " << area_min << endl;

  cerr << "vector grey_value: " << endl;
  cerr << "  ";
  cerr << "size: " << grey_value.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(grey_value) << endl;
  cerr << "grey_value_mean: " << grey_value_mean << endl;
  cerr << "grey_value_sum: " << grey_value_sum << endl;
  cerr << "grey_value_max: " << grey_value_max << endl;
  cerr << "grey_value_min: " << grey_value_min << endl;

  cerr << "vector length: " << endl;
  cerr << "  ";
  cerr << "size: " << length.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(length) << endl;
  cerr << "length_mean: " << length_mean << endl;
  cerr << "length_sum: " << length_sum << endl;
  cerr << "length_max: " << length_max << endl;
  cerr << "length_min: " << length_min << endl;

  cerr << "vector perimeter: " << endl;
  cerr << "  ";
  cerr << "size: " << perimeter.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(perimeter) << endl;
  cerr << "perimeter_mean: " << perimeter_mean << endl;
  cerr << "perimeter_sum: " << perimeter_sum << endl;
  cerr << "perimeter_max: " << perimeter_max << endl;
  cerr << "perimeter_min: " << perimeter_min << endl;

  cerr << "vector width: " << endl;
  cerr << "  ";
  cerr << "size: " << width.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(width) << endl;
  cerr << "width_mean: " << width_mean << endl;
  cerr << "width_sum: " << width_sum << endl;
  cerr << "width_max: " << width_max << endl;
  cerr << "width_min: " << width_min << endl;

  cerr << "vector headBodyAngle: " << endl;
  cerr << "  ";
  cerr << "size: " << headBodyAngle.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(headBodyAngle) << endl;

  cerr << "vector orientationAngle: " << endl;
  cerr << "  ";
  cerr << "size: " << orientationAngle.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(orientationAngle) << endl;

  cerr << "vector roundness: " << endl;
  cerr << "  ";
  cerr << "size: " << roundness.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(roundness) << endl;

  cerr << "vector angular_speed: " << endl;
  cerr << "  ";
  cerr << "size: " << angular_speed.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(angular_speed) << endl;

  cerr << "vector midpoint_speed_x: " << endl;
  cerr << "  ";
  cerr << "size: " << midpoint_speed_x.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(midpoint_speed_x) << endl;

  cerr << "vector midpoint_speed_y: " << endl;
  cerr << "  ";
  cerr << "size: " << midpoint_speed_y.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(midpoint_speed_y) << endl;

  cerr << "max_midpoint_speed: " << max_midpoint_speed << endl;
  cerr << "min_midpoint_speed: " << min_midpoint_speed << endl;

  cerr << "vector centroid_speed_x: " << endl;
  cerr << "  ";
  cerr << "size: " << centroid_speed_x.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(centroid_speed_x) << endl;

  cerr << "vector centroid_speed_y: " << endl;
  cerr << "  ";
  cerr << "size: " << centroid_speed_y.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(centroid_speed_y) << endl;

  cerr << "max_centroid_speed: " << max_centroid_speed << endl;
  cerr << "min_centroid_speed: " << min_centroid_speed << endl;

}

void larvaObject::csvLine(size_t CURRENT_FRAME, 
                          size_t VIDEO_FPS, 
                          cv::Point2f &cc, 
                          double ppm, 
                          string &csvline)
{
  size_t c_index=CURRENT_FRAME-start_frame;
  cv::Point2f bp=cv::Point2f(blobs[c_index].minx,
      blobs[c_index].miny);
  stringstream data;
  bool DATA_OK=true;
  //data << (float) CURRENT_FRAME/VIDEO_FPS << ",";
  data << CURRENT_FRAME << ",";

  bool rev=false;
  if(lrvDistances[c_index].Spine[0]==
      tails[c_index]+bp &&
      lrvDistances[c_index].Spine.back()!=tails[c_index]+bp)
  {
    auto &d=lrvDistances[c_index].Spine;
    /*if(c_index>0)
      {
      if(!round_flag[c_index] && !round_flag[c_index-1])
      {
      vector<double> spineAnglesPre;
      vector<double> spineAngles;
      vector<double> spineAngleDiff;
      getAnglesFromSpine(lrvDistances[c_index-1].Spine,spineAnglesPre,false);
      getAnglesFromSpine(d,spineAngles,false);
      if(spineAngles.size()==spineAnglesPre.size())
      {
      for(size_t i=0;i<spineAngles.size();i++)
      {
      spineAngleDiff.push_back(spineAngles[i]-spineAnglesPre[i]);
      }
      BOOST_LOG_TRIVIAL(debug) << "Angles: " << printVector(spineAngles);
      BOOST_LOG_TRIVIAL(debug) << "AngleDiff: " << printVector(spineAngleDiff);
      }
      }
      }*/
    for(auto &s : d)
    {
      data << (s.x-cc.x)*ppm << "," << (-s.y+cc.y)*ppm <<",";
    }
    rev=true;
  }
  else if(lrvDistances[c_index].Spine[0]!=
      tails[c_index]+bp &&
      lrvDistances[c_index].Spine.back()==tails[c_index]+bp)
  {
    auto &d=lrvDistances[c_index].Spine;
    /*if(c_index>0)
      {
      if(!round_flag[c_index] && !round_flag[c_index-1])
      {
      vector<double> spineAnglesPre;
      vector<double> spineAngles;
      vector<double> spineAngleDiff;
      getAnglesFromSpine(lrvDistances[c_index-1].Spine,spineAnglesPre,false);
      getAnglesFromSpine(d,spineAngles,false);
      if(spineAngles.size()==spineAnglesPre.size())
      {
      for(size_t i=0;i<spineAngles.size();i++)
      {
      spineAngleDiff.push_back(spineAngles[i]-spineAnglesPre[i]);
      }
      BOOST_LOG_TRIVIAL(debug) << "Angles: " << printVector(spineAngles);
      BOOST_LOG_TRIVIAL(debug) << "AngleDiff: " << printVector(spineAngleDiff);
      }
      }
      }*/
    for(auto s=d.rbegin();s!=d.rend();s++)
    {
      data << (s->x-cc.x)*ppm << ","
        << (-s->y+cc.y)*ppm <<",";
    }
    rev=false;
  }
  else
  {
    cout << "Head tail confusion detected: Larva: " << larva_ID 
      << " Frame: " << CURRENT_FRAME << endl;
    round_flag[c_index]=1;
    DATA_OK=false;
    //csvline="";
    //return;
    for(int di=0;di<12;di++)
    {
        data << "na,na,";
    }
  }

  if(rev & DATA_OK)
  {
    data << (lrvDistances[c_index].Spine[0].x-cc.x)*ppm << ","
      <<  (-lrvDistances[c_index].Spine[0].y+cc.y)*ppm << ",";
    for(auto &p: lrvDistances[c_index].spinePairs)
    {
      data  << (p.second.x-cc.x)*ppm << ","
        << (-p.second.y+cc.y)*ppm <<",";
    }
    data << (lrvDistances[c_index].Spine.back().x-cc.x)*ppm
      << ","
      <<   (-lrvDistances[c_index].Spine.back().y+cc.y)*ppm
      << ",";
    auto &sp=lrvDistances[c_index].spinePairs;
    for(auto p=sp.rbegin();p!=sp.rend();p++)
    {
      data << (p->first.x-cc.x)*ppm << ","
        << (-p->first.y+cc.y)*ppm <<",";
    }
  }
  else if(!rev & DATA_OK)
  {
    data << (lrvDistances[c_index].Spine.back().x-cc.x)*ppm
      << ","
      <<   (-lrvDistances[c_index].Spine.back().y+cc.y)*ppm
      << ",";
    auto &sp=lrvDistances[c_index].spinePairs;
    for(auto p=sp.rbegin();p!=sp.rend();p++)
    {
      data  << (p->first.x-cc.x)*ppm << ","
        << (-p->first.y+cc.y)*ppm <<",";
    }

    data << (lrvDistances[c_index].Spine[0].x-cc.x)*ppm << ","
      <<  (-lrvDistances[c_index].Spine[0].y+cc.y)*ppm << ",";
    for(auto &p: lrvDistances[c_index].spinePairs)
    {
      data << (p.second.x-cc.x)*ppm << ","
        << (-p.second.y+cc.y)*ppm <<",";
    }
  }
  else
  {
    for(int di=0;di<22;di++)
    {
        data << "na,na,";
    }
  }
  data << (blobs[c_index].centroid.x-cc.x)*ppm << ","
    << (blobs[c_index].centroid.y-cc.y)*ppm << ",";
  data << cvAngle(&blobs[c_index]) << ",";
  data << area[c_index] << ",";
  data << grey_value[c_index] << ",";
  data << length[c_index]*ppm << ",";
  data << width[c_index]*ppm << ",";
  data << perimeter[c_index]*ppm << ",";
  data << round_flag[c_index] << "," << ppm;
  data << endl;
  csvline=data.str();
}

int larvaObject::switchFaultyAssignment(
  std::map<size_t,std::vector<size_t> > &detected_clusters,
  std::map<size_t,larvaObject> &detected_larvae)
{
  int i=inCluster.size()-1;
  int collisionIndex;
  int exchangeDuration; //in Frames up to current moment
  size_t blob;

// Look for last collision
  for (; i>=0; i--)
    {
      if (inCluster[i]>0)
        {
          blob=inCluster[i];
          collisionIndex=i;
          break;
        }
    }
  if(i==0)
    {
      std::cerr << "No Collision found, do nothing" << std::endl;
      return(0);
    }

// Find range of values that need to be switched

  // exchangeDuration:
  // The duration of the exchange. In this simple case,
  // it also provides us with the range.
  // For all values:
  //    Range: array.size()-1-exchangeDuration  up to array.size()-1
  // Here it's -2 since the index is always shifted by
  // 1 (starts at 0 and we do not want to change the
  // information while the blob was still there (which is
  // indicated by the collision index)
  exchangeDuration=inCluster.size()-collisionIndex-2;

// Find the other larva involved
  size_t otherLarvaID=0;
  std::vector<size_t>::iterator otherLarvaIT=detected_clusters[blob].begin();
  for (; otherLarvaIT!=detected_clusters[blob].begin(); ++otherLarvaIT)
    {
      if (*otherLarvaIT!=larva_ID)
        {
          otherLarvaID=*otherLarvaIT;
          break;
        }
    }

// Exchange the values

  //larvaObject &otherLarva=detected_larvae[otherLarvaID];
  //int exchangeIdx=0;

  return 0;
}
