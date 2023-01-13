#ifndef __LRVTRACK_LARVAOBJECT_HPP__
#define __LRVTRACK_LARVAOBJECT_HPP__
#include "cvblob.h"
#include <vector>
#include "larvaDistanceMap.hpp"

/*
 * Class containing information about each blob
 * The class contains a bunch of synchronized vectors.
 */

/*class blobObject_OL{
  private:
    // ID of the larva
    //size_t ID;

    // ID of the originating blob of the larva.
    //  if the blob was there from the begining then
    //    originating_ID=ID
    //  if the blob is comming after dissapearing
    //    originating_ID=0
    size_t originating_ID;

    // ID where the blob entered before it dissappeared
    //  if the blob ends without a collision or dissapearing
    //  ending_ID=ID
    //  if the blob dissapears
    //  ending_ID=0
    
    size_t ending_ID;

    // Frame number of Nth data instance
    std::vector<size_t> FRAME;

    // Frames of existence until the Nth data instance
    std::vector<size_t> FRAMES_COUNT;

    //Vectors of length, width,size, perimeter, grey_value,various speeds etc.
    std::vector<double> length;
    std::vector<double> width;
    std::vector<double> size;
    std::vector<double> perimeter;
    std::vector<double> greyvalue;
    std::vector<double> spine0_speedx;
    std::vector<double> spine0_speedy;
    std::vector<double> spineM_speedx;
    std::vector<double> spineM_speedy;
    std::vector<double> spineL_speedx;
    std::vector<double> spineL_speedy;
    std::vector<double> centroid_speedx;
    std::vector<double> centroid_speedy;

    //Vector of centroids
    std::vector<double> centroids;

    //Vector of detailed Contour, spine reconstructions
    std::vector<larvaDistanceMap> larvaDistances;

    //Vector of original blobs
    std::vector<cvb::CvBlob> blobs;
    
    //Vector of original blobs
    //std::vector<larvaFit> fitLarva;

    bool collisionObject;
};*/

/*
 * Class containing information about each larva and its history on the plate.
 */
class larvaObject
{
public:
  size_t start_frame;
  size_t end_frame;
  size_t lifetimeWithStats;
  size_t lastBlobWithStats;
  size_t lastFrameWithStats;
  size_t larva_ID;
  size_t updated_ID;
  std::vector<size_t> diverged_to;
  size_t converged_to;
  bool isCluster;
  int blobSize;
  std::vector<double> capture_times;
  std::vector<cvb::CvBlob> blobs; //Blob for each frame for a given larva
  std::vector<double> area;
  double area_mean;
  double area_sum;
  double area_max;
  double area_min;

  std::vector<double> grey_value;
  double grey_value_mean;
  double grey_value_sum;
  double grey_value_max;
  double grey_value_min;

  std::vector<double> length;
  double length_mean;
  double length_sum;
  double length_max;
  double length_min;

  std::vector<double> perimeter;
  double perimeter_mean;
  double perimeter_sum;
  double perimeter_max;
  double perimeter_min;

  std::vector<double> width;
  double width_mean;
  double width_sum;
  double width_max;
  double width_min;

  std::vector<double> minHTDist;
  double minHTDist_mean;
  double minHTDist_sum;

  std::vector<double> headBodyAngle;
  std::vector<double> orientationAngle;

  std::vector<double> roundness;
  double roundness_sum;
  double roundness_mean;
  std::vector<double> angular_speed;

  std::vector<cv::Point2f> centroids;
  std::vector<cv::Point2f> centroids_full;
  std::vector<double> midpoint_speed_x;
  std::vector<double> midpoint_speed_y;
  double max_midpoint_speed;
  double min_midpoint_speed;

  std::vector<double> centroid_speed_x;
  std::vector<double> centroid_speed_y;
  double max_centroid_speed;
  double min_centroid_speed;

  std::vector<double> centroid_distance_x;
  std::vector<double> centroid_distance_y;
  double centroid_distance_x_sum;
  double centroid_distance_y_sum;

  std::vector<size_t> inCluster;

  std::vector<larvaDistanceMap> lrvDistances;

  std::vector<cv::Point2f> heads;
  std::vector<double> heads_brightness;
  std::vector<cv::Point2f> tails;
  std::vector<double> tails_brightness;

  std::vector<bool> round_flag;
  void csvLine(size_t CURRENT_FRAME, 
      size_t VIDEO_FPS, 
      cv::Point2f &cc, 
      double ppm, 
      std::string &csvline);

  larvaObject():
    start_frame(0),
    end_frame(0),
    lifetimeWithStats(0),
    lastBlobWithStats(0),
    larva_ID(0),
    converged_to(0),
    isCluster(false),
    blobSize(1),
    area_mean(0),
    area_sum(0),
    area_max(0),
    area_min(0),
    grey_value_mean(0),
    grey_value_sum(0),
    grey_value_max(0),
    grey_value_min(0),
    length_mean(0),
    length_sum(0),
    length_max(0),
    length_min(0),
    perimeter_mean(0),
    perimeter_sum(0),
    perimeter_max(0),
    perimeter_min(0),
    width_mean(0),
    width_sum(0),
    width_max(0),
    width_min(0),
    minHTDist_mean(0),
    minHTDist_sum(0),
    centroid_distance_x_sum(0),
    centroid_distance_y_sum(0),
    round_flag(false)
  {}
  
  void dump() const;
  int switchFaultyAssignment(
    std::map<size_t,std::vector<size_t> > &detected_clusters,
    std::map<size_t,larvaObject> &detected_larvae
  );
};
#endif
