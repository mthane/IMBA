#ifndef __LRVTRACK_LRVTRACKBASE_HPP__
#define __LRVTRACK_LRVTRACKBASE_HPP__
#include <opencv2/imgproc/imgproc.hpp>

#define MIN_DISCARD_DISTANCE 30
#define ROI_PADDING 0
#define MAX_HORIZ_RESOLUTION 22000 //Pixels
#define ROUNDNESS_THRESHOLD 230000

void lrvTrackNormalize(cv::Mat &src,
    cv::Mat &dst,
    double alpha,
    double beta,
    int norm_type);

void lrvTrackBilateral(cv::Mat &src,
    cv::Mat &dst,
    int ksize,
    float sigma,
    double sigmaSpace);
#endif
