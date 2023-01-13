#include "lrvTrackBase.hpp"
#ifdef LRVTRACK_WITH_CUDA
#include "opencv2/gpu/gpu.hpp"
#endif

#ifdef LRVTRACK_WITH_OPENCL
#include "opencv2/core/ocl.hpp"
#endif

#ifdef LRVTRACK_WITH_CUDA
void lrvTrackNormalize(cv::Mat &src,
                              cv::Mat &dst,
                              double alpha,
                              double beta,
                              int norm_type)
{
  cv::gpu::GpuMat dst_host, src_host;
  src_host.upload(src);
  cv::gpu::normalize(src_host,dst_host,alpha,beta,norm_type);
  dst_host.download(dst);
}

void lrvTrackBilateral(cv::Mat &src,
                              cv::Mat &dst,
                              int ksize,
                              float sigma,
                              double sigmaSpace)
{
  cv::gpu::GpuMat dst_host,src_host;
  src_host.upload(src);
  cv::gpu::bilateralFilter(src_host,dst_host,ksize,sigma,sigmaSpace);
}

/*#elif defined(LRVTRACK_WITH_OPENCL)
void lrvTrackNormalize(cv::Mat &_src , cv::Mat &_dst, double a, double b,
                              int norm_type)
{
  cv::ocl::oclMat src(_src);
  std::cout << "Using OPENCL normalization" << std::endl;
  int rtype;
  double scale = 1, shift = 0;
  if( norm_type == CV_MINMAX )
    {
      double smin = 0, smax = 0;
      double dmin = MIN( a, b ), dmax = MAX( a, b );
      cv::ocl::minMaxLoc( src, &smin, &smax, 0, 0);
      scale = (dmax - dmin)*(smax - smin > DBL_EPSILON ? 1./(smax - smin) : 0);
      shift = dmin - smin*scale;
    }
  else if( norm_type == CV_L2 || norm_type == CV_L1 || norm_type == CV_C )
    {
      scale = cv::ocl::norm( src, norm_type );
      scale = scale > DBL_EPSILON ? a/scale : 0.;
      shift = 0;
    }
  else
    CV_Error( CV_StsBadArg, "Unknown/unsupported norm type" );

  rtype =  src.depth();

  _dst.create(_src.dims , _src.size, CV_MAKETYPE(rtype, _src.channels()));
  cv::ocl::oclMat dst(_dst);

  src.convertTo( dst, rtype, scale, shift );
  dst.upload(_dst);
}

void lrvTrackBilateral(cv::Mat &_src,
                              cv::Mat &_dst,
                              int ksize,
                              float sigma,
                              double sigmaSpace)
{
  cv::ocl::oclMat src(_src);
  cv::ocl::oclMat dst(_dst);
  int rtype =  src.depth();
  _dst.create(_src.dims , _src.size, CV_MAKETYPE(rtype, _src.channels()));
  cv::ocl::bilateralFilter(src,dst,ksize,sigma,sigmaSpace);
  dst.upload(_dst);
}
*/
#else
void lrvTrackNormalize(cv::Mat &src,
                              cv::Mat &dst,
                              double alpha,
                              double beta,
                              int norm_type)
{
  cv::normalize(src,dst,alpha,beta,norm_type);
}

void lrvTrackBilateral(cv::Mat &src,
                              cv::Mat &dst,
                              int ksize,
                              float sigma,
                              double sigmaSpace)
{
  bilateralFilter(src,dst,ksize,sigma,sigmaSpace);
}
#endif
