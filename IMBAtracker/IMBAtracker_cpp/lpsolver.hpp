#ifndef __LRVTRACK_LPSOLVER_HPP__
#define __LRVTRACK_LPSOLVER_HPP__
#include <opencv2/core.hpp>
#include <vector>
#include <map>

bool TUmatrixTest(cv::Mat m);
void LPmaxFunction(cv::Mat &mf);
int LPconstr(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex);
int LPconstr(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex,
            double avg_size);
int LPconstrAvg(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex,
	    double avgSize);
//bool testPLQuery();
#endif
