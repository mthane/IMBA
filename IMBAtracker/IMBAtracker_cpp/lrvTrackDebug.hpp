#ifndef __LRVTRACK_DEBUG_HPP__
#define __LRVTRACK_DEBUG_HPP__
#include <string>
#include <iostream>
#include <sstream>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
//Debugging functions for output
extern int LRVTRACK_VERBOSE_LEVEL;

void log_init();
void verbosePrint(std::stringstream &toPrint);
void verbosePrint(const char * toPrint);
void verbosePrint(std::string &toPrint);
std::string printUIMap(std::map<size_t, size_t> &A);
namespace lrvTrack
{
  template <typename number >
    std::string printVector(std::vector<number> vec,int position=0)
    {
      if (vec.size()==0)
        return "";
      std::stringstream sstm;
      //bool const is_number= std::is_arithmetic<number>::value;
      //static_assert( is_number, "Provided type is not an arithmetic type");
      sstm << "[" ;
      typename std::vector<number>::const_iterator i=vec.begin()+position;
      sstm << *i ;
      ++i;
      for( ; i != vec.end(); ++i)
        sstm << ","<< *i;
      sstm << "]";
      return sstm.str();
    }
}

#endif
