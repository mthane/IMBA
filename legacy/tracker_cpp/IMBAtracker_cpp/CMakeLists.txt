cmake_minimum_required(VERSION 3.10)
set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBS OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project("Larvae Tracker")

add_executable(lrvTrack
	lrvTrackOL.cpp
	blobUtils.cpp
	lpsolver.cpp
	lrvTrackBase.cpp
	lrvTrackDebug.cpp
	larvaDistanceMap.cpp
	larvaObject.cpp
	lrvTrackFit.cpp)

configure_file(lrvTrackConfig.h.in lrvTrackConfig.h)

set(cvBlob_DIR ${CMAKE_SOURCE_DIR}/cvblob/)
find_package(cvBlob REQUIRED)
include_directories(${cvBlob_DIR}/cvBlob/)
link_directories(${cvBlob_DIR}/lib/)
target_link_libraries(lrvTrack ${cvBlob_DIR}/lib/libcvblob.so)

find_package(OpenCV REQUIRED COMPONENTS core highgui videoio video imgcodecs imgproc ml)
if(OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})
  link_directories(${OpenCV_LIBRARY_DIRS})
  target_link_libraries(lrvTrack ${OpenCV_LIBRARIES})
endif()

find_package(ALGLIB REQUIRED)
if (ALGLIB_FOUND)
  include_directories(${ALGLIB_INCLUDE_DIRS})
  #  link_directories(${ALGLIB_LIBRARY_DIRS})
  target_link_libraries(lrvTrack alglib)
endif()

include(FindLpSolve)
find_package(LpSolve REQUIRED)
if (LpSolve_FOUND)
	include_directories(${LPSOLVE_INCLUDE_DIRS})
	link_directories(${LPSOLVE_LIBRARY_DIRS})
	target_link_libraries(lrvTrack ${LPSOLVE_LIBRARIES})
endif()

find_package(Boost REQUIRED COMPONENTS log log_setup program_options filesystem system timer thread)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(lrvTrack ${Boost_LIBRARIES})
endif(Boost_FOUND)

include(FindTBB)
find_package(TBB REQUIRED)
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIRS})
  link_directories(${TBB_LIBRARY_DIRS})
  target_link_libraries (lrvTrack ${TBB_LIBRARIES})
endif(TBB_FOUND)

