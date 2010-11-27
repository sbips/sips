cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(FCM_FLAIR)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(FCM_FLAIR FCM_FLAIR.cxx FCM_FLAIR_histo.cxx )

TARGET_LINK_LIBRARIES(FCM_FLAIR ITKCommon ITKStatistics ITKNumerics  ITKIO ITKBasicFilters)
