cmake_minimum_required(VERSION 2.6)

# This project is designed to be built outside the Insight source tree.
PROJECT(FLEX_fcm)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(FLEX_fcm FLEX_fcm.cxx main.cxx )

TARGET_LINK_LIBRARIES(FLEX_fcm ITKCommon ITKIO ITKBasicFilters )
