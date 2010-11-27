# This project is designed to be built outside the Insight source tree.
PROJECT(keep_large_objects)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(keep_large_objects keep_large_objects.cxx )

TARGET_LINK_LIBRARIES(keep_large_objects ITKCommon ITKIO ITKBasicFilters ITKStatistics)