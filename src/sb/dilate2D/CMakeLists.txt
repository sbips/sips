PROJECT(dilate2D)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(dilate2D dilate2D.cxx)

TARGET_LINK_LIBRARIES (dilate2D
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




