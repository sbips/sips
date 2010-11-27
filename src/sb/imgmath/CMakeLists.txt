PROJECT(imgmath)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(imgmath imgmath.cxx)

TARGET_LINK_LIBRARIES (imgmath
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




