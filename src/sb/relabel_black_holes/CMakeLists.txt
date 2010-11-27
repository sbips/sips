PROJECT(relabel_black_holes)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(relabel_black_holes relabel_black_holes.cxx)

TARGET_LINK_LIBRARIES (relabel_black_holes
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




