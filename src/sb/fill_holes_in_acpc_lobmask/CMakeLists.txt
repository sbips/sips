PROJECT(fill_holes_in_acpc_lobmask)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(fill_holes_in_acpc_lobmask fill_holes_in_acpc_lobmask.cxx)

TARGET_LINK_LIBRARIES (fill_holes_in_acpc_lobmask
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




