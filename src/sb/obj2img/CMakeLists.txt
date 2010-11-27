PROJECT(obj2img)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(obj2img 
obj2img.cxx 
itkAnalyzeObjectEntry.cxx 
itkAnalyzeObjectLabelMapImageIO.cxx
itkAnalyzeObjectLabelMapImageIOFactory.cxx)

TARGET_LINK_LIBRARIES (obj2img
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




