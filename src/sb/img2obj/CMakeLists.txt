PROJECT(img2obj)

  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)



ADD_EXECUTABLE(img2obj 
img2obj.cxx 
itkAnalyzeObjectEntry.cxx 
itkAnalyzeObjectLabelMapImageIO.cxx
itkAnalyzeObjectLabelMapImageIOFactory.cxx)

TARGET_LINK_LIBRARIES (img2obj
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  )




