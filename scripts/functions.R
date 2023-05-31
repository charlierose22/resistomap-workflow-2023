# remove unsatisfactory flags
flag <- function(data) 
  {data[!grepl('CurveFitFail|MultipleMeltPeak|NoAmplification',
                                     data$flags),]}
