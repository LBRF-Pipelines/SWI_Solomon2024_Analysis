#### Define functions  used in the SWI analysis ####

#### Dependencies ####

library(gsignal)

#### Functions ####


# An adjustable filter function for force transducer. 
# It is a wrapper on gsignal functions

filtfiltnew = function(y, type, hz, samples_per_sec, 
                       keep_dc = T, filter_order = 1){
  # f is a 1-D input vector
  # type is the filter type either "low" "high", "stop", "pass" or "null"
  # hz is the frequency cutoff value(s) for the filter
  # keep_dc is a boolean that adds the median value back to the input vector 
  # after baseline correction
  # filter_order sets the filter order
  if(type=='null'){
    return(y)
  }
  
  y_median = median(y)
  
  butterfilt <- gsignal::butter(
    n = filter_order,
    w = hz / (samples_per_sec / 2),
    type = type, 
    plane = 'z',
    output = c("Sos")
  )
  
  newfilt <- gsignal::filtfilt(butterfilt,c(y - y_median))
  
  if(keep_dc){
    return(newfilt + y_median)
  }else{
    return(y)
  }
  
}

c_diff = function(f, h, n){
  # f is a 1-D input vector
  # h is the spacing between elements of f
  # n is the window size for central differentiation. 
  # It can only be a value from the set [3,5,7,9]
  midStartPoint = ceiling(n / 2)
  midEndPoint = length(f) - midStartPoint + 1
  if(n == 3){
    df = (f[midStartPoint:midEndPoint + 1] -
            f[midStartPoint:midEndPoint - 1]) /
      (2 * h)
  }
  if(n == 5){
    df = (f[midStartPoint:midEndPoint - 2] - 
            8 * (f[midStartPoint:midEndPoint - 1]) + 
            8 * (f[midStartPoint:midEndPoint + 1]) - 
            f[midStartPoint:midEndPoint + 2]) /
      (12 * h)
  }
  if(n == 7){
    df = (-f[midStartPoint:midEndPoint - 3] + 
            9 * (f[midStartPoint:midEndPoint - 2]) - 
            45 * (f[midStartPoint:midEndPoint - 1]) +
            45 * (f[midStartPoint:midEndPoint + 1]) -
            9 * (f[midStartPoint:midEndPoint + 2]) +
            f[midStartPoint:midEndPoint + 3]) /
      (60 * h)
  }
  if(n == 9){
    df = (3 * (f[midStartPoint:midEndPoint - 4])-
            32 * (f[midStartPoint:midEndPoint - 3]) +
            168 * (f[midStartPoint:midEndPoint - 2]) -
            672 * (f[midStartPoint:midEndPoint - 1]) +
            672 * (f[midStartPoint:midEndPoint + 1]) -
            168 * (f[midStartPoint:midEndPoint + 2]) +
            32 * (f[midStartPoint:midEndPoint + 3]) - 
            3 * (f[midStartPoint:midEndPoint + 4])) /
      (840 * h)
  }
  return(df)
}

n_point_cen_diff = function(f, h, n){
  # f is a 1-D input vector
  # h is the spacing between elements of f
  # n is the window size for central differentiation. 
  # It can only be a value from the set [3, 5, 7, 9]
  if(!n %in% c(3, 5, 7, 9)){
    stop("n is not a valid value of n point central differentiation.
         Please use one of values from the set [3,5,7,9].")
  }
  
  df_1 = diff(f[1:2]) / h
  df_End = diff(tail(f, 2)) / h
  # Calculate 3-point for all
  df_3pt = c_diff(f, h, 3)
  if (n==3){
    df = c(df_1, df_3pt, df_End)
  }
  if (n>=5){
    df_5pt = c_diff(f, h, 5)
    # For the 2nd and next-to-last grid point, use 3-point differencing.
    df_2 = head(df_3pt, 1)
    df_Endm1 = tail(df_3pt, 1)
    if (n==5){
      df = c(df_1, df_2, df_5pt, df_Endm1, df_End)
    }
  }
  if (n>=7) {
    df_7pt = c_diff(f, h, 7)
    # For the 3nd and 2nd from last grid point, use 5-point differencing.
    df_3 = head(df_5pt, 1)
    df_Endm2 = tail(df_5pt, 1)
    if (n==7){
      df = c(df_1, df_2, df_3, df_7pt, df_Endm2, df_Endm1, df_End)
    }
  }
  if (n>=9){
    df_9pt = c_diff(f,h,9)
    # For the 4nd and 3rd from last grid point, use 7-point differencing.
    df_4 = head(df_7pt,1)
    df_Endm3 = tail(df_7pt,1)
    df = c(df_1, df_2, df_3, df_4, df_9pt, df_Endm3, df_Endm2, df_Endm1, df_End)
  }
  return(df)
}