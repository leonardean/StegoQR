//
//  Inverse_Wavelet2d.h
//  QR_Code
//
//  Created by leonardo xu on 24/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#ifndef __QR_Code__Inverse_Wavelet2d__
#define __QR_Code__Inverse_Wavelet2d__

#include "Wavelet2d.h"

class Inverse_Wavelet2d {
public:
  void upsamp(vector<double> &, int, vector<double> &);
  //double op_sum(double, double);
  //int vector_sum(vector<double> &, vector<double> &, vector<double> &);
  void* idwt_1d(string, vector<double> &, vector<double> &, vector<double> &);
  void* idwt_2d(vector<double>  &,vector<double> &, string, vector<vector<double> > &, vector<int> &);
};

#endif /* defined(__QR_Code__Inverse_Wavelet2d__) */
