//
//  Wavelet2d.h
//  QR_Code
//
//  Created by leonardo xu on 23/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#ifndef __QR_Code__Wavelet2d__
#define __QR_Code__Wavelet2d__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <string>

using namespace cv;
using namespace std;

class Wavelet2d {
public:
  void* dwt_2d (vector<vector<double> > &, int, string, vector<double> &, vector<double> &, vector<int> &, string ,string &, bool &);
  void* dwt_2d_proc (string, vector<vector<double> > &, vector<vector<double> > &, vector<vector<double> > &, vector<vector<double> > &, vector<vector<double> > &);
  int filter_coef(string, vector<double> &, vector<double> &, vector<double> &, vector<double> &);
  void* dwt_1d(string, vector<double> &, vector<double> &, vector<double> &);
  void* symm_ext(vector<double> &s, int);
  double convfftm(vector<double> &, vector<double> &, vector<double> &);
  void down_sample(vector<double> &, int, vector<double> &);
  void* dwt_output_relen(vector<int> &, vector<int> &, int);
  void* embed(vector<vector<double> > &cH, vector<vector<double> > &cV, string message_bit, bool &);
  void* extract(vector<vector<double> > cH, vector<vector<double> > cV, string &message_bit);
};



#endif /* defined(__QR_Code__Wavelet2d__) */
