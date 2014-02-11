//
//  Embeder.h
//  QR_Code
//
//  Created by leonardo xu on 24/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#ifndef __QR_Code__Embeder__
#define __QR_Code__Embeder__

#include "Wavelet2d.h"
#include "Inverse_Wavelet2d.h"

class Embeder {
public:
  void* max_val_2d(vector<vector<double> > &, double &);
  void* max_val_1d(vector<double> &, double &);
  void* img_display(vector<double> &,vector<vector<double> > &, vector<int> &, vector<int> &, int);
  //void data_embed(char*, int);
  void* data_embed(char*, int, string &, vector<double> &, vector<double> &, vector<int> &, int &, int&, bool &, string, char*);
};

#endif /* defined(__QR_Code__Embeder__) */
