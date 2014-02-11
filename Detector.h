//
//  Detector.h
//  QR_Code
//
//  Created by leonardo xu on 21/08/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#ifndef __QR_Code__Detector__
#define __QR_Code__Detector__

#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"

#define MINHESSIAN 1100
#define NOCTAVES 2
#define NOCTAVELAYERS 4
#define R 32

using namespace cv;
using namespace std;

class Detector {
  
  Mat cover_img;
  vector<KeyPoint> keypoints;
  
public:
  
  Detector();
  void draw_region(KeyPoint point, Mat image);
  bool on_edge( KeyPoint point, Mat image);
  float calc_dist( KeyPoint pt1, KeyPoint pt2 );
  void detect (char* img_name);
  void eliminate ();
  void draw();
  void surf_detect ( char* img_name);
  vector<KeyPoint> get_keypoints();
  Mat get_cover_img();
};

#endif /* defined(__QR_Code__Detector__) */
