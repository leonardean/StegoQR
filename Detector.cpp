//
//  Detector.cpp
//  QR_Code
//
//  Created by leonardo xu on 21/08/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#include "Detector.h"

//--constructor
Detector::Detector(){}

//--draw the area of 2R*2R around key point
void Detector::draw_region (KeyPoint point, Mat image) {
  vector<Point2i> corner(4);
  corner[0] = cvPoint(point.pt.x - R, point.pt.y - R);
  corner[1] = cvPoint(point.pt.x + R, point.pt.y - R);
  corner[2] = cvPoint(point.pt.x - R, point.pt.y + R);
  corner[3] = cvPoint(point.pt.x + R, point.pt.y + R);
  line( image, corner[0], corner[1], Scalar(0, 255, 0), 4 );
  line( image, corner[0], corner[2], Scalar(0, 255, 0), 4 );
  line( image, corner[3], corner[1], Scalar(0, 255, 0), 4 );
  line( image, corner[3], corner[2], Scalar(0, 255, 0), 4 );
  imshow("Salient Areas", image);
}

//--if a keypoint is on the edge of the image
bool Detector::on_edge( KeyPoint point, Mat image) {
  bool result = false;
  if (point.pt.x < R)
    result = true;
  if (point.pt.y < R)
    result = true;
  if (image.cols - point.pt.x < R)
    result = true;
  if (image.rows - point.pt.y < R)
    result = true;
  return result;
}

//--calculate the distance between two key points
float Detector::calc_dist( KeyPoint pt1, KeyPoint pt2 ) {
  return sqrt((pt1.pt.x - pt2.pt.x) * (pt1.pt.x - pt2.pt.x) +
              (pt1.pt.y - pt2.pt.y) * (pt1.pt.y - pt2.pt.y));
}

//--detect the salient regions using SURF and obtain key points
void Detector::detect (char* img_name) {
  cover_img = imread (img_name, CV_LOAD_IMAGE_GRAYSCALE );
  
  int minHessian = MINHESSIAN;
  int noctaves = NOCTAVES;
  int noctaveLayers = NOCTAVELAYERS;
  
  SurfFeatureDetector detector( minHessian, noctaves, noctaveLayers);
  detector.detect(cover_img, keypoints);
}

//--to make points disjoint, remove closed points and those on edge
void Detector::eliminate () {
  //print number of points before eliminationg
  //cout<<"Number of points: "<<keypoints.size()<<endl;
  
  //eliminate closed points
  for (int i = 0; i < keypoints.size() - 1; i ++) {
    for (int j = i + 1; j < keypoints.size(); j ++) {
      if (calc_dist(keypoints.at(i), keypoints.at(j)) < 2 * R * sqrt(2.0)){
        keypoints.erase(keypoints.begin() + j);
        j --;
      }
    }
  }
  
  //remove points on edge
  for (int i = 0; i < keypoints.size(); i ++) {
    if (on_edge(keypoints.at(i), cover_img) == true) {//false) {
      keypoints.erase(keypoints.begin() + i);
      i --;
    }
  }
  
  //--print the number of points before eliminationg
  //cout<<"Number of points: "<<keypoints.size()<<endl;
  double str_sum = 0;
  for (int i = 0; i < keypoints.size(); i ++) {
    str_sum += keypoints.at(i).response;
  }
  //cout<<"Avg point str: "<<str_sum / keypoints.size()<<endl;
}

//--draw key points and salient regions
void Detector::draw() {
  //image of keypoints
  Mat img_keypoints;
  
  //draw detected (drawn) keypoints
  drawKeypoints( cover_img, keypoints, img_keypoints, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
  
  //draw salient regions
  for (int i = 0; i < keypoints.size(); i ++){
    draw_region(keypoints.at(i), img_keypoints);
  }
//  for (KeyPoint tmp : keypoints) {
//    draw_region(tmp, img_keypoints);
//  }
  
  waitKey(0);
}

//--surf detection
void Detector::surf_detect ( char* img_name) {
  detect(img_name);
  eliminate();
  draw();
}

//--get keypionts
vector<KeyPoint> Detector::get_keypoints() {
  return keypoints;
}

//--get cover image
Mat Detector::get_cover_img() {
  return cover_img;
}