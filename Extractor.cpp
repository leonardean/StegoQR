//
//  Extractor.cpp
//  QR_Code
//
//  Created by leonardo xu on 25/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#include "Extractor.h"
#include <string>

//--data embeder
void* Extractor::data_extract(char* file_name,
                          int J,
                          string &wavelet_name,
                          vector<double> &output,
                          vector<double> &flag,
                          vector<int> &length,
                          int &rows,
                          int &cols,
                          string &message) {
  
  //variables initialization
  Wavelet2d wavelet2d;
  Inverse_Wavelet2d idwt2d;
  int height, width, num_ch, pix_depth;
  CvSize size;
  
  //read img from file name
  IplImage* img_orig = cvLoadImage(file_name);
  
  //obtain properties of the input image
  height = img_orig->height;
  width = img_orig->width;
  num_ch = img_orig->nChannels;
  pix_depth = img_orig->depth;
  size.width =width;
  size.height=height;
  
  //split channels
  IplImage* R = cvCreateImage(size, pix_depth, 1);
  IplImage* G = cvCreateImage(size, pix_depth, 1);
  IplImage* B = cvCreateImage(size, pix_depth, 1);
  cvSplit(img_orig, R, G, B, NULL);
  
  Mat mat_G(G);
  Mat mat_B(B);
  
  //display original image as reference
  cvNamedWindow("Original Image");
  cvShowImage("Original Image", img_orig);
  Mat img_orig_mat(img_orig);
  //cv::imwrite("/Users/leonardo/Desktop/dadada.png", img_orig_mat);
  cvWaitKey();
  cvDestroyWindow("Original Image");
  cvSaveImage("orig.bmp",img_orig);
  
  //turn image to 2d signal for ease of use, use R channel for dwt
  rows =(int) height;
  cols =(int) width;
  Mat matimg(img_orig);
  vector<vector<double> > signal(rows, vector<double>(cols));
  int k =1;
  for (int i=0; i < rows; i++) {
    for (int j =0; j < cols; j++){
      unsigned char temp;
      temp = ((uchar*) matimg.data + i * matimg.step)[j  * matimg.elemSize() + k ];
      signal[i][j] = (double) temp;
    }
  }
  
  //obtain filter coefficients
  vector<double> lp,hp,lpr,hpr;
  wavelet2d.filter_coef(wavelet_name,lp,hp,lpr,hpr);
  
  //apply 2d dwt with symmetric extension
  bool full;
  string message_output;
  wavelet2d.dwt_2d(signal,J,wavelet_name,output,flag,length, "extract", message_output, full);
  cout<<"Hidden message is:"<<message_output<<endl;
  
  return 0;
}


