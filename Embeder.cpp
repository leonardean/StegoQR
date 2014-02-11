//
//  Embeder.cpp
//  QR_Code
//
//  Created by leonardo xu on 24/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#include "Embeder.h"

//--calculate the max value of a 2d matrix
void* Embeder::max_val_2d(vector<vector<double> > &signal, double &max) {
  max = 0;
  for (unsigned int i = 0; i < (unsigned int) signal.size(); i ++) {
    for (unsigned int j = 0; j < (unsigned int) signal[0].size(); j ++) {
      if (max <= signal[i][j]) {
        max = signal[i][j];
      }
    }
  }
  return 0;
}

//--calculate the max value of a 1d matrix
void* Embeder::max_val_1d(vector<double> &signal, double &max) {
  max = 0;
  for (unsigned int i = 0; i < (unsigned int) signal.size(); i ++) {
    if (max <= signal[i]) {
      max = signal[i];
    }
  }
  return 0;
}

//--calculate a matrix based on dwt output for display all the decomposed image
void* Embeder::img_display(vector<double> &output,vector<vector<double> > &dwtdisp, vector<int> &length , vector<int> &length2, int J) {
  int sum = 0;
  
  for (int iter =0; iter < J; iter++) {
    int d_rows=length[2*iter]-length2[2*iter];
    int d_cols=length[2*iter+1]-length2[2*iter + 1];
    
    int rows_n =length[2 * iter];
    int cols_n = length[2 * iter + 1];
    vector<vector<double> >  dwt_output(2 * rows_n, vector<double>(2 * cols_n));
    if (iter == 0) {
      for(int i =0; i < rows_n; i++){
        for (int j =0; j < cols_n; j++){
          dwt_output[i][j]=output[i*cols_n + j];
        }
      }
      
      for(int i =0; i < rows_n; i++){
        for (int j = cols_n; j < cols_n * 2; j++){
          dwt_output[i][j]= output[rows_n * cols_n + i * cols_n + (j - cols_n)];
        }
      }
      
      for(int i = rows_n; i < rows_n * 2; i++){
        for (int j =0; j < cols_n; j++){
          dwt_output[i][j]=output[2 * rows_n * cols_n+ (i - rows_n) * cols_n + j];
        }
      }
      
      for(int i = rows_n; i < rows_n * 2; i++){
        for (int j = cols_n; j < cols_n * 2; j++){
          dwt_output[i][j]=output[3 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
        }
      }
    } else {
      for(int i =0; i < rows_n; i++){
        for (int j = cols_n; j < cols_n * 2; j++){
          dwt_output[i][j]= output[sum + i * cols_n + (j - cols_n)];
        }
      }
      
      for(int i = rows_n; i < rows_n * 2; i++){
        for (int j =0; j < cols_n; j++){
          dwt_output[i][j]=output[sum + rows_n * cols_n+ (i - rows_n) * cols_n + j];
        }
      }
      
      for(int i = rows_n; i < rows_n * 2; i++){
        for (int j = cols_n; j < cols_n * 2; j++){
          dwt_output[i][j]=output[sum + 2 * rows_n * cols_n+ (i -rows_n) * cols_n + (j -cols_n)];
        }
      }
    }
    
    int rows_x = length2[2*iter];
    int cols_x =length2[2*iter +1];
    
    int d_cols2 = (int) ceil( (double) (d_cols - 1) / 2.0);
    int d_rows2 = (int) ceil( (double) (d_rows - 1) / 2.0);
    if (iter ==0) {
      for(int i =0; i < rows_x; i++){
        for (int j =0; j < cols_x; j++){
          if (i + d_rows -1 < 0){
            dwtdisp[i][j]=0;
          }
          else if (j + d_cols -1 < 0){
            dwtdisp[i][j]=0;
          } else {
            dwtdisp[i][j]=dwt_output[i+d_rows -1][j+d_cols -1];
          }
        }
      }
    }
    for(int i =0; i < rows_x; i++){
      for (int j = cols_x; j < cols_x * 2; j++){
        if (i + d_rows2 < 0){
          dwtdisp[i][j]=0;
        }
        else if (j + 2* (d_cols -1) +1 > (signed) dwt_output[0].size() - 1){
          dwtdisp[i][j]=0;
        } else {
          dwtdisp[i][j]= dwt_output[i+d_rows2 ][j + 2* (d_cols -1)+1 ];
        }
      }
    }
    
    for(int i = rows_x; i < rows_x * 2; i++){
      for (int j =0; j < cols_x; j++){
        if (i + 2* (d_rows -1) + 1 > (signed) dwt_output.size() - 1){
          dwtdisp[i][j]=0;
        }
        else if (j + d_cols2 < 0){
          dwtdisp[i][j]=0;
        } else {
          
          dwtdisp[i][j]=dwt_output[i+2 * (d_rows - 1) + 1 ][j+d_cols2 ];
        }
      }
    }
    
    for(int i = rows_x; i < rows_x * 2; i++){
      for (int j = cols_x; j < cols_x * 2; j++){
        if (i +  (d_rows -1) + 1 + d_rows2 > (signed) dwt_output.size() - 1){
          dwtdisp[i][j]=0;
        }
        else if (j + (d_cols -1) + 1 + d_cols2  > (signed) dwt_output[0].size() - 1){
          dwtdisp[i][j]=0;
        } else {
          dwtdisp[i][j]=dwt_output[i +  (d_rows -1) + 1 + d_rows2 ][j + (d_cols -1) + 1 + d_cols2 ];
        }
      }
    }
    if (iter == 0) {
      sum+= 4*rows_n*cols_n;
    } else {
      sum+= 3*rows_n * cols_n;
    }
  }
  
  return 0;
}

//--data embeder
void* Embeder::data_embed(char* file_name,
                          int J,
                          string &wavelet_name,
                          vector<double> &output,
                          vector<double> &flag,
                          vector<int> &length,
                          int &rows,
                          int &cols,
                          bool &full,
                          string message,
                          char* output_file) {
  
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
  string message_output;
  wavelet2d.dwt_2d(signal,J,wavelet_name,output,flag,length, "embed", message, full);
  cout<<message_output<<endl;
  double max;
  vector<int> length_recal;
  
  //recalculate length of coefficients due to the difference before and after dwt
  wavelet2d.dwt_output_relen(length,length_recal,J);
  int sz = (int)length_recal.size();
  int rows_n=length_recal[sz-2];
  int cols_n = length_recal[sz-1];
  
  //calculate a matrix based on dwt output for display
  vector<vector< double> > dwt_display(rows_n, vector<double>(cols_n));
  img_display(output,dwt_display, length , length_recal, J);
  vector<vector<double> >  dwt_output= dwt_display;
  max_val_2d(dwt_output,max);
  
  //Displaying Scaled Image
  IplImage *cvImg;
  CvSize imgSize;
  imgSize.width = cols_n;
  imgSize.height = rows_n;
  cvImg = cvCreateImage( imgSize, 8, 1 );
  vector<vector<double> > dwt_hold(rows_n, vector<double>( cols_n));
  dwt_hold = dwt_output; // a copy of dwt_output for later use
  
  // Setting coefficients of created image to the scaled DWT output values
  for (int i = 0; i < imgSize.height; i++ ) {
    for (int j = 0; j < imgSize.width; j++ ){
      if ( dwt_output[i][j] <= 0.0){
        dwt_output[i][j] = 0.0;
      }
      if ( i <= (length_recal[0]) && j <= (length_recal[1]) ) {
        ((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
        (char) ( (dwt_output[i][j] / max) * 255.0);
      } else {
        ((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
        (char) (dwt_output[i][j]) ;
      }
    }
  }
  
  //display image
  cvNamedWindow( "DWT Image", 1 );
  cvShowImage( "DWT Image", cvImg );
  cvWaitKey();
  cvDestroyWindow("DWT Image");
  cvSaveImage("dwt.bmp",cvImg);
  
  //reconstruct the image to RGB
  //calculate output of inverse discrete wavelet transform
  vector<vector<double> > idwt_output(rows, vector<double>(cols));
  idwt2d.idwt_2d(output,flag, wavelet_name, idwt_output,length);
  
  //Displaying Reconstructed Image
  IplImage *dvImg;
  CvSize dvSize; // size of output image
  dvSize.width = (int)idwt_output[0].size();
  dvSize.height = (int)idwt_output.size();
  dvImg = cvCreateImage( dvSize, 8, 1 );
  
  //put output inverse dwt into dvImg
  for (int i = 0; i < dvSize.height; i++ )
    for (int j = 0; j < dvSize.width; j++ )
      ((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
      (char) (idwt_output[i][j])  ;
  
  //merge 3 channels together
  IplImage* ipl_G = new IplImage(mat_G);
  IplImage* ipl_B = new IplImage(mat_B);
  IplImage* img_merged = cvCreateImage(dvSize, 8, 3);
  //cvMerge(dvImg, ipl_G, ipl_B, NULL, img_merged);
  
  //display reconstructed image
  cvNamedWindow( "Reconstructed Image", 1 ); // creation of a visualisation window
  cvShowImage( "Reconstructed Image", dvImg ); // image visualisation
  Mat mat_merged(dvImg);
  cv::imwrite(output_file, mat_merged);
  cout<<"Stego image produced."<<endl;
  cvWaitKey();
  cvDestroyWindow("Reconstructed Image");
  cvSaveImage("recon.bmp",dvImg);
  
  return 0;
}