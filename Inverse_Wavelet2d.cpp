//
//  Inverse_Wavelet2d.cpp
//  QR_Code
//
//  Created by leonardo xu on 24/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#include "Inverse_Wavelet2d.h"

//--upsample a 1d signal with factor of U
void Inverse_Wavelet2d::upsamp(vector<double> &signal_in, int U, vector<double> &signal_out) {
  int len = (int)signal_in.size();
  double len_n = ceil( (double) len * (double) U);
  
  for (int i = 0; i < (int) len_n; i++) {
    if ( i % U == 0) {
      double temp = signal_in[i / U];
      signal_out.push_back(temp);
    }
    else {
      signal_out.push_back(0);
    }
  }
}

//--binary operation for sum
double op_sum(double i, double j) {
  return (i + j);
}

//--output the sum of two vectors
int vector_sum(vector<double> &a, vector<double> &b, vector<double> &c){
  c.resize(a.size());
  transform (a.begin(), a.end(), b.begin(), b.begin(), op_sum);
  c = b;
  return 0;
}

//--1d inverse discrete wavelet transform with symmetric extension
void* Inverse_Wavelet2d::idwt_1d(string wavelet_name,    //wavelet name
              vector<double> &output, //output of inverse dwt
              vector<double> &cA,     //coefficient of approximation
              vector<double> &cD) {   //coefficient of details
  Wavelet2d wavelet2d;
  
  //lp: low pass decomposition filter
  //hp: high pass decomposition filter
  //lpr: low pass reconstruction filter
  //hpr: high pass reconstruction filter
  vector<double> lp, hp, lpr, hpr;
  wavelet2d.filter_coef(wavelet_name,lp,hp,lpr,hpr);
  
  int lf = (int)lpr.size(); //length of filter
  int U = 2; // Upsampling Factor
  
  // operations in the Low Frequency branch of the Synthesis Filter Bank
  vector<double> X_lp;
  vector<double> cA_up;
  upsamp(cA, U, cA_up );
  cA_up.pop_back();
  wavelet2d.convfftm(cA_up, lpr, X_lp);
  
  // operations in the High Frequency branch of the Synthesis Filter Bank
  vector<double> X_hp;
  vector<double> cD_up;
  upsamp(cD, U, cD_up);
  cD_up.pop_back();
  wavelet2d.convfftm(cD_up, hpr, X_hp);
  
  vector_sum(X_lp, X_hp, output);
  
  //remove symmetric extension
  output.erase(output.begin(), output.begin() + lf - 2);
  output.erase(output.end()-(lf - 2), output.end());
  
  return 0;
  
}

//--2d inverse discrete wavelet transform with symmetric extension
void* Inverse_Wavelet2d::idwt_2d(vector<double>  &dwt_output, //output of dwt
              vector<double> &flag,        // house keeping
              string wavelet_name,         //wavelet name
              vector<vector<double> > &idwt_output, //output inverse dwt
              vector<int> &length){        //length vector obtained from dwt
  //obtain basic information from dwt outputs
  int J = (int)flag[0];
  int rows = length[0];
  int cols = length[1];
  int sum_coef = 0;
  Wavelet2d wavelet2d;
  
  //obtain wavelet filter coefficients
  vector<double> lp, hp, lpr, hpr;
  wavelet2d.filter_coef(wavelet_name, lp, hp, lpr, hpr);
  vector<vector<double> > cLL(rows, vector<double>(cols));
  unsigned int lf = (unsigned int)lp.size();
  
  for (int iter = 0; iter < J; iter ++) {
    int rows_recal = length[2 * iter];
    int cols_recal = length[2 * iter + 1];
    
    //initialize wavelet reconstruction coefficients
    //cLH vertical details
    //cHL horizontal details
    //cHH diagonal details
    vector<vector<double> > cLH(rows_recal, vector<double>(cols_recal));
    vector<vector<double> > cHL(rows_recal, vector<double>(cols_recal));
    vector<vector<double> > cHH(rows_recal, vector<double>(cols_recal));
    
    //allocating wavelet reconstruction coefficients
    for (int i = 0 ; i < rows_recal; i++ ){
      for (int j = 0; j < cols_recal; j++){
        if (iter == 0) {
          cLL[i][j] = dwt_output[sum_coef+ i * cols_recal + j];
          cLH[i][j] = dwt_output[sum_coef+ rows_recal * cols_recal+ i * cols_recal + j];
          cHL[i][j] = dwt_output[sum_coef+ 2 * rows_recal * cols_recal + i * cols_recal + j];
          cHH[i][j] = dwt_output[sum_coef+ 3* rows_recal * cols_recal + i * cols_recal + j];
        } else {
          cLH[i][j] = dwt_output[sum_coef+  i * cols_recal + j];
          cHL[i][j] = dwt_output[sum_coef+ rows_recal * cols_recal + i * cols_recal + j];
          cHH[i][j] = dwt_output[sum_coef+ 2* rows_recal * cols_recal + i * cols_recal + j];
        }
      }
    }
    
    //obtain size of coefficients in both dimension
    unsigned int len_x = (unsigned int)cLH.size();
    unsigned int len_y = (unsigned int)cLH[0].size();
    
    //row upsampling and column filtering in first low pass branch
    vector<vector<double> > cL(2 * len_x - lf + 2, vector<double>(len_y ));
    vector<vector<double> > cH(2 * len_x - lf + 2, vector<double>(len_y ));
    
    if (iter == 0) {
      for (unsigned int j = 0; j < len_y; j ++) {
        vector<double> signal_LL, signal_LH, output;
        for (unsigned int i = 0; i < len_x; i ++) {
          double tmp1 = cLL[i][j];
          double tmp2 = cLH[i][j];
          signal_LL.push_back(tmp1);
          signal_LH.push_back(tmp2);
        }
        
        //apply 1d idwt
        idwt_1d(wavelet_name,output,signal_LL,signal_LH);
        for (int i = 0;i < (int)output.size(); i ++) {
          cL[i][j] = output[i];
        }
      }
    } else{
      //row upsampling and column filtering in low pass branch
      unsigned int rows1 = (unsigned int)cLH.size();
      unsigned int cols1 = (unsigned int)cLH[0].size();
      
      for (unsigned int j = 0; j < cols1; j++){
        vector<double> signal_LL, signal_LH, output;
        for (unsigned int i = 0; i < rows1; i++){
          double tmp1 = cLL[i][j];
          double tmp2 = cLH[i][j];
          signal_LL.push_back(tmp1);
          signal_LH.push_back(tmp2);
        }
        
        //apply 1d idwt
        idwt_1d(wavelet_name,output,signal_LL,signal_LH);
        for (unsigned int i = 0; i < output.size(); i ++){
          cL[i][j]=output[i];
        }
      }
    }
    
    //row upsampling and column filtering in high pass branch
    for (unsigned int j = 0; j < len_y; j ++) {
      vector<double> signal_HL, signal_HH, output;
      for (unsigned int i = 0;i < len_x; i ++) {
        double tmp1 = cHL[i][j];
        double tmp2 = cHH[i][j];
        signal_HL.push_back(tmp1);
        signal_HH.push_back(tmp2);
      }
      
      //apply 1d idwt
      idwt_1d(wavelet_name,output,signal_HL,signal_HH);
      for (int i = 0;i < (int)output.size(); i ++) {
        cH[i][j] = output[i];
      }
    }
    
    //column upsampling and row filtering in both branch
    vector<vector<double> > signal(2 * len_x-lf + 2, vector<double>(2 * len_y - lf +2));
    for (unsigned int i = 0; i < 2 * len_x - lf +2; i++) {
      vector<double> signal_L,signal_H,output;
      for (unsigned int j=0;j <  len_y;j++) {
        double tmp1 = cL[i][j];
        double tmp2 = cH[i][j];
        signal_L.push_back(tmp1);
        signal_H.push_back(tmp2);
      }
      
      //apply 1d idwt
      idwt_1d(wavelet_name,output,signal_L,signal_H);
      for (int j = 0;j < (int)output.size(); j ++) {
        signal[i][j] = output[j];
      }
    }
    
    idwt_output = signal;
    
    if (iter == 0) {
      sum_coef+= 4 *rows_recal * cols_recal;
    } else {
      sum_coef+= 3 *rows_recal * cols_recal;
    }
    cLL = signal;
  }
  return 0;
}