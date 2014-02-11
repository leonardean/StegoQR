//
//  Wavelet2d.cpp
//  QR_Code
//
//  Created by leonardo xu on 23/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#include "Wavelet2d.h"
#include "math.h"

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a < _b ? _a : _b; })

fftw_plan plan_forward_inp,plan_forward_filt, plan_backward;
static unsigned int transient_size_of_fft = 0;

//--emded payload bits
void* Wavelet2d::embed(vector<vector<double> > &cH, vector<vector<double> > &cV, string message_bit, bool & full){
  
  //check if the length of the message fill out the area
  int length_msg = (int)message_bit.length();
  if (length_msg >= cH.size() * cH[0].size()) {
    full = true;
  } else {
    full = false;
  }
  
  //store the length of message into the first vertical & horizontal coeff
  cH[0][0] = length_msg;
  cV[0][0] = length_msg;
 
  //modify vertical and horizontal coefficients according to payload bit
  for (int i = 1; i <= length_msg; i ++) {
    
    //coefficient location
    int m = floor((i) / cH.size());
    int n = (i) % cH.size();
    double T = 500;
    double D = cH[m][n] - cV[m][n];
    double D2 = cV[m][n] - cH[m][n];
    //cout<<i<<" "<<m<<" "<<n<<endl;
    //cout<<"i:"<<i<<endl;
    if ((message_bit[i - 1] == '1') && (D < T)) {
      cH[m][n] = cH[m][n] + (T - D) / 2;
      cV[m][n] = cV[m][n] - (T - D) / 2;
      //cout<<"1cH: "<<cH[m][n]<<" cV: "<<cV[m][n]<<" "<<i<<endl;
    } else if ((message_bit[i - 1] == '1') && (D >= T)) {
      //cout<<"11cH: "<<cH[m][n]<<" cV: "<<cV[m][n]<<" "<<i<<endl;
    }else if ((message_bit[i - 1] == '0') && (D2 < T)) {
      cH[m][n] = cH[m][n] - (T - D2) / 2;
      cV[m][n] = cV[m][n] + (T - D2) / 2;
      //cout<<"0cH: "<<cH[m][n]<<" cV: "<<cV[m][n]<<" "<<i<<endl;
    } else {
      //cout<<"00cH: "<<cH[m][n]<<" cV: "<<cV[m][n]<<" "<<i<<endl;
    }
  }
  return 0;
}

int round_up(int num, int multiple) {
  if (multiple == 0) {
    return num;
  }
  
  int remainder = num % multiple;
  if (remainder == 0)
    return num;
  return num + multiple - remainder;
}

int round_down(int num, int multiple) {
  if (multiple == 0) {
    return num;
  }
  
  int remainder = num % multiple;
  if (remainder == 0) {
    return num;
  }
  return num - remainder;
}

int round(int num, int multiple){
  int result = 0;
  
  int num_rup = round_up(num, multiple);
  int num_rdn = round_down(num, multiple);
  int abs_dist_up = abs(num_rup - num);
  int abs_dist_dn = abs(num - num_rdn);
  if (abs_dist_up > abs_dist_dn) {
    result = abs_dist_dn;
  } else if (abs_dist_up <= abs_dist_dn)
    result = abs_dist_up;
  return result;
}

//--extract payload bit
void* Wavelet2d::extract(vector<vector<double> > cH, vector<vector<double> > cV, string &message_bit) {
  //clear the string
  message_bit.clear();
  
  //set the length of message as fixed
  int length_msg = 16;
  
  //determine payload bit and append into string
  for (int i = 1; i <= length_msg; i ++) {
    int m = floor((i) / cH.size());
    int n = (i) % cH.size();
    //cout<<i<<" "<<m<<" "<<n<<endl;
    if (cH[m][n] > cV[m][n]) {
      message_bit.append("0");
      //cout<<"1cH: "<<cH[m][n]<<"cV: "<<cV[m][n]<<endl;
    } else if (cH[m][n] < cV[m][n]) {
      message_bit.append("1");
      //cout<<"0cH: "<<cH[m][n]<<"cV: "<<cV[m][n]<<endl;
    }
  }
  
  for (int i = 0; i < cH.size(); i ++) {
    for (int j = 0; j < cH[0].size(); j ++) {
      if (cH[i][j]  - cV[i][j] >=400) {
        //cout<<"dongxi: "<<i<<" "<<j<<" "<<cH[i][j]<<" "<<cV[i][j]<<endl;
      }
      //cout<<"cH: "<<cH[i][j]<<"cV: "<<cV[i][j]<<endl;
    }
  }
  return 0;
}

//--2d discrete wavelet transform with symmetric extension
void* Wavelet2d::dwt_2d (vector<vector<double> > &origsig, // original signal
                  int J, //number of decomposition levels
                  string nm, //wavelet name
                  vector<double> & dwt_output, //output coefficients
                  vector<double> &flag, //house keeping values for inverse
                  vector<int> &length,  //stores length of each output coefficients
                  string command,       //whether to embed or extract
                  String &message_bit,  //bits to be embedded
                  bool &full) {         //if the area is full of data
  //get rows and columns of a 2d signal and duplicate itself
  vector <vector<double> > sig = origsig;
  int rows = (int)sig.size();
  int cols = (int)sig[0].size();
  vector<vector<double> > orig_cp(rows, vector<double>(cols));
  orig_cp = sig;
  
  //check if the input level of decomposition can be proceeded
  int iteration_max = min((int) ceil(log(double(sig.size())))/log(2.0),
                           (int) ceil(log(double(sig[0].size()))/log(2.0)));
  if (iteration_max < J) {
    cout<<"Invalid decomposition level."<<endl;
    exit(1);
  }
  
  //set the size of flag and length, initialize filters
  vector<double> lp, hp, lpr, hpr;
  flag.push_back(double(J));
  length.insert(length.begin(), cols);
  length.insert(length.begin(), rows);
  
  int sum_coef = 0;
  for (int iter = 0; iter < J; iter ++) {
    //assign decomposition and restruction filter coefficients
    filter_coef(nm, lp, hp, lpr, hpr);
    
    //resize rows, cols, and length
    unsigned int lf = (unsigned int)lp.size();
    rows = (int) floor((double)(rows + lf - 1)/2);
    cols = (int) floor((double)(cols + lf - 1)/2);
    length.insert(length.begin(), cols);
    length.insert(length.begin(), rows);
    
    //initialize output coefficients:
    //cA-approximation
    //cH-horizontal
    //cV-vertical
    //cD-diagonal
    vector<vector<double> > cA(rows, vector<double>(cols));
    vector<vector<double> > cH(rows, vector<double>(cols));
    vector<vector<double> > cV(rows, vector<double>(cols));
    vector<vector<double> > cD(rows, vector<double>(cols));
    
    //apply dwt to get coefficients
    if (iter ==0) {
      dwt_2d_proc(nm, orig_cp, cA, cH, cV, cD);
    }else {
      dwt_2d_proc(nm, orig_cp, cA, cH, cV, cD);
    }
    
    //at the last level of dwt, embed the bits
    if (iter == J - 1) {
      if (command == "embed") {
        embed(cH, cV, message_bit, full);
      } else if (command == "extract") {
        extract(cH, cV, message_bit);
      }
    }
    
    //allocating coefficients
    vector<double> signal_tmp;
    orig_cp = cA;
    if (iter == J - 1) {
      for (int i = 0; i < rows; i ++) {
        for (int j = 0; j < cols; j ++) {
          double tmp = cA[i][j];
          //cout<<cA[i][j]<<endl;
          signal_tmp.push_back(tmp);
        }
      }
    }
    
    for (int i = 0; i < rows; i ++) {
      for (int j = cols; j < cols * 2; j ++) {
         double tmp = cH[i][j - cols];
         //cout<<tmp<<endl;
         signal_tmp.push_back(tmp);
       }
    }
      
    for (int i = rows; i < rows * 2; i ++) {
      for (int j = 0; j < cols; j ++) {
        double tmp = cV[i - rows][j];
        signal_tmp.push_back(tmp);
      }
    }
      
    for (int i = rows; i < rows * 2; i ++) {
      for (int j = cols; j < cols * 2 ; j ++) {
        double tmp = cD[i - rows][j - cols];
        signal_tmp.push_back(tmp);
      }
    }
    
    //appending output coefficients
    dwt_output.insert(dwt_output.begin(), signal_tmp.begin(), signal_tmp.end());
    sum_coef += 4 * rows * cols;
  }
  
  return 0;
}

//recompute the length of coefficients due to the resize of img
//after applying dwt
void* Wavelet2d::dwt_output_relen(vector<int> &length, vector<int> &length2, int J) {
  unsigned int sz=(unsigned int)length.size();
  int rows = length[sz-2];
  int cols = length[sz-1];
  for (int i =0; i < J; i++) {
    rows =(int) ceil((double) rows/ 2.0);
    cols =(int) ceil((double) cols/ 2.0);
  }
  for (int i =0; i < J + 1; i++) {
    length2.push_back(rows);
    length2.push_back(cols);
    rows = rows * 2;
    cols = cols*2;
  }
  return 0;
}

//--2d discrete wavelet transform detailed process
void* Wavelet2d::dwt_2d_proc(string name,
                              vector<vector<double> > &signal,
                              vector<vector<double> > &cLL, //approximation
                              vector<vector<double> > &cLH, //vertical details
                              vector<vector<double> > &cHL, //horizontal details
                              vector<vector<double> > &cHH) { //diagonal details
  int rows = (int)signal.size();
  int cols = (int)signal[0].size();
  int cols_lp = (int)cLL[0].size();
  int cols_hp = (int)cLL[0].size();
  
  vector<double> lp, hp, lpr, hpr;
  filter_coef(name, lp, hp, lpr, hpr);
  vector<vector<double> > lp_dn(rows, vector<double>(cols_lp));
  vector<vector<double> > hp_dn(rows, vector<double>(cols_hp));
  
  //row filtering and column downsampling in each branch
  for (int i = 0; i < rows; i ++) {
    vector<double> row_tmp, lp_output, hp_output;
    for (int j = 0; j < cols; j ++) {
      double tmp = signal[i][j];
      row_tmp.push_back(tmp);
    }
    dwt_1d(name, row_tmp, lp_output, hp_output);
    for (int j = 0; j < (int) lp_output.size(); j ++) {
      lp_dn[i][j] = lp_output[j];
      hp_dn[i][j] = hp_output[j];
    }
  }
  cols = cols_lp;
  
  //column filtering and row downsampling in low pass branch
  for (int j = 0; j < cols; j ++) {
    vector<double> row_tmp, lp_output, hp_output;
    for (int i = 0; i < rows; i ++) {
      double tmp = lp_dn[i][j];
      row_tmp.push_back(tmp);
    }
    dwt_1d(name, row_tmp, lp_output, hp_output);
    for (int i = 0; i < (int)lp_output.size(); i ++) {
      cLL[i][j] = lp_output[i];
      cLH[i][j] = hp_output[i];
    }
  }
  
  //column filtering and row downsampling in high pass branch
  for (int j = 0; j < cols; j ++) {
    vector<double> row_tmp, lp_output, hp_output;
    for (int i = 0; i < rows; i ++) {
      double tmp = hp_dn[i][j];
      row_tmp.push_back(tmp);
    }
    dwt_1d(name, row_tmp, lp_output, hp_output);
    for (int i = 0; i < (int)lp_output.size(); i ++) {
      cHL[i][j] = lp_output[i];
      cHH[i][j] = hp_output[i];
    }
  }
  return 0;
}

//--1d discrete wavelet transform
void* Wavelet2d::dwt_1d(string name, vector<double> &signal, vector<double> &cA, vector<double> &cD) {
  vector<double> lp, hp, lpr, hpr;
  filter_coef(name, lp, hp, lpr, hpr);
  int D = 2; //downsampling factor
  int lf = (int)lp.size();
  //make the signal symmetric
  symm_ext(signal, lf - 1);
  
  //low pass computation
  vector<double> cA_tmp;
  convfftm(signal, lp, cA_tmp);
  cA_tmp.erase(cA_tmp.begin(), cA_tmp.begin() + lf);
  cA_tmp.erase(cA_tmp.end() - lf + 1, cA_tmp.end());
  down_sample(cA_tmp, D, cA);
  
  //high pass computation
  vector<double> cD_tmp;
  convfftm(signal, hp, cD_tmp);
  cD_tmp.erase(cD_tmp.begin(), cD_tmp.begin() + lf);
  cD_tmp.erase(cD_tmp.end() -lf + 1, cD_tmp.end());
  down_sample(cD_tmp, D, cD);
  
  filter_coef(name, lp, hp, lpr, hpr);
  return 0;
}

//--signal down sampling
void Wavelet2d::down_sample(vector<double> &signal, int M, vector<double> &signal_d){
  int len = (int)signal.size();
  double len_n = ceil( (double) len / (double) M);
  for (int i = 0; i < (int) len_n; i++) {
    double temp = signal[i*M];
    signal_d.push_back(temp);
  }
}

//fast fourier transform convolution
double Wavelet2d::convfftm(vector<double> &a, vector<double> &b, vector<double> &c) {
  fftw_complex *inp_data, *filt_data, *inp_fft, *filt_fft, *temp_data, *temp_ifft;
  
  unsigned int sz = (unsigned int)(a.size() + b.size() - 1);
  inp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  filt_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  
  inp_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  filt_fft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  
  temp_data = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  temp_ifft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * sz );
  
  if (sz != transient_size_of_fft) {
    
    if (transient_size_of_fft != 0) {
      fftw_destroy_plan(plan_forward_inp);
      fftw_destroy_plan(plan_forward_filt);
      fftw_destroy_plan(plan_backward);
    }
    
    plan_forward_inp  = fftw_plan_dft_1d( sz, inp_data, inp_fft, FFTW_FORWARD, FFTW_MEASURE );
    plan_forward_filt  = fftw_plan_dft_1d( sz, filt_data, filt_fft, FFTW_FORWARD, FFTW_MEASURE );
    plan_backward = fftw_plan_dft_1d( sz, temp_data, temp_ifft, FFTW_BACKWARD, FFTW_MEASURE );
    transient_size_of_fft = sz;
    
  }
  
  for (unsigned int i =0; i < sz; i++) {
    if (i < a.size()) {
      inp_data[i][0] = a[i];
    } else {
      inp_data[i][0] = 0.0;
    }
    inp_data[i][1] = 0.0;
    if (i < b.size()) {
      filt_data[i][0] = b[i];
    } else {
      filt_data[i][0] = 0.0;
    }
    filt_data[i][1] = 0.0;
  }
  
  fftw_execute_dft( plan_forward_inp,inp_data, inp_fft);
  fftw_execute_dft( plan_forward_filt,filt_data, filt_fft);
  
  for (unsigned int i =0; i < sz; i++){
    temp_data[i][0] = inp_fft[i][0]*filt_fft[i][0] - inp_fft[i][1]*filt_fft[i][1];
    temp_data[i][1] = inp_fft[i][0]*filt_fft[i][1] + inp_fft[i][1]*filt_fft[i][0];
  }
  
  fftw_execute_dft( plan_backward, temp_data, temp_ifft);
  
  for (unsigned int i = 0; i < sz; i++) {
    double temp1;
    temp1 = temp_ifft[i][0] / (double) sz;
    c.push_back(temp1);
  }
  fftw_free(inp_data);
  fftw_free(filt_data);
  fftw_free(inp_fft);
  fftw_free(filt_fft);
  fftw_free(temp_data);
  fftw_free(temp_ifft);
  
  return 0;
}

//--make a signal symmetric
void* Wavelet2d::symm_ext(vector<double> &signal, int a) {
  unsigned int len = (unsigned int)signal.size();
  for (int i =0; i < a; i++) {
    double temp1= signal[i * 2];
    double temp2= signal[len - 1];
    signal.insert(signal.begin(),temp1);
    signal.insert(signal.end(),temp2);
  }
  return 0;
}

//--assign filter coefficients
int Wavelet2d::filter_coef(string name, vector<double> &lp, vector<double> &hp, vector<double> &lpr, vector<double> &hpr) {  
  if (name == "CDF") {
    double lp_a[] = {0.026748757411, -0.016864118443, -0.078223266529, 0.266864118443, 0.602949018236, 0.266864118443, -0.078223266529, -0.016864118443, 0.026748757411};
    double hp_a[] = {0.0, 0.091271763114, -0.057543526229, -0.591271763114, 1.11508705, -0.591271763114, -0.057543526229, 0.091271763114, 0.0};
    double lpr_a[] = {0.0, -0.091271763114, -0.057543526229, 0.591271763114, 1.11508705, 0.591271763114, -0.057543526229, -0.091271763114, 0.0};
    double hpr_a[] = {0.026748757411, 0.016864118443, -0.078223266529, -0.266864118443, 0.602949018236, -0.266864118443, -0.078223266529, 0.016864118443, 0.026748757411};
    lp.assign(lp_a, lp_a + sizeof(lp_a)/sizeof(double));
    hp.assign(hp_a, hp_a + sizeof(hp_a)/sizeof(double));
    lpr.assign(lpr_a, lpr_a + sizeof(lpr_a)/sizeof(double));
    hpr.assign(hpr_a, hpr_a + sizeof(hpr_a)/sizeof(double));
  }
  else if (name == "bior4.4") {
    double lp1_a[] = {0.0,
      0.03782845550726404,
      -0.023849465019556843,
      -0.11062440441843718,
      0.37740285561283066,
      0.85269867900889385,
      0.37740285561283066,
      -0.11062440441843718,
      -0.023849465019556843,
      0.03782845550726404
    };
    lp.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));
    
    double hp1_a[] = {0.0,
      -0.064538882628697058,
      0.040689417609164058,
      0.41809227322161724,
      -0.7884856164055829,
      0.41809227322161724,
      0.040689417609164058,
      -0.064538882628697058,
      0.0,
      0.0
    };
    hp.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));
    
    double lp2_a[] = {0.0,
      -0.064538882628697058,
      -0.040689417609164058,
      0.41809227322161724,
      0.7884856164055829,
      0.41809227322161724,
      -0.040689417609164058,
      -0.064538882628697058,
      0.0,
      0.0
    };
    lpr.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));
    
    double hp2_a[] = {0.0,
      -0.03782845550726404,
      -0.023849465019556843,
      0.11062440441843718,
      0.37740285561283066,
      -0.85269867900889385,
      0.37740285561283066,
      0.11062440441843718,
      -0.023849465019556843,
      -0.03782845550726404
    };
    hpr.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
  }
  else if ( name == "db3"){
    double lp1_a[] = {0.035226291882100656, -0.085441273882241486, -0.13501102001039084,
      0.45987750211933132, 0.80689150931333875, 0.33267055295095688};
    lp.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));
    double hp1_a[] = {-0.33267055295095688, 0.80689150931333875, -0.45987750211933132,
      -0.13501102001039084, 0.085441273882241486, 0.035226291882100656 };
    hp.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));
    double lp2_a[] = {0.33267055295095688, 0.80689150931333875, 0.45987750211933132,
      -0.13501102001039084, -0.085441273882241486, 0.035226291882100656 };
    lpr.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));    
    double hp2_a[] = {0.035226291882100656, 0.085441273882241486, -0.13501102001039084,
      -0.45987750211933132, 0.80689150931333875, -0.33267055295095688 };
    hpr.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
  }
  else if (name == "bior2.2") {
    double lp1_a[] = {0.0,
      -0.17677669529663689,
      0.35355339059327379,
      1.0606601717798214,
      0.35355339059327379,
      -0.17677669529663689};
    lp.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));
    double hp1_a[] = {0.0,
      -0.17677669529663689,
      0.35355339059327379,
      1.0606601717798214,
      0.35355339059327379,
      -0.17677669529663689};
    hp.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));
    double lp2_a[] = {0.0,
      -0.17677669529663689,
      0.35355339059327379,
      1.0606601717798214,
      0.35355339059327379,
      -0.17677669529663689};
    lpr.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));
    double hp2_a[] = {0.0,
      0.17677669529663689,
      0.35355339059327379,
      -1.0606601717798214,
      0.35355339059327379,
      0.17677669529663689};
    hpr.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
  }
  return 0;
}