//
//  main.cpp
//  QR_Code
//
//  Created by leonardo xu on 21/08/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//
#include "Detector.h"
#include "Wavelet2d.h"
#include "Embeder.h"
#include "Extractor.h"
#include <string.h>


 int main( int argc, char** argv )
 {
   if (argc != 5) {
     cout<<"Need 4 arguments:\nembed or extract\nPath to cover image\nPath to destination stego image\nPayload bits to be embedded"<<endl;
     exit(1);
   }
   Detector detector;
   detector.surf_detect(argv[2]);
   
   //initialize variables for embeder and extractor
   Embeder embeder;
   Extractor extractor;
   string wavelet_name = "bior4.4";
   vector<double> dwt_output, flag;
   vector<int> length;
   int rows, cols;
   Mat mat_G, mat_B;
   bool full;
   int J = 1;
   string message = argv[4];
   
   if (strncmp(argv[1],"embed", 2) == 0) {
     //data embeding
     embeder.data_embed(argv[2], J, wavelet_name, dwt_output, flag, length, rows, cols, full, message, argv[3]);
   }
   if (strncmp(argv[1],"extract", 2) == 0) {
     //data extracting
     extractor.data_extract(argv[3], J, wavelet_name, dwt_output, flag, length, rows, cols, message);
   }
   
   return 0;
}



