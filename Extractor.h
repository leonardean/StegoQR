//
//  Extractor.h
//  QR_Code
//
//  Created by leonardo xu on 25/09/2013.
//  Copyright (c) 2013 leonardo xu. All rights reserved.
//

#ifndef __QR_Code__Extractor__
#define __QR_Code__Extractor__

#include "Embeder.h"
class Extractor {
public:
  void* data_extract(char* file_name, int J, string &wavelet_name, vector<double> &output, vector<double> &flag, vector<int> &length, int &rows, int &cols, string &message);
};

#endif /* defined(__QR_Code__Extractor__) */
