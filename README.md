StegoQR
=======

New generation of QR code based on Steganography
## The Phenomena
QR code (Quick Response Code) is the trademark for a type of matrix barcode first designed for automotive industry in Japan. It is read by an imaging device such as a camera and interpreted algorithmically by underlying software based on patterns present in both horizontal and vertical components of the image. Due to the availability of “Smart phones”, QR code has been used as an approach to pass text information using still images in various areas such as promotion, events, and entertainment. However, current QR code has an obvious drawback of not being aesthetically pleasing because it is monochrome and does not make much sense to human beings. Inspired by steganography, this project is aimed to make beautiful QR code with an approach to embed text message into an arbitrary image without making noticeable change to the appearance of the image. 

## The Definition
Steganography provides an approach to hide information into media such as audio, image, and video. Combined with encryption, it proves to be even secured. Steganography is generally defined as the art and science of communicating in a covert fashion. Differs from cryptography, which changes plaintext to cipher text (not understandable to eavesdroppers), steganography techniques try to hide the very existence of the message itself. In digital steganography, the original image to be encoded is named cover-image and the image after encoding is named stego-image. One typical implementation described in is to embed information into a cover-image by modifying values of the Least Significant Bits (LSB) of pixels in pixel domain.

##The Goal
Information can be hidden into an image without being noticed by human eyes. However, the above scheme is vulnerable to the following cases:
* As the above approach involves color value modification in pixel domain, when decoding, it strongly relies on the accuracy of pixel sampling of the stego-images. When the points that data are embedded into may be within noise when being extracted, data extraction could be inaccurate.
* The above approach is color-altering intolerant. It requires lossless data preservation during image transfer from the encoder to the decoder. As a result, when the stego-image at the receiving end may be compressed, filtered or blurred, data extraction could be inaccurate.
* The steganography is not synchronized, that is, extracting data from the stego-image above needs the cover-image.

In order to overcome the above limitations, the proposed scheme will hide information into regions of a cover-image instead of each pixel. More specifically, information will be hidden in the cover-image in Discrete Wavelet Transform (DWF) domain in a content-based manner due to its advantage of being able to capture both frequency and location information (prior to Fourier transform) and clearly partition the high-frequency and low-frequency information on a pixel-by-pixel basis. In addition, the technique of Speed-Up Robust Features (SURF) will be applied to obtain the regions of interests generically. Therefore, steganography synchronization will be achieved so that data extraction will allow lossy compression and noise; the data extraction process will also be completely blind, since only the stega-image is required.

##The Added values
Generally, it would be evolutionary to redefine the current QR code by hugely improving its appearance while keeping its featuring functionalities. Technically, it is challenging to encode hidden information into an image without noticeably changing its appearance. In addition, enabling decoding data-loss tolerant is unprecedented. As a continuing project, I would also like to transplant the decoder to a mobile device platform such as Android and iOS once this project has been completed.

##Algorithm keywords:
* SURF Detection
  * Keypoints elimination

![SURF](http://www.mftp.info/20140202/1392084968x1927178161.png)
* Pre-DWT Transformation
  * Extension Techniques

![extension](http://www.mftp.info/20140202/1392085116x1927178161.png)
  * Biorthogonal and orthogonal
  * Choice of wavelet (BIOR 4.4, Cohen-Daubechies-Feauveau 9/7 wavelet)
* Data Embedding
  * DWT Transform

![dwtG](http://www.mftp.info/20140202/1392085269x1927178161.png)
![dwt](http://www.mftp.info/20140202/1392085321x1927178161.png)
![dwt](http://www.mftp.info/20140202/1392085346x1927178161.png)

```CPP
//--2d discrete wavelet transform with symmetric extension
void* Wavelet2d::dwt_2d (vector<vector<double> > &signal, // original signal
                  int J, //number of decomposition levels
                  string name, //wavelet name
                  vector<double> & dwt_output, //output coefficients 
                  string command, //whether to embed or extract 
                  String &message_bit) { //bits to be embedded
  // initialize filters
  vector<double> lp, hp, lpr, hpr;
  for (int iter = 0; iter < J; iter ++) {
    //assign decomposition and restruction filter coefficients
    filter_coef(name, lp, hp, lpr, hpr); 
    //initialize output coefficients: 
    //cA-approximation
    //cH-horizontal
    //cV-vertical
    //cD-diagonal
    vector<vector<double> > cA(rows, vector<double>(cols)); 
    vector<vector<double> > cH(rows, vector<double>(cols)); 
    vector<vector<double> > cV(rows, vector<double>(cols)); 
    vector<vector<double> > cD(rows, vector<double>(cols)); 
    //convert 2d signal to 1d signal
    for (int i = 0; i < rows; i ++) { 
      vector<double> signal_1d;
      for (int j = 0; j < cols; j ++) {
        double tmp = signal[i][j];
        signal_1d.push_back(tmp); 
      }
    }
    //apply dwt to get coefficients
    vector<vector<double> > lp_dn(rows, vector<double>(cols_lp)); 
    vector<vector<double> > hp_dn(rows, vector<double>(cols_hp)); 
    //row filtering and column down sampling in each branch/
    dwt_1d(name, signal_1d, lp_output, hp_output); //1-D DWT 
    for (int j = 0; j < (int) lp_output.size(); j ++) {
      lp_dn[i][j] = lp_output[j];
      hp_dn[i][j] = hp_output[j]; 
    }
    //column filtering and row down sampling in low pass branch/ 
    dwt_1d(name, col_lp_dn, lp_output, hp_output); //1-D DWT
    for (int j = 0; j < (int) lp_output.size(); j ++) { 
      cA[i][j] = lp_output[j];
      cV[i][j] = hp_output[j];
    }
    //column filtering and row down sampling in low pass branch/ 
    dwt_1d(name, col_hp_dn, lp_output, hp_output); //1-D DWT
    for (int j = 0; j < (int) lp_output.size(); j ++) { 
      cH[i][j] = lp_output[j];
      cD[i][j] = hp_output[j];
    }
    //at the last level of dwt, embed the bits
    if (iter == J - 1) {
      if (command == "embed") {
        embed(cH, cV, message_bit); //pay load bit embedding as discussed before
      } else if (command == "extract") {
        extract(cH, cV, message_bit); //pay load bit extraction to be explained
      } 
    }
    //allocating coefficients into the output signal, which is the frequency domain of the 
    //original image or region
    dwt_output = alloc_coeff(cA, cH, cV, cD); 
  }
}
```
  * Inverse DWT Transform
![idwt](http://www.mftp.info/20140202/1392085386x1927178161.png)
* Data Extraction
