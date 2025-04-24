//#include <vector>
//#include <fstream>
//#include <iostream>
//#include "dwt.h"
//
//using namespace std;
//
//// Haar 小波变换
//void haar_wavelet_transform(const vector<double>& signal, vector<double>& approx, vector<double>& detail) {
//    int n = signal.size();
//    int half = n / 2;
//
//    approx.resize(half);
//    detail.resize(half);
//
//    for (int i = 0; i < half; ++i) {
//        approx[i] = (signal[2 * i] + signal[2 * i + 1]) / sqrt(2.0); // 近似系数
//        detail[i] = (signal[2 * i] - signal[2 * i + 1]) / sqrt(2.0); // 细节系数
//    }
//}