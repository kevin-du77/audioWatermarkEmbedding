#pragma once	
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// 注意：需要安装 libsndfile 库来实现 WAV 读写
#include <sndfile.h>
#include <sndfile.hh> // libsndfile 头文件

// 定义音频信息结构体
struct AudioInfo {
    int sampleRate;     // 采样率
    int channels;       // 声道数
    int frames;         // 总帧数
    int bitsPerSample;  // 每个样本的位数
};

pair<vector<double>, AudioInfo> readAudio(const string& filePath);

void writeAudioFile(const std::string& filePath, const std::vector<double>& audioData, AudioInfo info);