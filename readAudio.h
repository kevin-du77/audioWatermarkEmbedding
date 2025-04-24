#pragma once	
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// ע�⣺��Ҫ��װ libsndfile ����ʵ�� WAV ��д
#include <sndfile.h>
#include <sndfile.hh> // libsndfile ͷ�ļ�

// ������Ƶ��Ϣ�ṹ��
struct AudioInfo {
    int sampleRate;     // ������
    int channels;       // ������
    int frames;         // ��֡��
    int bitsPerSample;  // ÿ��������λ��
};

pair<vector<double>, AudioInfo> readAudio(const string& filePath);

void writeAudioFile(const std::string& filePath, const std::vector<double>& audioData, AudioInfo info);