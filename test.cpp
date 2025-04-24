#include <iostream>
#include <vector>
#include "find_continuous_place.h"
#include "readAudio.h"

using namespace std;

int main()
{
	// ���� find_continuous_place ����
	vector<int> inputIndex = { 1, 2, 3, 4, 5, 6, 7, 201, 202, 303 };
	vector<double> inputValue = { 11, 12, 13, -14, -15, -16, 17, -18, -19, 20 };
	vector<int> place;
	vector<double> value;

	find_continuous_place(inputIndex, inputValue, place, value);

	// ������
	cout << "Place: ";
	for (const auto& p : place) {
		cout << p << " ";
	}
	cout << endl;

	cout << "Value: ";
	for (const auto& v : value) {
		cout << v << " ";
	}
	cout << endl;

	// --------------------------2. ��ȡ��Ƶ�ļ�---------------------------
	// 
	string audioFile = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1.wav";
	//vector<double> data;
	int Fs, bits;
	pair<vector<double>, AudioInfo> result = readAudio(audioFile); // ��ȡ��Ƶ�ļ�
	vector<double> data = result.first; // ��ȡ��Ƶ����
	AudioInfo info = result.second;     // ��ȡ��Ƶ��Ϣ
	// ��ӡ��Ƶ��Ϣ
	cout << "Sample Rate: " << info.sampleRate << " Hz\n";
	cout << "Channels: " << info.channels << "\n";
	cout << "Frames: " << info.frames << "\n";
	cout << "Bits Per Sample: " << info.bitsPerSample << "\n";

	// ��ӡ������Ƶ���ݣ�ǰ 10 ��������
	cout << "Audio Data (first 10 samples): ";
	int printBegin = 300000-1;//��ӡ��ʼ��
	for (size_t i = printBegin; i < printBegin + 10 && i < data.size(); ++i) {
		cout << data[i] << " ";
	}
	cout << endl;

	vector<vector<int>> joint(100); // �洢����λ��
	// ��ӡ������Ƶ���ݣ�ǰ 10 ��������
	cout << "Audio Data (first 10 samples): ";
	int printBegin = 1 - 1;//��ӡ��ʼ��
	for (size_t i = printBegin; i < printBegin + 10 && i < joint[1].size(); ++i) {
		cout << joint[1][i] << " ";
	}
	cout << endl;

	//��ӡa2_bin��d2_bin�ĵ�һ�е�ȫ������
	//for (size_t i = 0; i < a2_bin.size(); ++i) {
	//	cout << a2_bin[i].size() << " ";
	//}
	//cout << endl;
	for (size_t i = 0; i < a2_bin[0].size(); ++i) {
		cout << a2_bin[0][i] << " ";
	}
	cout << endl;

	// ��ӡ d2_bin �ĵ�һ�е�ȫ������
	//for (size_t i = 0; i < d2_bin.size(); ++i) {
	//	cout << d2_bin[i].size() << " ";
	//}
	//cout << endl;
	for (size_t i = 0; i < d2_bin[0].size(); ++i) {
		cout << d2_bin[0][i] << " ";
	}
	cout << endl;

	return 0;
}