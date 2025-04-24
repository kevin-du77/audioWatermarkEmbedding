#include <iostream>
#include <vector>
#include "find_continuous_place.h"
#include "readAudio.h"

using namespace std;

int main()
{
	// 测试 find_continuous_place 函数
	vector<int> inputIndex = { 1, 2, 3, 4, 5, 6, 7, 201, 202, 303 };
	vector<double> inputValue = { 11, 12, 13, -14, -15, -16, 17, -18, -19, 20 };
	vector<int> place;
	vector<double> value;

	find_continuous_place(inputIndex, inputValue, place, value);

	// 输出结果
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

	// --------------------------2. 读取音频文件---------------------------
	// 
	string audioFile = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1.wav";
	//vector<double> data;
	int Fs, bits;
	pair<vector<double>, AudioInfo> result = readAudio(audioFile); // 读取音频文件
	vector<double> data = result.first; // 获取音频数据
	AudioInfo info = result.second;     // 获取音频信息
	// 打印音频信息
	cout << "Sample Rate: " << info.sampleRate << " Hz\n";
	cout << "Channels: " << info.channels << "\n";
	cout << "Frames: " << info.frames << "\n";
	cout << "Bits Per Sample: " << info.bitsPerSample << "\n";

	// 打印部分音频数据（前 10 个样本）
	cout << "Audio Data (first 10 samples): ";
	int printBegin = 300000-1;//打印起始点
	for (size_t i = printBegin; i < printBegin + 10 && i < data.size(); ++i) {
		cout << data[i] << " ";
	}
	cout << endl;

	vector<vector<int>> joint(100); // 存储连续位置
	// 打印部分音频数据（前 10 个样本）
	cout << "Audio Data (first 10 samples): ";
	int printBegin = 1 - 1;//打印起始点
	for (size_t i = printBegin; i < printBegin + 10 && i < joint[1].size(); ++i) {
		cout << joint[1][i] << " ";
	}
	cout << endl;

	//打印a2_bin和d2_bin的第一行的全部数据
	//for (size_t i = 0; i < a2_bin.size(); ++i) {
	//	cout << a2_bin[i].size() << " ";
	//}
	//cout << endl;
	for (size_t i = 0; i < a2_bin[0].size(); ++i) {
		cout << a2_bin[0][i] << " ";
	}
	cout << endl;

	// 打印 d2_bin 的第一行的全部数据
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