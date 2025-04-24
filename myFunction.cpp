#include "myFunction.h"

// 辅助函数：将数字四舍五入到指定位数（例如：roundn(lamda, -2)）
double roundn(double value, int n) {
	double factor = pow(10.0, n);
	return round(value * factor) / factor;
}

//打印数组到excel文件
void printToExcel(const vector<double>& data) {
	//ofstream file(filename);  , const string& filename
	ofstream outfile("output.xls");
	//if (!file.is_open()) {
	//	cerr << "无法打开文件: " << filename << endl;
	//	return;
	//}
	for (size_t i = 0; i < data.size(); i++) {
		outfile << data[i] << endl;
	}
	outfile.close();
}

void printAudioData(const vector<double>& data, int start)
{
	cout << "Audio Data (10 samples) from : " << start << " to " << start+10;
	int printBegin = start;//打印起始点
	for (size_t i = printBegin; i < printBegin + 10 && i < data.size(); ++i) {
		cout << data[i] << " ";
	}
	cout << endl;
}

// 数组初始化
void initArrayDouble(vector<double>& data, int size, double value)
{
	data.resize(size);
	for (int i = 0; i < size; ++i) {
		data[i] = value;
	}
}

void initArrayInt(vector<int>& data, int size, int value)
{
	data.resize(size);
	for (int i = 0; i < size; ++i) {
		data[i] = value;
	}
}
// 二维数组初始化
void initArrayInt2D(vector<vector<int>>& data, int rows, int cols, int value) 
{
	data.resize(rows);
	for (int i = 0; i < rows; ++i) {
		data[i].resize(cols, value);
	}
}

void initArrayDouble2D(vector<vector<double>>& data, int rows, int cols, double value)
{
	data.resize(rows);
	for (int i = 0; i < rows; ++i) {
		data[i].resize(cols, value);
	}
}
