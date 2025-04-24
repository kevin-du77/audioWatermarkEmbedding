#pragma once
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// 辅助函数：将数字四舍五入到指定位数（例如：roundn(lamda, -2)）
double roundn(double value, int n);

//打印数组到excel文件
void printToExcel(const vector<double>& data);

// 打印部分音频数据（前 10 个样本）
void printAudioData(const vector<double>& data, int start);

// 数组初始化
void initArrayDouble(vector<double>& data, int size, double value);

void initArrayInt(vector<int>& data, int size, int value);

// 二维数组初始化
void initArrayInt2D(vector<vector<int>>& data, int rows, int cols, int value);

void initArrayDouble2D(vector<vector<double>>& data, int rows, int cols, double value);
