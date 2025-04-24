#pragma once
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// �����������������������뵽ָ��λ�������磺roundn(lamda, -2)��
double roundn(double value, int n);

//��ӡ���鵽excel�ļ�
void printToExcel(const vector<double>& data);

// ��ӡ������Ƶ���ݣ�ǰ 10 ��������
void printAudioData(const vector<double>& data, int start);

// �����ʼ��
void initArrayDouble(vector<double>& data, int size, double value);

void initArrayInt(vector<int>& data, int size, int value);

// ��ά�����ʼ��
void initArrayInt2D(vector<vector<int>>& data, int rows, int cols, int value);

void initArrayDouble2D(vector<vector<double>>& data, int rows, int cols, double value);
