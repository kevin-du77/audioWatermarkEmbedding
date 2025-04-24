#include <vector>
#include <iostream>
#include "find_continuous_place.h"


// ��������
void find_continuous_place(const vector<int>& inputIndex, const vector<double>& inputValue, vector<int>& place, vector<double>& value) {
    // ��������Ƿ�Ϊ��
    if (inputIndex.empty() || inputValue.empty() || inputIndex.size() != inputValue.size()) {
        printf("������������ʧ��\n");
        place.clear();
        value.clear();
        return;
    }

    // �����߼�ֵ������Ԫ�ز�ֵ�Ƿ�Ϊ1
    vector<int> logicValue(inputIndex.size() - 1);
    for (size_t i = 0; i < inputIndex.size() - 1; ++i) {
        logicValue[i] = (inputIndex[i + 1] - inputIndex[i] == 1) ? 1 : 0;
    }

    // ������������1��λ��
    vector<int> continueList;
    for (size_t i = 0; i < logicValue.size() - 2; ++i) {
        if (logicValue[i] == 1 && logicValue[i + 1] == 1 && logicValue[i + 2] == 1) {
            continueList.push_back(i);//����ʹ�� push_back ������ vector �����Ԫ�أ�
        }
    }

    // ���û���ҵ�����ģʽ�����ؿս��
    if (continueList.empty()) {
        place.clear();
        value.clear();
        return;
    }

    // ��ʼ�����
    place.clear();
    value.clear();
    size_t temp = 0;

    // ��ȡ��һ��������
    for (int j = 0; j < 4; ++j) {
        place.push_back(inputIndex[continueList[0] + j]);
        value.push_back(inputValue[continueList[0] + j]);
    }

    // ���� continueList����ȡ����������
    for (size_t i = 0; i < continueList.size() - 1; ++i) {
        if (continueList[i + 1] - continueList[temp] >= 4) {
            temp = i + 1;
            for (int j = 0; j < 4; ++j) {
                place.push_back(inputIndex[continueList[temp] + j]);
                value.push_back(inputValue[continueList[temp] + j]);
            }
        }
    }
	return;
}