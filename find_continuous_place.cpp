#include <vector>
#include <iostream>
#include "find_continuous_place.h"


// 函数定义
void find_continuous_place(const vector<int>& inputIndex, const vector<double>& inputValue, vector<int>& place, vector<double>& value) {
    // 检查输入是否为空
    if (inputIndex.empty() || inputValue.empty() || inputIndex.size() != inputValue.size()) {
        printf("查找连续索引失败\n");
        place.clear();
        value.clear();
        return;
    }

    // 计算逻辑值：相邻元素差值是否为1
    vector<int> logicValue(inputIndex.size() - 1);
    for (size_t i = 0; i < inputIndex.size() - 1; ++i) {
        logicValue[i] = (inputIndex[i + 1] - inputIndex[i] == 1) ? 1 : 0;
    }

    // 查找连续三个1的位置
    vector<int> continueList;
    for (size_t i = 0; i < logicValue.size() - 2; ++i) {
        if (logicValue[i] == 1 && logicValue[i + 1] == 1 && logicValue[i + 2] == 1) {
            continueList.push_back(i);//可以使用 push_back 方法向 vector 中添加元素：
        }
    }

    // 如果没有找到连续模式，返回空结果
    if (continueList.empty()) {
        place.clear();
        value.clear();
        return;
    }

    // 初始化结果
    place.clear();
    value.clear();
    size_t temp = 0;

    // 提取第一个连续段
    for (int j = 0; j < 4; ++j) {
        place.push_back(inputIndex[continueList[0] + j]);
        value.push_back(inputValue[continueList[0] + j]);
    }

    // 遍历 continueList，提取其他连续段
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