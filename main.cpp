#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "readAudio.h"  // 读取音频文件的头文件
#include <cmath>        // 用于 abs 和 floor
#include <numeric>      // 用于 accumulate
#include "myFunction.h" // 辅助函数的头文件
#include "find_continuous_place.h" // 查找连续位置的头文件
#include <bitset>       // 用于 bitset
//#include "dwt.h"        // Haar 小波变换的头文件

using namespace std;

// --------------------------1. 参数定义和水印定义--------------------------
// 定义水印向量 w（取值为 1 或 -1）
vector<int> w = { 1, 1, 1, 1, 1, 1, -1, -1, -1, 1,
                -1, 1 ,1, 1, -1, 1, 1, -1, 1, 1,
                1, -1, 1, 1, -1, -1, 1, 1,1, 1,
                1, -1, -1, -1,  1,  1, -1,  1, 1,  1,
                -1,  1,  1, -1,  1, -1, -1, -1, -1,  1,
                -1, -1,  1, -1,  1,  1, -1, -1,  1,  1 };
int Lw = w.size();

double lamda_unround = 4.0;  // lamda 越大 bin_width 越大
double lamda = roundn(lamda_unround, 2);
double T = 4.0;
int nBegin = 1;
int nBins = 3 * Lw + 3;
int nEnd = nBins;
int nBins_number = 10;
int nBinsFinal = nBins + nBins_number;// 统计 bin 的数量，实验中数值为 183 + 10 = 193

// Haar 小波变换
void haar_wavelet_transform(const vector<double>& signal, vector<double>& approx, vector<double>& detail) {
    int n = signal.size();
    int half = n / 2;

    approx.resize(half);
    detail.resize(half);

    for (int i = 0; i < half; ++i) {
        approx[i] = (signal[2 * i] + signal[2 * i + 1]) / sqrt(2.0); // 近似系数
        detail[i] = (signal[2 * i] - signal[2 * i + 1]) / sqrt(2.0); // 细节系数
    }
}

// Haar 小波变换的二阶变换
//void haar_wavelet_transform2D(int nBinsFinal, vector <vector< double >> & value_num_joint, vector <vector< double >>& dwt_a2)
//{
//    vector<vector<double>> ca1_bin(nBinsFinal); // 初始化 a1_bin
//    //vector<vector<double>> ca2_bin(nBinsFinal); // 初始化 a2_bin
//    vector<vector<double>> d1_bin(nBinsFinal);  // 初始化 d1_bin
//    vector<vector<double>> d2_bin(nBinsFinal);  // 初始化 d2_bin
//    for (size_t i = 0; i < value_num_joint.size(); ++i) {
//        if (!value_num_joint[i].empty()) {
//            vector<double> dwt_value = value_num_joint[i];
//            //vector<double> a1, d1, a2, d2;
//            haar_wavelet_transform(dwt_value, ca1_bin[i], d1_bin[i]); // 近似系数
//            haar_wavelet_transform(ca1_bin[i], dwt_a2[i], d2_bin[i]); // 2阶
//
//        }
//    }
//}

//// 统计各个bins中的样本数量
//pair< vector<vector<int>>, vector<vector<double>> > count_bins(const vector<double>& data, double dRange, vector<double> B1, vector<int> index)
//{
//    double bin_width = floor(2 * dRange / nBins);
//    size_t L_B1 = B1.size();
//    // 偏移 data 和 B1 到正数范围
//    vector<double> data1(data.size());
//    for (size_t i = 0; i < data.size(); ++i) {
//        data1[i] = data[i] + dRange + 1;
//    }
//
//    vector<double> S1(B1.size());
//    for (size_t i = 0; i < B1.size(); ++i) {
//        S1[i] = B1[i] + dRange + 1;
//    }
//
//    // 将 nBins_number 转成二进制 6 位
//    bitset<6> nBins_number_bin(nBins_number);
//
//    // 初始化每个 bin 的样本数量
//    vector<int> nNum(nBinsFinal, 0);
//
//    // 表示是第几个 bin 里面，进行了直方图统计?
//    for (int i = 0; i < L_B1; ++i) {
//        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
//        if (binIndex >= 0 && binIndex < nNum.size()) {
//            nNum[binIndex] += 1;
//        }
//    }
//
//    vector<int> lo(nNum.size(), 1); // 初始化 lo 为全 1
//    // 初始化 value_num 和 index_num
//    vector<vector<double>> value_num(nBinsFinal);
//    vector<vector<int>> index_num(nBinsFinal);
//
//    for (int i = 0; i < nBins; ++i) {
//        value_num[i] = vector<double>(nNum[i], 0.0); // 初始化为 0.0
//        index_num[i] = vector<int>(nNum[i], 0);      // 初始化为 0
//    }
//
//    // 这里的 nNum[i] 是每个 bin 中的样本数量
//    // value_num记录样本值，存在第(S1(i) / bin_width)个bin的第k个位置，k为lo对应bin记录的数，就是这个bin有了多少个数据，一个萝卜一个坑
//    // 填充 value_num 和 index_num
//    for (int i = 0; i < L_B1; ++i) {
//        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
//        if (binIndex >= 0 && binIndex < nBins) {
//            value_num[binIndex][lo[binIndex] - 1] = S1[i];
//            index_num[binIndex][lo[binIndex] - 1] = index[i];
//            lo[binIndex] += 1; // 更新 lo
//        }
//    }
//
//    //记录每个bin里面连续的有哪些变量
//    vector<vector<int>> joint(nBinsFinal);              // 存储连续位置
//    vector<vector<double>> value_num_joint(nBinsFinal); // 存储连续值
//
//    // 遍历 index_num
//    for (int i = 0; i < index_num.size(); ++i) {
//        if (index_num[i].size() >= 4) {
//            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
//        }
//    }
//    return { joint, value_num_joint };
//}

// 实现逆小波变换 (idwt) 使用 'db1' 小波基
vector<double> idwt(const vector<double>& a, const vector<double>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Input vectors 'a' and 'b' must have the same length.");
    }

    size_t n = a.size();
    vector<double> result(2 * n);

    // db1 小波基的逆变换滤波器
    const double h0 = 1.0 / sqrt(2.0); // 低通滤波器系数
    const double h1 = 1.0 / sqrt(2.0); // 高通滤波器系数

    for (size_t i = 0; i < n; ++i) {
        result[2 * i] = h0 * a[i] + h1 * b[i]; // 偶数索引
        result[2 * i + 1] = h0 * a[i] - h1 * b[i]; // 奇数索引
    }

    return result;
}

// 提取水印
void extractWatermark(const vector<vector<double>>& ca2_binw,  vector<int>& w1) 
{
    vector<double> dRw(Lw); // 存储重构后的比例
    w1.resize(Lw);          // 存储提取的水印

    for (int j = 0; j < Lw; ++j) {
        int nIndex;
        if (j % 2 == 0) {
			nIndex = nBegin + 3 * ((j + 1 + 1) / 2 - 1) - 1;// MATLAB 索引从 1 开始
        }
        else {
            nIndex = nEnd - 3 * ((j + 1) / 2) + 1 - 1;
        }

        int a = ca2_binw[nIndex].size();
        int b = ca2_binw[nIndex + 1].size();
        int c = ca2_binw[nIndex + 2].size();

        dRw[j] = (2.0 * b) / (a + c);

        if (dRw[j] >= 1) {
            w1[j] = 1;
        }
        else {
            w1[j] = -1;
        }
    }
    //// 输出 dRw 和 w1
    //cout << "dRw: ";
    //for (double val : dRw) {
    //    cout << val << " ";
    //}
    //cout << endl;

    //cout << "Extracted watermark (w1): ";
    //for (int val : w1) {
    //    cout << val << " ";
    //}
    //cout << endl;

    //ofstream outfile("output.xls");
    ////if (!file.is_open()) {
    ////	cerr << "无法打开文件: " << filename << endl;
    ////	return;
    ////}
    //for (size_t i = 0; i < w1.size(); i++) {
    //    outfile << w1[i] << endl;
    //}
    //outfile.close();
    //cout << "输出表格成功！" << endl;

    cout << endl;
}

double calculateBER(const vector<int>& w1, const vector<int>& w) {
    
    double dBer = 0;

    // 遍历每个水印位，计算错误位数
    for (int j = 0; j < Lw; ++j) {
        if (w1[j] != w[j]) {
            dBer++;
        }
    }

    // 计算 BER (比特错误率)
    return dBer / Lw;
}

int main()
{

    // --------------------------2. 读取音频文件--------------------------
	//
    string audioFile = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1.wav";
    pair<vector<double>, AudioInfo> result = readAudio(audioFile); // 读取音频文件
    vector<double> data = result.first; // 获取音频数据
    AudioInfo info = result.second;     // 获取音频信息

    int NA = 1 << (info.bitsPerSample - 1);  // 2^(bits-1)
    //将数据放大到整数范围
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= NA;
		//printf("%f ", data[i]);
    }
    vector<double> data_original = data;
    int L_data_original = data.size();

    //--------------------------3.计算音频的幅值均值A--------------------------
	double A = 0;
	for (size_t i = 0; i < L_data_original; i++) {
		A += abs(data[i]);
	}
	A /= L_data_original;
    // 计算分箱的覆盖范围和宽度
    double dRange = floor(lamda * A);
    double bin_width = floor(2 * dRange / nBins);

    // --------------------------4. 计算分箱的索引--------------------------

    // 筛选数据并记录索引
    vector<double> B1; // 存储符合条件的音频数据
    vector<int> index; // 存储对应的索引
    for (size_t i = 0; i < data.size(); ++i) {
        if (round(data[i]) >= round(-dRange) && round(data[i]) <= round(dRange)) {
            B1.push_back(data[i]);
            index.push_back(i + 1); // MATLAB 索引从 1 开始
        }
    }
    // 计算 B1 的长度
    

    //--------------------------5.统计各个bins中的样本数量--------------------------
	//pair< vector<vector<int>>, vector<vector<double>> > count_bins_result = count_bins(data, dRange, B1, index);
	//vector<vector<int>> joint = count_bins_result.first;              // 存储连续位置
	//vector<vector<double>> value_num_joint = count_bins_result.second; // 存储连续值

    size_t L_B1 = B1.size();
    // 偏移 data 和 B1 到正数范围
    vector<double> data1(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        data1[i] = data[i] + dRange + 1;
    }

    vector<double> S1(B1.size());
    for (size_t i = 0; i < B1.size(); ++i) {
        S1[i] = B1[i] + dRange + 1;
    }

    // 将 nBins_number 转成二进制 6 位
    bitset<6> nBins_number_bin(nBins_number);

    // 初始化每个 bin 的样本数量
    vector<int> nNum(nBinsFinal, 0);

    // 表示是第几个 bin 里面，进行了直方图统计?
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
        if (binIndex >= 0 && binIndex < nNum.size()) {
            nNum[binIndex] += 1;
        }
    }

    vector<int> lo(nNum.size(), 1); // 初始化 lo 为全 1
    // 初始化 value_num 和 index_num
    vector<vector<double>> value_num(nBinsFinal);
    vector<vector<int>> index_num(nBinsFinal);

    for (int i = 0; i < nBins; ++i) {
        value_num[i] = vector<double>(nNum[i], 0.0); // 初始化为 0.0
        index_num[i] = vector<int>(nNum[i], 0);      // 初始化为 0
    }

    // 这里的 nNum[i] 是每个 bin 中的样本数量
    // value_num记录样本值，存在第(S1(i) / bin_width)个bin的第k个位置，k为lo对应bin记录的数，就是这个bin有了多少个数据，一个萝卜一个坑
    // 填充 value_num 和 index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
        if (binIndex >= 0 && binIndex < nBins) {
            value_num[binIndex][lo[binIndex] - 1] = S1[i];
            index_num[binIndex][lo[binIndex] - 1] = index[i];
            lo[binIndex] += 1; // 更新 lo
        }
    }

    //记录每个bin里面连续的有哪些变量
    vector<vector<int>> joint(nBinsFinal);              // 存储连续位置
    vector<vector<double>> value_num_joint(nBinsFinal); // 存储连续值

    // 遍历 index_num
    for (int i = 0; i < index_num.size(); ++i) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
        }
    }


	// --------------------------6. 计算 Haar 小波变换--------------------------
	//vector<vector<double>> dwt_a2(nBinsFinal);
	//haar_wavelet_transform2D(nBinsFinal, value_num_joint, dwt_a2);

    vector<vector<double>> ca1_bin(nBinsFinal); // 初始化 a1_bin
    vector<vector<double>> dwt_a2(nBinsFinal); // 初始化 a2_bin
    vector<vector<double>> d1_bin(nBinsFinal);  // 初始化 d1_bin
    vector<vector<double>> d2_bin(nBinsFinal);  // 初始化 d2_bin
    for (size_t i = 0; i < value_num_joint.size(); ++i) {
        if (!value_num_joint[i].empty()) {
            vector<double> dwt_value = value_num_joint[i];
            //vector<double> a1, d1, a2, d2;
            haar_wavelet_transform(dwt_value, ca1_bin[i], d1_bin[i]); // 近似系数
            haar_wavelet_transform(ca1_bin[i], dwt_a2[i], d2_bin[i]); // 2阶

        }
    }
    //  dwt_a2四舍五入
    for (size_t i = 0; i < dwt_a2.size(); ++i)
    {
        for (size_t j = 0; j < dwt_a2[i].size(); ++j)
        {
            dwt_a2[i][j] = roundn(dwt_a2[i][j], 1);
        }
    }

	//--------------------------7.根据阈值 T 进行水印嵌入--------------------------
    //嵌入前的数据，用于比较水印嵌入前后的差异
	vector<double> dR;//这里不要声明大小，否则会从最后一个开始赋值
    int a_w[60], b_w[60], c_w[60];
    for (int j = 0; j < Lw; ++j) {
        static int nIndex = 0;
        
        if (j % 2 == 0) {
            nIndex = nBegin + 3 * ((j+1 + 1) / 2 - 1);
        }
        else {
            nIndex = nEnd - 3 * ((j + 1) / 2) + 1;//逆天j+1忘记加括号
        }
		nIndex -= 1; // MATLAB 索引从 1 开始
		//cout << "nIndex of " << j << " : " << nIndex << endl;

        int a = dwt_a2[nIndex].size();
        int b = dwt_a2[nIndex + 1].size();
        int c = dwt_a2[nIndex + 2].size();
        dR.push_back((2.0 * b) / (a + c));//dR[j-1] = (2.0 * b) / (a + c);

		a_w[j] = a;
		b_w[j] = b;
		c_w[j] = c;

        if (w[j] == 1) {
            if ((2.0 * b) / (a + c) <= T) {
                int T0 = ceil(((a + c) * T - 2.0 * b) / (2.0 + T));
                int Ta = T0 * a / (a + c);
                int Tc = T0 * c / (a + c);
                for (int i = 0; i < Ta; ++i) {
                    dwt_a2[nIndex][i] += 2 * bin_width;
                }
                for (int i = 0; i < Tc; ++i) {
                    dwt_a2[nIndex + 2][i] -= 2 * bin_width;
                }
            }
        }
        else {
            if ((2.0 * b) / (a + c) >= 1.0 / T) {
                int T0 = ceil((2.0 * b * T - (a + c)) / (1.0 + 2.0 * T));
                int Ta = ceil(T0 * a / (a + c));
                int Tc = ceil(T0 * c / (a + c));
                for (int i = 0; i < Ta; ++i) {
                    dwt_a2[nIndex + 1][i] -= 2 * bin_width;
                }
                for (int i = Ta; i < Tc; ++i) {
                    dwt_a2[nIndex + 1][i] += 2 * bin_width;
                }
            }
        }
    }

 //   //输出dwt_a2每行的列数到excel
	//ofstream outfile("dwt_a2.xls");
	////if (!file.is_open()) {
	////	cerr << "无法打开文件: " << filename << endl;
	////	return;
	////}
	//for (size_t i = 0; i < dwt_a2.size(); i++) {
	//	outfile << dwt_a2[i].size() << endl;
	//}
	//outfile.close();
	//cout << "输出dwt_a2表格成功！" << endl;


	// --------------------------8. 计算 Haar 小波逆变换--------------------------
    vector<vector<double>> value_num_recover;
    value_num_recover.resize(dwt_a2.size());

    for (size_t i = 0; i < dwt_a2.size(); ++i) {
        if (!dwt_a2[i].empty()) {
            vector<double> idwt_value_ca2 = dwt_a2[i];
            vector<double> idwt_value_ca1 = ca1_bin[i];
            vector<double> idwt_value_cd2 = d2_bin[i];
            vector<double> idwt_value_cd1 = d1_bin[i];

            // Perform inverse wavelet transform (placeholder for actual implementation)
            vector<double> ca1r = idwt(idwt_value_ca2, idwt_value_cd2);

            if (idwt_value_ca2.size() == 1) {
                reverse(ca1r.begin(), ca1r.end());
            }

            vector<double> dwt_recover = idwt(vector<double>(ca1r.begin(), ca1r.begin() + idwt_value_cd1.size()), idwt_value_cd1);

            if (value_num[i].size() % 2 == 1 && dwt_recover.size() == value_num[i].size() + 1) {
                dwt_recover.pop_back();
            }

            value_num_recover[i] = dwt_recover;
        }
    }
    //将value_num_recover四舍五入
	for (size_t i = 0; i < value_num_recover.size(); ++i) {
		if (!value_num_recover[i].empty()) {
			for (size_t j = 0; j < value_num_recover[i].size(); ++j) {
				value_num_recover[i][j] = round(value_num_recover[i][j]);
			}
		}
	}

	// --------------------------9. 计算嵌入后的音频(恢复S1）--------------------------
    // 偏移S1_recover到正数范围
    vector<double> S1_recover(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        S1_recover[i] = data[i] + dRange + 1;
    }
	// 将嵌入后的值恢复到原始范围
	for (size_t i = 0; i < value_num_recover.size(); ++i) {
		if (!value_num_recover[i].empty()) {
			for (size_t j = 0; j < value_num_recover[i].size(); ++j) {
                //cout << S1_recover[50359-1] << endl;
                //cout << joint[i][j] << endl;
                //cout << value_num_recover[i][j] << endl;
				S1_recover[joint[i][j]-1] = value_num_recover[i][j];
                //cout << S1_recover[50359-1] << endl;
			}
		}
	}
    //Iw = S1_recover - dRange - 1;
    //Sw_embedding=round(S1_recover);
	// 将 S1_recover 四舍五入
    vector <double> Sw_embedding(S1_recover.size());
    vector <double> embeddingData;   // 将整数信号归一化为小数 [-1, 1] 区间,/NA
	for (size_t i = 0; i < S1_recover.size(); ++i) {
        S1_recover[i] = S1_recover[i] - dRange - 1;
        Sw_embedding[i] = round(S1_recover[i]);
		embeddingData.push_back(Sw_embedding[i]/NA);
	}

    //--------------------------10.保存和比较重建音频--------------------------
    string audioFile_Rebuild = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1_watermarked.wav";

    writeAudioFile(audioFile_Rebuild, embeddingData, info);

    //计算信噪比
    double nSum = 0.0;
    
    double nModify = 0.0;
    for (int i = 0; i < data_original.size(); ++i)
    {
        nModify = nModify + pow(double(Sw_embedding[i] - data_original[i]), 2);
        nSum = nSum + pow(double(data_original[i]), 2);
    }
    double SNR = -10 * log10(double(nModify) / double(nSum));

	cout << "SNR: " << SNR << endl;

    //--------------------------11.计算音频的幅值均值A--------------------------
    //把Sw_embedding赋值给data
	data.clear();
	for (size_t i = 0; i < Sw_embedding.size(); ++i) {
        data.push_back(Sw_embedding[i]);
	}
    A = 0;
    for (size_t i = 0; i < data.size(); i++) {
        A += abs(data[i]);
    }
    A /= L_data_original;
    // 计算分箱的覆盖范围和宽度
    dRange = floor(lamda * A);
    bin_width = floor(2 * dRange / nBins);

    // --------------------------12. 计算分箱的索引--------------------------
    // 筛选数据并记录索引
    B1.clear();
	index.clear();
    for (size_t i = 0; i < data.size(); ++i) {
        if (round(data[i]) >= round(-dRange) && round(data[i]) <= round(dRange)) {
            B1.push_back(data[i]);
            index.push_back(i + 1); // MATLAB 索引从 1 开始
        }
    }
    L_B1 = B1.size();

    //--------------------------13.统计各个bins中的样本数量--------------------------
    nBins_number = 10;
    nBinsFinal = nBins + nBins_number;// 统计 bin 的数量，实验中数值为 183 + 10 = 193

    // 偏移 data 和 B1 到正数范围
	data1.clear();
    for (size_t i = 0; i < data.size(); ++i) {
        data1.push_back(data[i] + dRange + 1);
    }

    S1.clear();
    for (size_t i = 0; i < B1.size(); ++i) {
        S1.push_back(B1[i] + dRange + 1);
    }

    // 将 nBins_number 转成二进制 6 位
    /*bitset<6> nBins_number_bin(nBins_number);*/

    // 初始化每个 bin 的样本数量
    //vector<double> nNum(nBinsFinal, 0);
	initArrayInt(nNum, nBinsFinal, 0); 

    // nNum表示是第几个 bin 里面，进行了直方图统计?
    // 为了初始化value_num和index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
        if (binIndex >= 0 && binIndex < nNum.size()) {
            nNum[binIndex] += 1;
        }
    }

	initArrayInt(lo, nNum.size(), 1); // 初始化 lo 为全 1
    // 初始化 value_num 和 index_num
    for (int i = 0; i < nBins; ++i) {
        value_num[i] = vector<double>(nNum[i], 0.0); // 初始化为 0.0
        index_num[i] = vector<int>(nNum[i], 0);      // 初始化为 0
    }

    // 这里的 nNum[i] 是每个 bin 中的样本数量
    // value_num记录样本值，存在第(S1(i) / bin_width)个bin的第k个位置，k为lo对应bin记录的数，就是这个bin有了多少个数据，一个萝卜一个坑
    // 填充 value_num 和 index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // 转换为 0-based 索引
        if (binIndex >= 0 && binIndex < nBins) {
            value_num[binIndex][lo[binIndex] - 1] = S1[i];
            index_num[binIndex][lo[binIndex] - 1] = index[i];
            lo[binIndex] += 1; // 更新 lo
        }
    }

    //记录每个bin里面连续的有哪些变量
    vector<vector<int>> joint_Rebuild(nBinsFinal);              // 存储连续位置
    vector<vector<double>> value_num_joint_Rebuild(nBinsFinal); // 存储连续值

    // 遍历 index_num
    for (int i = 0; i < index_num.size(); ++i) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint_Rebuild[i], value_num_joint_Rebuild[i]);
        }
    }

    // --------------------------14. 计算 Haar 小波变换--------------------------
    vector<vector<double>> ca1_bin_Rebuild(nBinsFinal); // 初始化 a1_bin
    vector<vector<double>> ca2_bin_Rebuild(nBinsFinal); // 初始化 a2_bin
    vector<vector<double>> d1_bin_Rebuild(nBinsFinal);  // 初始化 d1_bin
    vector<vector<double>> d2_bin_Rebuild(nBinsFinal);  // 初始化 d2_bin
    for (size_t i = 0; i < value_num_joint_Rebuild.size(); ++i) {
        if (!value_num_joint_Rebuild[i].empty()) {
            vector<double> dwt_value = value_num_joint_Rebuild[i];
            haar_wavelet_transform(dwt_value, ca1_bin_Rebuild[i], d1_bin_Rebuild[i]); // 近似系数
            haar_wavelet_transform(ca1_bin_Rebuild[i], ca2_bin_Rebuild[i], d2_bin_Rebuild[i]); // 2阶
        }
    }

    //水印提取
    vector<int> w1;
    for (size_t i = 0; i < ca2_bin_Rebuild.size(); ++i)
    {
        for (size_t j = 0; j < ca2_bin_Rebuild[i].size(); ++j)
        {
            ca2_bin_Rebuild[i][j] = roundn(ca2_bin_Rebuild[i][j], 1);
        }
    }
    extractWatermark(ca2_bin_Rebuild, w1);
    
    // 计算 BER
    double dBer = calculateBER(w1, w);

    // 输出结果
    cout << "Bit Error Rate (BER): " << dBer << endl;
	return 0;
}