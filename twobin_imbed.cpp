// twobin_imbed.cpp
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// 注意：需要安装 libsndfile 库来实现 WAV 读写
#include <sndfile.h>

using namespace std;

// 辅助函数：将数字四舍五入到指定位数（例如：roundn(lamda, -2)）
double roundn(double value, int n) {
    double factor = pow(10.0, -n);
    return round(value * factor) / factor;
}

// 简单的 dec2bin_double 实现（返回一个表示该数字的二进制位数组）
vector<int> dec2bin_double(int num, int bits) {
    vector<int> bin(bits, 0);
    for (int i = bits - 1; i >= 0; i--) {
        bin[i] = num & 1;
        num = num >> 1;
    }
    return bin;
}

// 实现 Haar 小波变换 (db1)
// 输入数据长度最好为偶数
void dwt(const vector<double>& data, vector<double>& approx, vector<double>& detail) {
    int n = data.size();
    int half = n / 2;
    approx.resize(half);
    detail.resize(half);
    const double invSqrt2 = 1.0 / sqrt(2.0);
    for (int i = 0; i < half; i++) {
        approx[i] = (data[2 * i] + data[2 * i + 1]) * invSqrt2;
        detail[i] = (data[2 * i] - data[2 * i + 1]) * invSqrt2;
    }
}

// 逆 Haar 小波变换
vector<double> idwt(const vector<double>& approx, const vector<double>& detail) {
    int n = approx.size();
    vector<double> data(2 * n);
    const double invSqrt2 = 1.0 / sqrt(2.0);
    for (int i = 0; i < n; i++) {
        data[2 * i] = (approx[i] + detail[i]) * invSqrt2;
        data[2 * i + 1] = (approx[i] - detail[i]) * invSqrt2;
    }
    return data;
}

// find_continuous_place：查找索引数组中连续的部分
// 这里我们简单实现：返回连续序列（差值为1）的第一个段落，要求至少长度为4
// joint: 存储原始索引； value: 对应的信号值
void find_continuous_place(const vector<int>& indices, const vector<double>& values,
    vector<int>& joint, vector<double>& value_joint) {
    if (indices.empty()) return;
    vector<int> tempJoint;
    vector<double> tempValues;
    tempJoint.push_back(indices[0]);
    tempValues.push_back(values[0]);
    for (size_t i = 1; i < indices.size(); i++) {
        if (indices[i] == indices[i - 1] + 1) {
            tempJoint.push_back(indices[i]);
            tempValues.push_back(values[i]);
        }
        else {
            if (tempJoint.size() >= 4) {
                joint = tempJoint;
                value_joint = tempValues;
                return;
            }
            tempJoint.clear();
            tempValues.clear();
            tempJoint.push_back(indices[i]);
            tempValues.push_back(values[i]);
        }
    }
    // 检查最后一段
    if (tempJoint.size() >= 4) {
        joint = tempJoint;
        value_joint = tempValues;
    }
}

// 读取 WAV 文件（使用 libsndfile）
bool readWav(const string& filename, vector<double>& data, int& sampleRate, int& bitsPerSample) {
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!infile) {
        cerr << "无法打开文件: " << filename << endl;
        return false;
    }
    sampleRate = sfinfo.samplerate;
    bitsPerSample = 8 * sizeof(short); // 假定16位
    int num_items = sfinfo.frames * sfinfo.channels;
    vector<short> buffer(num_items);
    sf_read_short(infile, buffer.data(), num_items);
    sf_close(infile);
    // 若为多通道，取第一个通道
    data.resize(sfinfo.frames);
    for (int i = 0; i < sfinfo.frames; i++) {
        data[i] = buffer[i * sfinfo.channels];
    }
    return true;
}

// 写入 WAV 文件（使用 libsndfile）
bool writeWav(const string& filename, const vector<double>& data, int sampleRate) {
    SF_INFO sfinfo;
    sfinfo.samplerate = sampleRate;
    sfinfo.channels = 1;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    SNDFILE* outfile = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
    if (!outfile) {
        cerr << "无法写入文件: " << filename << endl;
        return false;
    }
    // 将 double 转换为 short
    vector<short> buffer(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        buffer[i] = static_cast<short>(round(data[i] * 32768)); // 假定归一化范围[-1,1]
    }
    sf_write_short(outfile, buffer.data(), buffer.size());
    sf_close(outfile);
    return true;
}

int main() {
    // --------------------------
    // 1. 参数定义和水印定义
    // 定义水印向量 w（取值为 1 或 -1）
    vector<int> w = { 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1,
                      1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1,
                      -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, -1,
                      -1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1 };
    int Lw = w.size();

    double lamda = 4.0;  // lamda 越大 bin_width 越大
    lamda = roundn(lamda, -2);
    double T = 4.0;
    int nBegin = 1;
    int nBins = 3 * Lw + 3;
    int nEnd = nBins;

    // --------------------------
    // 2. 读取音频文件
    string audioFile = "D:\\xy\\code\\zhiftu\\track1.wav";
    vector<double> data;
    int Fs, bits;
    if (!readWav(audioFile, data, Fs, bits)) {
        return -1;
    }
    int NA = 1 << (bits - 1);  // 2^(bits-1)
    // 将数据放大到整数范围
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= NA;
    }
    vector<double> data_original = data;
    int L_data_original = data.size();

    // --------------------------
    // 3. 计算幅值均值 A 及 dRange, bin_width
    double sumA = 0;
    for (double d : data) {
        sumA += fabs(d);
    }
    double A = sumA / data.size();
    int dRange = floor(lamda * A);
    int bin_width = static_cast<int>(floor(2.0 * dRange / nBins));

    // --------------------------
    // 4. 选择满足条件的样本 B1 及其索引
    vector<double> B1;
    vector<int> indexB;
    for (size_t i = 0; i < data.size(); i++) {
        if (round(data[i]) >= round(-dRange) && round(data[i]) <= round(dRange)) {
            B1.push_back(data[i]);
            indexB.push_back(i);
        }
    }
    int L_B1 = B1.size();

    // --------------------------
    // 5. 统计各个 bin 中的样本数量
    vector<double> S1(L_B1);
    for (int i = 0; i < L_B1; i++) {
        S1[i] = B1[i] + dRange + 1;
    }
    int nBins_number = 10;
    vector<int> nBins_number_bin = dec2bin_double(nBins_number, 6);
    vector<int> nNum(nBins + nBins_number, 0);
    for (int i = 0; i < L_B1; i++) {
        int binIdx = static_cast<int>(ceil(S1[i] / bin_width));
        if (binIdx >= 0 && binIdx < nNum.size())
            nNum[binIdx]++;
    }

    // --------------------------
    // 6. 为每个 bin 分配存储空间，并保存样本值与索引
    vector< vector<double> > value_num(nNum.size());
    vector< vector<int> > index_num(nNum.size());
    // 初始化每个 bin 的数组
    for (size_t i = 0; i < nNum.size(); i++) {
        value_num[i].resize(nNum[i], 0.0);
        index_num[i].resize(nNum[i], 0);
    }
    vector<int> lo(nNum.size(), 0);
    for (int i = 0; i < L_B1; i++) {
        int binIdx = static_cast<int>(ceil(S1[i] / bin_width));
        if (binIdx >= 0 && binIdx < nNum.size()) {
            value_num[binIdx][lo[binIdx]] = S1[i];
            index_num[binIdx][lo[binIdx]] = indexB[i];
            lo[binIdx]++;
        }
    }

    // --------------------------
    // 7. 对每个 bin 找到连续区段（至少长度为4）
    int totalBins = nBins + nBins_number;
    vector< vector<int> > joint(totalBins);
    vector< vector<double> > value_num_joint(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
        }
    }

    // --------------------------
    // 8. 对每个 bin 的连续区段做两级 Haar 小波变换
    // 保存小波变换结果
    vector< vector<double> > ca1_bin(totalBins), ca2_bin(totalBins),
        d1_bin(totalBins), d2_bin(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (!value_num_joint[i].empty()) {
            vector<double> approx1, detail1;
            dwt(value_num_joint[i], approx1, detail1);
            vector<double> approx2, detail2;
            if (approx1.size() >= 2) { // 要求长度足够
                dwt(approx1, approx2, detail2);
            }
            ca1_bin[i] = approx1;
            ca2_bin[i] = approx2;
            d1_bin[i] = detail1;
            d2_bin[i] = detail2;
        }
    }
    // 保存原始 ca2_bin 用于后续比对
    vector< vector<double> > ca2_bin_original = ca2_bin;

    // --------------------------
    // 9. 水印嵌入处理
    vector<double> dR(Lw, 0.0);
    for (int j = 0; j < Lw; j++) {
        int nIndex = 0;
        if (j % 2 == 0)  // j 为偶数时（MATLAB从1开始，这里 j从0 开始）
            nIndex = nBegin - 1 + 3 * (j / 2);
        else
            nIndex = nEnd - 3 * ((j + 1) / 2);
        // 检查下标范围
        if (nIndex + 2 >= totalBins) continue;
        int a = ca2_bin[nIndex].size();
        int b = ca2_bin[nIndex + 1].size();
        int c = ca2_bin[nIndex + 2].size();
        if (a + c == 0) continue;
        dR[j] = (2.0 * b) / (a + c);
        // 根据水印位进行调控
        if (w[j] == 1) {
            if (dR[j] <= T) {
                int T0 = ceil(((a + c) * T - 2 * b) / (2 + T));
                int Ta = ceil(T0 * a / double(a + c));
                int Tc = ceil(T0 * c / double(a + c));
                // 修改 ca2_bin[nIndex] 和 ca2_bin[nIndex+2]
                for (int i = 0; i < Ta && i < ca2_bin[nIndex].size(); i++) {
                    ca2_bin[nIndex][i] += 2 * bin_width;
                }
                for (int i = 0; i < Tc && i < ca2_bin[nIndex + 2].size(); i++) {
                    ca2_bin[nIndex + 2][i] -= 2 * bin_width;
                }
            }
        }
        else { // w[j] == -1
            if (dR[j] >= 1.0 / T) {
                int T0 = ceil((2 * b * T - (a + c)) / double(1 + 2 * T));
                int Ta = ceil(T0 * a / double(a + c));
                int Tc = ceil(T0 * c / double(a + c));
                for (int i = 0; i < Ta && i < ca2_bin[nIndex + 1].size(); i++) {
                    ca2_bin[nIndex + 1][i] -= 2 * bin_width;
                }
                for (int i = Ta; i < Tc && i < ca2_bin[nIndex + 1].size(); i++) {
                    ca2_bin[nIndex + 1][i] += 2 * bin_width;
                }
            }
        }
    }

    // --------------------------
    // 10. 逆小波变换恢复信号
    vector< vector<double> > value_num_recover(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (!ca2_bin[i].empty()) {
            // 先对第二级逆变换
            vector<double> ca1r = idwt(ca2_bin[i], d2_bin[i]);
            // 调整长度
            if (ca1r.size() > d1_bin[i].size())
                ca1r.resize(d1_bin[i].size());
            vector<double> dwt_recover = idwt(ca1r, d1_bin[i]);
            // 若长度不匹配，可做调整（这里简化处理）
            if (value_num[i].size() % 2 == 1 && dwt_recover.size() == value_num[i].size() + 1) {
                dwt_recover.pop_back();
            }
            value_num_recover[i] = dwt_recover;
        }
    }
    // 利用恢复后的 value_num_recover 和 joint 重构 S1_recover
    vector<double> S1_recover = S1;
    for (int i = 0; i < totalBins; i++) {
        if (!value_num_recover[i].empty() && !joint[i].empty()) {
            for (size_t j = 0; j < value_num_recover[i].size() && j < joint[i].size(); j++) {
                int idx = joint[i][j];
                if (idx >= 0 && idx < S1_recover.size()) {
                    S1_recover[idx] = value_num_recover[i][j];
                }
            }
        }
    }
    // 恢复整数信号
    vector<double> Iw(S1_recover.size());
    for (size_t i = 0; i < S1_recover.size(); i++) {
        Iw[i] = S1_recover[i] - dRange - 1;
    }
    // 四舍五入得到嵌入后信号
    vector<double> Sw_embedding(Iw.size());
    for (size_t i = 0; i < Iw.size(); i++) {
        Sw_embedding[i] = round(Iw[i]);
    }

    // --------------------------
    // 11. 保存重构音频
    vector<double> y_reconstructed(Sw_embedding.size());
    for (size_t i = 0; i < Sw_embedding.size(); i++) {
        y_reconstructed[i] = Sw_embedding[i] / double(NA);
    }
    string outputFile = "D:\\xy\\code\\zhiftu\\track1_water.wav";
    if (!writeWav(outputFile, y_reconstructed, Fs)) {
        return -1;
    }

    // --------------------------
    // 12. 计算 SNR
    double nModify = 0.0, nSum = 0.0;
    for (size_t i = 0; i < data_original.size(); i++) {
        double diff = Sw_embedding[i] - data_original[i];
        nModify += diff * diff;
        nSum += data_original[i] * data_original[i];
    }
    double SNR = -10.0 * log10(nModify / nSum);
    cout << "SNR = " << SNR << " dB" << endl;

    // --------------------------
    // 13. 以下部分为水印提取过程，与前面类似：
    // 对重构音频再次处理，提取水印
    // 重新计算 B1、S1、分箱等步骤……
    // 这里略去部分重复过程，实际可将公共代码封装为函数

    // 示例：对每个 bin 重新做小波变换（使用 Haar 变换），存入 ca2_binw
    vector< vector<double> > ca1_binw(totalBins), ca2_binw(totalBins),
        d1_binw(totalBins), d2_binw(totalBins);
    // 这里假定对 Sw_embedding 做相同的分箱及连续区段查找，
    // 实际应重复步骤 4~7，以下仅示意
    for (int i = 0; i < totalBins; i++) {
        if (!value_num_joint[i].empty()) {
            vector<double> approx1, detail1;
            dwt(value_num_joint[i], approx1, detail1);
            vector<double> approx2, detail2;
            if (approx1.size() >= 2)
                dwt(approx1, approx2, detail2);
            ca1_binw[i] = approx1;
            ca2_binw[i] = approx2;
            d1_binw[i] = detail1;
            d2_binw[i] = detail2;
        }
    }
    // 提取水印 w1
    vector<int> w1(Lw, 0);
    vector<double> dRw(Lw, 0.0);
    for (int j = 0; j < Lw; j++) {
        int nIndex = 0;
        if (j % 2 == 0)
            nIndex = nBegin - 1 + 3 * (j / 2);
        else
            nIndex = nEnd - 3 * ((j + 1) / 2);
        if (nIndex + 2 >= totalBins) continue;
        int a = ca2_binw[nIndex].size();
        int b = ca2_binw[nIndex + 1].size();
        int c = ca2_binw[nIndex + 2].size();
        if (a + c == 0) continue;
        dRw[j] = (2.0 * b) / (a + c);
        w1[j] = (dRw[j] >= 1.0) ? 1 : -1;
    }
    // 计算 BER
    int dBer = 0;
    for (int j = 0; j < Lw; j++) {
        if (w1[j] != w[j])
            dBer++;
    }
    double ber = double(dBer) / Lw;
    cout << "BER = " << ber << endl;

    // --------------------------
    // 14. （可选）绘图：C++ 中可使用第三方库（如 matplotlib-cpp 或 gnuplot）进行结果绘图

    return 0;
}
