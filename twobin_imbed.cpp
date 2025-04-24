// twobin_imbed.cpp
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

// ע�⣺��Ҫ��װ libsndfile ����ʵ�� WAV ��д
#include <sndfile.h>

using namespace std;

// �����������������������뵽ָ��λ�������磺roundn(lamda, -2)��
double roundn(double value, int n) {
    double factor = pow(10.0, -n);
    return round(value * factor) / factor;
}

// �򵥵� dec2bin_double ʵ�֣�����һ����ʾ�����ֵĶ�����λ���飩
vector<int> dec2bin_double(int num, int bits) {
    vector<int> bin(bits, 0);
    for (int i = bits - 1; i >= 0; i--) {
        bin[i] = num & 1;
        num = num >> 1;
    }
    return bin;
}

// ʵ�� Haar С���任 (db1)
// �������ݳ������Ϊż��
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

// �� Haar С���任
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

// find_continuous_place���������������������Ĳ���
// �������Ǽ�ʵ�֣������������У���ֵΪ1���ĵ�һ�����䣬Ҫ�����ٳ���Ϊ4
// joint: �洢ԭʼ������ value: ��Ӧ���ź�ֵ
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
    // ������һ��
    if (tempJoint.size() >= 4) {
        joint = tempJoint;
        value_joint = tempValues;
    }
}

// ��ȡ WAV �ļ���ʹ�� libsndfile��
bool readWav(const string& filename, vector<double>& data, int& sampleRate, int& bitsPerSample) {
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!infile) {
        cerr << "�޷����ļ�: " << filename << endl;
        return false;
    }
    sampleRate = sfinfo.samplerate;
    bitsPerSample = 8 * sizeof(short); // �ٶ�16λ
    int num_items = sfinfo.frames * sfinfo.channels;
    vector<short> buffer(num_items);
    sf_read_short(infile, buffer.data(), num_items);
    sf_close(infile);
    // ��Ϊ��ͨ����ȡ��һ��ͨ��
    data.resize(sfinfo.frames);
    for (int i = 0; i < sfinfo.frames; i++) {
        data[i] = buffer[i * sfinfo.channels];
    }
    return true;
}

// д�� WAV �ļ���ʹ�� libsndfile��
bool writeWav(const string& filename, const vector<double>& data, int sampleRate) {
    SF_INFO sfinfo;
    sfinfo.samplerate = sampleRate;
    sfinfo.channels = 1;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    SNDFILE* outfile = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
    if (!outfile) {
        cerr << "�޷�д���ļ�: " << filename << endl;
        return false;
    }
    // �� double ת��Ϊ short
    vector<short> buffer(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        buffer[i] = static_cast<short>(round(data[i] * 32768)); // �ٶ���һ����Χ[-1,1]
    }
    sf_write_short(outfile, buffer.data(), buffer.size());
    sf_close(outfile);
    return true;
}

int main() {
    // --------------------------
    // 1. ���������ˮӡ����
    // ����ˮӡ���� w��ȡֵΪ 1 �� -1��
    vector<int> w = { 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1,
                      1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1,
                      -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, -1,
                      -1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1 };
    int Lw = w.size();

    double lamda = 4.0;  // lamda Խ�� bin_width Խ��
    lamda = roundn(lamda, -2);
    double T = 4.0;
    int nBegin = 1;
    int nBins = 3 * Lw + 3;
    int nEnd = nBins;

    // --------------------------
    // 2. ��ȡ��Ƶ�ļ�
    string audioFile = "D:\\xy\\code\\zhiftu\\track1.wav";
    vector<double> data;
    int Fs, bits;
    if (!readWav(audioFile, data, Fs, bits)) {
        return -1;
    }
    int NA = 1 << (bits - 1);  // 2^(bits-1)
    // �����ݷŴ�������Χ
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= NA;
    }
    vector<double> data_original = data;
    int L_data_original = data.size();

    // --------------------------
    // 3. �����ֵ��ֵ A �� dRange, bin_width
    double sumA = 0;
    for (double d : data) {
        sumA += fabs(d);
    }
    double A = sumA / data.size();
    int dRange = floor(lamda * A);
    int bin_width = static_cast<int>(floor(2.0 * dRange / nBins));

    // --------------------------
    // 4. ѡ���������������� B1 ��������
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
    // 5. ͳ�Ƹ��� bin �е���������
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
    // 6. Ϊÿ�� bin ����洢�ռ䣬����������ֵ������
    vector< vector<double> > value_num(nNum.size());
    vector< vector<int> > index_num(nNum.size());
    // ��ʼ��ÿ�� bin ������
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
    // 7. ��ÿ�� bin �ҵ��������Σ����ٳ���Ϊ4��
    int totalBins = nBins + nBins_number;
    vector< vector<int> > joint(totalBins);
    vector< vector<double> > value_num_joint(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
        }
    }

    // --------------------------
    // 8. ��ÿ�� bin ���������������� Haar С���任
    // ����С���任���
    vector< vector<double> > ca1_bin(totalBins), ca2_bin(totalBins),
        d1_bin(totalBins), d2_bin(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (!value_num_joint[i].empty()) {
            vector<double> approx1, detail1;
            dwt(value_num_joint[i], approx1, detail1);
            vector<double> approx2, detail2;
            if (approx1.size() >= 2) { // Ҫ�󳤶��㹻
                dwt(approx1, approx2, detail2);
            }
            ca1_bin[i] = approx1;
            ca2_bin[i] = approx2;
            d1_bin[i] = detail1;
            d2_bin[i] = detail2;
        }
    }
    // ����ԭʼ ca2_bin ���ں����ȶ�
    vector< vector<double> > ca2_bin_original = ca2_bin;

    // --------------------------
    // 9. ˮӡǶ�봦��
    vector<double> dR(Lw, 0.0);
    for (int j = 0; j < Lw; j++) {
        int nIndex = 0;
        if (j % 2 == 0)  // j Ϊż��ʱ��MATLAB��1��ʼ������ j��0 ��ʼ��
            nIndex = nBegin - 1 + 3 * (j / 2);
        else
            nIndex = nEnd - 3 * ((j + 1) / 2);
        // ����±귶Χ
        if (nIndex + 2 >= totalBins) continue;
        int a = ca2_bin[nIndex].size();
        int b = ca2_bin[nIndex + 1].size();
        int c = ca2_bin[nIndex + 2].size();
        if (a + c == 0) continue;
        dR[j] = (2.0 * b) / (a + c);
        // ����ˮӡλ���е���
        if (w[j] == 1) {
            if (dR[j] <= T) {
                int T0 = ceil(((a + c) * T - 2 * b) / (2 + T));
                int Ta = ceil(T0 * a / double(a + c));
                int Tc = ceil(T0 * c / double(a + c));
                // �޸� ca2_bin[nIndex] �� ca2_bin[nIndex+2]
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
    // 10. ��С���任�ָ��ź�
    vector< vector<double> > value_num_recover(totalBins);
    for (int i = 0; i < totalBins; i++) {
        if (!ca2_bin[i].empty()) {
            // �ȶԵڶ�����任
            vector<double> ca1r = idwt(ca2_bin[i], d2_bin[i]);
            // ��������
            if (ca1r.size() > d1_bin[i].size())
                ca1r.resize(d1_bin[i].size());
            vector<double> dwt_recover = idwt(ca1r, d1_bin[i]);
            // �����Ȳ�ƥ�䣬��������������򻯴���
            if (value_num[i].size() % 2 == 1 && dwt_recover.size() == value_num[i].size() + 1) {
                dwt_recover.pop_back();
            }
            value_num_recover[i] = dwt_recover;
        }
    }
    // ���ûָ���� value_num_recover �� joint �ع� S1_recover
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
    // �ָ������ź�
    vector<double> Iw(S1_recover.size());
    for (size_t i = 0; i < S1_recover.size(); i++) {
        Iw[i] = S1_recover[i] - dRange - 1;
    }
    // ��������õ�Ƕ����ź�
    vector<double> Sw_embedding(Iw.size());
    for (size_t i = 0; i < Iw.size(); i++) {
        Sw_embedding[i] = round(Iw[i]);
    }

    // --------------------------
    // 11. �����ع���Ƶ
    vector<double> y_reconstructed(Sw_embedding.size());
    for (size_t i = 0; i < Sw_embedding.size(); i++) {
        y_reconstructed[i] = Sw_embedding[i] / double(NA);
    }
    string outputFile = "D:\\xy\\code\\zhiftu\\track1_water.wav";
    if (!writeWav(outputFile, y_reconstructed, Fs)) {
        return -1;
    }

    // --------------------------
    // 12. ���� SNR
    double nModify = 0.0, nSum = 0.0;
    for (size_t i = 0; i < data_original.size(); i++) {
        double diff = Sw_embedding[i] - data_original[i];
        nModify += diff * diff;
        nSum += data_original[i] * data_original[i];
    }
    double SNR = -10.0 * log10(nModify / nSum);
    cout << "SNR = " << SNR << " dB" << endl;

    // --------------------------
    // 13. ���²���Ϊˮӡ��ȡ���̣���ǰ�����ƣ�
    // ���ع���Ƶ�ٴδ�����ȡˮӡ
    // ���¼��� B1��S1������Ȳ��衭��
    // ������ȥ�����ظ����̣�ʵ�ʿɽ����������װΪ����

    // ʾ������ÿ�� bin ������С���任��ʹ�� Haar �任�������� ca2_binw
    vector< vector<double> > ca1_binw(totalBins), ca2_binw(totalBins),
        d1_binw(totalBins), d2_binw(totalBins);
    // ����ٶ��� Sw_embedding ����ͬ�ķ��估�������β��ң�
    // ʵ��Ӧ�ظ����� 4~7�����½�ʾ��
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
    // ��ȡˮӡ w1
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
    // ���� BER
    int dBer = 0;
    for (int j = 0; j < Lw; j++) {
        if (w1[j] != w[j])
            dBer++;
    }
    double ber = double(dBer) / Lw;
    cout << "BER = " << ber << endl;

    // --------------------------
    // 14. ����ѡ����ͼ��C++ �п�ʹ�õ������⣨�� matplotlib-cpp �� gnuplot�����н����ͼ

    return 0;
}
