#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "readAudio.h"  // ��ȡ��Ƶ�ļ���ͷ�ļ�
#include <cmath>        // ���� abs �� floor
#include <numeric>      // ���� accumulate
#include "myFunction.h" // ����������ͷ�ļ�
#include "find_continuous_place.h" // ��������λ�õ�ͷ�ļ�
#include <bitset>       // ���� bitset
//#include "dwt.h"        // Haar С���任��ͷ�ļ�

using namespace std;

// --------------------------1. ���������ˮӡ����--------------------------
// ����ˮӡ���� w��ȡֵΪ 1 �� -1��
vector<int> w = { 1, 1, 1, 1, 1, 1, -1, -1, -1, 1,
                -1, 1 ,1, 1, -1, 1, 1, -1, 1, 1,
                1, -1, 1, 1, -1, -1, 1, 1,1, 1,
                1, -1, -1, -1,  1,  1, -1,  1, 1,  1,
                -1,  1,  1, -1,  1, -1, -1, -1, -1,  1,
                -1, -1,  1, -1,  1,  1, -1, -1,  1,  1 };
int Lw = w.size();

double lamda_unround = 4.0;  // lamda Խ�� bin_width Խ��
double lamda = roundn(lamda_unround, 2);
double T = 4.0;
int nBegin = 1;
int nBins = 3 * Lw + 3;
int nEnd = nBins;
int nBins_number = 10;
int nBinsFinal = nBins + nBins_number;// ͳ�� bin ��������ʵ������ֵΪ 183 + 10 = 193

// Haar С���任
void haar_wavelet_transform(const vector<double>& signal, vector<double>& approx, vector<double>& detail) {
    int n = signal.size();
    int half = n / 2;

    approx.resize(half);
    detail.resize(half);

    for (int i = 0; i < half; ++i) {
        approx[i] = (signal[2 * i] + signal[2 * i + 1]) / sqrt(2.0); // ����ϵ��
        detail[i] = (signal[2 * i] - signal[2 * i + 1]) / sqrt(2.0); // ϸ��ϵ��
    }
}

// Haar С���任�Ķ��ױ任
//void haar_wavelet_transform2D(int nBinsFinal, vector <vector< double >> & value_num_joint, vector <vector< double >>& dwt_a2)
//{
//    vector<vector<double>> ca1_bin(nBinsFinal); // ��ʼ�� a1_bin
//    //vector<vector<double>> ca2_bin(nBinsFinal); // ��ʼ�� a2_bin
//    vector<vector<double>> d1_bin(nBinsFinal);  // ��ʼ�� d1_bin
//    vector<vector<double>> d2_bin(nBinsFinal);  // ��ʼ�� d2_bin
//    for (size_t i = 0; i < value_num_joint.size(); ++i) {
//        if (!value_num_joint[i].empty()) {
//            vector<double> dwt_value = value_num_joint[i];
//            //vector<double> a1, d1, a2, d2;
//            haar_wavelet_transform(dwt_value, ca1_bin[i], d1_bin[i]); // ����ϵ��
//            haar_wavelet_transform(ca1_bin[i], dwt_a2[i], d2_bin[i]); // 2��
//
//        }
//    }
//}

//// ͳ�Ƹ���bins�е���������
//pair< vector<vector<int>>, vector<vector<double>> > count_bins(const vector<double>& data, double dRange, vector<double> B1, vector<int> index)
//{
//    double bin_width = floor(2 * dRange / nBins);
//    size_t L_B1 = B1.size();
//    // ƫ�� data �� B1 ��������Χ
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
//    // �� nBins_number ת�ɶ����� 6 λ
//    bitset<6> nBins_number_bin(nBins_number);
//
//    // ��ʼ��ÿ�� bin ����������
//    vector<int> nNum(nBinsFinal, 0);
//
//    // ��ʾ�ǵڼ��� bin ���棬������ֱ��ͼͳ��?
//    for (int i = 0; i < L_B1; ++i) {
//        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
//        if (binIndex >= 0 && binIndex < nNum.size()) {
//            nNum[binIndex] += 1;
//        }
//    }
//
//    vector<int> lo(nNum.size(), 1); // ��ʼ�� lo Ϊȫ 1
//    // ��ʼ�� value_num �� index_num
//    vector<vector<double>> value_num(nBinsFinal);
//    vector<vector<int>> index_num(nBinsFinal);
//
//    for (int i = 0; i < nBins; ++i) {
//        value_num[i] = vector<double>(nNum[i], 0.0); // ��ʼ��Ϊ 0.0
//        index_num[i] = vector<int>(nNum[i], 0);      // ��ʼ��Ϊ 0
//    }
//
//    // ����� nNum[i] ��ÿ�� bin �е���������
//    // value_num��¼����ֵ�����ڵ�(S1(i) / bin_width)��bin�ĵ�k��λ�ã�kΪlo��Ӧbin��¼�������������bin���˶��ٸ����ݣ�һ���ܲ�һ����
//    // ��� value_num �� index_num
//    for (int i = 0; i < L_B1; ++i) {
//        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
//        if (binIndex >= 0 && binIndex < nBins) {
//            value_num[binIndex][lo[binIndex] - 1] = S1[i];
//            index_num[binIndex][lo[binIndex] - 1] = index[i];
//            lo[binIndex] += 1; // ���� lo
//        }
//    }
//
//    //��¼ÿ��bin��������������Щ����
//    vector<vector<int>> joint(nBinsFinal);              // �洢����λ��
//    vector<vector<double>> value_num_joint(nBinsFinal); // �洢����ֵ
//
//    // ���� index_num
//    for (int i = 0; i < index_num.size(); ++i) {
//        if (index_num[i].size() >= 4) {
//            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
//        }
//    }
//    return { joint, value_num_joint };
//}

// ʵ����С���任 (idwt) ʹ�� 'db1' С����
vector<double> idwt(const vector<double>& a, const vector<double>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Input vectors 'a' and 'b' must have the same length.");
    }

    size_t n = a.size();
    vector<double> result(2 * n);

    // db1 С��������任�˲���
    const double h0 = 1.0 / sqrt(2.0); // ��ͨ�˲���ϵ��
    const double h1 = 1.0 / sqrt(2.0); // ��ͨ�˲���ϵ��

    for (size_t i = 0; i < n; ++i) {
        result[2 * i] = h0 * a[i] + h1 * b[i]; // ż������
        result[2 * i + 1] = h0 * a[i] - h1 * b[i]; // ��������
    }

    return result;
}

// ��ȡˮӡ
void extractWatermark(const vector<vector<double>>& ca2_binw,  vector<int>& w1) 
{
    vector<double> dRw(Lw); // �洢�ع���ı���
    w1.resize(Lw);          // �洢��ȡ��ˮӡ

    for (int j = 0; j < Lw; ++j) {
        int nIndex;
        if (j % 2 == 0) {
			nIndex = nBegin + 3 * ((j + 1 + 1) / 2 - 1) - 1;// MATLAB ������ 1 ��ʼ
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
    //// ��� dRw �� w1
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
    ////	cerr << "�޷����ļ�: " << filename << endl;
    ////	return;
    ////}
    //for (size_t i = 0; i < w1.size(); i++) {
    //    outfile << w1[i] << endl;
    //}
    //outfile.close();
    //cout << "������ɹ���" << endl;

    cout << endl;
}

double calculateBER(const vector<int>& w1, const vector<int>& w) {
    
    double dBer = 0;

    // ����ÿ��ˮӡλ���������λ��
    for (int j = 0; j < Lw; ++j) {
        if (w1[j] != w[j]) {
            dBer++;
        }
    }

    // ���� BER (���ش�����)
    return dBer / Lw;
}

int main()
{

    // --------------------------2. ��ȡ��Ƶ�ļ�--------------------------
	//
    string audioFile = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1.wav";
    pair<vector<double>, AudioInfo> result = readAudio(audioFile); // ��ȡ��Ƶ�ļ�
    vector<double> data = result.first; // ��ȡ��Ƶ����
    AudioInfo info = result.second;     // ��ȡ��Ƶ��Ϣ

    int NA = 1 << (info.bitsPerSample - 1);  // 2^(bits-1)
    //�����ݷŴ�������Χ
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= NA;
		//printf("%f ", data[i]);
    }
    vector<double> data_original = data;
    int L_data_original = data.size();

    //--------------------------3.������Ƶ�ķ�ֵ��ֵA--------------------------
	double A = 0;
	for (size_t i = 0; i < L_data_original; i++) {
		A += abs(data[i]);
	}
	A /= L_data_original;
    // �������ĸ��Ƿ�Χ�Ϳ��
    double dRange = floor(lamda * A);
    double bin_width = floor(2 * dRange / nBins);

    // --------------------------4. ������������--------------------------

    // ɸѡ���ݲ���¼����
    vector<double> B1; // �洢������������Ƶ����
    vector<int> index; // �洢��Ӧ������
    for (size_t i = 0; i < data.size(); ++i) {
        if (round(data[i]) >= round(-dRange) && round(data[i]) <= round(dRange)) {
            B1.push_back(data[i]);
            index.push_back(i + 1); // MATLAB ������ 1 ��ʼ
        }
    }
    // ���� B1 �ĳ���
    

    //--------------------------5.ͳ�Ƹ���bins�е���������--------------------------
	//pair< vector<vector<int>>, vector<vector<double>> > count_bins_result = count_bins(data, dRange, B1, index);
	//vector<vector<int>> joint = count_bins_result.first;              // �洢����λ��
	//vector<vector<double>> value_num_joint = count_bins_result.second; // �洢����ֵ

    size_t L_B1 = B1.size();
    // ƫ�� data �� B1 ��������Χ
    vector<double> data1(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        data1[i] = data[i] + dRange + 1;
    }

    vector<double> S1(B1.size());
    for (size_t i = 0; i < B1.size(); ++i) {
        S1[i] = B1[i] + dRange + 1;
    }

    // �� nBins_number ת�ɶ����� 6 λ
    bitset<6> nBins_number_bin(nBins_number);

    // ��ʼ��ÿ�� bin ����������
    vector<int> nNum(nBinsFinal, 0);

    // ��ʾ�ǵڼ��� bin ���棬������ֱ��ͼͳ��?
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
        if (binIndex >= 0 && binIndex < nNum.size()) {
            nNum[binIndex] += 1;
        }
    }

    vector<int> lo(nNum.size(), 1); // ��ʼ�� lo Ϊȫ 1
    // ��ʼ�� value_num �� index_num
    vector<vector<double>> value_num(nBinsFinal);
    vector<vector<int>> index_num(nBinsFinal);

    for (int i = 0; i < nBins; ++i) {
        value_num[i] = vector<double>(nNum[i], 0.0); // ��ʼ��Ϊ 0.0
        index_num[i] = vector<int>(nNum[i], 0);      // ��ʼ��Ϊ 0
    }

    // ����� nNum[i] ��ÿ�� bin �е���������
    // value_num��¼����ֵ�����ڵ�(S1(i) / bin_width)��bin�ĵ�k��λ�ã�kΪlo��Ӧbin��¼�������������bin���˶��ٸ����ݣ�һ���ܲ�һ����
    // ��� value_num �� index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
        if (binIndex >= 0 && binIndex < nBins) {
            value_num[binIndex][lo[binIndex] - 1] = S1[i];
            index_num[binIndex][lo[binIndex] - 1] = index[i];
            lo[binIndex] += 1; // ���� lo
        }
    }

    //��¼ÿ��bin��������������Щ����
    vector<vector<int>> joint(nBinsFinal);              // �洢����λ��
    vector<vector<double>> value_num_joint(nBinsFinal); // �洢����ֵ

    // ���� index_num
    for (int i = 0; i < index_num.size(); ++i) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint[i], value_num_joint[i]);
        }
    }


	// --------------------------6. ���� Haar С���任--------------------------
	//vector<vector<double>> dwt_a2(nBinsFinal);
	//haar_wavelet_transform2D(nBinsFinal, value_num_joint, dwt_a2);

    vector<vector<double>> ca1_bin(nBinsFinal); // ��ʼ�� a1_bin
    vector<vector<double>> dwt_a2(nBinsFinal); // ��ʼ�� a2_bin
    vector<vector<double>> d1_bin(nBinsFinal);  // ��ʼ�� d1_bin
    vector<vector<double>> d2_bin(nBinsFinal);  // ��ʼ�� d2_bin
    for (size_t i = 0; i < value_num_joint.size(); ++i) {
        if (!value_num_joint[i].empty()) {
            vector<double> dwt_value = value_num_joint[i];
            //vector<double> a1, d1, a2, d2;
            haar_wavelet_transform(dwt_value, ca1_bin[i], d1_bin[i]); // ����ϵ��
            haar_wavelet_transform(ca1_bin[i], dwt_a2[i], d2_bin[i]); // 2��

        }
    }
    //  dwt_a2��������
    for (size_t i = 0; i < dwt_a2.size(); ++i)
    {
        for (size_t j = 0; j < dwt_a2[i].size(); ++j)
        {
            dwt_a2[i][j] = roundn(dwt_a2[i][j], 1);
        }
    }

	//--------------------------7.������ֵ T ����ˮӡǶ��--------------------------
    //Ƕ��ǰ�����ݣ����ڱȽ�ˮӡǶ��ǰ��Ĳ���
	vector<double> dR;//���ﲻҪ������С�����������һ����ʼ��ֵ
    int a_w[60], b_w[60], c_w[60];
    for (int j = 0; j < Lw; ++j) {
        static int nIndex = 0;
        
        if (j % 2 == 0) {
            nIndex = nBegin + 3 * ((j+1 + 1) / 2 - 1);
        }
        else {
            nIndex = nEnd - 3 * ((j + 1) / 2) + 1;//����j+1���Ǽ�����
        }
		nIndex -= 1; // MATLAB ������ 1 ��ʼ
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

 //   //���dwt_a2ÿ�е�������excel
	//ofstream outfile("dwt_a2.xls");
	////if (!file.is_open()) {
	////	cerr << "�޷����ļ�: " << filename << endl;
	////	return;
	////}
	//for (size_t i = 0; i < dwt_a2.size(); i++) {
	//	outfile << dwt_a2[i].size() << endl;
	//}
	//outfile.close();
	//cout << "���dwt_a2���ɹ���" << endl;


	// --------------------------8. ���� Haar С����任--------------------------
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
    //��value_num_recover��������
	for (size_t i = 0; i < value_num_recover.size(); ++i) {
		if (!value_num_recover[i].empty()) {
			for (size_t j = 0; j < value_num_recover[i].size(); ++j) {
				value_num_recover[i][j] = round(value_num_recover[i][j]);
			}
		}
	}

	// --------------------------9. ����Ƕ������Ƶ(�ָ�S1��--------------------------
    // ƫ��S1_recover��������Χ
    vector<double> S1_recover(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        S1_recover[i] = data[i] + dRange + 1;
    }
	// ��Ƕ����ֵ�ָ���ԭʼ��Χ
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
	// �� S1_recover ��������
    vector <double> Sw_embedding(S1_recover.size());
    vector <double> embeddingData;   // �������źŹ�һ��ΪС�� [-1, 1] ����,/NA
	for (size_t i = 0; i < S1_recover.size(); ++i) {
        S1_recover[i] = S1_recover[i] - dRange - 1;
        Sw_embedding[i] = round(S1_recover[i]);
		embeddingData.push_back(Sw_embedding[i]/NA);
	}

    //--------------------------10.����ͱȽ��ؽ���Ƶ--------------------------
    string audioFile_Rebuild = "D:\\Programming\\Code\\VsRepos\\WatermarkEmbedding\\track1_watermarked.wav";

    writeAudioFile(audioFile_Rebuild, embeddingData, info);

    //���������
    double nSum = 0.0;
    
    double nModify = 0.0;
    for (int i = 0; i < data_original.size(); ++i)
    {
        nModify = nModify + pow(double(Sw_embedding[i] - data_original[i]), 2);
        nSum = nSum + pow(double(data_original[i]), 2);
    }
    double SNR = -10 * log10(double(nModify) / double(nSum));

	cout << "SNR: " << SNR << endl;

    //--------------------------11.������Ƶ�ķ�ֵ��ֵA--------------------------
    //��Sw_embedding��ֵ��data
	data.clear();
	for (size_t i = 0; i < Sw_embedding.size(); ++i) {
        data.push_back(Sw_embedding[i]);
	}
    A = 0;
    for (size_t i = 0; i < data.size(); i++) {
        A += abs(data[i]);
    }
    A /= L_data_original;
    // �������ĸ��Ƿ�Χ�Ϳ��
    dRange = floor(lamda * A);
    bin_width = floor(2 * dRange / nBins);

    // --------------------------12. ������������--------------------------
    // ɸѡ���ݲ���¼����
    B1.clear();
	index.clear();
    for (size_t i = 0; i < data.size(); ++i) {
        if (round(data[i]) >= round(-dRange) && round(data[i]) <= round(dRange)) {
            B1.push_back(data[i]);
            index.push_back(i + 1); // MATLAB ������ 1 ��ʼ
        }
    }
    L_B1 = B1.size();

    //--------------------------13.ͳ�Ƹ���bins�е���������--------------------------
    nBins_number = 10;
    nBinsFinal = nBins + nBins_number;// ͳ�� bin ��������ʵ������ֵΪ 183 + 10 = 193

    // ƫ�� data �� B1 ��������Χ
	data1.clear();
    for (size_t i = 0; i < data.size(); ++i) {
        data1.push_back(data[i] + dRange + 1);
    }

    S1.clear();
    for (size_t i = 0; i < B1.size(); ++i) {
        S1.push_back(B1[i] + dRange + 1);
    }

    // �� nBins_number ת�ɶ����� 6 λ
    /*bitset<6> nBins_number_bin(nBins_number);*/

    // ��ʼ��ÿ�� bin ����������
    //vector<double> nNum(nBinsFinal, 0);
	initArrayInt(nNum, nBinsFinal, 0); 

    // nNum��ʾ�ǵڼ��� bin ���棬������ֱ��ͼͳ��?
    // Ϊ�˳�ʼ��value_num��index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
        if (binIndex >= 0 && binIndex < nNum.size()) {
            nNum[binIndex] += 1;
        }
    }

	initArrayInt(lo, nNum.size(), 1); // ��ʼ�� lo Ϊȫ 1
    // ��ʼ�� value_num �� index_num
    for (int i = 0; i < nBins; ++i) {
        value_num[i] = vector<double>(nNum[i], 0.0); // ��ʼ��Ϊ 0.0
        index_num[i] = vector<int>(nNum[i], 0);      // ��ʼ��Ϊ 0
    }

    // ����� nNum[i] ��ÿ�� bin �е���������
    // value_num��¼����ֵ�����ڵ�(S1(i) / bin_width)��bin�ĵ�k��λ�ã�kΪlo��Ӧbin��¼�������������bin���˶��ٸ����ݣ�һ���ܲ�һ����
    // ��� value_num �� index_num
    for (int i = 0; i < L_B1; ++i) {
        int binIndex = static_cast<int>(ceil(S1[i] / bin_width)) - 1; // ת��Ϊ 0-based ����
        if (binIndex >= 0 && binIndex < nBins) {
            value_num[binIndex][lo[binIndex] - 1] = S1[i];
            index_num[binIndex][lo[binIndex] - 1] = index[i];
            lo[binIndex] += 1; // ���� lo
        }
    }

    //��¼ÿ��bin��������������Щ����
    vector<vector<int>> joint_Rebuild(nBinsFinal);              // �洢����λ��
    vector<vector<double>> value_num_joint_Rebuild(nBinsFinal); // �洢����ֵ

    // ���� index_num
    for (int i = 0; i < index_num.size(); ++i) {
        if (index_num[i].size() >= 4) {
            find_continuous_place(index_num[i], value_num[i], joint_Rebuild[i], value_num_joint_Rebuild[i]);
        }
    }

    // --------------------------14. ���� Haar С���任--------------------------
    vector<vector<double>> ca1_bin_Rebuild(nBinsFinal); // ��ʼ�� a1_bin
    vector<vector<double>> ca2_bin_Rebuild(nBinsFinal); // ��ʼ�� a2_bin
    vector<vector<double>> d1_bin_Rebuild(nBinsFinal);  // ��ʼ�� d1_bin
    vector<vector<double>> d2_bin_Rebuild(nBinsFinal);  // ��ʼ�� d2_bin
    for (size_t i = 0; i < value_num_joint_Rebuild.size(); ++i) {
        if (!value_num_joint_Rebuild[i].empty()) {
            vector<double> dwt_value = value_num_joint_Rebuild[i];
            haar_wavelet_transform(dwt_value, ca1_bin_Rebuild[i], d1_bin_Rebuild[i]); // ����ϵ��
            haar_wavelet_transform(ca1_bin_Rebuild[i], ca2_bin_Rebuild[i], d2_bin_Rebuild[i]); // 2��
        }
    }

    //ˮӡ��ȡ
    vector<int> w1;
    for (size_t i = 0; i < ca2_bin_Rebuild.size(); ++i)
    {
        for (size_t j = 0; j < ca2_bin_Rebuild[i].size(); ++j)
        {
            ca2_bin_Rebuild[i][j] = roundn(ca2_bin_Rebuild[i][j], 1);
        }
    }
    extractWatermark(ca2_bin_Rebuild, w1);
    
    // ���� BER
    double dBer = calculateBER(w1, w);

    // ������
    cout << "Bit Error Rate (BER): " << dBer << endl;
	return 0;
}