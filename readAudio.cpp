#include "readAudio.h"

using namespace std;

/*
 * @brief ��ȡ��Ƶ�ļ���������Ƶ���ݺ���Ϣ
 * @param filePath ��Ƶ�ļ�·��
 * @return ������Ƶ���ݺ���Ϣ�� pair
            data: ��Ƶ���ݣ�double ���ͣ�
            info: ��Ƶ��Ϣ�������ʡ�����������֡����ÿ��������λ����
 * @throws runtime_error ����޷�����Ƶ�ļ�
 *
 * ע�⣺��Ҫ��װ libsndfile ����ʵ�� WAV ��д
 * �÷���
        pair<vector<double>, AudioInfo> result = readAudio(audioFile); // ��ȡ��Ƶ�ļ�
        vector<double> data = result.first; // ��ȡ��Ƶ����
        AudioInfo info = result.second;     // ��ȡ��Ƶ��Ϣ
 */

// ��ȡ��Ƶ�ļ���������Ƶ���ݺ���Ϣ
pair<vector<double>, AudioInfo> readAudio(const string& filePath) {
    // ����Ƶ�ļ�
    SndfileHandle fileHandle(filePath);
    printf("�ļ���: %s\n", filePath.c_str());
    // ����ļ��Ƿ�ɹ���
    if (fileHandle.error()) {
        throw runtime_error("Error reading audio file: " + filePath);
    }

    // ��ȡ��Ƶ��Ϣ
    AudioInfo info;
    info.sampleRate = fileHandle.samplerate();
    info.channels = fileHandle.channels();
    info.frames = fileHandle.frames();
    info.bitsPerSample = 16; // ������Ƶ�ļ��� 16 λ���

    // ��ȡ��Ƶ����
    vector<double> data(fileHandle.frames() * fileHandle.channels());
	vector<double> buffer(fileHandle.frames() * fileHandle.channels());//����������һ��Ҫ��double����
    fileHandle.read(buffer.data(), buffer.size());

    // ����Ƶ����ת��Ϊ double ����
    for (size_t i = 0; i < buffer.size(); ++i) {
        data[i] = static_cast<double>(buffer[i]);
    }

    return { data, info };
}

// д����Ƶ�ļ��ĺ���
void writeAudioFile(const string& filePath, const vector<double>& audioData, AudioInfo info) {
    // ������Ƶ�ļ���ʽ
    SF_INFO sfInfo;
    sfInfo.samplerate = info.sampleRate; // ������
    sfInfo.channels = info.channels;    // ������
    sfInfo.format = SF_FORMAT_WAV;      // �ļ���ʽ��WAV

    // ����λ������� PCM ��ʽ
    if (info.bitsPerSample == 16) {
        sfInfo.format |= SF_FORMAT_PCM_16;
    }
    else if (info.bitsPerSample == 24) {
        sfInfo.format |= SF_FORMAT_PCM_24;
    }
    else if (info.bitsPerSample == 32) {
        sfInfo.format |= SF_FORMAT_PCM_32;
    }
    else {
        cerr << "Unsupported bit depth: " << info.bitsPerSample << endl;
        return;
    }

    // ���ļ���д��
    SNDFILE* outFile = sf_open(filePath.c_str(), SFM_WRITE, &sfInfo);
    if (!outFile) {
        cerr << "Error: Unable to open file for writing: " << sf_strerror(outFile) << endl;
        return;
    }

    //// �� double ����ת��Ϊ short �� int������λ��ȣ�
    //vector<int> pcmData(audioData.size());
    //int maxAmplitude = (1 << (info.bitsPerSample - 1)) - 1; // ������ֵ������ 16 λΪ 32767
    //for (size_t i = 0; i < audioData.size(); ++i) {
    //    pcmData[i] = static_cast<int>(audioData[i] * maxAmplitude); // ��һ���� [-maxAmplitude, maxAmplitude]
    //}

    // д����Ƶ����
    sf_count_t framesWritten = sf_write_double(outFile, audioData.data(), audioData.size());
    if (framesWritten != static_cast<sf_count_t>(audioData.size())) {
        cerr << "Error: Failed to write all frames to file." << endl;
    }

    // �ر��ļ�
    sf_close(outFile);
    cout << "Audio file written successfully to: " << filePath << endl;
}