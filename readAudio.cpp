#include "readAudio.h"

using namespace std;

/*
 * @brief 读取音频文件并返回音频数据和信息
 * @param filePath 音频文件路径
 * @return 包含音频数据和信息的 pair
            data: 音频数据（double 类型）
            info: 音频信息（采样率、声道数、总帧数、每个样本的位数）
 * @throws runtime_error 如果无法打开音频文件
 *
 * 注意：需要安装 libsndfile 库来实现 WAV 读写
 * 用法：
        pair<vector<double>, AudioInfo> result = readAudio(audioFile); // 读取音频文件
        vector<double> data = result.first; // 获取音频数据
        AudioInfo info = result.second;     // 获取音频信息
 */

// 读取音频文件并返回音频数据和信息
pair<vector<double>, AudioInfo> readAudio(const string& filePath) {
    // 打开音频文件
    SndfileHandle fileHandle(filePath);
    printf("文件名: %s\n", filePath.c_str());
    // 检查文件是否成功打开
    if (fileHandle.error()) {
        throw runtime_error("Error reading audio file: " + filePath);
    }

    // 获取音频信息
    AudioInfo info;
    info.sampleRate = fileHandle.samplerate();
    info.channels = fileHandle.channels();
    info.frames = fileHandle.frames();
    info.bitsPerSample = 16; // 假设音频文件是 16 位深度

    // 读取音频数据
    vector<double> data(fileHandle.frames() * fileHandle.channels());
	vector<double> buffer(fileHandle.frames() * fileHandle.channels());//！！！这里一定要用double类型
    fileHandle.read(buffer.data(), buffer.size());

    // 将音频数据转换为 double 类型
    for (size_t i = 0; i < buffer.size(); ++i) {
        data[i] = static_cast<double>(buffer[i]);
    }

    return { data, info };
}

// 写入音频文件的函数
void writeAudioFile(const string& filePath, const vector<double>& audioData, AudioInfo info) {
    // 配置音频文件格式
    SF_INFO sfInfo;
    sfInfo.samplerate = info.sampleRate; // 采样率
    sfInfo.channels = info.channels;    // 声道数
    sfInfo.format = SF_FORMAT_WAV;      // 文件格式：WAV

    // 根据位深度设置 PCM 格式
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

    // 打开文件以写入
    SNDFILE* outFile = sf_open(filePath.c_str(), SFM_WRITE, &sfInfo);
    if (!outFile) {
        cerr << "Error: Unable to open file for writing: " << sf_strerror(outFile) << endl;
        return;
    }

    //// 将 double 数据转换为 short 或 int（根据位深度）
    //vector<int> pcmData(audioData.size());
    //int maxAmplitude = (1 << (info.bitsPerSample - 1)) - 1; // 最大振幅值，例如 16 位为 32767
    //for (size_t i = 0; i < audioData.size(); ++i) {
    //    pcmData[i] = static_cast<int>(audioData[i] * maxAmplitude); // 归一化到 [-maxAmplitude, maxAmplitude]
    //}

    // 写入音频数据
    sf_count_t framesWritten = sf_write_double(outFile, audioData.data(), audioData.size());
    if (framesWritten != static_cast<sf_count_t>(audioData.size())) {
        cerr << "Error: Failed to write all frames to file." << endl;
    }

    // 关闭文件
    sf_close(outFile);
    cout << "Audio file written successfully to: " << filePath << endl;
}