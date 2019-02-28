
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"

#define DR_MP3_IMPLEMENTATION


#include "dr_mp3.h"

#include "timing.h"
#include "loudnorm.h"

#ifndef MAX
#define MAX(a, b)                    (((a) > (b)) ? (a) : (b))
#endif

void wavWrite_f32(char *filename, float *buffer, int sampleRate, uint32_t totalSampleCount, uint32_t channels) {
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = channels;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 32;
    drwav *pWav = drwav_open_file_write(filename, &format);
    if (pWav) {
        drwav_uint64 samplesWritten = drwav_write(pWav, totalSampleCount, buffer);
        drwav_uninit(pWav);
        if (samplesWritten != totalSampleCount) {
            fprintf(stderr, "write file [%s] error.\n", filename);
            exit(1);
        }
    }
}

float *wavRead_f32(const char *filename, uint32_t *sampleRate, uint64_t *sampleCount, uint32_t *channels) {
    drwav_uint64 totalSampleCount = 0;
    float *input = drwav_open_file_and_read_pcm_frames_f32(filename, channels, sampleRate, &totalSampleCount);
    if (input == NULL) {
        drmp3_config pConfig;
        input = drmp3_open_file_and_read_f32(filename, &pConfig, &totalSampleCount);
        if (input != NULL) {
            *channels = pConfig.outputChannels;
            *sampleRate = pConfig.outputSampleRate;
        }
    }
    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    *sampleCount = totalSampleCount * (*channels);
    return input;
}


void splitpath(const char *path, char *drv, char *dir, char *name, char *ext) {
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    } else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}

void printUsage() {
    printf("usage:\n");
    printf("./Loudness input.wav\n");
    printf("./Loudness input.mp3\n");
    printf("or\n");
    printf("./Loudness input.wav output.wav\n");
    printf("./Loudness input.mp3 output.wav\n");
    printf("press any key to exit.\n");
    getchar();
}

void LoudnessProcess(char *in_file, char *out_file) {
    uint32_t sampleRate = 0;
    uint64_t sampleCount = 0;
    uint32_t channels = 0;
    float *input = wavRead_f32(in_file, &sampleRate, &sampleCount, &channels);
    if (input) {
        float *output = (float *) malloc(sampleCount * sizeof(float));
        if (output) {
            double startTime = now();
            float I = -16;
            float LRA = 11;
            float TP = -1.5f;
            float offset = 0;
            float measured_I = 0;
            float measured_TP = 99;
            float measured_LRA = 0;
            float measured_thresh = -70;
            LoudNormContext *analysisCtx = initLoudNormContext(sampleRate, channels, I, LRA, TP, measured_I,
                                                               measured_TP, measured_LRA, measured_thresh, offset);
            if (analysisCtx != NULL) {
                LoudNormFilter(analysisCtx, sampleCount, input, output);
                getAnalysisInfo(analysisCtx, &measured_I, &measured_LRA, &measured_TP, &measured_thresh,
                                &offset);
                uninitLoudNormContext(analysisCtx);
            }
            LoudNormContext *normContext = initLoudNormContext(sampleRate, channels, I, LRA, TP, measured_I,
                                                               measured_TP,
                                                               measured_LRA, measured_thresh, offset);
            if (normContext != NULL) {
                LoudNormFilter(normContext, sampleCount, input, output);
                uninitLoudNormContext(normContext);
            }
            double time_interval = calcElapsed(startTime, now());
            printf("time interval: %f ms\n ", (time_interval * 1000));
            wavWrite_f32(out_file, output, sampleRate, (uint32_t) sampleCount, channels);
            free(output);
        }
        free(input);
    }
}


int main(int argc, char *argv[]) {
    printf("Audio Processing\n");
    printf("blog:http://cpuimage.cnblogs.com/\n");
    printf("Audio Loudness Normalization Filter\n");
    if (argc < 2) {
        printUsage();
        return -1;
    }
    char *in_file = argv[1];
    if (argc > 2) {
        char *out_file = argv[2];
        LoudnessProcess(in_file, out_file);
    } else {
        char drive[3];
        char dir[256];
        char fname[256];
        char ext[256];
        char out_file[1024];
        splitpath(in_file, drive, dir, fname, ext);
        sprintf(out_file, "%s%s%s_out.wav", drive, dir, fname);
        LoudnessProcess(in_file, out_file);
    }

    printf("press any key to exit. \n");
    getchar();
    return 0;
}

#ifdef __cplusplus
}
#endif