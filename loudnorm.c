/*
 * Copyright (c) 2016 Kyle Swanson <k@ylo.ph>.
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/* http://k.ylo.ph/2016/04/04/loudnorm.html */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>              /* You may have to define _USE_MATH_DEFINES if you use MSVC */

#include "loudnorm.h"

#define CHECK_ERROR(condition, errorcode, goto_point)                          \
    if ((condition)) {                                                         \
        errcode = (errorcode);                                                 \
        goto goto_point;                                                       \
    }

#define ALMOST_ZERO 0.000001

#define RELATIVE_GATE         (-10.0)
#define RELATIVE_GATE_FACTOR  pow(10.0, RELATIVE_GATE / 10.0)
#define MINUS_20DB            pow(10.0, -20.0 / 10.0)

struct FFEBUR128StateInternal {
    /** Filtered audio data (used as ring buffer). */
    double *audio_data;
    /** Size of audio_data array. */
    size_t audio_data_frames;
    /** Current index for audio_data. */
    size_t audio_data_index;
    /** How many frames are needed for a gating block. Will correspond to 400ms
     *  of audio at initialization, and 100ms after the first block (75% overlap
     *  as specified in the 2011 revision of BS1770). */
    unsigned long needed_frames;
    /** The channel map. Has as many elements as there are channels. */
    int *channel_map;
    /** How many samples fit in 100ms (rounded). */
    unsigned long samples_in_100ms;
    /** BS.1770 filter coefficients (nominator). */
    double b[5];
    /** BS.1770 filter coefficients (denominator). */
    double a[5];
    /** BS.1770 filter state. */
    double v[5][5];
    /** Histograms, used to calculate LRA. */
    unsigned long *block_energy_histogram;
    unsigned long *short_term_block_energy_histogram;
    /** Keeps track of when a new short term block is needed. */
    size_t short_term_frame_counter;
    /** Maximum sample peak, one per channel */
    double *sample_peak;
    /** The maximum window duration in ms. */
    unsigned long window;
    /** Data pointer array for interleaved data */
    void **data_ptrs;
};

static char histogram_init = 0;
static DECLARE_ALIGNED(32, double, histogram_energies)[1000];
static DECLARE_ALIGNED(32, double, histogram_energy_boundaries)[1001];

static void ebur128_init_filter(FFEBUR128State *st) {
    int i, j;

    double f0 = 1681.974450955533;
    double G = 3.999843853973347;
    double Q = 0.7071752369554196;

    double K = tan(M_PI * f0 / (double) st->samplerate);
    double Vh = pow(10.0, G / 20.0);
    double Vb = pow(Vh, 0.4996667741545416);

    double pb[3] = {0.0, 0.0, 0.0};
    double pa[3] = {1.0, 0.0, 0.0};
    double rb[3] = {1.0, -2.0, 1.0};
    double ra[3] = {1.0, 0.0, 0.0};

    double a0 = 1.0 + K / Q + K * K;
    pb[0] = (Vh + Vb * K / Q + K * K) / a0;
    pb[1] = 2.0 * (K * K - Vh) / a0;
    pb[2] = (Vh - Vb * K / Q + K * K) / a0;
    pa[1] = 2.0 * (K * K - 1.0) / a0;
    pa[2] = (1.0 - K / Q + K * K) / a0;

    f0 = 38.13547087602444;
    Q = 0.5003270373238773;
    K = tan(M_PI * f0 / (double) st->samplerate);

    ra[1] = 2.0 * (K * K - 1.0) / (1.0 + K / Q + K * K);
    ra[2] = (1.0 - K / Q + K * K) / (1.0 + K / Q + K * K);

    st->d->b[0] = pb[0] * rb[0];
    st->d->b[1] = pb[0] * rb[1] + pb[1] * rb[0];
    st->d->b[2] = pb[0] * rb[2] + pb[1] * rb[1] + pb[2] * rb[0];
    st->d->b[3] = pb[1] * rb[2] + pb[2] * rb[1];
    st->d->b[4] = pb[2] * rb[2];

    st->d->a[0] = pa[0] * ra[0];
    st->d->a[1] = pa[0] * ra[1] + pa[1] * ra[0];
    st->d->a[2] = pa[0] * ra[2] + pa[1] * ra[1] + pa[2] * ra[0];
    st->d->a[3] = pa[1] * ra[2] + pa[2] * ra[1];
    st->d->a[4] = pa[2] * ra[2];

    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) {
            st->d->v[i][j] = 0.0;
        }
    }
}

static int ebur128_init_channel_map(FFEBUR128State *st) {
    size_t i;
    st->d->channel_map =
            (int *) calloc(st->channels, sizeof(int));
    if (!st->d->channel_map)
        return -1;
    if (st->channels == 4) {
        st->d->channel_map[0] = FF_EBUR128_LEFT;
        st->d->channel_map[1] = FF_EBUR128_RIGHT;
        st->d->channel_map[2] = FF_EBUR128_LEFT_SURROUND;
        st->d->channel_map[3] = FF_EBUR128_RIGHT_SURROUND;
    } else if (st->channels == 5) {
        st->d->channel_map[0] = FF_EBUR128_LEFT;
        st->d->channel_map[1] = FF_EBUR128_RIGHT;
        st->d->channel_map[2] = FF_EBUR128_CENTER;
        st->d->channel_map[3] = FF_EBUR128_LEFT_SURROUND;
        st->d->channel_map[4] = FF_EBUR128_RIGHT_SURROUND;
    } else {
        for (i = 0; i < st->channels; ++i) {
            switch (i) {
                case 0:
                    st->d->channel_map[i] = FF_EBUR128_LEFT;
                    break;
                case 1:
                    st->d->channel_map[i] = FF_EBUR128_RIGHT;
                    break;
                case 2:
                    st->d->channel_map[i] = FF_EBUR128_CENTER;
                    break;
                case 3:
                    st->d->channel_map[i] = FF_EBUR128_UNUSED;
                    break;
                case 4:
                    st->d->channel_map[i] = FF_EBUR128_LEFT_SURROUND;
                    break;
                case 5:
                    st->d->channel_map[i] = FF_EBUR128_RIGHT_SURROUND;
                    break;
                default:
                    st->d->channel_map[i] = FF_EBUR128_UNUSED;
                    break;
            }
        }
    }
    return 0;
}

static inline void init_histogram(void) {
    int i;
    /* initialize static constants */
    histogram_energy_boundaries[0] = pow(10.0, (-70.0 + 0.691) / 10.0);
    for (i = 0; i < 1000; ++i) {
        histogram_energies[i] =
                pow(10.0, ((double) i / 10.0 - 69.95 + 0.691) / 10.0);
    }
    for (i = 1; i < 1001; ++i) {
        histogram_energy_boundaries[i] =
                pow(10.0, ((double) i / 10.0 - 70.0 + 0.691) / 10.0);
    }
}

static inline int ff_thread_once(char *control, void (*routine)(void)) {
    if (!*control) {
        routine();
        *control = 1;
    }
    return 0;
}

FFEBUR128State *ff_ebur128_init(size_t channels,
                                size_t samplerate,
                                unsigned long window, int mode) {
    int errcode;
    FFEBUR128State *st;

    st = (FFEBUR128State *) malloc(sizeof(FFEBUR128State));
    CHECK_ERROR(!st, 0, exit)
    st->d = (struct FFEBUR128StateInternal *)
            malloc(sizeof(struct FFEBUR128StateInternal));
    CHECK_ERROR(!st->d, 0, free_state)
    st->channels = channels;
    errcode = ebur128_init_channel_map(st);
    CHECK_ERROR(errcode, 0, free_internal)

    st->d->sample_peak =
            (double *) calloc(channels, sizeof(double));
    CHECK_ERROR(!st->d->sample_peak, 0, free_channel_map)

    st->samplerate = samplerate;
    st->d->samples_in_100ms = (st->samplerate + 5) / 10;
    st->mode = mode;
    if ((mode & FF_EBUR128_MODE_S) == FF_EBUR128_MODE_S) {
        st->d->window = FFMAX(window, 3000);
    } else if ((mode & FF_EBUR128_MODE_M) == FF_EBUR128_MODE_M) {
        st->d->window = FFMAX(window, 400);
    } else {
        goto free_sample_peak;
    }
    st->d->audio_data_frames = st->samplerate * st->d->window / 1000;
    if (st->d->audio_data_frames % st->d->samples_in_100ms) {
        /* round up to multiple of samples_in_100ms */
        st->d->audio_data_frames = st->d->audio_data_frames
                                   + st->d->samples_in_100ms
                                   - (st->d->audio_data_frames % st->d->samples_in_100ms);
    }
    st->d->audio_data =
            (double *) calloc(st->d->audio_data_frames,
                              st->channels * sizeof(double));
    CHECK_ERROR(!st->d->audio_data, 0, free_sample_peak)

    ebur128_init_filter(st);

    st->d->block_energy_histogram =
            calloc(1000, sizeof(unsigned long));
    CHECK_ERROR(!st->d->block_energy_histogram, 0, free_audio_data)
    st->d->short_term_block_energy_histogram =
            calloc(1000, sizeof(unsigned long));
    CHECK_ERROR(!st->d->short_term_block_energy_histogram, 0,
                free_block_energy_histogram)
    st->d->short_term_frame_counter = 0;

    /* the first block needs 400ms of audio data */
    st->d->needed_frames = st->d->samples_in_100ms * 4;
    /* start at the beginning of the buffer */
    st->d->audio_data_index = 0;

    if (ff_thread_once(&histogram_init, &init_histogram) != 0)
        goto free_short_term_block_energy_histogram;

    st->d->data_ptrs = calloc(channels, sizeof(void *));
    CHECK_ERROR(!st->d->data_ptrs, 0,
                free_short_term_block_energy_histogram);

    return st;

    free_short_term_block_energy_histogram:
    free(st->d->short_term_block_energy_histogram);
    free_block_energy_histogram:
    free(st->d->block_energy_histogram);
    free_audio_data:
    free(st->d->audio_data);
    free_sample_peak:
    free(st->d->sample_peak);
    free_channel_map:
    free(st->d->channel_map);
    free_internal:
    free(st->d);
    free_state:
    free(st);
    exit:
    return NULL;
}

void ff_ebur128_destroy(FFEBUR128State **st) {
    free((*st)->d->block_energy_histogram);
    free((*st)->d->short_term_block_energy_histogram);
    free((*st)->d->audio_data);
    free((*st)->d->channel_map);
    free((*st)->d->sample_peak);
    free((*st)->d->data_ptrs);
    free((*st)->d);
    free(*st);
    *st = NULL;
}

#define EBUR128_FILTER(type, scaling_factor)                                       \
static void ebur128_filter_##type(FFEBUR128State* st, const type** srcs,           \
                                  size_t src_index, size_t frames,                 \
                                  int stride) {                                    \
    double* audio_data = st->d->audio_data + st->d->audio_data_index;              \
    size_t i, c;                                                                   \
                                                                                   \
    if ((st->mode & FF_EBUR128_MODE_SAMPLE_PEAK) == FF_EBUR128_MODE_SAMPLE_PEAK) { \
        for (c = 0; c < st->channels; ++c) {                                       \
            double max = 0.0;                                                      \
            for (i = 0; i < frames; ++i) {                                         \
                type v = srcs[c][src_index + i * stride];                          \
                if (v > max) {                                                     \
                    max =        v;                                                \
                } else if (-v > max) {                                             \
                    max = -1.0 * v;                                                \
                }                                                                  \
            }                                                                      \
            max /= scaling_factor;                                                 \
            if (max > st->d->sample_peak[c]) st->d->sample_peak[c] = max;          \
        }                                                                          \
    }                                                                              \
    for (c = 0; c < st->channels; ++c) {                                           \
        int ci = st->d->channel_map[c] - 1;                                        \
        if (ci < 0) continue;                                                      \
        else if (ci == FF_EBUR128_DUAL_MONO - 1) ci = 0; /*dual mono */            \
        for (i = 0; i < frames; ++i) {                                             \
            st->d->v[ci][0] = (double) (srcs[c][src_index + i * stride] / scaling_factor) \
                         - st->d->a[1] * st->d->v[ci][1]                           \
                         - st->d->a[2] * st->d->v[ci][2]                           \
                         - st->d->a[3] * st->d->v[ci][3]                           \
                         - st->d->a[4] * st->d->v[ci][4];                          \
            audio_data[i * st->channels + c] =                                     \
                           st->d->b[0] * st->d->v[ci][0]                           \
                         + st->d->b[1] * st->d->v[ci][1]                           \
                         + st->d->b[2] * st->d->v[ci][2]                           \
                         + st->d->b[3] * st->d->v[ci][3]                           \
                         + st->d->b[4] * st->d->v[ci][4];                          \
            st->d->v[ci][4] = st->d->v[ci][3];                                     \
            st->d->v[ci][3] = st->d->v[ci][2];                                     \
            st->d->v[ci][2] = st->d->v[ci][1];                                     \
            st->d->v[ci][1] = st->d->v[ci][0];                                     \
        }                                                                          \
        st->d->v[ci][4] = fabs(st->d->v[ci][4]) < DBL_MIN ? 0.0 : st->d->v[ci][4]; \
        st->d->v[ci][3] = fabs(st->d->v[ci][3]) < DBL_MIN ? 0.0 : st->d->v[ci][3]; \
        st->d->v[ci][2] = fabs(st->d->v[ci][2]) < DBL_MIN ? 0.0 : st->d->v[ci][2]; \
        st->d->v[ci][1] = fabs(st->d->v[ci][1]) < DBL_MIN ? 0.0 : st->d->v[ci][1]; \
    }                                                                              \
}

EBUR128_FILTER(short, -((double) SHRT_MIN))

EBUR128_FILTER(int, -((double) INT_MIN))

EBUR128_FILTER(float, 1.0)

EBUR128_FILTER(double, 1.0)

static double ebur128_energy_to_loudness(double energy) {
    return 10 * (log(energy) / log(10.0)) - 0.691;
}

static size_t find_histogram_index(double energy) {
    size_t index_min = 0;
    size_t index_max = 1000;
    size_t index_mid;

    do {
        index_mid = (index_min + index_max) / 2;
        if (energy >= histogram_energy_boundaries[index_mid]) {
            index_min = index_mid;
        } else {
            index_max = index_mid;
        }
    } while (index_max - index_min != 1);

    return index_min;
}

static void ebur128_calc_gating_block(FFEBUR128State *st,
                                      size_t frames_per_block,
                                      double *optional_output) {
    size_t i, c;
    double sum = 0.0;
    double channel_sum;
    for (c = 0; c < st->channels; ++c) {
        if (st->d->channel_map[c] == FF_EBUR128_UNUSED)
            continue;
        channel_sum = 0.0;
        if (st->d->audio_data_index < frames_per_block * st->channels) {
            for (i = 0; i < st->d->audio_data_index / st->channels; ++i) {
                channel_sum += st->d->audio_data[i * st->channels + c] *
                               st->d->audio_data[i * st->channels + c];
            }
            for (i = st->d->audio_data_frames -
                     (frames_per_block -
                      st->d->audio_data_index / st->channels);
                 i < st->d->audio_data_frames; ++i) {
                channel_sum += st->d->audio_data[i * st->channels + c] *
                               st->d->audio_data[i * st->channels + c];
            }
        } else {
            for (i =
                         st->d->audio_data_index / st->channels - frames_per_block;
                 i < st->d->audio_data_index / st->channels; ++i) {
                channel_sum +=
                        st->d->audio_data[i * st->channels +
                                          c] * st->d->audio_data[i *
                                                                 st->channels +
                                                                 c];
            }
        }
        if (st->d->channel_map[c] == FF_EBUR128_Mp110 ||
            st->d->channel_map[c] == FF_EBUR128_Mm110 ||
            st->d->channel_map[c] == FF_EBUR128_Mp060 ||
            st->d->channel_map[c] == FF_EBUR128_Mm060 ||
            st->d->channel_map[c] == FF_EBUR128_Mp090 ||
            st->d->channel_map[c] == FF_EBUR128_Mm090) {
            channel_sum *= 1.41;
        } else if (st->d->channel_map[c] == FF_EBUR128_DUAL_MONO) {
            channel_sum *= 2.0;
        }
        sum += channel_sum;
    }
    sum /= (double) frames_per_block;
    if (optional_output) {
        *optional_output = sum;
    } else if (sum >= histogram_energy_boundaries[0]) {
        ++st->d->block_energy_histogram[find_histogram_index(sum)];
    }
}

int ff_ebur128_set_channel(FFEBUR128State *st,
                           unsigned int channel_number, int value) {
    if (channel_number >= st->channels) {
        return 1;
    }
    if (value == FF_EBUR128_DUAL_MONO &&
        (st->channels != 1 || channel_number != 0)) {
        return 1;
    }
    st->d->channel_map[channel_number] = value;
    return 0;
}

static int ebur128_energy_shortterm(FFEBUR128State *st, double *out);

#define FF_EBUR128_ADD_FRAMES_PLANAR(type)                                             \
void ff_ebur128_add_frames_planar_##type(FFEBUR128State* st, const type** srcs,        \
                                 size_t frames, int stride) {                          \
    size_t src_index = 0;                                                              \
    while (frames > 0) {                                                               \
        if (frames >= st->d->needed_frames) {                                          \
            ebur128_filter_##type(st, srcs, src_index, st->d->needed_frames, stride);  \
            src_index += st->d->needed_frames * stride;                                \
            frames -= st->d->needed_frames;                                            \
            st->d->audio_data_index += st->d->needed_frames * st->channels;            \
            /* calculate the new gating block */                                       \
            if ((st->mode & FF_EBUR128_MODE_I) == FF_EBUR128_MODE_I) {                 \
                ebur128_calc_gating_block(st, st->d->samples_in_100ms * 4, NULL);      \
            }                                                                          \
            if ((st->mode & FF_EBUR128_MODE_LRA) == FF_EBUR128_MODE_LRA) {             \
                st->d->short_term_frame_counter += st->d->needed_frames;               \
                if (st->d->short_term_frame_counter == st->d->samples_in_100ms * 30) { \
                    double st_energy;                                                  \
                    ebur128_energy_shortterm(st, &st_energy);                          \
                    if (st_energy >= histogram_energy_boundaries[0]) {                 \
                        ++st->d->short_term_block_energy_histogram[                    \
                                                    find_histogram_index(st_energy)];  \
                    }                                                                  \
                    st->d->short_term_frame_counter = st->d->samples_in_100ms * 20;    \
                }                                                                      \
            }                                                                          \
            /* 100ms are needed for all blocks besides the first one */                \
            st->d->needed_frames = st->d->samples_in_100ms;                            \
            /* reset audio_data_index when buffer full */                              \
            if (st->d->audio_data_index == st->d->audio_data_frames * st->channels) {  \
                st->d->audio_data_index = 0;                                           \
            }                                                                          \
        } else {                                                                       \
            ebur128_filter_##type(st, srcs, src_index, frames, stride);                \
            st->d->audio_data_index += frames * st->channels;                          \
            if ((st->mode & FF_EBUR128_MODE_LRA) == FF_EBUR128_MODE_LRA) {             \
                st->d->short_term_frame_counter += frames;                             \
            }                                                                          \
            st->d->needed_frames -= frames;                                            \
            frames = 0;                                                                \
        }                                                                              \
    }                                                                                  \
}

FF_EBUR128_ADD_FRAMES_PLANAR(short)

FF_EBUR128_ADD_FRAMES_PLANAR(int)

FF_EBUR128_ADD_FRAMES_PLANAR(float)

FF_EBUR128_ADD_FRAMES_PLANAR(double)

#define FF_EBUR128_ADD_FRAMES(type)                                            \
void ff_ebur128_add_frames_##type(FFEBUR128State* st, const type* src,         \
                                    size_t frames) {                           \
  int i;                                                                       \
  const type **buf = (const type**)st->d->data_ptrs;                           \
  for (i = 0; i < st->channels; i++)                                           \
    buf[i] = src + i;                                                          \
  ff_ebur128_add_frames_planar_##type(st, buf, frames, st->channels);          \
}

FF_EBUR128_ADD_FRAMES(short)

FF_EBUR128_ADD_FRAMES(int)

FF_EBUR128_ADD_FRAMES(float)

FF_EBUR128_ADD_FRAMES(double)

static int ebur128_calc_relative_threshold(FFEBUR128State **sts, size_t size,
                                           double *relative_threshold) {
    size_t i, j;
    int above_thresh_counter = 0;
    *relative_threshold = 0.0;

    for (i = 0; i < size; i++) {
        unsigned long *block_energy_histogram = sts[i]->d->block_energy_histogram;
        for (j = 0; j < 1000; ++j) {
            *relative_threshold += block_energy_histogram[j] * histogram_energies[j];
            above_thresh_counter += block_energy_histogram[j];
        }
    }

    if (above_thresh_counter != 0) {
        *relative_threshold /= (double) above_thresh_counter;
        *relative_threshold *= RELATIVE_GATE_FACTOR;
    }

    return above_thresh_counter;
}

static int ebur128_gated_loudness(FFEBUR128State **sts, size_t size,
                                  double *out) {
    double gated_loudness = 0.0;
    double relative_threshold;
    size_t above_thresh_counter;
    size_t i, j, start_index;

    for (i = 0; i < size; i++)
        if ((sts[i]->mode & FF_EBUR128_MODE_I) != FF_EBUR128_MODE_I)
            return -1;

    if (!ebur128_calc_relative_threshold(sts, size, &relative_threshold)) {
        *out = -HUGE_VAL;
        return 0;
    }

    above_thresh_counter = 0;
    if (relative_threshold < histogram_energy_boundaries[0]) {
        start_index = 0;
    } else {
        start_index = find_histogram_index(relative_threshold);
        if (relative_threshold > histogram_energies[start_index]) {
            ++start_index;
        }
    }
    for (i = 0; i < size; i++) {
        for (j = start_index; j < 1000; ++j) {
            gated_loudness += sts[i]->d->block_energy_histogram[j] *
                              histogram_energies[j];
            above_thresh_counter += sts[i]->d->block_energy_histogram[j];
        }
    }
    if (!above_thresh_counter) {
        *out = -HUGE_VAL;
        return 0;
    }
    gated_loudness /= (double) above_thresh_counter;
    *out = ebur128_energy_to_loudness(gated_loudness);
    return 0;
}

int ff_ebur128_relative_threshold(FFEBUR128State *st, double *out) {
    double relative_threshold;

    if ((st->mode & FF_EBUR128_MODE_I) != FF_EBUR128_MODE_I)
        return -1;

    if (!ebur128_calc_relative_threshold(&st, 1, &relative_threshold)) {
        *out = -70.0;
        return 0;
    }

    *out = ebur128_energy_to_loudness(relative_threshold);
    return 0;
}

int ff_ebur128_loudness_global(FFEBUR128State *st, double *out) {
    return ebur128_gated_loudness(&st, 1, out);
}

int ff_ebur128_loudness_global_multiple(FFEBUR128State **sts, size_t size,
                                        double *out) {
    return ebur128_gated_loudness(sts, size, out);
}

static int ebur128_energy_in_interval(FFEBUR128State *st,
                                      size_t interval_frames, double *out) {
    if (interval_frames > st->d->audio_data_frames) {
        return -1;
    }
    ebur128_calc_gating_block(st, interval_frames, out);
    return 0;
}

static int ebur128_energy_shortterm(FFEBUR128State *st, double *out) {
    return ebur128_energy_in_interval(st, st->d->samples_in_100ms * 30,
                                      out);
}

int ff_ebur128_loudness_momentary(FFEBUR128State *st, double *out) {
    double energy;
    int error = ebur128_energy_in_interval(st, st->d->samples_in_100ms * 4,
                                           &energy);
    if (error) {
        return error;
    } else if (energy <= 0.0) {
        *out = -HUGE_VAL;
        return 0;
    }
    *out = ebur128_energy_to_loudness(energy);
    return 0;
}

int ff_ebur128_loudness_shortterm(FFEBUR128State *st, double *out) {
    double energy;
    int error = ebur128_energy_shortterm(st, &energy);
    if (error) {
        return error;
    } else if (energy <= 0.0) {
        *out = -HUGE_VAL;
        return 0;
    }
    *out = ebur128_energy_to_loudness(energy);
    return 0;
}

int ff_ebur128_loudness_window(FFEBUR128State *st,
                               unsigned long window, double *out) {
    double energy;
    size_t interval_frames = st->samplerate * window / 1000;
    int error = ebur128_energy_in_interval(st, interval_frames, &energy);
    if (error) {
        return error;
    } else if (energy <= 0.0) {
        *out = -HUGE_VAL;
        return 0;
    }
    *out = ebur128_energy_to_loudness(energy);
    return 0;
}

/* EBU - TECH 3342 */
int ff_ebur128_loudness_range_multiple(FFEBUR128State **sts, size_t size,
                                       double *out) {
    size_t i, j;
    size_t stl_size;
    double stl_power, stl_integrated;
    /* High and low percentile energy */
    double h_en, l_en;
    unsigned long hist[1000] = {0};
    size_t percentile_low, percentile_high;
    size_t index;

    for (i = 0; i < size; ++i) {
        if (sts[i]) {
            if ((sts[i]->mode & FF_EBUR128_MODE_LRA) !=
                FF_EBUR128_MODE_LRA) {
                return -1;
            }
        }
    }

    stl_size = 0;
    stl_power = 0.0;
    for (i = 0; i < size; ++i) {
        if (!sts[i])
            continue;
        for (j = 0; j < 1000; ++j) {
            hist[j] += sts[i]->d->short_term_block_energy_histogram[j];
            stl_size += sts[i]->d->short_term_block_energy_histogram[j];
            stl_power += sts[i]->d->short_term_block_energy_histogram[j]
                         * histogram_energies[j];
        }
    }
    if (!stl_size) {
        *out = 0.0;
        return 0;
    }

    stl_power /= stl_size;
    stl_integrated = MINUS_20DB * stl_power;

    if (stl_integrated < histogram_energy_boundaries[0]) {
        index = 0;
    } else {
        index = find_histogram_index(stl_integrated);
        if (stl_integrated > histogram_energies[index]) {
            ++index;
        }
    }
    stl_size = 0;
    for (j = index; j < 1000; ++j) {
        stl_size += hist[j];
    }
    if (!stl_size) {
        *out = 0.0;
        return 0;
    }

    percentile_low = (size_t) ((stl_size - 1) * 0.1 + 0.5);
    percentile_high = (size_t) ((stl_size - 1) * 0.95 + 0.5);

    stl_size = 0;
    j = index;
    while (stl_size <= percentile_low) {
        stl_size += hist[j++];
    }
    l_en = histogram_energies[j - 1];
    while (stl_size <= percentile_high) {
        stl_size += hist[j++];
    }
    h_en = histogram_energies[j - 1];
    *out =
            ebur128_energy_to_loudness(h_en) -
            ebur128_energy_to_loudness(l_en);
    return 0;
}

int ff_ebur128_loudness_range(FFEBUR128State *st, double *out) {
    return ff_ebur128_loudness_range_multiple(&st, 1, out);
}

int ff_ebur128_sample_peak(FFEBUR128State *st,
                           unsigned int channel_number, double *out) {
    if ((st->mode & FF_EBUR128_MODE_SAMPLE_PEAK) !=
        FF_EBUR128_MODE_SAMPLE_PEAK) {
        return -1;
    } else if (channel_number >= st->channels) {
        return -1;
    }
    *out = st->d->sample_peak[channel_number];
    return 0;
}

static inline int frame_size(int sample_rate, int frame_len_msec) {
    const int frame_size = round((double) sample_rate * (frame_len_msec / 1000.0));
    return frame_size + (frame_size % 2);
}

static void init_gaussian_filter(LoudNormContext *s) {
    double total_weight = 0.0;
    const double sigma = 3.5;
    double adjust;
    int i;

    const int offset = 21 / 2;
    const double c1 = 1.0 / (sigma * sqrt(2.0 * M_PI));
    const double c2 = 2.0 * pow(sigma, 2.0);

    for (i = 0; i < 21; i++) {
        const int x = i - offset;
        s->weights[i] = c1 * exp(-(pow(x, 2.0) / c2));
        total_weight += s->weights[i];
    }

    adjust = 1.0 / total_weight;
    for (i = 0; i < 21; i++)
        s->weights[i] *= adjust;
}

static double gaussian_filter(LoudNormContext *s, int index) {
    double result = 0.;
    int i;

    index = index - 10 > 0 ? index - 10 : index + 20;
    for (i = 0; i < 21; i++)
        result += s->delta[((index + i) < 30) ? (index + i) : (index + i - 30)] * s->weights[i];

    return result;
}

static void
detect_peak(LoudNormContext *s, int offset, int nb_samples, int channels, int *peak_delta, float *peak_value) {
    int n, c, i, index;
    double ceiling;
    float *buf;

    *peak_delta = -1;
    buf = s->limiter_buf;
    ceiling = s->target_tp;
    int samples = s->samplerate / 100;
    index = s->limiter_buf_index + (offset * channels) + (samples * channels);
    if (index >= s->limiter_buf_size)
        index -= s->limiter_buf_size;

    if (s->frame_type == FIRST_FRAME) {
        for (c = 0; c < channels; c++)
            s->prev_smp[c] = fabsf(buf[index + c - channels]);
    }

    for (n = 0; n < nb_samples; n++) {
        for (c = 0; c < channels; c++) {
            double
                    this, next, max_peak;

            this = fabsf(buf[(index + c) < s->limiter_buf_size ? (index + c) : (index + c - s->limiter_buf_size)]);
            next = fabsf(
                    buf[(index + c + channels) < s->limiter_buf_size ? (index + c + channels) : (index + c + channels -
                                                                                                 s->limiter_buf_size)]);

            if ((s->prev_smp[c] <= this) && (next <= this) && (this > ceiling) && (n > 0)) {
                int detected;

                detected = 1;
                for (i = 2; i < 12; i++) {
                    next = fabsf(
                            buf[(index + c + (i * channels)) < s->limiter_buf_size ? (index + c + (i * channels)) : (
                                    index + c + (i * channels) - s->limiter_buf_size)]);
                    if (next > this) {
                        detected = 0;
                        break;
                    }
                }

                if (!detected)
                    continue;

                for (c = 0; c < channels; c++) {
                    if (c == 0 || fabsf(buf[index + c]) > max_peak)
                        max_peak = fabsf(buf[index + c]);

                    s->prev_smp[c] = fabsf(
                            buf[(index + c) < s->limiter_buf_size ? (index + c) : (index + c - s->limiter_buf_size)]);
                }

                *peak_delta = n;
                s->peak_index = index;
                *peak_value = max_peak;
                return;
            }

            s->prev_smp[c] = this;
        }

        index += channels;
        if (index >= s->limiter_buf_size)
            index -= s->limiter_buf_size;
    }
}

static void true_peak_limiter(LoudNormContext *s, float *out, int nb_samples, int channels) {
    int n, c, index, peak_delta, smp_cnt;
    float ceiling, peak_value;
    float *buf;
    int samples = s->samplerate / 100;
    buf = s->limiter_buf;
    ceiling = s->target_tp;
    index = s->limiter_buf_index;
    smp_cnt = 0;

    if (s->frame_type == FIRST_FRAME) {
        double max;

        max = 0.;
        for (n = 0; n < samples; n++) {
            for (c = 0; c < channels; c++) {
                max = fabsf(buf[c]) > max ? fabsf(buf[c]) : max;
            }
            buf += channels;
        }

        if (max > ceiling) {
            s->gain_reduction[1] = ceiling / max;
            s->limiter_state = LS_SUSTAIN;
            buf = s->limiter_buf;

            for (n = 0; n < samples; n++) {
                for (c = 0; c < channels; c++) {
                    double env;
                    env = s->gain_reduction[1];
                    buf[c] *= env;
                }
                buf += channels;
            }
        }

        buf = s->limiter_buf;
    }

    do {

        switch (s->limiter_state) {
            case LS_OUT:
                detect_peak(s, smp_cnt, nb_samples - smp_cnt, channels, &peak_delta, &peak_value);
                if (peak_delta != -1) {
                    s->env_cnt = 0;
                    smp_cnt += (peak_delta - s->attack_length);
                    s->gain_reduction[0] = 1.;
                    s->gain_reduction[1] = ceiling / peak_value;
                    s->limiter_state = LS_ATTACK;

                    s->env_index = s->peak_index - (s->attack_length * channels);
                    if (s->env_index < 0)
                        s->env_index += s->limiter_buf_size;

                    s->env_index += (s->env_cnt * channels);
                    if (s->env_index > s->limiter_buf_size)
                        s->env_index -= s->limiter_buf_size;

                } else {
                    smp_cnt = nb_samples;
                }
                break;

            case LS_ATTACK:
                for (; s->env_cnt < s->attack_length; s->env_cnt++) {
                    for (c = 0; c < channels; c++) {
                        double env;
                        env = s->gain_reduction[0] - ((double) s->env_cnt / (s->attack_length - 1) *
                                                      (s->gain_reduction[0] - s->gain_reduction[1]));
                        buf[s->env_index + c] *= env;
                    }

                    s->env_index += channels;
                    if (s->env_index >= s->limiter_buf_size)
                        s->env_index -= s->limiter_buf_size;

                    smp_cnt++;
                    if (smp_cnt >= nb_samples) {
                        s->env_cnt++;
                        break;
                    }
                }

                if (smp_cnt < nb_samples) {
                    s->env_cnt = 0;
                    s->attack_length = samples;
                    s->limiter_state = LS_SUSTAIN;
                }
                break;

            case LS_SUSTAIN:
                detect_peak(s, smp_cnt, nb_samples, channels, &peak_delta, &peak_value);
                if (peak_delta == -1) {
                    s->limiter_state = LS_RELEASE;
                    s->gain_reduction[0] = s->gain_reduction[1];
                    s->gain_reduction[1] = 1.;
                    s->env_cnt = 0;
                    break;
                } else {
                    double gain_reduction;
                    gain_reduction = ceiling / peak_value;

                    if (gain_reduction < s->gain_reduction[1]) {
                        s->limiter_state = LS_ATTACK;

                        s->attack_length = peak_delta;
                        if (s->attack_length <= 1)
                            s->attack_length = 2;

                        s->gain_reduction[0] = s->gain_reduction[1];
                        s->gain_reduction[1] = gain_reduction;
                        s->env_cnt = 0;
                        break;
                    }

                    for (s->env_cnt = 0; s->env_cnt < peak_delta; s->env_cnt++) {
                        for (c = 0; c < channels; c++) {
                            double env;
                            env = s->gain_reduction[1];
                            buf[s->env_index + c] *= env;
                        }

                        s->env_index += channels;
                        if (s->env_index >= s->limiter_buf_size)
                            s->env_index -= s->limiter_buf_size;

                        smp_cnt++;
                        if (smp_cnt >= nb_samples) {
                            s->env_cnt++;
                            break;
                        }
                    }
                }
                break;

            case LS_RELEASE:
                for (; s->env_cnt < s->release_length; s->env_cnt++) {
                    for (c = 0; c < channels; c++) {
                        double env;
                        env = s->gain_reduction[0] + (((double) s->env_cnt / (s->release_length - 1)) *
                                                      (s->gain_reduction[1] - s->gain_reduction[0]));
                        buf[s->env_index + c] *= env;
                    }

                    s->env_index += channels;
                    if (s->env_index >= s->limiter_buf_size)
                        s->env_index -= s->limiter_buf_size;

                    smp_cnt++;
                    if (smp_cnt >= nb_samples) {
                        s->env_cnt++;
                        break;
                    }
                }

                if (smp_cnt < nb_samples) {
                    s->env_cnt = 0;
                    s->limiter_state = LS_OUT;
                }

                break;
        }

    } while (smp_cnt < nb_samples);

    for (n = 0; n < nb_samples; n++) {
        for (c = 0; c < channels; c++) {
            out[c] = buf[index + c];
            if (fabsf(out[c]) > ceiling) {
                out[c] = ceiling * (out[c] < 0 ? -1 : 1);
            }
        }
        out += channels;
        index += channels;
        if (index >= s->limiter_buf_size)
            index -= s->limiter_buf_size;
    }
}

static int filter_frame(LoudNormContext *s, const float *in_data, float *out_data, size_t in_nb_samples
) {
    const float *src;
    float *dst;
    float *buf;
    float *limiter_buf;
    int i, n, c, subframe_length, src_index;
    double gain, gain_next, env_global, env_shortterm,
            global, shortterm, lra, relative_threshold;

    if (in_data != out_data)
        memcpy(out_data, in_data, sizeof(float) * in_nb_samples);
    if (s->pts == NOPTS_VALUE) {
        s->pts = 0;
    }

    src = in_data;
    dst = out_data;
    buf = s->buf;
    limiter_buf = s->limiter_buf;

    ff_ebur128_add_frames_float(s->r128_in, src, in_nb_samples);

    if (s->frame_type == FIRST_FRAME && in_nb_samples < frame_size(s->samplerate, 3000)) {
        double offset, offset_tp, true_peak;

        ff_ebur128_loudness_global(s->r128_in, &global);
        for (c = 0; c < s->channels; c++) {
            double tmp;
            ff_ebur128_sample_peak(s->r128_in, c, &tmp);
            if (c == 0 || tmp > true_peak)
                true_peak = tmp;
        }

        offset = s->target_i - global;
        offset_tp = true_peak + offset;
        s->offset = offset_tp < s->target_tp ? offset : s->target_tp - true_peak;
        s->offset = pow(10., s->offset / 20.);
        s->frame_type = LINEAR_MODE;
    }

    switch (s->frame_type) {
        case FIRST_FRAME:
            for (n = 0; n < in_nb_samples; n++) {
                for (c = 0; c < s->channels; c++) {
                    buf[s->buf_index + c] = src[c];
                }
                src += s->channels;
                s->buf_index += s->channels;
            }

            ff_ebur128_loudness_shortterm(s->r128_in, &shortterm);

            if (shortterm < s->measured_thresh) {
                s->above_threshold = 0;
                env_shortterm = shortterm <= -70. ? 0. : s->target_i - s->measured_i;
            } else {
                s->above_threshold = 1;
                env_shortterm = shortterm <= -70. ? 0. : s->target_i - shortterm;
            }

            for (n = 0; n < 30; n++)
                s->delta[n] = pow(10., env_shortterm / 20.);
            s->prev_delta = s->delta[s->index];

            s->buf_index =
            s->limiter_buf_index = 0;

            for (n = 0; n < (s->limiter_buf_size / s->channels); n++) {
                for (c = 0; c < s->channels; c++) {
                    limiter_buf[s->limiter_buf_index + c] = buf[s->buf_index + c] * s->delta[s->index] * s->offset;
                }
                s->limiter_buf_index += s->channels;
                if (s->limiter_buf_index >= s->limiter_buf_size)
                    s->limiter_buf_index -= s->limiter_buf_size;

                s->buf_index += s->channels;
            }

            subframe_length = frame_size(s->samplerate, 100);
            true_peak_limiter(s, dst, subframe_length, s->channels);
            ff_ebur128_add_frames_float(s->r128_out, dst, subframe_length);

            s->pts += subframe_length;

            s->frame_type = INNER_FRAME;
            break;

        case INNER_FRAME:
            gain = gaussian_filter(s, s->index + 10 < 30 ? s->index + 10 : s->index + 10 - 30);
            gain_next = gaussian_filter(s, s->index + 11 < 30 ? s->index + 11 : s->index + 11 - 30);

            for (n = 0; n < in_nb_samples; n++) {
                for (c = 0; c < s->channels; c++) {
                    buf[s->prev_buf_index + c] = src[c];
                    limiter_buf[s->limiter_buf_index + c] =
                            buf[s->buf_index + c] * (gain + (((double) n / in_nb_samples) * (gain_next - gain))) *
                            s->offset;
                }
                src += s->channels;

                s->limiter_buf_index += s->channels;
                if (s->limiter_buf_index >= s->limiter_buf_size)
                    s->limiter_buf_index -= s->limiter_buf_size;

                s->prev_buf_index += s->channels;
                if (s->prev_buf_index >= s->buf_size)
                    s->prev_buf_index -= s->buf_size;

                s->buf_index += s->channels;
                if (s->buf_index >= s->buf_size)
                    s->buf_index -= s->buf_size;
            }

            subframe_length = (frame_size(s->samplerate, 100) - in_nb_samples) * s->channels;
            s->limiter_buf_index = s->limiter_buf_index + subframe_length < s->limiter_buf_size ? s->limiter_buf_index +
                                                                                                  subframe_length :
                                   s->limiter_buf_index + subframe_length - s->limiter_buf_size;

            true_peak_limiter(s, dst, in_nb_samples, s->channels);
            ff_ebur128_add_frames_float(s->r128_out, dst, in_nb_samples);

            ff_ebur128_loudness_range(s->r128_in, &lra);
            ff_ebur128_loudness_global(s->r128_in, &global);
            ff_ebur128_loudness_shortterm(s->r128_in, &shortterm);
            ff_ebur128_relative_threshold(s->r128_in, &relative_threshold);

            if (s->above_threshold == 0) {
                double shortterm_out;

                if (shortterm > s->measured_thresh)
                    s->prev_delta *= 1.0058;

                ff_ebur128_loudness_shortterm(s->r128_out, &shortterm_out);
                if (shortterm_out >= s->target_i)
                    s->above_threshold = 1;
            }

            if (shortterm < relative_threshold || shortterm <= -70. || s->above_threshold == 0) {
                s->delta[s->index] = s->prev_delta;
            } else {
                env_global =
                        fabsf(shortterm - global) < (s->target_lra / 2.) ? shortterm - global : (s->target_lra / 2.) *
                                                                                                ((shortterm - global) <
                                                                                                 0
                                                                                                 ? -1 : 1);
                env_shortterm = s->target_i - shortterm;
                s->delta[s->index] = pow(10., (env_global + env_shortterm) / 20.);
            }

            s->prev_delta = s->delta[s->index];
            s->index++;
            if (s->index >= 30)
                s->index -= 30;
            s->prev_nb_samples = in_nb_samples;
            s->pts += in_nb_samples;
            break;

        case FINAL_FRAME:
            gain = gaussian_filter(s, s->index + 10 < 30 ? s->index + 10 : s->index + 10 - 30);
            s->limiter_buf_index = 0;
            src_index = 0;

            for (n = 0; n < s->limiter_buf_size / s->channels; n++) {
                for (c = 0; c < s->channels; c++) {
                    s->limiter_buf[s->limiter_buf_index + c] = src[src_index + c] * gain * s->offset;
                }
                src_index += s->channels;

                s->limiter_buf_index += s->channels;
                if (s->limiter_buf_index >= s->limiter_buf_size)
                    s->limiter_buf_index -= s->limiter_buf_size;
            }

            subframe_length = frame_size(s->samplerate, 100);
            for (i = 0; i < in_nb_samples / subframe_length; i++) {
                true_peak_limiter(s, dst, subframe_length, s->channels);

                for (n = 0; n < subframe_length; n++) {
                    for (c = 0; c < s->channels; c++) {
                        if (src_index < (in_nb_samples * s->channels)) {
                            limiter_buf[s->limiter_buf_index + c] = src[src_index + c] * gain * s->offset;
                        } else {
                            limiter_buf[s->limiter_buf_index + c] = 0.;
                        }
                    }

                    if (src_index < (in_nb_samples * s->channels))
                        src_index += s->channels;

                    s->limiter_buf_index += s->channels;
                    if (s->limiter_buf_index >= s->limiter_buf_size)
                        s->limiter_buf_index -= s->limiter_buf_size;
                }

                dst += (subframe_length * s->channels);
            }

            dst = out_data;
            ff_ebur128_add_frames_float(s->r128_out, dst, in_nb_samples);
            break;

        case LINEAR_MODE:
            for (n = 0; n < in_nb_samples; n++) {
                for (c = 0; c < s->channels; c++) {
                    dst[c] = src[c] * s->offset;
                }
                src += s->channels;
                dst += s->channels;
            }

            dst = out_data;
            ff_ebur128_add_frames_float(s->r128_out, dst, in_nb_samples);
            s->pts += in_nb_samples;
            break;
    }
    return 1;
}


static int request_frame(LoudNormContext *s, float *in_data, float *frame_data, size_t frameCount) {
    int ret = 0;
    if (s->frame_type == INNER_FRAME) {
        float *src;
        float *buf;
        int nb_samples, n, c, offset;

        nb_samples = (s->buf_size / s->channels) - s->prev_nb_samples;
        nb_samples -= (frame_size(s->samplerate, 100) - s->prev_nb_samples);

        buf = s->buf;
        src = frame_data;
        offset = ((s->limiter_buf_size / s->channels) - s->prev_nb_samples) * s->channels;
        offset -= (frame_size(s->samplerate, 100) - s->prev_nb_samples) * s->channels;
        s->buf_index = s->buf_index - offset < 0 ? s->buf_index - offset + s->buf_size : s->buf_index - offset;
        for (n = 0; n < nb_samples; n++) {
            for (c = 0; c < s->channels; c++) {
                src[c] = buf[s->buf_index + c];
            }
            src += s->channels;
            s->buf_index += s->channels;
            if (s->buf_index >= s->buf_size)
                s->buf_index -= s->buf_size;
        }
        s->frame_type = FINAL_FRAME;
        ret = filter_frame(s, in_data, frame_data, frameCount);
    }
    return ret;
}

static int config_input(LoudNormContext *s, size_t channels, size_t sample_rate) {
    s->r128_in = ff_ebur128_init(channels, sample_rate, 0,
                                 FF_EBUR128_MODE_I | FF_EBUR128_MODE_S | FF_EBUR128_MODE_LRA |
                                 FF_EBUR128_MODE_SAMPLE_PEAK);
    if (!s->r128_in)
        return -1;

    s->r128_out = ff_ebur128_init(channels, sample_rate, 0,
                                  FF_EBUR128_MODE_I | FF_EBUR128_MODE_S | FF_EBUR128_MODE_LRA |
                                  FF_EBUR128_MODE_SAMPLE_PEAK);
    if (!s->r128_out)
        return -1;

    if (channels == 1 && s->dual_mono) {
        ff_ebur128_set_channel(s->r128_in, 0, FF_EBUR128_DUAL_MONO);
        ff_ebur128_set_channel(s->r128_out, 0, FF_EBUR128_DUAL_MONO);
    }

    s->buf_size = frame_size(sample_rate, 3000) * channels;
    s->buf = calloc(s->buf_size, sizeof(*s->buf));
    if (!s->buf)
        return -1;

    s->limiter_buf_size = frame_size(sample_rate, 210) * channels;
    s->limiter_buf = calloc(s->buf_size, sizeof(*s->limiter_buf));
    if (!s->limiter_buf)
        return -1;

    s->prev_smp = calloc(channels, sizeof(*s->prev_smp));
    if (!s->prev_smp)
        return -1;

    init_gaussian_filter(s);

    s->pts = NOPTS_VALUE;
    s->buf_index =
    s->prev_buf_index =
    s->limiter_buf_index = 0;
    s->channels = channels;
    s->index = 1;
    s->samplerate = sample_rate;
    s->limiter_state = LS_OUT;
    s->offset = pow(10., s->offset / 20.);
    s->target_tp = pow(10., s->target_tp / 20.);
    s->attack_length = frame_size(sample_rate, 10);
    s->release_length = frame_size(sample_rate, 100);

    return 0;
}


static int init(LoudNormContext *s) {
    s->frame_type = FIRST_FRAME;

    if (s->linear) {
        double offset, offset_tp;
        offset = s->target_i - s->measured_i;
        offset_tp = s->measured_tp + offset;

        if (s->measured_tp != 99 && s->measured_thresh != -70 && s->measured_lra != 0 && s->measured_i != 0) {
            if ((offset_tp <= s->target_tp) && (s->measured_lra <= s->target_lra)) {
                s->frame_type = LINEAR_MODE;
                s->offset = offset;
            }
        }
    }

    return 0;
}

void freep(void *arg) {
    void *val;

    memcpy(&val, arg, sizeof(val));
    memcpy(arg, &(void *) {NULL}, sizeof(val));
    free(val);
}

void getAnalysisInfo(const LoudNormContext *s, float *I,
                     float *LRA,
                     float *TP,
                     float *thresh,
                     float *offset) {
    double i_in, i_out, lra_in, lra_out, thresh_in, thresh_out, tp_in = 0, tp_out = 0;
    uint32_t c;

    ff_ebur128_loudness_range(s->r128_in, &lra_in);
    ff_ebur128_loudness_global(s->r128_in, &i_in);
    ff_ebur128_relative_threshold(s->r128_in, &thresh_in);
    for (c = 0; c < s->channels; c++) {
        double tmp;
        ff_ebur128_sample_peak(s->r128_in, c, &tmp);
        if ((c == 0) || (tmp > tp_in))
            tp_in = tmp;
    }

    ff_ebur128_loudness_range(s->r128_out, &lra_out);
    ff_ebur128_loudness_global(s->r128_out, &i_out);
    ff_ebur128_relative_threshold(s->r128_out, &thresh_out);
    for (c = 0; c < s->channels; c++) {
        double tmp;
        ff_ebur128_sample_peak(s->r128_out, c, &tmp);
        if ((c == 0) || (tmp > tp_out))
            tp_out = tmp;
    }
    *I = (float) i_in;
    *LRA = (float) lra_in;
    *TP = (float) (20. * log10(tp_in));
    *offset = (float) (s->target_i - i_out);
    *thresh = (float) thresh_in;
}

void uninitLoudNormContext(LoudNormContext *s) {
    if (s->r128_in)
        ff_ebur128_destroy(&s->r128_in);
    if (s->r128_out)
        ff_ebur128_destroy(&s->r128_out);
    freep(&s->limiter_buf);
    freep(&s->prev_smp);
    freep(&s->buf);
    free(s);
}


void LoudNormFilter(LoudNormContext *ctx, uint64_t sampleCount, float *input, float *output) {
    size_t sec = 3;
    size_t nr_frames = sec * ctx->samplerate;
    size_t step_size = (nr_frames * ctx->channels);
    size_t nSample = sampleCount / step_size;
    float *in = input;
    float *out = output;

    for (int i = 0; i < nSample; i++) {
        filter_frame(ctx, in, out, nr_frames);
        request_frame(ctx, in, out, nr_frames);
        in += step_size;
        out += step_size;
    }
    size_t last_frames = (sampleCount / ctx->channels) % nr_frames;
    if (last_frames != 0) {
        float *pad_buffer = (float *) calloc(step_size * 2, sizeof(float));
        if (pad_buffer != NULL) {
            float *input_buffer = pad_buffer;
            float *output_buffer = pad_buffer + step_size;
            size_t size = last_frames * ctx->channels;
            memcpy(input_buffer, in, size * sizeof(float));
            filter_frame(ctx, input_buffer, output_buffer, nr_frames);
            request_frame(ctx, input_buffer, output_buffer, nr_frames);
            memcpy(out, output_buffer, size * sizeof(float));
            free(pad_buffer);
        }
    }
}

LoudNormContext *initLoudNormContext(size_t sampleRate, size_t channels, float I, float LRA, float TP, float measured_I,
                                     float measured_TP,
                                     float measured_LRA, float measured_thresh, float Offset) {
    LoudNormContext *ctx = (LoudNormContext *) calloc(sizeof(LoudNormContext), 1);
    if (ctx == NULL)
        return NULL;
    // set integrated loudness target   [ -70., -5. ] def:-24
    ctx->target_i = I;
    // set loudness range target [ 1., 20. ] def: 7
    ctx->target_lra = LRA;
    // set maximum true peak  [ -9., 0. ]  def: -2
    ctx->target_tp = TP;//
    //measured IL of input file   [ -99., 0. ] def: 0
    ctx->measured_i = measured_I;// [ -99., 0. ]
    //measured LRA of input file [ 0., 99. ]  def: 0
    ctx->measured_lra = measured_LRA;
    //measured true peak of input file  [ -99., 99. ] def: 99
    ctx->measured_tp = measured_TP;
    //measured threshold of input file [ -99., 0. ] def:-70.
    ctx->measured_thresh = measured_thresh;
    //set offset gain [ -99., 99. ] def: 0
    ctx->offset = Offset;
    //normalize linearly if possible" [ 0, 1 ] def: 1
    ctx->linear = 1;
    //treat mono input as dual-mono [ 0, 1 ] def: 0
    ctx->dual_mono = 0;
    init(ctx);
    int ret = config_input(ctx, channels, sampleRate);
    if (ret == -1) {
        free(ctx);
        return NULL;
    }
    return ctx;
}

#ifdef __cplusplus
}
#endif
