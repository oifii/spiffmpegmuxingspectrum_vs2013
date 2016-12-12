/*
 * Copyright (c) 2016 Stephane Poirier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
/*
 * 2016jan15, modifying ffmpeg muxing.c so it reads audio from file relying upon
 *            libsndfile, the lenght of the video will now be determined by audio
 *            file.
 * 2016jan16, modifying ffmpeg muxing.c further so it creates graphics from input
 *            audio file provided. the drawing code integrated from spispectrumlive_pa
 * 2016jan18, problem integrating mode 0, fft() returns all zeros.
 * 2016jan18, problem solved, mode 0 is working swell, in c fabs() must be used instead of abs()
 *            when manipulating double
 *
 * stephane.poirier@oifii.org or spi@oifii.org
 */
/*
 * Copyright (c) 2003 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * @file
 * libavformat API example.
 *
 * Output a media file in any supported libavformat format. The default
 * codecs are used.
 * @example muxing.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//spi, begin
//#define snprintf(buf,len, format,...) _snprintf_s(buf, len,len, format, __VA_ARGS__)
#ifdef __cplusplus
extern "C" {
#endif
//spi, end
#include <libavutil/avassert.h>
#include <libavutil/channel_layout.h>
#include <libavutil/opt.h>
#include <libavutil/mathematics.h>
#include <libavutil/timestamp.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
#include <libswresample/swresample.h>
//spi, begin
#ifdef __cplusplus
}
#endif
#include <sndfile.h>
#include <assert.h>
//fourier.h is the addon from audio programming book to the rfftw library,
//see audio programming book page 536 for details. in short, fourier.cpp
//wraps the rfftw by providing 2 functions: fft() and ifft().
//fourier.h also depends on libsndfile so it makes rfftw.lib depends on
//libsndfile for compiling (no need for linking it if you don't use it
//elsewhere).
//fft() function can only be called always with the same sample size N,
//this because within fft() implementation rfftw_create_plan() is called
//only once (the first time fft() is called).
#include <fourier.h> //in rfftw.lib (static library)

//spi, end


#define STREAM_FRAME_RATE 25 /* 25 images/s */
#define STREAM_PIX_FMT    AV_PIX_FMT_YUV420P /* default pix_fmt */

#define SCALE_FLAGS SWS_BICUBIC

//spi, begin
//#define STREAM_DURATION   10.0
int64_t global_STREAM_DURATION;
SNDFILE* global_pSNDFILE=NULL;
SF_INFO global_mySF_INFO;
int global_NUM_CHANNELS=2;
int16_t global_intbuffer[44100*2*10];
float global_floatbuffer[44100*2*10];
int specmode;
int RGBtoYCbCr_Y(int red, int green, int blue)
{
    float R = red;
    float G = green;
    float B = blue;
    int Y = (int)(0.257*R + 0.504*G + 0.098*B + 16);
    return Y;
}
int RGBtoYCbCr_Cb(int red, int green, int blue)
{
    float R = red;
    float G = green;
    float B = blue;
    int Cb = (int)(-0.148*R - 0.291*G + 0.439*B + 128);
    return Cb;
}
int RGBtoYCbCr_Cr(int red, int green, int blue)
{
    float R = red;
    float G = green;
    float B = blue;
    int Cr = (int)(0.439*R - 0.368*G - 0.071*B + 128);
    return Cr;
}
//spi, end



// a wrapper around a single output AVStream
typedef struct OutputStream {
    AVStream *st;

    /* pts of the next frame that will be generated */
    int64_t next_pts;
    int samples_count;

    AVFrame *frame;
    AVFrame *tmp_frame;

    float t, tincr, tincr2;

    struct SwsContext *sws_ctx;
    struct SwrContext *swr_ctx;
} OutputStream;

static void log_packet(const AVFormatContext *fmt_ctx, const AVPacket *pkt)
{
    AVRational *time_base = &fmt_ctx->streams[pkt->stream_index]->time_base;

    printf("pts:%s pts_time:%s dts:%s dts_time:%s duration:%s duration_time:%s stream_index:%d\n",
           av_ts2str(pkt->pts), av_ts2timestr(pkt->pts, time_base),
           av_ts2str(pkt->dts), av_ts2timestr(pkt->dts, time_base),
           av_ts2str(pkt->duration), av_ts2timestr(pkt->duration, time_base),
           pkt->stream_index);
}

static int write_frame(AVFormatContext *fmt_ctx, const AVRational *time_base, AVStream *st, AVPacket *pkt)
{
    /* rescale output packet timestamp values from codec to stream timebase */
    av_packet_rescale_ts(pkt, *time_base, st->time_base);
    pkt->stream_index = st->index;

    /* Write the compressed frame to the media file. */
    log_packet(fmt_ctx, pkt);
    return av_interleaved_write_frame(fmt_ctx, pkt);
}

/* Add an output stream. */
static void add_stream(OutputStream *ost, AVFormatContext *oc,
                       AVCodec **codec,
                       enum AVCodecID codec_id)
{
    AVCodecContext *c;
    int i;

    /* find the encoder */
    *codec = avcodec_find_encoder(codec_id);
    if (!(*codec)) {
        fprintf(stderr, "Could not find encoder for '%s'\n",
                avcodec_get_name(codec_id));
        exit(1);
    }

    ost->st = avformat_new_stream(oc, *codec);
    if (!ost->st) {
        fprintf(stderr, "Could not allocate stream\n");
        exit(1);
    }
    ost->st->id = oc->nb_streams-1;
    c = ost->st->codec;

    switch ((*codec)->type) {
    case AVMEDIA_TYPE_AUDIO:
        c->sample_fmt  = (*codec)->sample_fmts ?
            (*codec)->sample_fmts[0] : AV_SAMPLE_FMT_FLTP;
        c->bit_rate    = 64000; //original
        //c->bit_rate    = 128000; //spi
        c->sample_rate = 44100;
        if ((*codec)->supported_samplerates) {
            c->sample_rate = (*codec)->supported_samplerates[0];
            for (i = 0; (*codec)->supported_samplerates[i]; i++) {
                if ((*codec)->supported_samplerates[i] == 44100)
                    c->sample_rate = 44100;
            }
        }
        c->channels        = av_get_channel_layout_nb_channels(c->channel_layout);
        c->channel_layout = AV_CH_LAYOUT_STEREO;
        if ((*codec)->channel_layouts) {
            c->channel_layout = (*codec)->channel_layouts[0];
            for (i = 0; (*codec)->channel_layouts[i]; i++) {
                if ((*codec)->channel_layouts[i] == AV_CH_LAYOUT_STEREO)
                    c->channel_layout = AV_CH_LAYOUT_STEREO;
            }
        }
        c->channels        = av_get_channel_layout_nb_channels(c->channel_layout);
        ost->st->time_base = (AVRational){ 1, c->sample_rate };
        break;

    case AVMEDIA_TYPE_VIDEO:
        c->codec_id = codec_id;

        c->bit_rate = 400000;
        /* Resolution must be a multiple of two. */
        //spi, begin
        //c->width    = 352; //original
        //c->height   = 288; //original
        //c->width    = 1920; //spi
        //c->height   = 1080; //spi
        c->width    = 1024; //spi
        c->height   = 576; //spi
        //spi, end
        /* timebase: This is the fundamental unit of time (in seconds) in terms
         * of which frame timestamps are represented. For fixed-fps content,
         * timebase should be 1/framerate and timestamp increments should be
         * identical to 1. */
        ost->st->time_base = (AVRational){ 1, STREAM_FRAME_RATE };
        c->time_base       = ost->st->time_base;

        c->gop_size      = 12; /* emit one intra frame every twelve frames at most */
        c->pix_fmt       = STREAM_PIX_FMT;
        if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
            /* just for testing, we also add B frames */
            c->max_b_frames = 2;
        }
        if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
            /* Needed to avoid using macroblocks in which some coeffs overflow.
             * This does not happen with normal video, it just happens here as
             * the motion of the chroma plane does not match the luma plane. */
            c->mb_decision = 2;
        }
    break;

    default:
        break;
    }

    /* Some formats want stream headers to be separate. */
    if (oc->oformat->flags & AVFMT_GLOBALHEADER)
        c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}

/**************************************************************/
/* audio output */

static AVFrame *alloc_audio_frame(enum AVSampleFormat sample_fmt,
                                  uint64_t channel_layout,
                                  int sample_rate, int nb_samples)
{
    AVFrame *frame = av_frame_alloc();
    int ret;

    if (!frame) {
        fprintf(stderr, "Error allocating an audio frame\n");
        exit(1);
    }

    frame->format = sample_fmt;
    frame->channel_layout = channel_layout;
    frame->sample_rate = sample_rate;
    frame->nb_samples = nb_samples;

    if (nb_samples) {
        ret = av_frame_get_buffer(frame, 0);
        if (ret < 0) {
            fprintf(stderr, "Error allocating an audio buffer\n");
            exit(1);
        }
    }

    return frame;
}

static void open_audio(AVFormatContext *oc, AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg)
{
    AVCodecContext *c;
    int nb_samples;
    int ret;
    AVDictionary *opt = NULL;

    c = ost->st->codec;

    /* open it */
    av_dict_copy(&opt, opt_arg, 0);
    ret = avcodec_open2(c, codec, &opt);
    av_dict_free(&opt);
    if (ret < 0) {
        fprintf(stderr, "Could not open audio codec: %s\n", av_err2str(ret));
        exit(1);
    }

    /* init signal generator */
    ost->t     = 0;
    ost->tincr = 2 * M_PI * 110.0 / c->sample_rate;
    /* increment frequency by 110 Hz per second */
    ost->tincr2 = 2 * M_PI * 110.0 / c->sample_rate / c->sample_rate;

    if (c->codec->capabilities & AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
        nb_samples = 10000;
    else
        nb_samples = c->frame_size; //original
        //nb_samples = 1920; //spi, doesn't work

    ost->frame     = alloc_audio_frame(c->sample_fmt, c->channel_layout,
                                       c->sample_rate, nb_samples);
    ost->tmp_frame = alloc_audio_frame(AV_SAMPLE_FMT_S16, c->channel_layout,
                                       c->sample_rate, nb_samples);

    /* create resampler context */
        ost->swr_ctx = swr_alloc();
        if (!ost->swr_ctx) {
            fprintf(stderr, "Could not allocate resampler context\n");
            exit(1);
        }

        /* set options */
        av_opt_set_int       (ost->swr_ctx, "in_channel_count",   c->channels,       0);
        av_opt_set_int       (ost->swr_ctx, "in_sample_rate",     c->sample_rate,    0);
        av_opt_set_sample_fmt(ost->swr_ctx, "in_sample_fmt",      AV_SAMPLE_FMT_S16, 0);
        av_opt_set_int       (ost->swr_ctx, "out_channel_count",  c->channels,       0);
        av_opt_set_int       (ost->swr_ctx, "out_sample_rate",    c->sample_rate,    0);
        av_opt_set_sample_fmt(ost->swr_ctx, "out_sample_fmt",     c->sample_fmt,     0);

        /* initialize the resampling context */
        if ((ret = swr_init(ost->swr_ctx)) < 0) {
            fprintf(stderr, "Failed to initialize the resampling context\n");
            exit(1);
        }
}

/* Prepare a 16 bit dummy audio frame of 'frame_size' samples and
 * 'nb_channels' channels. */
static AVFrame *get_audio_frame(OutputStream *ost)
{
    AVFrame *frame = ost->tmp_frame;
    int j, i, v;
    int16_t *q = (int16_t*)frame->data[0];

    /* check if we want to generate more frames */
    if (av_compare_ts(ost->next_pts, ost->st->codec->time_base,
                      global_STREAM_DURATION, (AVRational){ 1, 1 }) >= 0)
        return NULL;

    //spi, begin
    /*
    for (j = 0; j <frame->nb_samples; j++) {
        v = (int)(sin(ost->t) * 10000);
        for (i = 0; i < ost->st->codec->channels; i++)
            *q++ = v;
        ost->t     += ost->tincr;
        ost->tincr += ost->tincr2;
    }
    */
    /*
    int count;
    count = sf_readf_short(global_pSNDFILE, global_intbuffer, frame->nb_samples);
    assert(count==(frame->nb_samples));
    global_NUM_CHANNELS = ost->st->codec->channels;
    for (j = 0; j <frame->nb_samples; j++)
    {
        // *q++ = global_intbuffer[2*j];
        // *q++ = global_intbuffer[2*j+1];
        for (i = 0; i < global_NUM_CHANNELS; i++)
        {
            *q++ = global_intbuffer[global_NUM_CHANNELS*j+i];
            global_floatbuffer[global_NUM_CHANNELS*j+i]=global_intbuffer[global_NUM_CHANNELS*j+i]/32768.0f;
        }

    }
    */
    int count;
    count = sf_readf_float(global_pSNDFILE, global_floatbuffer, frame->nb_samples);
    assert(count==(frame->nb_samples));
    global_NUM_CHANNELS = ost->st->codec->channels;
    for (j = 0; j <frame->nb_samples; j++)
    {
        // *q++ = global_intbuffer[2*j];
        // *q++ = global_intbuffer[2*j+1];
        for (i = 0; i < global_NUM_CHANNELS; i++)
        {
            global_intbuffer[global_NUM_CHANNELS*j+i] = global_floatbuffer[global_NUM_CHANNELS*j+i]*32767;
            *q++ = global_intbuffer[global_NUM_CHANNELS*j+i];
        }

    }
    //spi, end

    frame->pts = ost->next_pts;
    ost->next_pts  += frame->nb_samples;

    return frame;
}

/*
 * encode one audio frame and send it to the muxer
 * return 1 when encoding is finished, 0 otherwise
 */
static int write_audio_frame(AVFormatContext *oc, OutputStream *ost)
{
    AVCodecContext *c;
    AVPacket pkt = { 0 }; // data and size must be 0;
    AVFrame *frame;
    int ret;
    int got_packet;
    int dst_nb_samples;

    av_init_packet(&pkt);
    c = ost->st->codec;

    frame = get_audio_frame(ost);

    if (frame) {
        /* convert samples from native format to destination codec format, using the resampler */
            /* compute destination number of samples */
            dst_nb_samples = av_rescale_rnd(swr_get_delay(ost->swr_ctx, c->sample_rate) + frame->nb_samples,
                                            c->sample_rate, c->sample_rate, AV_ROUND_UP);
            av_assert0(dst_nb_samples == frame->nb_samples);

        /* when we pass a frame to the encoder, it may keep a reference to it
         * internally;
         * make sure we do not overwrite it here
         */
        ret = av_frame_make_writable(ost->frame);
        if (ret < 0)
            exit(1);

            /* convert to destination format */
            ret = swr_convert(ost->swr_ctx,
                              ost->frame->data, dst_nb_samples,
                              (const uint8_t **)frame->data, frame->nb_samples);
            if (ret < 0) {
                fprintf(stderr, "Error while converting\n");
                exit(1);
            }
            frame = ost->frame;

        frame->pts = av_rescale_q(ost->samples_count, (AVRational){1, c->sample_rate}, c->time_base);
        ost->samples_count += dst_nb_samples;
    }

    ret = avcodec_encode_audio2(c, &pkt, frame, &got_packet);
    if (ret < 0) {
        fprintf(stderr, "Error encoding audio frame: %s\n", av_err2str(ret));
        exit(1);
    }

    if (got_packet) {
        ret = write_frame(oc, &c->time_base, ost->st, &pkt);
        if (ret < 0) {
            fprintf(stderr, "Error while writing audio frame: %s\n",
                    av_err2str(ret));
            exit(1);
        }
    }

    return (frame || got_packet) ? 0 : 1;
}

/**************************************************************/
/* video output */

static AVFrame *alloc_picture(enum AVPixelFormat pix_fmt, int width, int height)
{
    AVFrame *picture;
    int ret;

    picture = av_frame_alloc();
    if (!picture)
        return NULL;

    picture->format = pix_fmt;
    picture->width  = width;
    picture->height = height;

    /* allocate the buffers for the frame data */
    ret = av_frame_get_buffer(picture, 32);
    if (ret < 0) {
        fprintf(stderr, "Could not allocate frame data.\n");
        exit(1);
    }

    return picture;
}

static void open_video(AVFormatContext *oc, AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg)
{
    int ret;
    AVCodecContext *c = ost->st->codec;
    AVDictionary *opt = NULL;

    av_dict_copy(&opt, opt_arg, 0);

    /* open the codec */
    ret = avcodec_open2(c, codec, &opt);
    av_dict_free(&opt);
    if (ret < 0) {
        fprintf(stderr, "Could not open video codec: %s\n", av_err2str(ret));
        exit(1);
    }

    /* allocate and init a re-usable frame */
    ost->frame = alloc_picture(c->pix_fmt, c->width, c->height);
    if (!ost->frame) {
        fprintf(stderr, "Could not allocate video frame\n");
        exit(1);
    }

    /* If the output format is not YUV420P, then a temporary YUV420P
     * picture is needed too. It is then converted to the required
     * output format. */
    ost->tmp_frame = NULL;
    if (c->pix_fmt != AV_PIX_FMT_YUV420P) {
        ost->tmp_frame = alloc_picture(AV_PIX_FMT_YUV420P, c->width, c->height);
        if (!ost->tmp_frame) {
            fprintf(stderr, "Could not allocate temporary picture\n");
            exit(1);
        }
    }
}

// Prepare a dummy image.
static void fill_yuv_image(AVFrame *pict, int frame_index,
                           int SPECWIDTH, int SPECHEIGHT)
                           //int width, int height)
{
    /*
    int x, y, ii, ret;

    // when we pass a frame to the encoder, it may keep a reference to it
    // internally;
    // make sure we do not overwrite it here
    //
    ret = av_frame_make_writable(pict);
    if (ret < 0)
        exit(1);

    ii = frame_index;

    // Y
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            pict->data[0][y * pict->linesize[0] + x] = x + y + ii * 3; //original
            //pict->data[0][y * pict->linesize[0] + x] = x + y + ii ; //spi
            //pict->data[0][y * pict->linesize[0] + x] = 128 ; //spi

    // Cb and Cr
    for (y = 0; y < height / 2; y++) //original
    //for (y = 0; y < height; y++) //spi
    {
        for (x = 0; x < width / 2; x++) //original
        //for (x = 0; x < width ; x++) //spi
        {
            pict->data[1][y * pict->linesize[1] + x] = 128 + y + ii * 2; //original
            pict->data[2][y * pict->linesize[2] + x] = 64 + x + ii * 5; //original
            //pict->data[1][y * pict->linesize[1] + x] = y + ii; //spi
            //pict->data[2][y * pict->linesize[2] + x] = x + ii; //spi
            //pict->data[1][y * pict->linesize[1] + x] = 128; //spi
            //pict->data[2][y * pict->linesize[2] + x] = 128; //spi
        }
    }
    */
    int i, j, ret;
    int x,y,y1;
    // when we pass a frame to the encoder, it may keep a reference to it
    // internally;
    // make sure we do not overwrite it here
    //
    ret = av_frame_make_writable(pict);
    if (ret < 0)
        exit(1);

    if (specmode==3 || specmode==4 || specmode==5 || specmode==6 ||
        specmode==7 || specmode==8 || specmode==9 || specmode==10 ||
        specmode==11 || specmode==12 || specmode==13 || specmode==14 ||
        specmode==15 || specmode==16 || specmode==17 || specmode==18)
    { // waveform and filled waveform

        int c, x, y;
        //float *buf;

        if(specmode==3 || specmode==4 || specmode==5 || specmode==6)
        {
            //solid background
            //int Y = RGBtoYCbCr_Y(0, 0, 0); //black
            //int Y = RGBtoYCbCr_Y(255, 0, 0); //red
            int Y = RGBtoYCbCr_Y(0, 0, 255); //blue
            //Y
            for(i=0; i<SPECWIDTH; i++)
            {
                for (j=0; j<SPECHEIGHT; j++)
                {
                    pict->data[0][j*pict->linesize[0]+i]=Y;
                }
            }
            //Cb and Cr
            //int Cb = RGBtoYCbCr_Cb(0, 0, 0); //black
            //int Cr = RGBtoYCbCr_Cr(0, 0, 0); //black
            //int Cb = RGBtoYCbCr_Cb(255, 0, 0); //red
            //int Cr = RGBtoYCbCr_Cr(255, 0, 0); //red
            int Cb = RGBtoYCbCr_Cb(0, 0, 255); //blue
            int Cr = RGBtoYCbCr_Cr(0, 0, 255); //blue
            for(i=0; i<SPECWIDTH/2; i++)
            {
                for (j=0; j<SPECHEIGHT/2; j++)
                {
                    pict->data[1][j*pict->linesize[1]+i]=Cb;
                    pict->data[2][j*pict->linesize[2]+i]=Cr;
                }
            }
        }
        //noisy background
        else if(specmode==7 || specmode==8 || specmode==9 || specmode==10)
        {
            //Y
            for(i=0; i<SPECWIDTH; i++)
            {
                for (j=0; j<SPECHEIGHT; j++)
                {
                    int random_integer;
                    //int lowest=1, highest=127;
                    //int lowest=64, highest=96; //good
                    int lowest=1, highest=255; //good
                    int range=(highest-lowest)+1;
                    random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[0][j*pict->linesize[0]+i]=RGBtoYCbCr_Y(random_integer, random_integer, random_integer);
                }
            }
            //Cb and Cr
            for(i=0; i<SPECWIDTH/2; i++)
            {
                for (j=0; j<SPECHEIGHT/2; j++)
                {
                    int random_integer;
                    //int lowest=1, highest=127;
                    //int lowest=64, highest=96; //good
                    int lowest=1, highest=255; //good
                    int range=(highest-lowest)+1;
                    random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[1][j*pict->linesize[1]+i]=RGBtoYCbCr_Cb(random_integer, random_integer, random_integer);
                    random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));
                    pict->data[2][j*pict->linesize[2]+i]=RGBtoYCbCr_Cr(random_integer, random_integer, random_integer);
                }
            }
        }
        //solid slightly shifting background
        else if(specmode==11 || specmode==12 || specmode==13 || specmode==14)
        {
            int random_integer;
            //int lowest=1, highest=127;
            int lowest=64, highest=96; //good
            //int lowest=1, highest=255; //good
            int range=(highest-lowest)+1;
            random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));

            int Y = RGBtoYCbCr_Y(random_integer, random_integer, random_integer);
            int Cb = RGBtoYCbCr_Cb(random_integer, random_integer, random_integer);
            int Cr = RGBtoYCbCr_Cr(random_integer, random_integer, random_integer);
            //Y
            for(i=0; i<SPECWIDTH; i++)
            {
                for (j=0; j<SPECHEIGHT; j++)
                {
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[0][j*pict->linesize[0]+i]=Y;
                }
            }
            //Cb and Cr
            for(i=0; i<SPECWIDTH/2; i++)
            {
                for (j=0; j<SPECHEIGHT/2; j++)
                {
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[1][j*pict->linesize[1]+i]=Cb;
                    random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));
                    pict->data[2][j*pict->linesize[2]+i]=Cr;
                }
            }

        }
        //solid radically shifting background
        else //specmode==15 || specmode==16 || specmode==17 || specmode==18
        {
            int random_integer;
            //int lowest=1, highest=127;
            //int lowest=64, highest=96; //good
            int lowest=1, highest=255; //good
            int range=(highest-lowest)+1;
            random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));

            int Y = RGBtoYCbCr_Y(random_integer, random_integer, random_integer);
            int Cb = RGBtoYCbCr_Cb(random_integer, random_integer, random_integer);
            int Cr = RGBtoYCbCr_Cr(random_integer, random_integer, random_integer);
            //Y
            for(i=0; i<SPECWIDTH; i++)
            {
                for (j=0; j<SPECHEIGHT; j++)
                {
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[0][j*pict->linesize[0]+i]=Y;
                }
            }
            //Cb and Cr
            for(i=0; i<SPECWIDTH/2; i++)
            {
                for (j=0; j<SPECHEIGHT/2; j++)
                {
                    //specbuf[j*SPECWIDTH+i]=random_integer;
                    pict->data[1][j*pict->linesize[1]+i]=Cb;
                    random_integer = lowest+(int)(range*rand()/(RAND_MAX + 1.0));
                    pict->data[2][j*pict->linesize[2]+i]=Cr;
                }
            }
        }

        /*
        buf=(float*)alloca(NUM_CHANNELS*SPECWIDTH*sizeof(float)); // allocate buffer for data
        global_err = Pa_ReadStream( global_stream, buf, NUM_CHANNELS*SPECWIDTH );
        if( global_err != paNoError )
        {
            //char errorbuf[2048];
            //sprintf(errorbuf, "Error reading stream: %s\n", Pa_GetErrorText(global_err));
            //MessageBox(0,errorbuf,0,MB_ICONERROR);
            return;
        }
        */

        //under waveform filled down to bottom
        if (specmode==6 || specmode==10 || specmode==14 || specmode==18)
        {
            for (c=0;c<global_NUM_CHANNELS;c++)
            {
                for (x=0;x<SPECWIDTH;x++)
                {
                    int v=(1-global_floatbuffer[x*global_NUM_CHANNELS+c])*SPECHEIGHT/2; // invert and scale to fit display
                    if (v<0) v=0;
                    else if (v>=SPECHEIGHT) v=SPECHEIGHT-1;
                    //under waveform filled down to bottom
                    y=v;
                    while(--y>=0) pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                }
            }
        }
        //waveform filled towards center
        else if (specmode==5 || specmode==9 || specmode==13 || specmode==17)
        {
            for (c=0;c<global_NUM_CHANNELS;c++)
            {
                for (x=0;x<SPECWIDTH;x++)
                {
                    int v=(1-global_floatbuffer[x*global_NUM_CHANNELS+c])*SPECHEIGHT/2; // invert and scale to fit display
                    if (v<0) v=0;
                    else if (v>=SPECHEIGHT) v=SPECHEIGHT-1;
                    //waveform filled towards center
                    y=v;
                    if(y>(SPECHEIGHT/2))
                        while(--y>=(SPECHEIGHT/2)) pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                    else if(y<(SPECHEIGHT/2))
                        while(++y<=(SPECHEIGHT/2)) pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                    else pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                }
            }
        }
        //waveform filled towards opposite
        else if (specmode==4 || specmode==8 || specmode==12 || specmode==16)
        {
            for (c=0;c<global_NUM_CHANNELS;c++)
            {
                for (x=0;x<SPECWIDTH;x++)
                {
                    int v=(1-global_floatbuffer[x*global_NUM_CHANNELS+c])*SPECHEIGHT/2; // invert and scale to fit display
                    if (v<0) v=0;
                    else if (v>=SPECHEIGHT) v=SPECHEIGHT-1;
                    //waveform filled towards opposite
                    y=v;
                    if(y>(SPECHEIGHT/2))
                        while(--y>=(SPECHEIGHT/2-(v-(SPECHEIGHT/2)))) pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                    else if(y<(SPECHEIGHT/2))
                        while(++y<=(SPECHEIGHT/2+((SPECHEIGHT/2)-v))) pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                    else pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                }
            }
        }
        //waveform (original)
        else if(specmode==3 || specmode==7 || specmode==11 || specmode==15)
        {
            for (c=0;c<global_NUM_CHANNELS;c++)
            {
                for (x=0;x<SPECWIDTH;x++)
                {
                    int v=(1-global_floatbuffer[x*global_NUM_CHANNELS+c])*SPECHEIGHT/2; // invert and scale to fit display
                    if (v<0) v=0;
                    else if (v>=SPECHEIGHT) v=SPECHEIGHT-1;
                    if (!x) y=v;
                    do
                    { // draw line from previous sample...
                        if (y<v) y++;
                        else if (y>v) y--;
                        pict->data[0][y*SPECWIDTH+x]=c&1?127:1;
                    } while (y!=v);
                }
            }
        }

    }
    else
    {
        //specmode = 0, 1 or 2

        /*
        assert(0);
        printf("modes 0, 1 and 2 not defined!\n");
        exit(1);
        */

        //specmode = 0, 1 or 2

        assert(SPECWIDTH==1024);
        float fftbuf[1024];

        //calling fft() that is defined in rfftw\\fourier.cpp does not work from here, may be the calling convention is a problem
        //calling fft() that is defined in rfftw-c\\fourier.c does not work from here, it returns all zeros
        //maybe rfftw-c will have to be compiled with qt creator mingw compiler
        //calling fft() that is defined in rfftw-c_qt2016\\fourier.c does not work from here, it returns all zeros
        fft(global_floatbuffer, fftbuf, 1024);

        //FILE* pFILE = fopen("debug-dat.txt", "w");
        //FILE* pFILE2 = fopen("debug-fft.txt", "w");
        for(i=0; i<1024; i++)
        {
            fftbuf[i]=fabs(fftbuf[i]);
            //if(pFILE2) fprintf(pFILE2, "%f\n", fftbuf[i]);
        }
        //if(pFILE) fclose(pFILE);
        //if(pFILE2) fclose(pFILE2);

        if (!specmode) //when specmode==0
        { // "normal" FFT


            //solid background
            //int Y = RGBtoYCbCr_Y(0, 0, 0); //black
            int Y = RGBtoYCbCr_Y(255, 0, 0); //red
            //int Y = RGBtoYCbCr_Y(0, 0, 255); //blue
            //Y
            for(i=0; i<SPECWIDTH; i++)
            {
                for (j=0; j<SPECHEIGHT; j++)
                {
                    pict->data[0][j*pict->linesize[0]+i]=Y;
                }
            }
            //Cb and Cr
            //int Cb = RGBtoYCbCr_Cb(0, 0, 0); //black
            //int Cr = RGBtoYCbCr_Cr(0, 0, 0); //black
            int Cb = RGBtoYCbCr_Cb(255, 0, 0); //red
            int Cr = RGBtoYCbCr_Cr(255, 0, 0); //red
            //int Cb = RGBtoYCbCr_Cb(0, 0, 255); //blue
            //int Cr = RGBtoYCbCr_Cr(0, 0, 255); //blue
            for(i=0; i<SPECWIDTH/2; i++)
            {
                for (j=0; j<SPECHEIGHT/2; j++)
                {
                    pict->data[1][j*pict->linesize[1]+i]=Cb;
                    pict->data[2][j*pict->linesize[2]+i]=Cr;
                }
            }

            //draw spectrum
            Y = RGBtoYCbCr_Y(0, 0, 0); //black
            for (x=0;x<SPECWIDTH/2;x++)
            {
#if 1
                y=sqrt(fftbuf[x+1])*3*SPECHEIGHT-4; // scale it (sqrt to make low values more visible)
#else
                y=fftbuf[x+1]*10*SPECHEIGHT; // scale it (linearly)
#endif
                if (y>SPECHEIGHT) y=SPECHEIGHT; // cap it
                if (x && (y1=(y+y1)/2)) // interpolate from previous to make the display smoother
                    //while (--y1>=0) specbuf[y1*SPECWIDTH+x*2-1]=y1+1;
                    //while (--y1>=0) pict->data[0][y1*SPECWIDTH+x*2-1]=(127*y1/SPECHEIGHT)+1;
                    //while (--y1>=0) pict->data[0][y1*SPECWIDTH+x*2-1]=Y;
                    while (--y1>=0) pict->data[0][(SPECHEIGHT-y1)*SPECWIDTH+x*2-1]=Y;
                y1=y;
                //while (--y>=0) specbuf[y*SPECWIDTH+x*2]=y+1; // draw level
                //while (--y>=0) pict->data[0][y*SPECWIDTH+x*2]=(127*y/SPECHEIGHT)+1; // draw level
                //while (--y>=0) pict->data[0][y*SPECWIDTH+x*2]=Y; // draw level
                while (--y>=0) pict->data[0][(SPECHEIGHT-y)*SPECWIDTH+x*2]=Y; // draw level
            }

        }
        else if (specmode==1)
        {
            //todo:
            assert(0);
        }
        else if (specmode==2)
        {
            //todo:
            assert(0);
        }
        else
        {
            //no other mode
            assert(0);
        }
    }
}

static AVFrame *get_video_frame(OutputStream *ost)
{
    AVCodecContext *c = ost->st->codec;

    /* check if we want to generate more frames */
    if (av_compare_ts(ost->next_pts, ost->st->codec->time_base,
                      global_STREAM_DURATION, (AVRational){ 1, 1 }) >= 0)
        return NULL;

    if (c->pix_fmt != AV_PIX_FMT_YUV420P) {
        /* as we only generate a YUV420P picture, we must convert it
         * to the codec pixel format if needed */
        if (!ost->sws_ctx) {
            ost->sws_ctx = sws_getContext(c->width, c->height,
                                          AV_PIX_FMT_YUV420P,
                                          c->width, c->height,
                                          c->pix_fmt,
                                          SCALE_FLAGS, NULL, NULL, NULL);
            if (!ost->sws_ctx) {
                fprintf(stderr,
                        "Could not initialize the conversion context\n");
                exit(1);
            }
        }
        fill_yuv_image(ost->tmp_frame, ost->next_pts, c->width, c->height);
        sws_scale(ost->sws_ctx,
                  (const uint8_t * const *)ost->tmp_frame->data, ost->tmp_frame->linesize,
                  0, c->height, ost->frame->data, ost->frame->linesize);
    } else {
        fill_yuv_image(ost->frame, ost->next_pts, c->width, c->height);
    }

    ost->frame->pts = ost->next_pts++;

    return ost->frame;
}

/*
 * encode one video frame and send it to the muxer
 * return 1 when encoding is finished, 0 otherwise
 */
static int write_video_frame(AVFormatContext *oc, OutputStream *ost)
{
    int ret;
    AVCodecContext *c;
    AVFrame *frame;
    int got_packet = 0;
    AVPacket pkt = { 0 };

    c = ost->st->codec;

    frame = get_video_frame(ost);

    av_init_packet(&pkt);

    /* encode the image */
    ret = avcodec_encode_video2(c, &pkt, frame, &got_packet);
    if (ret < 0) {
        fprintf(stderr, "Error encoding video frame: %s\n", av_err2str(ret));
        exit(1);
    }

    if (got_packet) {
        ret = write_frame(oc, &c->time_base, ost->st, &pkt);
    } else {
        ret = 0;
    }

    if (ret < 0) {
        fprintf(stderr, "Error while writing video frame: %s\n", av_err2str(ret));
        exit(1);
    }

    return (frame || got_packet) ? 0 : 1;
}

static void close_stream(AVFormatContext *oc, OutputStream *ost)
{
    avcodec_close(ost->st->codec);
    av_frame_free(&ost->frame);
    av_frame_free(&ost->tmp_frame);
    sws_freeContext(ost->sws_ctx);
    swr_free(&ost->swr_ctx);
}

/**************************************************************/
/* media file output */

int main(int argc, char **argv)
{
    OutputStream video_st = { 0 }, audio_st = { 0 };
    const char *filename;
    AVOutputFormat *fmt;
    AVFormatContext *oc;
    AVCodec *audio_codec, *video_codec;
    int ret;
    int have_video = 0, have_audio = 0;
    int encode_video = 0, encode_audio = 0;
    AVDictionary *opt = NULL;

    /* Initialize libavcodec, and register all codecs and formats. */
    av_register_all();

    if (argc < 2) {
        printf("usage: %s output_file\n"
               "API example program to output a media file with libavformat.\n"
               "This program generates a synthetic audio and video stream, encodes and\n"
               "muxes them into a file named output_file.\n"
               "The output format is automatically guessed according to the file extension.\n"
               "Raw images can also be output by using '%%d' in the filename.\n"
               "\n", argv[0]);
        return 1;
    }

    //spi, begin
    specmode = 4;
    //specmode = 0; //does not work for now

    //////////////////////////
    //initialize random number
    //////////////////////////
    srand((unsigned)time(0));
    //////////////
    //open sndfile
    //////////////
    global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\bookaudio_wav\\Nietzsche-Biography and Analysis\\Friedrich Nietzsche - Audio Book Friedrich Nietzsche II(time-stretched-120pc)\\nietzsche-is-a-hard-critic-of-both-philosophy-and-religion.wav", SFM_READ, &global_mySF_INFO);
    //global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\techno_wav\\techno5(dj-oifii-mix_underworld-rez_cymbals).wav", SFM_READ, &global_mySF_INFO);
    //global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\techno_wav\\techno6(dj-oifii-mix_underworld-rowla_cymbals).wav", SFM_READ, &global_mySF_INFO);
    //global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\trance_wav\\dj-oifii_trance-euphoric-uplifting_2015feb06.wav", SFM_READ, &global_mySF_INFO);
    //global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\worldaudio_wav\\00min30sec-and-less\\GF - Subramanian - track 03(intro)_18sec.wav", SFM_READ, &global_mySF_INFO);
    //global_pSNDFILE = sf_open("D:\\oifii-org\\httpdocs\\ha-org\\had\\dj-oifii\\worldaudio_wav\\00min15sec-and-less\\Geoffrey Oryema - TAO -  mara(introlater)_9sec.wav", SFM_READ, &global_mySF_INFO);
    assert(global_mySF_INFO.samplerate == 44100);
    assert(global_mySF_INFO.channels == 2);
    float file1duration_s;
    file1duration_s = ((float)global_mySF_INFO.frames) / ((float)global_mySF_INFO.samplerate);
    global_STREAM_DURATION = file1duration_s; //will round flooring the value, so it should be fine
    //spi, end


    filename = argv[1];
    if (argc > 3 && !strcmp(argv[2], "-flags")) {
        av_dict_set(&opt, argv[2]+1, argv[3], 0);
    }

    /* allocate the output media context */
    avformat_alloc_output_context2(&oc, NULL, NULL, filename);
    if (!oc) {
        printf("Could not deduce output format from file extension: using MPEG.\n");
        avformat_alloc_output_context2(&oc, NULL, "mpeg", filename);
    }
    if (!oc)
        return 1;

    fmt = oc->oformat;

    /* Add the audio and video streams using the default format codecs
     * and initialize the codecs. */
    if (fmt->video_codec != AV_CODEC_ID_NONE) {
        add_stream(&video_st, oc, &video_codec, fmt->video_codec);
        have_video = 1;
        encode_video = 1;
    }
    if (fmt->audio_codec != AV_CODEC_ID_NONE) {
        add_stream(&audio_st, oc, &audio_codec, fmt->audio_codec);
        have_audio = 1;
        encode_audio = 1;
    }

    /* Now that all the parameters are set, we can open the audio and
     * video codecs and allocate the necessary encode buffers. */
    if (have_video)
        open_video(oc, video_codec, &video_st, opt);

    if (have_audio)
        open_audio(oc, audio_codec, &audio_st, opt);

    av_dump_format(oc, 0, filename, 1);

    /* open the output file, if needed */
    if (!(fmt->flags & AVFMT_NOFILE)) {
        ret = avio_open(&oc->pb, filename, AVIO_FLAG_WRITE);
        if (ret < 0) {
            fprintf(stderr, "Could not open '%s': %s\n", filename,
                    av_err2str(ret));
            return 1;
        }
    }

    /* Write the stream header, if any. */
    ret = avformat_write_header(oc, &opt);
    if (ret < 0) {
        fprintf(stderr, "Error occurred when opening output file: %s\n",
                av_err2str(ret));
        return 1;
    }

    while (encode_video || encode_audio) {
        /* select the stream to encode */
        if (encode_video &&
            (!encode_audio || av_compare_ts(video_st.next_pts, video_st.st->codec->time_base,
                                            audio_st.next_pts, audio_st.st->codec->time_base) <= 0)) {
            encode_video = !write_video_frame(oc, &video_st);
        } else {
            encode_audio = !write_audio_frame(oc, &audio_st);
        }
    }

    /* Write the trailer, if any. The trailer must be written before you
     * close the CodecContexts open when you wrote the header; otherwise
     * av_write_trailer() may try to use memory that was freed on
     * av_codec_close(). */
    av_write_trailer(oc);

    /* Close each codec. */
    if (have_video)
        close_stream(oc, &video_st);
    if (have_audio)
        close_stream(oc, &audio_st);

    if (!(fmt->flags & AVFMT_NOFILE))
        /* Close the output file. */
        avio_closep(&oc->pb);

    /* free the stream */
    avformat_free_context(oc);

    //spi, begin
    if(global_pSNDFILE) sf_close(global_pSNDFILE);
    //spi, end
    return 0;
}
