# FFmpeg_Loudness

Audio Loudness Normalization Filter Port From FFmpeg

The most well known audio loudness normalization method currently is ebu R128.
This standard was defined by the European brodcasting union the sepcification of the standard can be found [here](https://tech.ebu.ch/docs/r/r128.pdf).

FFmpeg implements the [loudnorm](http://k.ylo.ph/2016/04/04/loudnorm.html) filter that supports ebu R128 standard.

The aim of this standard is to make audio equally loud everywhere.

## Parameters

When using FFmpeg to normalize audio with the loudnorm filter you can use the following parameters.

I : The target loudness in LUFs ( LUF stands for Loudness Units Full Scale and one unit of LUF equals one dB ).

LRA: This stands for Loudness Range. This describes the overall loduness range, from the softest part to the loudest part.

TP: This stands for true peak.

## Overview
When using FFmpeg to normalize audio with ebu R128 you have two options.
Single is ideal for live normalization, produces worse results but is faster.
With a single pass you can set parameters for the loudness and apply them directly.

Double pass is ideal for postprocessing, produces better results but slower.
With double pass first you scan the media file you want to normalize and then apply target loudness parameters with the measured values.

## Links:
[loudness explained](https://www.tcelectronic.com/brand/tcelectronic/loudness-explained)

# Donating

If you found this project useful, consider buying me a coffee

<a href="https://www.buymeacoffee.com/gaozhihan" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/black_img.png" alt="Buy Me A Coffee" style="height: auto !important;width: auto !important;" ></a>
  
 
