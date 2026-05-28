import math
from pathlib import Path

import numpy as np
from pydub import AudioSegment
from scipy import signal

from audio_blindmark.utils.random import get_rng

from .utils import read_wave, write_wave


def white_noise(audio: str | Path, audio_with_noise: str | Path, factor: float) -> None:
    channels, sampwidth, framerate = read_wave(audio)
    channels_with_noise = []
    for channel in channels:
        white_noise = np.random.Generator(get_rng()).normal(0, 1, channel.shape)
        channels_with_noise.append(channel + white_noise * np.sqrt(np.dot(channel, channel) / channel.shape[0]) * factor)
    write_wave(audio_with_noise, channels_with_noise, sampwidth, framerate)

def pink_noise(audio: str | Path, audio_with_noise: str | Path, factor: float) -> None:
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
    A = [1, -2.494956002, 2.017265875, -0.522189400]

    channels, sampwidth, framerate = read_wave(audio)
    channels_with_noise = []
    for channel in channels:
        white_noise = np.random.Generator(get_rng()).normal(0, 1, channel.shape)
        pink_noise = signal.lfilter(B, A, white_noise)
        pink_noise /= np.sqrt(np.dot(pink_noise, pink_noise) / pink_noise.shape[0])
        channels_with_noise.append(channel + pink_noise * np.sqrt(np.dot(channel, channel) / channel.shape[0]) * factor)
    write_wave(audio_with_noise, channels_with_noise, sampwidth, framerate)

def lossy_compression(audio: str | Path, compressed_audio: str | Path, form: str) -> None:
    wave = AudioSegment.from_wav(audio)
    wave.export(compressed_audio, format=form)
    wave = AudioSegment.from_file(compressed_audio, format=form)
    wave.export(compressed_audio, format='wav')

def resample(audio: str | Path, resampled_audio: str | Path, frame_rate: int) -> None:
    wave = AudioSegment.from_wav(audio)
    original_frame_rate = wave.frame_rate
    resampled_wave = wave.set_frame_rate(frame_rate).set_frame_rate(original_frame_rate)
    resampled_wave.export(resampled_audio, format='wav')

def zoom(audio: str | Path, zoomed_audio: str | Path, factor: float) -> None:
    wave = AudioSegment.from_wav(audio)
    zoomed_wave = wave.apply_gain(math.log10(factor) * 20)
    zoomed_wave.export(zoomed_audio, format='wav')
