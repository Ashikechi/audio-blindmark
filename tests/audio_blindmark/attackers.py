import math
from pathlib import Path

import numpy as np
from pydub import AudioSegment
from audio_blindmark.utils.random import get_rng

from .utils import read_wave, write_wave


def white_noise(audio: str | Path, audio_with_noise: str | Path, factor: float) -> None:
    channels, sampwidth, framerate = read_wave(audio)
    channels_with_noise = [channel + np.random.Generator(get_rng()).uniform(-1, 1, channel.shape) * np.sqrt(np.dot(channel, channel) / channel.shape[0]) * factor for channel in channels]
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
