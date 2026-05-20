from pathlib import Path

from pydub import AudioSegment

from .utils import read_wave, write_wave


def add_compression_noise(audio: str | Path, compressed_audio: str | Path, form: str) -> None:
    wave = AudioSegment.from_wav(audio)
    wave.export(compressed_audio, format=form)
    wave = AudioSegment.from_file(compressed_audio, format=form)
    wave.export(compressed_audio, format='wav')

def zoom(audio: str | Path, zoomed_audio: str | Path, factor: float) -> None:
    wave, width, framerate = read_wave(audio)
    wave = [each * factor for each in wave]
    write_wave(zoomed_audio, wave, width, framerate)
