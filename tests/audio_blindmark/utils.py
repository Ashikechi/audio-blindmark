import wave as w
from pathlib import Path

import numpy as np
from aquatk.metrics.PEAQ import PEAQ

from .audios import output_audio_path, raw_audio_path, report_path


def read_wave(path: str | Path) -> tuple[list[np.ndarray[tuple[int], np.dtype[np.float64]]], int, int]:
    with open(path, 'rb') as f:
        with w.open(f, 'rb') as wave:
            channel_num = wave.getnchannels()
            width = wave.getsampwidth()
            data = wave.readframes(-1)
            return [*np.array(np.frombuffer(data, f'<i{width}').astype(np.float64).reshape([-1, channel_num]).T, order='C')], wave.getsampwidth(), wave.getframerate()

def write_wave(path: str | Path, channels: list[np.ndarray[tuple[int], np.dtype[np.float64]]], width: int = 2, framerate: int = 48000) -> None:
    assert len(channels) > 0
    frame_num = channels[0].shape[0]
    for each in channels:
        assert each.shape[0] == frame_num

    with open(path, 'wb') as f:
        with w.open(f, 'wb') as wave:
            wave.setnchannels(len(channels))
            wave.setsampwidth(width)
            wave.setframerate(framerate)
            wave.setnframes(frame_num)
            wave.setcomptype('NONE', 'not compressed')
            wave.writeframes(np.stack(channels, axis=-1).reshape(-1).astype(f'<i{width}').tobytes())

def generate_PEAQ_report(audio_id: str, embedder: str) -> None:  # pylint: disable=C0103
    peaq = PEAQ(version="advanced")
    result = peaq.analyze_files(raw_audio_path(audio_id), output_audio_path(audio_id, embedder))
    with open(report_path(audio_id, embedder), 'w', encoding='utf-8') as f:
        f.write(str(result))
