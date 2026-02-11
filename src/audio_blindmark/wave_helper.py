import wave

import numpy as np


class WaveReadHelper:
    def __init__(self, f: wave.Wave_read) -> None:
        self.f = f
        self.channel_num = f.getnchannels()
        self.width = f.getsampwidth()

    def read(self, frame_num: int) -> tuple[list[np.ndarray], int]:
        data = self.f.readframes(frame_num)
        frame_num = len(data) // (self.channel_num * self.width)
        channels = [*np.frombuffer(data, f'<i{self.width}').astype(np.float64).reshape([-1, self.channel_num]).T]
        return channels, frame_num

class WaveWriteHelper:
    def __init__(self, f: wave.Wave_write) -> None:
        self.f = f
        self.channel_num = f.getnchannels()
        self.width = f.getsampwidth()

    def write(self, channels: list[np.ndarray]) -> None:
        assert(len(channels) == self.channel_num)
        for each in channels:
            assert(each.ndim == 1)
        frame_num = channels[0].shape[0]
        for each in channels:
            assert(each.shape[0] == frame_num)

        int_stack = np.stack(channels, axis = -1).reshape(-1).astype(f'<i{self.width}')
        self.f.writeframes(int_stack.tobytes())
