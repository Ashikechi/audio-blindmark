import wave

import numpy as np
from numpy.random import Generator

from . import base
from .coder import Encoder
from .utils.random import get_rng
from .utils.random_tools import get_random_bytes  # pylint: disable = E0611, E0401
from .wave_helper import WaveReadHelper, WaveWriteHelper


class EmbedError(Exception):
    pass

def embed(raw_audio: wave.Wave_read, data: bytes, output_audio: wave.Wave_write, encoder: Encoder, embedder: base.BaseEmbedder) -> None:
    read_helper = WaveReadHelper(raw_audio)
    raw_channels, frame_num = read_helper.read(-1)
    output_channels = [np.empty(frame_num) for _ in raw_channels]

    full_data = len(data).to_bytes(4, 'little') + data

    wave_length = embedder.wave_length()
    data_length = encoder.data_length(embedder.data_length())

    packets = [full_data[i * data_length:(i + 1) * data_length] for i in range(len(full_data) // data_length)]
    if len(full_data) > len(packets) * data_length:
        packets.append(full_data[len(packets) * data_length:] + get_random_bytes(data_length - (len(full_data) - len(packets) * data_length)))

    channel_pointers = [0] * len(raw_channels)
    i = 0
    while i < len(packets):
        earliest_channel = 0
        for j, each in enumerate(channel_pointers):
            if each < channel_pointers[earliest_channel]:
                earliest_channel = j
        p = channel_pointers[earliest_channel]
        if frame_num - p < wave_length:
            raise EmbedError(f'Only {i * data_length - 4} bytes can be embedded')

        try:
            output_channels[earliest_channel][p:p + wave_length] = embedder.embed(
                raw_channels[earliest_channel][p:p + wave_length],
                encoder.encode(i, packets[i])
            )
            print(f'Embed {i} at {p} of {earliest_channel}')
            channel_pointers[earliest_channel] += wave_length
            i += 1
        except base.EmbedError:
            output_channels[earliest_channel][p] = raw_channels[earliest_channel][p]
            channel_pointers[earliest_channel] += 1

    while True:
        earliest_channel = 0
        for j, each in enumerate(channel_pointers):
            if each < channel_pointers[earliest_channel]:
                earliest_channel = j
        p = channel_pointers[earliest_channel]
        if frame_num - p < wave_length:
            break

        seq = int(Generator(get_rng()).integers(0, len(packets)))
        try:
            output_channels[earliest_channel][p:p + wave_length] = embedder.embed(
                raw_channels[earliest_channel][p:p + wave_length],
                encoder.encode(seq, packets[seq])
            )
            print(f'Embed {seq} at {p} of {earliest_channel}')
            channel_pointers[earliest_channel] += wave_length
        except base.EmbedError:
            output_channels[earliest_channel][p] = raw_channels[earliest_channel][p]
            channel_pointers[earliest_channel] += 1

    for i, each in enumerate(channel_pointers):
        output_channels[i][each:] = raw_channels[i][each:]

    write_helper = WaveWriteHelper(output_audio)
    write_helper.write(output_channels)
