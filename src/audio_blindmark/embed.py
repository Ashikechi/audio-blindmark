import numpy as np
from numpy.random import Generator

from . import base
from .coder import Encoder
from .utils.random import get_rng
from .utils.random_tools import get_random_bytes


class EmbedError(Exception):
    pass

def embed(raw_channels: list[np.ndarray[tuple[int], np.dtype[np.float64]]], data: bytes, encoder: Encoder, embedder: base.BaseEmbedder) -> list[np.ndarray[tuple[int], np.dtype[np.float64]]]:
    assert len(raw_channels) > 0
    frame_num = raw_channels[0].shape[0]
    for channel in raw_channels:
        assert channel.shape[0] == frame_num

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
        for j, pointer in enumerate(channel_pointers):
            if pointer < channel_pointers[earliest_channel]:
                earliest_channel = j
        p = channel_pointers[earliest_channel]
        if frame_num - p < wave_length:
            raise EmbedError(f'Only {i * data_length - 4} bytes can be embedded')

        try:
            output_channels[earliest_channel][p:p + wave_length] = embedder.embed(
                raw_channels[earliest_channel][p:p + wave_length],
                encoder.encode(i, packets[i])
            )
            channel_pointers[earliest_channel] += wave_length
            i += 1
        except base.EmbedError:
            output_channels[earliest_channel][p] = raw_channels[earliest_channel][p]
            channel_pointers[earliest_channel] += 1

    while True:
        earliest_channel = 0
        for j, pointer in enumerate(channel_pointers):
            if pointer < channel_pointers[earliest_channel]:
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
            channel_pointers[earliest_channel] += wave_length
        except base.EmbedError:
            output_channels[earliest_channel][p] = raw_channels[earliest_channel][p]
            channel_pointers[earliest_channel] += 1

    for i, pointer in enumerate(channel_pointers):
        output_channels[i][pointer:] = raw_channels[i][pointer:]

    return output_channels
