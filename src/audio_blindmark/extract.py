import numpy as np

from . import base
from .coder import DecodeError, Decoder


class ExtractError(Exception):
    pass

def extract(channels: list[np.ndarray[tuple[int], np.dtype[np.float64]]], decoder: Decoder, extractor: base.BaseExtractor) -> bytes:
    assert len(channels) > 0
    frame_num = channels[0].shape[0]
    for each in channels:
        assert each.shape[0] == frame_num

    packets = {}

    wave_length = extractor.wave_length()
    data_length = decoder.data_length(extractor.data_length())

    for each in channels:
        p = 0
        while frame_num - p >= wave_length:
            try:
                seq, data = decoder.decode(extractor.extract(each[p:p + wave_length]))
                if seq not in packets:
                    packets[seq] = data
                p += wave_length
            except (base.ExtractError, DecodeError):
                p += 1

    if 0 not in packets:
        raise ExtractError('The header is lost')

    l = int.from_bytes(packets[0][:4], 'little')
    full_data = [packets[0][4:l + 4]]

    for i in range(1, (l + 4 - 1) // data_length + 1):
        if i not in packets:
            raise ExtractError(f'{i * data_length - 4}-{min((i + 1) * data_length - 4, l)} is lost')
        full_data.append(packets[i][:l - (i * data_length - 4)])

    return b''.join(full_data)
