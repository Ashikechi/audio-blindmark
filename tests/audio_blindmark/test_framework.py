from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.utils.random import seed

from .audios import SHORT_AUDIOS, output_audio_path, raw_audio_path
from .mock_steganographier import MockEmbedder, MockExtractor
from .utils import read_wave, write_wave

KEY = b'Kei'
DATA = b'Aris -- Princesses of Unnamed Gods' * 64

def test_framework():
    ECC_LENGTH = 2  # pylint: disable=C0103

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = MockEmbedder()
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'mock'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        extractor = MockExtractor()
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(output_audio_path(audio, 'mock'))[0], decoder, extractor) == DATA
