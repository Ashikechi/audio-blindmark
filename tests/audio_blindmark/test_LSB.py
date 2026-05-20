# pylint: disable=C0103
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.LSB import LSBEmbedder, LSBExtractor
from audio_blindmark.utils.random import seed

from .attackers import add_compression_noise, zoom
from .audios import LONG_AUDIOS, SHORT_AUDIOS, attacked_audio_path, output_audio_path, raw_audio_path
from .utils import generate_PEAQ_report, read_wave, write_wave

KEY = b'Kei'
DATA = b'Aris -- Princesses of Unnamed Gods' * 32

DATA_LENGTH = 128

def test_LSB():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'LSB'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(output_audio_path(audio, 'LSB'))[0], decoder, extractor) == DATA

        generate_PEAQ_report(audio, 'LSB')

@pytest.mark.skip
def test_LSB_with_mp3():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'mp3'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        add_compression_noise(attacked_audio_path(audio, 'LSB', 'mp3'), attacked_audio_path(audio, 'LSB', 'mp3'), 'mp3')

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'mp3'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_ogg():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'ogg'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        add_compression_noise(attacked_audio_path(audio, 'LSB', 'ogg'), attacked_audio_path(audio, 'LSB', 'ogg'), 'ogg')

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'ogg'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_DCT_with_zoom():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'zoom'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        zoom(attacked_audio_path(audio, 'LSB', 'zoom'), attacked_audio_path(audio, 'LSB', 'zoom'), 1.2)

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'zoom'))[0], decoder, extractor) == DATA

@pytest.mark.benchmark(group='LSB-embed')
def test_benchmark_LSB_embed(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            seed(42)

            embedder = LSBEmbedder(DATA_LENGTH)
            encoder = Encoder(KEY, ECC_LENGTH)
            raw_channels, _, _ = read_wave(raw_audio_path(audio))
            embed(raw_channels, DATA, encoder, embedder)

    benchmark(do)

@pytest.mark.benchmark(group='LSB-extract')
def test_benchmark_LSB_extract(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            extractor = LSBExtractor(DATA_LENGTH)
            decoder = Decoder(KEY, ECC_LENGTH)
            extract(read_wave(output_audio_path(audio, 'LSB'))[0], decoder, extractor)

    for audio in LONG_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'LSB'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

    benchmark(do)
