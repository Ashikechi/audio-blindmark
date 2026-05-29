# pylint: disable=C0103
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.LSB import LSBEmbedder, LSBExtractor
from audio_blindmark.utils.random import seed

from .attackers import lossy_compression, pink_noise, resample, white_noise, zoom
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
def test_LSB_with_white_noise():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'white_noise'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        white_noise(attacked_audio_path(audio, 'LSB', 'white_noise'), attacked_audio_path(audio, 'LSB', 'white_noise'), 0.05)

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'white_noise'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_pink_noise():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'pink_noise'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        pink_noise(attacked_audio_path(audio, 'LSB', 'pink_noise'), attacked_audio_path(audio, 'LSB', 'pink_noise'), 0.05)

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'pink_noise'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_mp3():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'mp3'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        lossy_compression(attacked_audio_path(audio, 'LSB', 'mp3'), attacked_audio_path(audio, 'LSB', 'mp3'), 'mp3')

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

        lossy_compression(attacked_audio_path(audio, 'LSB', 'ogg'), attacked_audio_path(audio, 'LSB', 'ogg'), 'ogg')

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'ogg'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_resample():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'resample'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        resample(attacked_audio_path(audio, 'LSB', 'resample'), attacked_audio_path(audio, 'LSB', 'resample'), 24000)

        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'LSB', 'resample'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_zoom():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'LSB', 'zoom'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        zoom(attacked_audio_path(audio, 'LSB', 'zoom'), attacked_audio_path(audio, 'LSB', 'zoom'), 0.8)

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
