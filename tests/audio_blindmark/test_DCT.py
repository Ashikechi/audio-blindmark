# pylint: disable=C0103
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.DCT import DCTEmbedder, DCTExtractor
from audio_blindmark.utils.random import seed

from .attackers import add_compression_noise, zoom
from .audios import LONG_AUDIOS, SHORT_AUDIOS, attacked_audio_path, output_audio_path, raw_audio_path
from .utils import generate_PEAQ_report, read_wave, write_wave

KEY = b'Ema'
DATA = b'Sakuraba Ema daisuki' * 4

WAVE_LENGTH = 4096
DATA_LENGTH = 128

def test_DCT():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'DCT'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(output_audio_path(audio, 'DCT'))[0], decoder, extractor) == DATA

        generate_PEAQ_report(audio, 'DCT')

@pytest.mark.skip
def test_DCT_with_mp3():
    ECC_LENGTH = 2
    QUANTUM = 1 / 16

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'DCT', 'mp3'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        add_compression_noise(attacked_audio_path(audio, 'DCT', 'mp3'), attacked_audio_path(audio, 'DCT', 'mp3'), 'mp3')

        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'DCT', 'mp3'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_DCT_with_ogg():
    ECC_LENGTH = 2
    QUANTUM = 1 / 16

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'DCT', 'ogg'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        add_compression_noise(attacked_audio_path(audio, 'DCT', 'ogg'), attacked_audio_path(audio, 'DCT', 'ogg'), 'ogg')

        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'DCT', 'ogg'))[0], decoder, extractor) == DATA

def test_DCT_with_zoom():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'DCT', 'zoom'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        zoom(attacked_audio_path(audio, 'DCT', 'zoom'), attacked_audio_path(audio, 'DCT', 'zoom'), 1.2)

        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'DCT', 'zoom'))[0], decoder, extractor) == DATA

@pytest.mark.benchmark(group='DCT-embed')
def test_benchmark_DCT_embed(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            seed(42)

            embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
            encoder = Encoder(KEY, ECC_LENGTH)
            raw_channels, _, _ = read_wave(raw_audio_path(audio))
            embed(raw_channels, DATA, encoder, embedder)

    benchmark(do)

@pytest.mark.benchmark(group='DCT-extract')
def test_benchmark_DCT_extract(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
            decoder = Decoder(KEY, ECC_LENGTH)
            extract(read_wave(output_audio_path(audio, 'DCT'))[0], decoder, extractor)

    for audio in LONG_AUDIOS:
        seed(42)

        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'DCT'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

    benchmark(do)
