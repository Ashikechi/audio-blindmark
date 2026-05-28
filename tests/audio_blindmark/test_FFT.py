# pylint: disable=C0103
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.FFT import FFTEmbedder, FFTExtractor
from audio_blindmark.utils.random import seed

from .attackers import lossy_compression, resample, white_noise, zoom
from .audios import LONG_AUDIOS, SHORT_AUDIOS, attacked_audio_path, output_audio_path, raw_audio_path
from .utils import generate_PEAQ_report, read_wave, write_wave

KEY = b'Madoka'
DATA = b'Saa, kanaete yo! Inkyubeitaa!' * 4

WAVE_LENGTH = 4096
DATA_LENGTH = 128

def test_FFT():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'FFT'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(output_audio_path(audio, 'FFT'))[0], decoder, extractor) == DATA

        generate_PEAQ_report(audio, 'FFT')

def test_FFT_with_white_noise():
    ECC_LENGTH = 16
    QUANTUM = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'FFT', 'white_noise'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        white_noise(attacked_audio_path(audio, 'FFT', 'white_noise'), attacked_audio_path(audio, 'FFT', 'white_noise'), 0.01)

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'FFT', 'white_noise'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_FFT_with_mp3():
    ECC_LENGTH = 2
    QUANTUM = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'FFT', 'mp3'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        lossy_compression(attacked_audio_path(audio, 'FFT', 'mp3'), attacked_audio_path(audio, 'FFT', 'mp3'), 'mp3')

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'FFT', 'mp3'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_FFT_with_ogg():
    ECC_LENGTH = 2
    QUANTUM = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'FFT', 'ogg'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        lossy_compression(attacked_audio_path(audio, 'FFT', 'ogg'), attacked_audio_path(audio, 'FFT', 'ogg'), 'ogg')

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'FFT', 'ogg'))[0], decoder, extractor) == DATA

@pytest.mark.skip
def test_FFT_with_resample():
    ECC_LENGTH = 16
    QUANTUM = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'FFT', 'resample'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        resample(attacked_audio_path(audio, 'FFT', 'resample'), attacked_audio_path(audio, 'FFT', 'resample'), 36000)

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=QUANTUM)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'FFT', 'resample'))[0], decoder, extractor) == DATA

def test_FFT_with_zoom():
    ECC_LENGTH = 2

    for audio in SHORT_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(attacked_audio_path(audio, 'FFT', 'zoom'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

        zoom(attacked_audio_path(audio, 'FFT', 'zoom'), attacked_audio_path(audio, 'FFT', 'zoom'), 1.2)

        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        assert extract(read_wave(attacked_audio_path(audio, 'FFT', 'zoom'))[0], decoder, extractor) == DATA

@pytest.mark.benchmark(group='FFT-embed')
def test_benchmark_FFT_embed(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            seed(42)

            embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
            encoder = Encoder(KEY, ECC_LENGTH)
            raw_channels, _, _ = read_wave(raw_audio_path(audio))
            embed(raw_channels, DATA, encoder, embedder)

    benchmark(do)

@pytest.mark.benchmark(group='FFT-extract')
def test_benchmark_FFT_extract(benchmark: BenchmarkFixture):
    ECC_LENGTH = 2

    def do():
        for audio in LONG_AUDIOS:
            extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
            decoder = Decoder(KEY, ECC_LENGTH)
            extract(read_wave(output_audio_path(audio, 'FFT'))[0], decoder, extractor)

    for audio in LONG_AUDIOS:
        seed(42)

        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        raw_channels, width, framerate = read_wave(raw_audio_path(audio))
        write_wave(output_audio_path(audio, 'FFT'), embed(raw_channels, DATA, encoder, embedder), width, framerate)

    benchmark(do)
