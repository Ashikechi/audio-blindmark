# pylint: disable=C0103
import wave

import pytest
from aquatk.metrics.PEAQ import PEAQ
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.FFT import FFTEmbedder, FFTExtractor
from audio_blindmark.utils.random import seed

from .utils import add_compression_noise, zoom

KEY = b'Madoka'
ECC_LENGTH = 2
DATA = b'Saa, kanaete yo! Inkyubeitaa!' * 4

WAVE_LENGTH = 4096
DATA_LENGTH = 128

MEDIA_DIR = 'assets/audio/'

REF_AUDIO = MEDIA_DIR + 'sample.wav'
OUTPUT_AUDIO = MEDIA_DIR + 'output_FFT.wav'

PEAQ_REPORT_FILE = 'reports/FFT.txt'

OUTPUT_AUDIO_MP3 = MEDIA_DIR + 'output_FFT.mp3.wav'
OUTPUT_AUDIO_OGG = MEDIA_DIR + 'output_FFT.ogg.wav'
OUTPUT_AUDIO_ZOOM = MEDIA_DIR + 'output_FFT.zoom.wav'

def test_FFT():
    seed(42)
    embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

    peaq = PEAQ(version="advanced")
    result = peaq.analyze_files(REF_AUDIO, OUTPUT_AUDIO)
    with open(PEAQ_REPORT_FILE, 'w', encoding='utf-8') as f:
        f.write(str(result))

@pytest.mark.skip
def test_FFT_with_mp3():
    quantum = 2

    seed(42)
    embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=quantum)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO_MP3, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    add_compression_noise(OUTPUT_AUDIO_MP3, OUTPUT_AUDIO_MP3, 'mp3')

    extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=quantum)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO_MP3, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.skip
def test_FFT_with_ogg():
    quantum = 2

    seed(42)
    embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH, quantum=quantum)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO_OGG, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    add_compression_noise(OUTPUT_AUDIO_OGG, OUTPUT_AUDIO_OGG, 'ogg')

    extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH, quantum=quantum)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO_OGG, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

def test_DCT_with_zoom():
    seed(42)
    embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO_ZOOM, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    zoom(OUTPUT_AUDIO_ZOOM, OUTPUT_AUDIO_ZOOM, 1.2)

    extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO_ZOOM, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.benchmark(group='FFT-embed')
def test_benchmark_FFT_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        with wave.open(REF_AUDIO, 'r') as raw_audio:
            with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)

@pytest.mark.benchmark(group='FFT-extract')
def test_benchmark_FFT_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = FFTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        with wave.open(OUTPUT_AUDIO, 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = FFTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)
