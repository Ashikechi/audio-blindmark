# pylint: disable=C0103
import wave

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.DCT import DCTEmbedder, DCTExtractor
from audio_blindmark.utils.random import seed

MEDIA_DIR = 'assets/audio/'

KEY = b'Ema'
DATA = b'Sakuraba Ema daisuki' * 4

WAVE_LENGTH = 4096
DATA_LENGTH = 64

def test_DCT():
    seed(42)
    embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
        with wave.open(MEDIA_DIR + 'output_DCT.wav', 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
    decoder = Decoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'output_DCT.wav', 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.benchmark(group='DCT-embed')
def test_benchmark_DCT_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, 2)
        with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
            with wave.open(MEDIA_DIR + 'output_DCT.wav', 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)

@pytest.mark.benchmark(group='DCT-extract')
def test_benchmark_DCT_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, 2)
        with wave.open(MEDIA_DIR + 'output_DCT.wav', 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
        with wave.open(MEDIA_DIR + 'output_DCT.wav', 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    benchmark(do)
