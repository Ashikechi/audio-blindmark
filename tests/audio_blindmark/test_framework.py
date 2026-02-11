import wave

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable = W0401
from audio_blindmark.utils.random import seed

from .fake_steganographier import FakeEmbedder, FakeExtractor

MEDIA_DIR = 'assets/audio/'

KEY = b'Kei'
DATA = b'Aris -- Princesses of Unnamed Gods' * 100

def test_framework():
    seed(42)
    embedder = FakeEmbedder()
    encoder = Encoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
        with wave.open(MEDIA_DIR + 'output.wav', 'w') as output_audio:
            output_audio: wave.Wave_write
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = FakeExtractor()
    decoder = Decoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'output.wav', 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.benchmark(group = 'framework-embed')
def test_benchmark_framework_embed(benchmark: BenchmarkFixture):
    def foo():
        seed(42)
        embedder = FakeEmbedder()
        encoder = Encoder(KEY, 2)
        with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
            with wave.open(MEDIA_DIR + 'output.wav', 'w') as output_audio:
                output_audio: wave.Wave_write
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(foo)

@pytest.mark.benchmark(group = 'framework-extract')
def test_benchmark_framework_extract(benchmark: BenchmarkFixture):
    def foo():
        extractor = FakeExtractor()
        decoder = Decoder(KEY, 2)
        with wave.open(MEDIA_DIR + 'output.wav', 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = FakeEmbedder()
    encoder = Encoder(KEY, 2)
    with wave.open(MEDIA_DIR + 'test_short.wav', 'r') as raw_audio:
        with wave.open(MEDIA_DIR + 'output.wav', 'w') as output_audio:
            output_audio: wave.Wave_write
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    benchmark(foo)
