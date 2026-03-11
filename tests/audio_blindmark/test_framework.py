import wave

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.utils.random import seed

from .fake_steganographier import FakeEmbedder, FakeExtractor

KEY = b'Kei'
DATA = b'Aris -- Princesses of Unnamed Gods' * 64

MEDIA_DIR = 'assets/audio/'

REF_AUDIO = MEDIA_DIR + 'sample.wav'
OUTPUT_AUDIO = MEDIA_DIR + 'output_fake.wav'

def test_framework():
    seed(42)
    embedder = FakeEmbedder()
    encoder = Encoder(KEY, 2)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = FakeExtractor()
    decoder = Decoder(KEY, 2)
    with wave.open(OUTPUT_AUDIO, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.benchmark(group='framework-embed')
def test_benchmark_framework_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = FakeEmbedder()
        encoder = Encoder(KEY, 2)
        with wave.open(REF_AUDIO, 'r') as raw_audio:
            with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)

@pytest.mark.benchmark(group='framework-extract')
def test_benchmark_framework_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = FakeExtractor()
        decoder = Decoder(KEY, 2)
        with wave.open(OUTPUT_AUDIO, 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = FakeEmbedder()
    encoder = Encoder(KEY, 2)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    benchmark(do)
