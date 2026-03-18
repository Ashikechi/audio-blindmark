# pylint: disable=C0103
import wave

import pytest
from aquatk.metrics.PEAQ import PEAQ
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.LSB import LSBEmbedder, LSBExtractor
from audio_blindmark.utils.random import seed

from .noise import add_compression_noise

KEY = b'Kei'
ECC_LENGTH = 64
DATA = b'Aris -- Princesses of Unnamed Gods' * 32

DATA_LENGTH = 128

MEDIA_DIR = 'assets/audio/'

REF_AUDIO = MEDIA_DIR + 'sample.wav'
OUTPUT_AUDIO = MEDIA_DIR + 'output_LSB.wav'

PEAQ_REPORT_FILE = 'reports/LSB.txt'

OUTPUT_AUDIO_MP3 = MEDIA_DIR + 'output_LSB.mp3.wav'
OUTPUT_AUDIO_OGG = MEDIA_DIR + 'output_LSB.ogg.wav'

def test_LSB():
    seed(42)
    embedder = LSBEmbedder(DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = LSBExtractor(DATA_LENGTH)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

    peaq = PEAQ(version="advanced")
    result = peaq.analyze_files(REF_AUDIO, OUTPUT_AUDIO)
    with open(PEAQ_REPORT_FILE, 'w', encoding='utf-8') as f:
        f.write(str(result))

@pytest.mark.skip
def test_LSB_with_mp3():
    seed(42)
    embedder = LSBEmbedder(DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO_MP3, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    add_compression_noise(OUTPUT_AUDIO_MP3, OUTPUT_AUDIO_MP3, 'mp3')

    extractor = LSBExtractor(DATA_LENGTH)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO_MP3, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.skip
def test_LSB_with_ogg():
    seed(42)
    embedder = LSBEmbedder(DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO_OGG, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    add_compression_noise(OUTPUT_AUDIO_OGG, OUTPUT_AUDIO_OGG, 'ogg')

    extractor = LSBExtractor(DATA_LENGTH)
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO_OGG, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

@pytest.mark.benchmark(group='LSB-embed')
def test_benchmark_LSB_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = LSBEmbedder(DATA_LENGTH)
        encoder = Encoder(KEY, ECC_LENGTH)
        with wave.open(REF_AUDIO, 'r') as raw_audio:
            with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)

@pytest.mark.benchmark(group='LSB-extract')
def test_benchmark_LSB_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = LSBExtractor(DATA_LENGTH)
        decoder = Decoder(KEY, ECC_LENGTH)
        with wave.open(OUTPUT_AUDIO, 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = LSBEmbedder(DATA_LENGTH)
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    benchmark(do)
