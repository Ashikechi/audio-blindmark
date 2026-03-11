# pylint: disable=C0103
import wave

import pytest
from aquatk.metrics.PEAQ import PEAQ
from pytest_benchmark.fixture import BenchmarkFixture

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.steganographier.DCT import DCTEmbedder, DCTExtractor
from audio_blindmark.utils.random import seed

KEY = b'Ema'
DATA = b'Sakuraba Ema daisuki' * 4

WAVE_LENGTH = 4096
DATA_LENGTH = 64

MEDIA_DIR = 'assets/audio/'

REF_AUDIO = MEDIA_DIR + 'sample.wav'
OUTPUT_AUDIO = MEDIA_DIR + 'output_DCT.wav'

PEAQ_REPORT_FILE = 'reports/DCT.txt'

def test_DCT():
    seed(42)
    embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, 2)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
    decoder = Decoder(KEY, 2)
    with wave.open(OUTPUT_AUDIO, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA

    peaq = PEAQ(version="advanced")
    result = peaq.analyze_files(REF_AUDIO, OUTPUT_AUDIO)
    with open(PEAQ_REPORT_FILE, 'w', encoding='utf-8') as f:
        f.write(str(result))

@pytest.mark.benchmark(group='DCT-embed')
def test_benchmark_DCT_embed(benchmark: BenchmarkFixture):
    def do():
        seed(42)
        embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
        encoder = Encoder(KEY, 2)
        with wave.open(REF_AUDIO, 'r') as raw_audio:
            with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
                output_audio.setparams(raw_audio.getparams())
                embed(raw_audio, DATA, output_audio, encoder, embedder)
    benchmark(do)

@pytest.mark.benchmark(group='DCT-extract')
def test_benchmark_DCT_extract(benchmark: BenchmarkFixture):
    def do():
        extractor = DCTExtractor(WAVE_LENGTH, DATA_LENGTH)
        decoder = Decoder(KEY, 2)
        with wave.open(OUTPUT_AUDIO, 'r') as audio:
            extract(audio, decoder, extractor)

    seed(42)
    embedder = DCTEmbedder(WAVE_LENGTH, DATA_LENGTH)
    encoder = Encoder(KEY, 2)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    benchmark(do)
