import wave

from audio_blindmark import *  # pylint: disable=W0401
from audio_blindmark.utils.random import seed

from .fake_steganographier import FakeEmbedder, FakeExtractor

KEY = b'Kei'
ECC_LENGTH = 2
DATA = b'Aris -- Princesses of Unnamed Gods' * 64

MEDIA_DIR = 'assets/audio/'

REF_AUDIO = MEDIA_DIR + 'sample.wav'
OUTPUT_AUDIO = MEDIA_DIR + 'output_fake.wav'

def test_framework():
    seed(42)
    embedder = FakeEmbedder()
    encoder = Encoder(KEY, ECC_LENGTH)
    with wave.open(REF_AUDIO, 'r') as raw_audio:
        with wave.open(OUTPUT_AUDIO, 'w') as output_audio:
            output_audio.setparams(raw_audio.getparams())
            embed(raw_audio, DATA, output_audio, encoder, embedder)

    extractor = FakeExtractor()
    decoder = Decoder(KEY, ECC_LENGTH)
    with wave.open(OUTPUT_AUDIO, 'r') as audio:
        assert extract(audio, decoder, extractor) == DATA
