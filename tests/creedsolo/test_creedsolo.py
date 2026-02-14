import random

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from creedsolo import ReedSolomonError, RSCodec


def test_creedsolo_sample():
    rsc = RSCodec(10)
    msg = b'hello world ' * 10
    enc = bytes(iter(rsc.encode(iter(msg))))
    dec, dec_enc, errata_pos = rsc.decode(iter(enc))
    assert bytes(iter(dec)) == msg
    assert bytes(iter(dec_enc)) == enc
    assert list(errata_pos) == []  # pylint: disable=C1803

def test_creedsolo_correction():
    rsc = RSCodec(10)
    msg = b'hello world ' * 10
    enc = bytes(iter(rsc.encode(iter(msg))))
    wrong_enc = bytearray(enc)
    for each in [0, 1, 2, 3, 4]:
        wrong_enc[each] = 0
        dec, dec_enc, errata_pos = rsc.decode(iter(wrong_enc))
        assert bytes(iter(dec)) == msg
        assert bytes(iter(dec_enc)) == enc
        assert list(errata_pos) == list(range(each + 1))
    wrong_enc[5] = 0
    with pytest.raises(ReedSolomonError):
        rsc.decode(iter(wrong_enc))

def test_creedsolo_long():
    rsc = RSCodec(2)
    msg = b'hello world ' * 100
    enc = bytes(iter(rsc.encode(iter(msg))))
    wrong_enc = bytearray(enc)
    wrong_enc[0] = wrong_enc[-1] = 0
    dec, dec_enc, errata_pos = rsc.decode(iter(wrong_enc))
    assert bytes(iter(dec)) == msg
    assert bytes(iter(dec_enc)) == enc
    assert list(errata_pos) == [0, len(enc) - 1]

def test_creedsolo_erase_pos():
    rsc = RSCodec(4)
    msg = b'hello world '
    enc = bytes(iter(rsc.encode(iter(msg))))
    wrong_enc = bytearray(enc)
    wrong_enc[0] = wrong_enc[1] = wrong_enc[2] = 0
    dec, dec_enc, errata_pos = rsc.decode(iter(wrong_enc), erase_pos = [0, 1])
    assert bytes(iter(dec)) == msg
    assert bytes(iter(dec_enc)) == enc
    assert list(errata_pos) == [0, 1, 2]

@pytest.mark.benchmark(group = 'reedsolo-encode')
def test_benchmark_creedsolo_encode(benchmark: BenchmarkFixture):
    random.seed(42)
    rsc = RSCodec(32)
    msg = list(random.randbytes(1024))
    benchmark(rsc.encode, msg)

@pytest.mark.benchmark(group = 'reedsolo-decode')
def test_benchmark_creedsolo_decode(benchmark: BenchmarkFixture):
    random.seed(42)
    rsc = RSCodec(32)
    msg = random.randbytes(1024)
    enc = rsc.encode(iter(msg))
    for i in random.sample(range(0, len(enc)), 16):
        enc[i] = random.randint(0, 255)
    benchmark(rsc.decode, enc)
