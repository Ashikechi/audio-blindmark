from typing import Optional, Union

from numpy.random import PCG64
from numpy.random.bit_generator import BitGenerator, SeedSequence

rng: BitGenerator = PCG64()

def to_bit_generator(random_state: Optional[Union[BitGenerator, SeedSequence, int]]) -> BitGenerator:
    if random_state is None:
        return PCG64()
    if isinstance(random_state, BitGenerator):
        return random_state
    return PCG64(random_state)

def get_rng() -> BitGenerator:
    return rng

def seed(random_state: Optional[Union[BitGenerator, SeedSequence, int]]) -> None:
    global rng  # pylint: disable=W0603
    rng = to_bit_generator(random_state)
