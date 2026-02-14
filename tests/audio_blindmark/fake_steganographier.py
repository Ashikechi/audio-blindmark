import numpy as np

from audio_blindmark.base import BaseEmbedder, BaseExtractor, EmbedError, ExtractError

data_list: list[bytes] = []
data_to_index: dict[bytes, int] = {}

class FakeEmbedder(BaseEmbedder):
    def wave_length(self) -> int:
        return self.data_length() << 3

    def data_length(self) -> int:
        return 255

    def embed(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]], data: bytes) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        assert wave.shape[0] == self.wave_length()

        if round(wave[0]) == round(wave[1]) == round(wave[2]) == round(wave[-3]) == round(wave[-2]) == round(wave[-1]) == 0:
            raise EmbedError

        if data in data_to_index:
            index = data_to_index[data]
        else:
            data_list.append(data)
            index = data_to_index[data] = len(data_list)

        r = wave.copy()
        r[0] = r[1] = r[2] = r[-3] = r[-2] = r[-1] = index
        return r

class FakeExtractor(BaseExtractor):
    def wave_length(self) -> int:
        return self.data_length() << 3

    def data_length(self) -> int:
        return 255

    def extract(self, wave: np.ndarray[tuple[int], np.dtype[np.float64]]) -> bytes:
        assert wave.shape[0] == self.wave_length()

        if not round(wave[0]) == round(wave[1]) == round(wave[2]) == round(wave[-3]) == round(wave[-2]) == round(wave[-1]):
            raise ExtractError

        index = round(wave[0])
        if not 0 < index <= len(data_list):
            raise ExtractError

        return data_list[index - 1]
