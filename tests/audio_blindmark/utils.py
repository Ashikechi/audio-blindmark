import wave as w

from pydub import AudioSegment

from audio_blindmark.wave_helper import WaveReadHelper, WaveWriteHelper


def add_compression_noise(audio: str, audio_with_noise: str, form: str) -> None:
    wave = AudioSegment.from_wav(audio)
    wave.export(audio_with_noise, format=form)
    wave = AudioSegment.from_file(audio_with_noise, format=form)
    wave.export(audio_with_noise, format='wav')

def zoom(audio: str, zoomed_audio: str, factor: float) -> None:
    with w.open(audio, 'r') as raw_audio:
        params = raw_audio.getparams()
        helper = WaveReadHelper(raw_audio)
        wave, _ = helper.read(-1)
    wave = [each * factor for each in wave]
    with w.open(audio, 'w') as output_audio:
        output_audio.setparams(params)
        helper = WaveWriteHelper(output_audio)
        helper.write(wave)
