from pydub import AudioSegment


def add_compression_noise(audio: str, audio_with_noise: str, form: str) -> None:
    wave = AudioSegment.from_wav(audio)
    wave.export(audio_with_noise, format=form)
    wave = AudioSegment.from_file(audio_with_noise, format=form)
    wave.export(audio_with_noise, format='wav')
