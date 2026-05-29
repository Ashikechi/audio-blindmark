from pathlib import Path

AUDIO_DIR = Path('assets/audio/')
REPORT_DIR = Path('reports/')

SHORT_AUDIOS = [
    '1_short',
    '2_short',
    '3_short',
    # '4_short',
    '5_short',
]

LONG_AUDIOS = [
    '1_long',
    '2_long',
    '3_long',
    # '4_long',
    '5_long',
]

ALL_AUDIOS = SHORT_AUDIOS + LONG_AUDIOS

def raw_audio_path(audio_id: str) -> Path:
    return AUDIO_DIR.joinpath(f'{audio_id}.wav')

def output_audio_path(audio_id: str, embedder: str) -> Path:
    return AUDIO_DIR.joinpath(f'{audio_id}-{embedder}.wav')

def report_path(audio_id: str, embedder: str) -> Path:
    return REPORT_DIR.joinpath(f'{audio_id}-{embedder}.txt')

def attacked_audio_path(audio_id: str, embedder: str, attacker: str) -> Path:
    return AUDIO_DIR.joinpath(f'{audio_id}-{embedder}-{attacker}.wav')
