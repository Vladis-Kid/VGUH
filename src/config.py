"""
Конфигурация модуля интеграции данных.
"""

import os
from pathlib import Path

# Базовые пути
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Создаём директории, если их нет
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Настройки API
PUBCHEM_TIMEOUT = 30
NIST_USER_AGENT = "HPLC-Data-Integration/1.0"
DRUGBANK_API_KEY = os.getenv("DRUGBANK_API_KEY", "")
DRUGBANK_API_URL = "https://go.drugbank.com"

# Список целевых соединений для примера (можно расширять)
TARGET_COMPOUNDS = [
    {"name": "aspirin", "cid": 2244},
    {"name": "caffeine", "cid": 2519},
    {"name": "paracetamol", "cid": 1983},
    {"name": "ibuprofen", "cid": 3672},
    {"name": "naproxen", "cid": 156391},
]

# Параметры хроматографии для моделирования (пример)
CHROMATOGRAPHY_PARAMS = {
    "column_length_mm": 150,
    "column_diameter_mm": 4.6,
    "particle_size_um": 5,
    "flow_rate_ml_min": 1.0,
    "temperature_c": 25,
}
