

# 1. Клонируем репозиторий
git clone https://github.com/your-username/hplc-data-integration.git
cd hplc-data-integration

# 2. Создаём виртуальное окружение
python -m venv venv
source venv/bin/activate  # Linux/Mac
# или venv\Scripts\activate  # Windows

# 3. Устанавливаем зависимости
pip install -r requirements.txt

# 4. Настраиваем переменные окружения (опционально)
export DRUGBANK_API_KEY="your_api_key_here"  # Linux/Mac
# или set DRUGBANK_API_KEY=your_api_key_here  # Windows

# 5. Запускаем все интеграции
python src/orchestrator.py

# 6. Запускаем конкретную интеграцию
python src/orchestrator.py pubchem


Настройка GitHub Secrets
Для работы с DrugBank API в GitHub Actions:
Перейдите в Settings → Secrets and variables → Actions
Нажмите New repository secret
Имя: DRUGBANK_API_KEY
Значение: ваш API ключ от DrugBank


Структура репозитория:
hplc-data-integration/
├── .github/
│   └── workflows/
│       └── daily_update.yml          # GitHub Actions для автоматического запуска
├── src/
│   ├── __init__.py
│   ├── config.py                     # Конфигурация (API ключи, настройки)
│   ├── orchestrator.py               # Главный оркестратор, запускает все интеграции
│   ├── integrators/
│   │   ├── __init__.py
│   │   ├── pubchem_integrator.py     # PubChem API
│   │   ├── nist_integrator.py        # NIST Chemistry WebBook
│   │   ├── hmdb_integrator.py        # Human Metabolome Database
│   │   ├── chebi_integrator.py       # ChEBI
│   │   ├── drugbank_integrator.py    # DrugBank
│   │   └── openchrom_integrator.py   # OpenChrom (обработка файлов)
│   └── utils/
│       ├── __init__.py
│       └── file_utils.py             # Утилиты для работы с файлами
├── data/                             # Сюда сохраняются собранные данные
│   ├── raw/                          # Сырые JSON-ответы от API
│   └── processed/                    # Обработанные, структурированные данные
├── requirements.txt
├── setup.py
└── README.md




# HPLC Simulation & Modeling Module

Программный модуль для симуляции и математического моделирования ВЭЖХ:
задание условий (колонка, подвижная фаза, температура, давление), live-правки
без перезапуска, изократический/градиентный/ступенчатый режимы, расчёт
tR, k', α, Rs, N, HETP, формы пика (EMG), статистика по повторам,
Monte-Carlo/GA оптимизация условий, и открытый веб-дашборд в стиле
Waters Empower / Agilent OpenLab.

## Структура репозитория

```
hplc_sim/
  models.py          # аналитические формулы: LSS-удержание, ван’т Гофф, градиент (Snyder-Dolan), EMG-пик
  components.py       # Column, MobilePhase, Compound, Method, GradientStep, Detector (dataclasses)
  metrics.py           # k', α, N, HETP, Rs, asymmetry, доверительные интервалы
  simulator.py         # HPLCSimulator: live-edit API + run() + run_replicates()
  optimization.py      # Monte-Carlo и генетический алгоритм подбора условий
  api.py                # FastAPI-сервер (сессии, PATCH для live-правок, /run, /optimize)
  integrations/
    pubchem.py          # PubChem PUG-REST -> logP/MW -> QSAR-seed для k_w/S
    nist_webbook.py      # NIST Chemistry WebBook + локальная NIST RI Database
    retip_client.py       # мост к Retip/pyRetip (ML-предсказание tR по SMILES)
    chromrs_bridge.py     # мост к chrom-rs (Rust, Ленгмюр-изотермы, RK4) для препаративных расчётов
    github_datasets.py     # загрузка референсных CSV с GitHub + сравнение sim vs эксперимент
    lims.py                 # экспорт в LIMS (REST push, AnIML-lite XML)
tests/
  test_models.py                 # модульные тесты (retention/metrics/live-edit/gradient)
  reference_mixture_data.csv      # шаблон референсных данных для валидации (ЗАМЕНИТЬ на реальные)
  validate_against_reference.py   # пример: sim vs эксперимент, RMSE/MAE
dashboard/
  index.html            # интерактивный веб-дашборд (Plotly.js), стиль Empower/OpenLab
requirements.txt
```

## Быстрый старт

```bash
pip install -r requirements.txt
pytest tests/ -q
uvicorn hplc_sim.api:app --reload --port 8000
# затем открыть dashboard/index.html в браузере (по умолчанию бьёт в http://localhost:8000)
```

### Программное использование

```python
from hplc_sim import Column, MobilePhase, Method, Compound, Detector, HPLCSimulator

col = Column(length_mm=150, id_mm=4.6, particle_um=3.5)
mp = MobilePhase(flow_ml_min=1.0, temperature_C=25)
method = Method(mode="gradient", run_time_min=20, dwell_time_min=0.4)
compounds = [Compound(name="Caffeine", k_w=6, S=3.4), Compound(name="Phenol", k_w=4, S=3.2)]

sim = HPLCSimulator(col, mp, method, Detector(), compounds)
sim.set_gradient([{"time_min": 0, "phi_B": 0.05}, {"time_min": 15, "phi_B": 0.9}])

result = sim.run()          # -> {time_min, signal_mAU, peaks:[{tR,k,N,HETP_um,Rs_vs_prev,alpha_vs_prev,...}], ...}

# LIVE-ПРАВКА без пересоздания объекта:
sim.patch_flow(1.5)
sim.patch_temperature(35)
sim.add_component(Compound(name="NewAnalyte", k_w=9, S=3.6))
result2 = sim.run()
```

## Математические модели (реализовано)

| Модель | Формула / метод | Файл |
|---|---|---|
| Изократическое удержание (LSS) | `log k = log k_w - S·φ` | `models.k_isocratic` |
| Температурная поправка | ван’т Гофф: `ln k_T = ln k_ref - (ΔH/R)(1/T - 1/T_ref)` | `models.van_t_hoff_correction` |
| Градиентное элюирование | численное интегрирование фундаментального уравнения Snyder–Dolan LSS-градиента | `models.gradient_retention_time` |
| Форма пика | Exponentially Modified Gaussian (Foley–Dorsey), численно устойчивая реализация через `erfcx` | `models.emg_peak` |
| k', α, N, HETP, Rs, As | классические хроматографические соотношения (USP-совместимые) | `metrics.py` |
| Доверительные интервалы | нормальное приближение по повторным инжекциям (RSD%, CI 95%) | `metrics.confidence_interval` |
| Оптимизация условий | Monte-Carlo random search / real-valued генетический алгоритм | `optimization.py` |

Архитектура моделей расширяема: для случаев, требующих высокой точности
(перегруженные колонки, нелинейная конкурентная адсорбция), ядро может
делегировать расчёт внешнему PDE-решателю **chrom-rs** (см. `chromrs_bridge.py`)
— выбор модели "по требуемой точности" реализован как отдельный опциональный бэкенд,
а не встроен в основной путь расчёта (чтобы live-правки оставались мгновенными).

## Каталог интеграций (бесплатные/открытые ресурсы)

| Категория | Инструмент | Лицензия/доступ | Файл-обёртка в этом репо |
|---|---|---|---|
| Симуляция (референс) | HPLC Simulator (hplcsimulator.org) | CC, бесплатно | — используется как источник валидации/сверки формул |
| Симуляция (PDE, высокая точность) | chrom-rs (github.com/biface/chromatography) | открытый, Rust | `integrations/chromrs_bridge.py` |
| Синтетика данных | mzrtsim (R) | открытый | см. `integrations/github_datasets.py` (шаблон загрузки) |
| Предсказание tR (ML) | Retip / pyRetip (github.com/oloBion/Retip) | открытый | `integrations/retip_client.py` |
| Предсказание tR (ML, вебсервис) | RT-Pred (rtpred.ca) | бесплатный веб-сервис | вызывается как внешний REST (см. пример ниже) |
| Предсказание tR (deep learning) | RT-Transformer | открытый (PyTorch) | точка расширения — заменить `predict_rt` в `retip_client.py` |
| Обработка/визуализация | OpenChrom, Appia | открытые | конвейер импорта вендорских форматов (внешний шаг) |
| Хим. свойства (QSAR-seed) | PubChem PUG-REST | бесплатно, без ключа | `integrations/pubchem.py` |
| Хим. свойства / удержание (референс) | NIST Chemistry WebBook, NIST RI Database | бесплатно | `integrations/nist_webbook.py` |
| Доп. базы | HMDB, ChEBI, DrugBank, KNApSAcK | открытые API | добавляются по аналогии с `pubchem.py` |
| Датасеты для валидации | icredd-cheminfo/chromatography-modeling, mpho-mafata/Chromatographic-data | открытые GitHub-репо | `integrations/github_datasets.py::KNOWN_REPOS` |
| Формат обмена | GC2ASM (Allotrope Simple Model) | открытый | точка расширения для стандартизации выгрузки |
| LIMS | generic REST / AnIML-lite | — | `integrations/lims.py` |

### Пример: подтянуть QSAR-затравку параметров из PubChem

```python
from hplc_sim.integrations.pubchem import fetch_properties, estimate_rp_params_from_logp

props = fetch_properties("ibuprofen")
seed = estimate_rp_params_from_logp(props["XLogP"])
# {'k_w': ..., 'S': ...}  -> используем как стартовые Compound(k_w=..., S=...)
```

### Пример: RT-Pred как внешний REST (без обёртки в этом репо, т.к. это сторонний сервис)

```python
import requests
r = requests.post("https://rtpred.ca/api/predict", json={"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"})
```

## Валидация

`tests/validate_against_reference.py` считает симуляцию для набора соединений
и сравнивает с CSV-таблицей экспериментальных tR (RMSE/MAE, относительная
ошибка на компонент). Поставляемый `reference_mixture_data.csv` —
**иллюстративный шаблон** (для демонстрации пайплайна offline); для реальной
валидации на «нескольких десятках тестовых смесей» замените его на:
- собственные лабораторные данные,
- выгрузку из NIST Retention Index Database (`nist_webbook.load_local_ri_database`),
- набор из открытого GitHub-репозитория (`github_datasets.fetch_raw_csv`).

## Возможные расширения (реализуются подключением модулей выше, без переписывания ядра)

- **LIMS**: `integrations/lims.py::push_result_rest` — REST-пуш результатов; или `export_animl_lite` для файлового обмена.
- **QSAR по структуре**: замените seed-эвристику в `pubchem.py` на вызов `retip_client.predict_rt(smiles)` или RT-Transformer.
- **Оптимизация условий**: `optimization.monte_carlo_optimize` / `genetic_optimize` — минимизация времени анализа при заданном целевом Rs.
- **Другие методы (ГХ, ИЭФ)**: `Column.ctype` уже поддерживает произвольные типы; для ГХ замените `k_isocratic`/градиентную модель на модель линейного индекса удерживания (Kovats RI) — модели изолированы в `models.py`, остальной код (metrics/simulator/api) переиспользуется без изменений.
- **Авто-интеграции**: `api.py` спроектирован как сессионный REST-сервис — любой внешний сервис (n8n, Zapier, корпоративный LIMS-коннектор) может дергать `/sessions`, `PATCH /sessions/{id}`, `/run` без изменений в коде модуля.
