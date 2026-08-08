

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
