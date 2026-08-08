"""
Интегратор с DrugBank.
Использует REST API (требуется API ключ).
"""

import requests
from typing import Dict, Optional, List
import json
from pathlib import Path

from ..config import RAW_DIR, DRUGBANK_API_KEY, DRUGBANK_API_URL, TARGET_COMPOUNDS


class DrugBankIntegrator:
    """Класс для работы с DrugBank API."""

    def __init__(self):
        self.source_name = "drugbank"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.api_key = DRUGBANK_API_KEY
        self.base_url = DRUGBANK_API_URL
        self.session = requests.Session()
        if self.api_key:
            self.session.headers.update({
                "Authorization": f"Bearer {self.api_key}",
            })

    def get_drug_by_id(self, drugbank_id: str) -> Optional[Dict]:
        """
        Получить информацию о препарате по DrugBank ID.

        Args:
            drugbank_id: DrugBank идентификатор (например, DB00001)

        Returns:
            Словарь с данными
        """
        if not self.api_key:
            print("⚠️ DrugBank API ключ не установлен. Пропускаем.")
            return None

        url = f"{self.base_url}/drugs/{drugbank_id}.json"
        try:
            response = self.session.get(url, timeout=30)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"Ошибка {response.status_code} для {drugbank_id}")
                return None
        except Exception as e:
            print(f"Ошибка при получении {drugbank_id}: {e}")
            return None

    def fetch_all_targets(self) -> List[Dict]:
        """
        Собрать данные по всем целевым соединениям.

        Returns:
            Список словарей с данными
        """
        results = []
        # Примерные DrugBank ID для целевых соединений
        drugbank_ids = {
            "aspirin": "DB00945",
            "caffeine": "DB00201",
            "paracetamol": "DB00316",
            "ibuprofen": "DB01050",
            "naproxen": "DB00788",
        }

        for name, drugbank_id in drugbank_ids.items():
            data = self.get_drug_by_id(drugbank_id)
            if data:
                data["name"] = name
                results.append(data)
                output_file = self.output_dir / f"{name.replace(' ', '_')}.json"
                with open(output_file, 'w') as f:
                    json.dump(data, f, indent=2, default=str)

        return results

    def run(self) -> Dict:
        """
        Запустить интеграцию.

        Returns:
            Словарь с результатами выполнения
        """
        print("🔄 Запуск интеграции с DrugBank...")
        if not self.api_key:
            print("⚠️ DrugBank API ключ не найден. Интеграция пропущена.")
            return {
                "source": self.source_name,
                "status": "skipped",
                "message": "API key not configured",
            }

        results = self.fetch_all_targets()
        print(f"✅ Собрано данных по {len(results)} препаратам")
        return {
            "source": self.source_name,
            "total_compounds": len(results),
            "compounds": results,
        }


if __name__ == "__main__":
    integrator = DrugBankIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
