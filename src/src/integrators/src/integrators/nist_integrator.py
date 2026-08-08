"""
Интегратор с NIST Chemistry WebBook.
Использует неофициальную библиотеку nistchempy.
"""

import requests
from typing import Dict, Optional, List
import json
from pathlib import Path
import time

from ..config import RAW_DIR, NIST_USER_AGENT, TARGET_COMPOUNDS


class NISTIntegrator:
    """Класс для работы с NIST Chemistry WebBook."""

    def __init__(self):
        self.source_name = "nist"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": NIST_USER_AGENT,
        })
        # Базовая URL для NIST Webbook
        self.base_url = "https://webbook.nist.gov/cgi/cbook.cgi"

    def get_compound_properties(self, name: str) -> Optional[Dict]:
        """
        Получить свойства соединения по названию.

        Args:
            name: Название соединения

        Returns:
            Словарь с данными
        """
        params = {
            "Name": name,
            "Mask": 200,  # Маска для получения всех свойств
            "Units": "SI",
        }

        try:
            response = self.session.get(self.base_url, params=params, timeout=30)
            response.raise_for_status()

            # Парсинг HTML — в реальном проекте используйте BeautifulSoup
            # Здесь упрощённая заглушка для примера
            return {
                "name": name,
                "source": "NIST Webbook",
                "url": response.url,
                "status": "success",
                "html_length": len(response.text),
            }
        except Exception as e:
            print(f"Ошибка при получении данных для {name}: {e}")
            return None

    def fetch_all_targets(self) -> List[Dict]:
        """
        Собрать данные по всем целевым соединениям.

        Returns:
            Список словарей с данными
        """
        results = []
        for compound_info in TARGET_COMPOUNDS:
            name = compound_info.get("name")
            if not name:
                continue

            data = self.get_compound_properties(name)
            if data:
                results.append(data)
                # Сохраняем в файл
                output_file = self.output_dir / f"{name.replace(' ', '_')}.json"
                with open(output_file, 'w') as f:
                    json.dump(data, f, indent=2, default=str)

            time.sleep(1)  # Вежливость к серверу

        return results

    def run(self) -> Dict:
        """
        Запустить интеграцию.

        Returns:
            Словарь с результатами выполнения
        """
        print("🔄 Запуск интеграции с NIST Chemistry WebBook...")
        results = self.fetch_all_targets()
        print(f"✅ Собрано данных по {len(results)} соединениям")
        return {
            "source": self.source_name,
            "total_compounds": len(results),
            "compounds": results,
        }


if __name__ == "__main__":
    integrator = NISTIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
