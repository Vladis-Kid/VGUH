"""
Интегратор с Human Metabolome Database (HMDB).
Использует SOAP API или библиотеку pypath.
"""

import requests
from typing import Dict, Optional, List
import json
from pathlib import Path
import time

from ..config import RAW_DIR, TARGET_COMPOUNDS


class HMDBIntegrator:
    """Класс для работы с HMDB API."""

    def __init__(self):
        self.source_name = "hmdb"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # SOAP API endpoint для HMDB
        self.soap_url = "http://www.hmdb.ca/services/ws/hmdbws.php"

    def search_metabolite(self, query: str) -> Optional[Dict]:
        """
        Поиск метаболита по названию через SOAP API.

        Args:
            query: Название метаболита

        Returns:
            Словарь с данными
        """
        # SOAP запрос для поиска
        soap_body = f"""<?xml version="1.0" encoding="UTF-8"?>
        <SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/">
            <SOAP-ENV:Body>
                <search xmlns="http://www.hmdb.ca">
                    <query>{query}</query>
                </search>
            </SOAP-ENV:Body>
        </SOAP-ENV:Envelope>
        """

        headers = {
            "Content-Type": "text/xml; charset=utf-8",
            "SOAPAction": "http://www.hmdb.ca/search",
        }

        try:
            response = requests.post(
                self.soap_url,
                data=soap_body,
                headers=headers,
                timeout=30
            )
            response.raise_for_status()

            # В реальном проекте нужно парсить SOAP ответ
            return {
                "query": query,
                "status": "success",
                "response_length": len(response.text),
            }
        except Exception as e:
            print(f"Ошибка при поиске {query}: {e}")
            return None

    def get_metabolite_by_id(self, hmdb_id: str) -> Optional[Dict]:
        """
        Получить метаболит по HMDB ID.

        Args:
            hmdb_id: HMDB идентификатор (например, HMDB0000122)

        Returns:
            Словарь с данными
        """
        soap_body = f"""<?xml version="1.0" encoding="UTF-8"?>
        <SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/">
            <SOAP-ENV:Body>
                <getMetabolite xmlns="http://www.hmdb.ca">
                    <accession>{hmdb_id}</accession>
                </getMetabolite>
            </SOAP-ENV:Body>
        </SOAP-ENV:Envelope>
        """

        headers = {
            "Content-Type": "text/xml; charset=utf-8",
            "SOAPAction": "http://www.hmdb.ca/getMetabolite",
        }

        try:
            response = requests.post(
                self.soap_url,
                data=soap_body,
                headers=headers,
                timeout=30
            )
            response.raise_for_status()
            return {
                "hmdb_id": hmdb_id,
                "status": "success",
                "response_length": len(response.text),
            }
        except Exception as e:
            print(f"Ошибка при получении {hmdb_id}: {e}")
            return None

    def fetch_all_targets(self) -> List[Dict]:
        """
        Собрать данные по всем целевым соединениям.

        Returns:
            Список словарей с данными
        """
        results = []
        # Примерные HMDB ID для целевых соединений
        hmdb_ids = {
            "aspirin": "HMDB0001879",
            "caffeine": "HMDB0001847",
            "paracetamol": "HMDB0001859",
            "ibuprofen": "HMDB0001925",
            "naproxen": "HMDB0001926",
        }

        for name, hmdb_id in hmdb_ids.items():
            data = self.get_metabolite_by_id(hmdb_id)
            if data:
                data["name"] = name
                results.append(data)
                output_file = self.output_dir / f"{name.replace(' ', '_')}.json"
                with open(output_file, 'w') as f:
                    json.dump(data, f, indent=2, default=str)

            time.sleep(1)

        return results

    def run(self) -> Dict:
        """
        Запустить интеграцию.

        Returns:
            Словарь с результатами выполнения
        """
        print("🔄 Запуск интеграции с HMDB...")
        results = self.fetch_all_targets()
        print(f"✅ Собрано данных по {len(results)} метаболитам")
        return {
            "source": self.source_name,
            "total_compounds": len(results),
            "compounds": results,
        }


if __name__ == "__main__":
    integrator = HMDBIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
