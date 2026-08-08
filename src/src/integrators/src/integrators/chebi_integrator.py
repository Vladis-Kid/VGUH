"""
Интегратор с ChEBI (Chemical Entities of Biological Interest).
Использует библиотеку bioservices.
"""

from bioservices import ChEBI
from typing import Dict, Optional, List
import json
from pathlib import Path

from ..config import RAW_DIR, TARGET_COMPOUNDS


class ChEBIIntegrator:
    """Класс для работы с ChEBI API."""

    def __init__(self):
        self.source_name = "chebi"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.chebi = ChEBI(verbose=False)

    def get_complete_entity(self, chebi_id: str) -> Optional[Dict]:
        """
        Получить полную информацию о сущности по ChEBI ID.

        Args:
            chebi_id: ChEBI идентификатор (например, CHEBI:15365)

        Returns:
            Словарь с данными сущности
        """
        try:
            entity = self.chebi.getCompleteEntity(chebi_id)
            return {
                "chebi_id": chebi_id,
                "entity": entity,
            }
        except Exception as e:
            print(f"Ошибка при получении {chebi_id}: {e}")
            return None

    def get_entity_by_name(self, name: str) -> Optional[Dict]:
        """
        Получить сущность по названию.

        Args:
            name: Название соединения

        Returns:
            Словарь с данными
        """
        try:
            # Поиск через LITE API
            results = self.chebi.findLiteByName(name)
            if results:
                # Берём первый результат
                chebi_id = results[0]
                return self.get_complete_entity(chebi_id)
            return None
        except Exception as e:
            print(f"Ошибка при поиске {name}: {e}")
            return None

    def fetch_all_targets(self) -> List[Dict]:
        """
        Собрать данные по всем целевым соединениям.

        Returns:
            Список словарей с данными
        """
        results = []
        # Примерные ChEBI ID для целевых соединений
        chebi_ids = {
            "aspirin": "CHEBI:15365",
            "caffeine": "CHEBI:27732",
            "paracetamol": "CHEBI:46195",
            "ibuprofen": "CHEBI:5855",
            "naproxen": "CHEBI:7471",
        }

        for name, chebi_id in chebi_ids.items():
            data = self.get_complete_entity(chebi_id)
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
        print("🔄 Запуск интеграции с ChEBI...")
        results = self.fetch_all_targets()
        print(f"✅ Собрано данных по {len(results)} сущностям")
        return {
            "source": self.source_name,
            "total_compounds": len(results),
            "compounds": results,
        }


if __name__ == "__main__":
    integrator = ChEBIIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
