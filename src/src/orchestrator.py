"""
Главный оркестратор для запуска всех интеграций.
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List

# Добавляем путь к модулям
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROCESSED_DIR
from integrators.pubchem_integrator import PubChemIntegrator
from integrators.nist_integrator import NISTIntegrator
from integrators.hmdb_integrator import HMDBIntegrator
from integrators.chebi_integrator import ChEBIIntegrator
from integrators.drugbank_integrator import DrugBankIntegrator
from integrators.openchrom_integrator import OpenChromIntegrator


class Orchestrator:
    """
    Оркестратор для управления всеми интеграциями.
    """

    def __init__(self):
        self.integrators = [
            PubChemIntegrator(),
            NISTIntegrator(),
            HMDBIntegrator(),
            ChEBIIntegrator(),
            DrugBankIntegrator(),
            OpenChromIntegrator(),
        ]
        self.timestamp = datetime.now().isoformat()

    def run_all(self) -> Dict:
        """
        Запустить все интеграции последовательно.

        Returns:
            Словарь с результатами всех интеграций
        """
        results = {
            "timestamp": self.timestamp,
            "integrations": [],
        }

        for integrator in self.integrators:
            try:
                result = integrator.run()
                results["integrations"].append(result)
            except Exception as e:
                results["integrations"].append({
                    "source": integrator.source_name,
                    "status": "error",
                    "error": str(e),
                })

        # Сохраняем общий отчёт
        report_file = PROCESSED_DIR / f"integration_report_{self.timestamp.replace(':', '-')}.json"
        with open(report_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        # Сохраняем последний отчёт
        latest_file = PROCESSED_DIR / "latest_report.json"
        with open(latest_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        print(f"\n📊 Отчёт сохранён: {report_file}")
        return results

    def run_source(self, source_name: str) -> Dict:
        """
        Запустить конкретную интеграцию по имени источника.

        Args:
            source_name: Имя источника (pubchem, nist, hmdb, chebi, drugbank, openchrom)

        Returns:
            Словарь с результатами
        """
        for integrator in self.integrators:
            if integrator.source_name == source_name:
                return integrator.run()
        return {"error": f"Источник '{source_name}' не найден"}


if __name__ == "__main__":
    orchestrator = Orchestrator()

    # Если передан аргумент командной строки, запускаем конкретный источник
    if len(sys.argv) > 1:
        source = sys.argv[1]
        result = orchestrator.run_source(source)
        print(json.dumps(result, indent=2, default=str))
    else:
        # Иначе запускаем все
        result = orchestrator.run_all()
        print(json.dumps(result, indent=2, default=str))
