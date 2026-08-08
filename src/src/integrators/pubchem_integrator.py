"""
Интегратор с PubChem PUG REST API.
Использует библиотеку pubchempy.
"""

import pubchempy as pcp
from typing import Dict, List, Optional
import json
from pathlib import Path

from ..config import RAW_DIR, TARGET_COMPOUNDS


class PubChemIntegrator:
    """Класс для работы с PubChem API."""

    def __init__(self):
        self.source_name = "pubchem"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def get_compound_by_cid(self, cid: int) -> Optional[Dict]:
        """
        Получить данные о соединении по CID.

        Args:
            cid: PubChem Compound Identifier

        Returns:
            Словарь с данными соединения
        """
        try:
            compound = pcp.Compound.from_cid(cid)
            return {
                "cid": compound.cid,
                "name": compound.iupac_name or compound.synonyms[0] if compound.synonyms else None,
                "molecular_formula": compound.molecular_formula,
                "molecular_weight": compound.molecular_weight,
                "canonical_smiles": compound.canonical_smiles,
                "iupac_name": compound.iupac_name,
                "synonyms": compound.synonyms,
                "xlogp": compound.xlogp,
                "h_bond_donor_count": compound.h_bond_donor_count,
                "h_bond_acceptor_count": compound.h_bond_acceptor_count,
                "rotatable_bond_count": compound.rotatable_bond_count,
                "heavy_atom_count": compound.heavy_atom_count,
            }
        except Exception as e:
            print(f"Ошибка при получении CID {cid}: {e}")
            return None

    def search_compound(self, name: str) -> Optional[Dict]:
        """
        Найти соединение по названию и получить его свойства.

        Args:
            name: Название соединения

        Returns:
            Словарь с данными соединения
        """
        try:
            compounds = pcp.get_compounds(name, 'name')
            if not compounds:
                return None
            compound = compounds[0]
            return {
                "cid": compound.cid,
                "name": compound.iupac_name or name,
                "molecular_formula": compound.molecular_formula,
                "molecular_weight": compound.molecular_weight,
                "canonical_smiles": compound.canonical_smiles,
                "iupac_name": compound.iupac_name,
                "xlogp": compound.xlogp,
            }
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
        for compound_info in TARGET_COMPOUNDS:
            cid = compound_info.get("cid")
            name = compound_info.get("name")

            if cid:
                data = self.get_compound_by_cid(cid)
            else:
                data = self.search_compound(name)

            if data:
                results.append(data)
                # Сохраняем в файл
                output_file = self.output_dir / f"{data['cid']}.json"
                with open(output_file, 'w') as f:
                    json.dump(data, f, indent=2, default=str)

        return results

    def run(self) -> Dict:
        """
        Запустить интеграцию.

        Returns:
            Словарь с результатами выполнения
        """
        print("🔄 Запуск интеграции с PubChem...")
        results = self.fetch_all_targets()
        print(f"✅ Собрано данных по {len(results)} соединениям")
        return {
            "source": self.source_name,
            "total_compounds": len(results),
            "compounds": results,
        }


if __name__ == "__main__":
    integrator = PubChemIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
