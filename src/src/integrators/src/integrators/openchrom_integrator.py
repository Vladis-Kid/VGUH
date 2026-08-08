"""
Интегратор с OpenChrom для обработки хроматографических данных.
"""

import subprocess
import json
from typing import Dict, Optional, List
from pathlib import Path

from ..config import RAW_DIR, PROCESSED_DIR, CHROMATOGRAPHY_PARAMS


class OpenChromIntegrator:
    """
    Класс для работы с OpenChrom.
    Предполагает наличие установленного OpenChrom с CLI-интерфейсом.
    """

    def __init__(self):
        self.source_name = "openchrom"
        self.output_dir = RAW_DIR / self.source_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.processed_dir = PROCESSED_DIR / self.source_name
        self.processed_dir.mkdir(parents=True, exist_ok=True)

    def convert_chromatogram(self, input_file: Path, output_format: str = "csv") -> Optional[Path]:
        """
        Конвертировать хроматографический файл в другой формат.

        Args:
            input_file: Путь к входному файлу
            output_format: Выходной формат (csv, json, mzml)

        Returns:
            Путь к сконвертированному файлу или None
        """
        output_file = self.processed_dir / f"{input_file.stem}.{output_format}"

        # Пример команды для OpenChrom CLI (зависит от версии)
        # В реальном проекте нужно использовать официальный API
        cmd = [
            "openchrom",
            "convert",
            "-i", str(input_file),
            "-o", str(output_file),
            "-f", output_format,
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                return output_file
            else:
                print(f"Ошибка конвертации: {result.stderr}")
                return None
        except Exception as e:
            print(f"Ошибка при запуске OpenChrom: {e}")
            return None

    def detect_peaks(self, input_file: Path) -> Optional[List[Dict]]:
        """
        Обнаружить пики в хроматограмме.

        Args:
            input_file: Путь к файлу с хроматограммой

        Returns:
            Список обнаруженных пиков
        """
        # В реальном проекте здесь должен быть вызов OpenChrom API
        # для обнаружения и интеграции пиков

        # Заглушка для примера
        return [
            {
                "peak_number": 1,
                "retention_time_min": 3.42,
                "area": 1245678.9,
                "height": 98765.4,
                "width_at_half_height": 0.12,
                "symmetry": 1.05,
            },
            {
                "peak_number": 2,
                "retention_time_min": 5.78,
                "area": 876543.2,
                "height": 65432.1,
                "width_at_half_height": 0.15,
                "symmetry": 0.98,
            },
        ]

    def run(self) -> Dict:
        """
        Запустить интеграцию с OpenChrom.

        Returns:
            Словарь с результатами выполнения
        """
        print("🔄 Запуск интеграции с OpenChrom...")

        # Поиск файлов для обработки
        input_files = list(self.output_dir.glob("*.chrom")) + list(self.output_dir.glob("*.cdf"))

        if not input_files:
            print("⚠️ Файлы для обработки не найдены. Создаём тестовые данные.")
            # Создаём тестовый файл для примера
            test_file = self.output_dir / "test_sample.chrom"
            test_file.touch()
            input_files = [test_file]

        results = []
        for input_file in input_files:
            # Конвертация
            converted = self.convert_chromatogram(input_file, "json")
            if converted:
                # Обнаружение пиков
                peaks = self.detect_peaks(input_file)
                results.append({
                    "input_file": str(input_file),
                    "converted_file": str(converted) if converted else None,
                    "peaks": peaks,
                    "chromatography_params": CHROMATOGRAPHY_PARAMS,
                })

        print(f"✅ Обработано {len(results)} файлов")
        return {
            "source": self.source_name,
            "total_files": len(results),
            "results": results,
        }


if __name__ == "__main__":
    integrator = OpenChromIntegrator()
    result = integrator.run()
    print(json.dumps(result, indent=2, default=str))
