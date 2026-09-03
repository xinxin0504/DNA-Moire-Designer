#!/usr/bin/env python3
"""Fast scale-bar preflight for the Moiré-analysis upload step.

This intentionally reuses the canonical detector/OCR implementation from the
full image-analysis worker while avoiding FFT and real-space domain analysis.
"""

from __future__ import annotations

import argparse
import json

try:
    from .image_analysis_worker import detect_scale_bar, ocr_scale, read_pgm
except ImportError:  # Executed directly in a development checkout.
    from image_analysis_worker import detect_scale_bar, ocr_scale, read_pgm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pgm")
    parser.add_argument("--original")
    args = parser.parse_args()
    image = read_pgm(args.pgm)
    bar = detect_scale_bar(image)
    scale_nm, ocr_text = ocr_scale(args.original, bar)
    result = {
        "width": int(image.shape[1]),
        "height": int(image.shape[0]),
        "scale_bar": bar,
        "scale_value_nm": scale_nm,
        "ocr_text": ocr_text,
    }
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()
