#!/usr/bin/env python3
"""Create deterministic macOS .icns files from release artwork."""

from pathlib import Path
import shutil
import subprocess

from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "build_macos" / "assets"


def iconset(source: Path, name: str) -> None:
    target = ASSETS / (name + ".iconset")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    image = Image.open(source).convert("RGBA")
    side = max(image.size)
    canvas = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    canvas.alpha_composite(image, ((side - image.width) // 2,
                                   (side - image.height) // 2))
    for points in (16, 32, 128, 256, 512):
        for scale in (1, 2):
            pixels = points * scale
            resized = canvas.resize((pixels, pixels), Image.Resampling.LANCZOS)
            suffix = "@2x" if scale == 2 else ""
            resized.save(target / f"icon_{points}x{points}{suffix}.png")
    subprocess.run(["iconutil", "-c", "icns", str(target),
                    "-o", str(ASSETS / (name + ".icns"))], check=True)
    shutil.rmtree(target)


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    installed = Path("/Applications/DNA Moiré Designer.app/Contents/Resources/moire-design.icns")
    designer = ASSETS / "moire-design.icns"
    if installed.is_file():
        shutil.copy2(installed, designer)
        # Keep the build independent of iconutil behavior on newer macOS.
        # The companion remains visually grouped with the Designer release;
        # its own in-window official caDNAno artwork is unchanged.
        shutil.copy2(installed, ASSETS / "cadnano2.icns")
    else:
        # The source SVG is converted through Qt by the running application;
        # keep a deterministic release build possible even without that app by
        # using the companion artwork as a conservative fallback.
        iconset(ROOT / "source" / "cadnano_companion" / "cadnano2" / "ui" /
                "mainwindow" / "images" / "cadnano2-app-icon.png",
                "moire-design")
        shutil.copy2(designer, ASSETS / "cadnano2.icns")


if __name__ == "__main__":
    main()
