"""Create deterministic multi-resolution Windows icons."""

from pathlib import Path
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
ASSETS.mkdir(parents=True, exist_ok=True)
SIZES = [(16, 16), (24, 24), (32, 32), (48, 48), (64, 64),
         (128, 128), (256, 256)]


def designer_icon():
    size = 1024
    image = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle((40, 40, 984, 984), radius=190,
                           fill=(23, 40, 59, 255))
    for color, angle in (((85, 168, 255, 238), 0),
                         ((255, 113, 141, 218), 11)):
        layer = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        layer_draw = ImageDraw.Draw(layer)
        for position in (300, 470, 640):
            layer_draw.line((145, position, 879, position),
                            fill=color, width=30)
            layer_draw.line((position, 145, position, 879),
                            fill=color, width=30)
        if angle:
            layer = layer.rotate(angle, Image.Resampling.BICUBIC,
                                 center=(size // 2, size // 2))
        image.alpha_composite(layer)
    draw.ellipse((420, 420, 604, 604), fill=(198, 181, 232, 255),
                 outline=(255, 255, 255, 255), width=22)
    return image


def cadnano_icon():
    source = (ROOT.parent / "source" / "cadnano_companion" / "cadnano2" /
              "ui" / "mainwindow" / "images" / "cadnano2-app-icon.png")
    return Image.open(source).convert("RGBA")


for name, image in (("moire-designer.ico", designer_icon()),
                    ("cadnano2.ico", cadnano_icon())):
    image.save(ASSETS / name, format="ICO", sizes=SIZES)
    print(ASSETS / name)
