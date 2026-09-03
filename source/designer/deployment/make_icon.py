"""Generate a simple native-resolution application icon."""

from pathlib import Path

from PIL import Image, ImageDraw


size = 1024
image = Image.new("RGBA", (size, size), (0, 0, 0, 0))
draw = ImageDraw.Draw(image)
draw.rounded_rectangle((40, 40, size-40, size-40), radius=190,
                       fill=(23, 40, 59, 255))
for color, angle in (((85, 168, 255, 238), 0), ((255, 113, 141, 218), 11)):
    layer = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    for position in (300, 470, 640):
        ld.line((145, position, 879, position), fill=color, width=30)
        ld.line((position, 145, position, 879), fill=color, width=30)
    if angle:
        layer = layer.rotate(angle, resample=Image.Resampling.BICUBIC,
                             center=(size//2, size//2))
    image.alpha_composite(layer)
draw.ellipse((420, 420, 604, 604), fill=(156, 123, 217, 255),
             outline=(255, 255, 255, 255), width=22)
target = Path(__file__).with_name("moire-design-1024.png")
image.save(target)
print(target)
