"""Create an annotated GIF from the distributed movie binary fragments.

Run the solver from a separate working directory, then run this script from
that directory. It expects four movie folders produced by ``mpirun -np 4``.
The GIF is written directly; no PNG frames are retained on disk.
"""

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


SLICE_Y = 15
MPI_BOUNDARIES = (8, 16, 24)
SCALE = 12
MARGIN = 80
FRAME_COUNT = 80


def load_records(root):
    directories = sorted(root.glob("*_electric_field_movie_ME_*"))
    records = []
    for directory in directories:
        values = np.fromfile(directory / f"{directory.name}.bin", dtype="<f8")
        if values.size % 7:
            raise RuntimeError(f"Invalid movie binary layout: {directory}")
        records.append(values.reshape(-1, 7))
    if len(records) != 4:
        raise RuntimeError(f"Expected four MPI movie fragments, found {len(records)}")
    return np.concatenate(records)


def colourise(values, limit):
    normalised = np.clip(values / limit, -1.0, 1.0)
    image = np.empty(values.shape + (3,), dtype=np.uint8)
    positive = normalised >= 0.0
    image[..., 0] = np.where(positive, 255, 255 * (1.0 + normalised))
    image[..., 1] = 255 * (1.0 - np.abs(normalised))
    image[..., 2] = np.where(positive, 255 * (1.0 - normalised), 255)
    return image


def draw_legend(canvas, field_image, limit):
    draw = ImageDraw.Draw(canvas)
    gradient = colourise(np.linspace(limit, -limit, 256).reshape(256, 1), limit)
    bar = Image.fromarray(gradient, mode="RGB").resize((20, field_image.height), Image.Resampling.NEAREST)
    x = MARGIN + field_image.width + 18
    canvas.paste(bar, (x, MARGIN))
    draw.rectangle((x, MARGIN, x + 20, MARGIN + field_image.height), outline="black", width=1)
    draw.text((x, MARGIN - 20), "Ex (V/m)", fill="black")
    draw.text((x + 28, MARGIN), f"+{limit:.3e}", fill="black")
    draw.text((x + 28, MARGIN + field_image.height // 2), "0", fill="black")
    draw.text((x + 28, MARGIN + field_image.height - 10), f"-{limit:.3e}", fill="black")


def draw_axes_and_boundaries(canvas, field_image):
    draw = ImageDraw.Draw(canvas)
    left, top = MARGIN, MARGIN
    draw.rectangle((left, top, left + field_image.width, top + field_image.height), outline="black", width=1)
    for x_value in (0, 15, 29):
        x = left + x_value * SCALE
        draw.line((x, top + field_image.height, x, top + field_image.height + 5), fill="black", width=1)
        draw.text((x - 5, top + field_image.height + 8), str(x_value), fill="black")
    for z_value in (0, 8, 16, 24, 31):
        z = top + z_value * SCALE
        draw.line((left - 5, z, left, z), fill="black", width=1)
        draw.text((left - 25, z - 5), str(z_value), fill="black")
    for z_value in MPI_BOUNDARIES:
        z = top + z_value * SCALE
        for x in range(left, left + field_image.width, 8):
            draw.line((x, z, min(x + 4, left + field_image.width), z), fill="darkorange", width=2)
    draw.text((left, top - 20), "Orange dashed lines: MPI Z boundaries", fill="darkorange")
    draw.text((left, top + field_image.height + 32), "x cell index", fill="black")
    draw.text((10, top), "z", fill="black")


def make_frame(records, time, limit):
    values = records[(records[:, 0] == time) & (records[:, 2].astype(int) == SLICE_Y)]
    field = np.zeros((32, 30), dtype=float)
    field[values[:, 3].astype(int), values[:, 1].astype(int)] = values[:, 4]
    field_image = Image.fromarray(colourise(field, limit), mode="RGB").resize(
        (30 * SCALE, 32 * SCALE), Image.Resampling.NEAREST
    )
    canvas = Image.new("RGB", (field_image.width + 2 * MARGIN + 100, field_image.height + 2 * MARGIN), "white")
    canvas.paste(field_image, (MARGIN, MARGIN))
    ImageDraw.Draw(canvas).text((MARGIN, 20), f"Ex slice at y={SLICE_Y}, t={time * 1e9:.3f} ns", fill="black")
    draw_axes_and_boundaries(canvas, field_image)
    draw_legend(canvas, field_image, limit)
    return canvas


def main():
    root = Path.cwd()
    records = load_records(root)
    times = np.unique(records[:, 0])
    expected_coordinates = 30 * 30 * 32
    for time in times:
        coordinates = records[records[:, 0] == time, 1:4]
        if len(coordinates) != expected_coordinates or len(np.unique(coordinates, axis=0)) != expected_coordinates:
            raise RuntimeError(f"Incomplete or duplicated MPI coverage at t={time}")

    limit = np.max(np.abs(records[records[:, 2].astype(int) == SLICE_Y, 4]))
    frame_indices = np.linspace(0, len(times) - 1, FRAME_COUNT, dtype=int)
    frames = [make_frame(records, times[index], limit) for index in frame_indices]
    output_path = root / "mpi_layer_continuity.gif"
    frames[0].save(output_path, save_all=True, append_images=frames[1:], duration=120, loop=0)
    print(f"Created {output_path} from {len(times)} complete MPI frames")


if __name__ == "__main__":
    main()
