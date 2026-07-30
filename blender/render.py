import argparse
import json
import pathlib
import subprocess

import common


def render_blend(in_file, out_png):
    subprocess.run(
        [
            common.BLENDER,
            "--background",
            str(in_file),
            # Order matters here! Set up the render format and output path first, then pick the frame.
            "--render-output", out_png,
            "--render-format", "PNG",
            "--render-frame", "1",
        ],
        check=True,
    )


def render_blends(in_dir, out_dir):
    out_dir.mkdir()
    for in_blend in in_dir.glob("*.blend"):
        out_png = out_dir / (in_blend.stem + ".png")
        render_blend(in_blend, out_png)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in", type=str)
    parser.add_argument("out", type=str)
    args = parser.parse_args()

    in_path = pathlib.Path(getattr(args, "in"))
    out_path = pathlib.Path(args.out)

    if in_path.is_file():
        render_blend(in_path, out_path)
    else:
        render_blends(in_path, out_path)


if __name__ == "__main__":
    main()
