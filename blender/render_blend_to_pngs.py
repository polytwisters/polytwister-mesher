import argparse
import json
import pathlib
import subprocess

import common


def render_blend(in_file, out_dir):
    subprocess.run(
        [
            common.BLENDER,
            "--background",
            str(in_file),
            # Order matters here! Set up the render format and output path first, then --render-anim.
            "--render-output",
            out_dir.resolve() / "render_####.png",
            "--render-format",
            "PNG",
            "--render-anim",
        ],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in", type=str)
    parser.add_argument("out", type=str)
    args = parser.parse_args()

    in_path = pathlib.Path(getattr(args, "in"))
    out_path = pathlib.Path(args.out)
    out_path.mkdir()

    render_blend(in_path, out_path)


if __name__ == "__main__":
    main()
