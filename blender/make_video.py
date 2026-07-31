import argparse
from pathlib import Path
import subprocess

from common import FFMPEG, IMAGEMAGICK_IDENTIFY


def make_webm(in_dir, out_file):
    fps = 24

    command = [
        FFMPEG,
        "-r", f"{fps}",
        "-i", in_dir / "section_%04d_0001.png",
        # Set video codec. VP9 supports transparency.
        "-c:v", "libvpx-vp9",
        "-b:v", "2M",
        # Do not overwrite.
        "-n",
        out_file
    ]
    subprocess.run(command, check=True)


def make_gif(in_dir, out_file):
    subprocess.run([
        FFMPEG,
        "-i", in_dir / "render_%04d.png",
        "-vf", "scale=400:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse",
        "-loop", "0",
        out_file
    ], check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_dir")
    parser.add_argument("out_file")
    args = parser.parse_args()

    in_dir = Path(args.in_dir)
    out_file = Path(args.out_file)

    if out_file.suffix == ".webm":
        make_webm(in_dir, out_file)
    elif out_file.suffix == ".gif":
        make_gif(in_dir, out_file)
    else:
        raise ValueError("Unrecognized output formate, use .webm or .gif")


if __name__ == "__main__":
    main()