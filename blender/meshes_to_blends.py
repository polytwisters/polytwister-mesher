import argparse
import json
import pathlib

import common


def export_mesh_directory_as_blend(in_dir, out_file, config=None):
    args = [str(in_dir.resolve()), "-o", str(out_file)]
    if config is not None:
        args += ["-c", json.dumps(config)]
    common.run_blender_script(
        common.BLENDER_SCRIPT,
        blender_args=[],
        script_args=args,
        interactive=False,
    )


def export_mesh_directory_as_multiple_blends(in_dir, out_dir, config=None):
    out_dir.mkdir()
    mesh_dirs = [x for x in in_dir.iterdir() if x.is_dir()]
    for mesh_dir in mesh_dirs:
        out_file = out_dir / (mesh_dir.name + ".blend")
        export_mesh_directory_as_blend(mesh_dir, out_file, config=config)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_dir", type=str)
    parser.add_argument("out_blend", type=str)
    args = parser.parse_args()

    in_dir = pathlib.Path(args.in_dir)
    out_blend = pathlib.Path(args.out_blend)

    export_mesh_directory_as_multiple_blends(in_dir, out_blend)


if __name__ == "__main__":
    main()