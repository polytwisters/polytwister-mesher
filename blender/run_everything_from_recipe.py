import subprocess
import argparse
import json
import pathlib

import common
import meshes_to_blends
import render

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str)
    parser.add_argument("out_dir", type=str)
    args = parser.parse_args()

    input_file = pathlib.Path(args.input_file)
    out_dir = pathlib.Path(args.out_dir)

    with open(input_file) as f:
        config_root = json.load(f)

    out_dir.mkdir(parents=True)
    
    polytwister = config_root["polytwister"]
    if isinstance(polytwister, dict):
        # Convex polytwister.
        polytwister = json.dumps(polytwister)

    w = config_root.get("w", None)
    frames = config_root.get("frames", None)
    is_animation = frames is not None

    mesher_config_file = out_dir / "mesher_config.json"
    with open(mesher_config_file, "w") as f:
        json.dump(config_root["mesher_config"], f)
    
    blender_config_file = out_dir / "blender_config.json"
    with open(blender_config_file, "w") as f:
        json.dump(config_root["blender_config"], f)
    
    mesh_dir = out_dir / "meshes"
    blend_dir = out_dir / "blends"

    if is_animation:
        subprocess.run([
            "cargo", "run", "--release",
            "--",
            "--config", mesher_config_file,
            polytwister,
            "-n", str(frames),
            mesh_dir,
        ], check=True)
    else:
        subprocess.run([
            "cargo", "run", "--release",
            "--",
            "--config", mesher_config_file,
            polytwister,
            "-w", str(w),
            mesh_dir,
        ], check=True)

    meshes_to_blends.export_mesh_directory_as_multiple_blends(
        mesh_dir,
        blend,
        config=blender_config_file
    )

    if w is not None:
        subprocess.run([
            common.BLENDER,
            "--background",
            str(blend),
            "--render-output",
            out_dir / "render",
            "--render-format", "PNG",
            "--render-frame", "1",
        ])
    else:
        raise ValueError("too lazy to fix this")

if __name__ == "__main__":
    main()