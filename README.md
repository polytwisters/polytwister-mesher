# Polytwister Mesher

Tool for creating 3D meshes of cross sections of [polytwisters](https://polytwisters.com/). Output is in Stanford PLY files as a triangle mesh with vertex normals.

This is a research codebase, so it might be a little janky.

## Usage

`polytwisters.json` contains the 4D geometric and combinatorial of all uniform polytwisters. It is checked in for convenience but you can generate it yourself from the polytwisters.com repo using `npm run export-geometry all polytwisters.json`.

```
cargo build --release
```

The executable is now at `./target/release/polytwister_mesher`.

### Single cross section, merged mesh

Render a single cross section showing rings, strips, and twisters as a single PLY file:

```
./target/release/polytwister_mesher polytwisters.json tetter -w 0.1 --merged out.ply
```

In place of "tetter" you can use any Bowers acronym (`gaquapiditer`) or full name (`"cube twister"` or `cube-twister` or `cube_twister`) or index (`34`) or symbol `3.3`.

MeshLab is a good way to quickly preview and inspect such meshes.

### Single cross section, split meshes

For production rendering, you usually want to produce individual PLY meshes for rings, strips, and twisters:

```
./target/release/polytwister_mesher section polytwisters.json sadtadoditer \
    -w 0.1 \
    --split sadtadoditer_meshes/
```

This command will create the directory `sadtadoditer_meshes/` and the following four files:

```
sadtadoditer_meshes/rings.ply
sadtadoditer_meshes/strips.ply
sadtadoditer_meshes/twisters_1.ply
sadtadoditer_meshes/twisters_2.ply
```

`twisters_1` and `twisters_2` are the two orbits of the twisters. If the polytwister is regular, then it has only one twister orbit and `twisters_2` is an empty mesh.

The directory `sadtadoditer_meshes` is referred to as a "section directory."

### Animation

Render evenly spaced W values to a directory of meshes:

```
./target/release/polytwister_mesher animation polytwisters.json sadtadoditer -n 24 out_dir
```

This creates `out_dir` and the following structure:

```
out_dir/
    section_0000/
        rings.ply
        strips.ply
        twisters_1.ply
        twisters_2.ply
    section_0001/
        rings.ply
        strips.ply
        twisters_1.ply
        twisters_2.ply
    ...
```

where there is one subdirectory for each frame. Each subdirectory is a section directory. The directory `out_dir` is referred to as an animation directory.

## Rendering with Blender

A Blender Python script at `blender/blender_script.py` can take either a section directory or an animation directory and load it as a single Blender file. To invoke this script and view the resulting Blender file interactively, you can run

```
python3 blender/blender_wrapper.py <directory>
```

If the directory is a section directory, it imports the PLY files and sets up a camera, materials, and render settings for you. If the directory is an animation directory, it imports all the section directories and creates an animation out of them.

If you want to produce a .blend file non-interactively:

```
python3 blender/meshes_to_blend.py <directory> polytwister.blend
```

If the Blender file is an animation, the following scripts are also provided for rendering:

```
python3 blender/render_blend_pngs.py polytwister.blend out_pngs_dir/
python3 blender/make_video.py out_pngs_dir/ out.mp4  # requires ffmpeg
```

**NOTE:** these Python scripts are a bit janky, so you might want to modify them for your use.

## Development

To easily preview the results, install MeshLab and create a file at the root of this directory called `local_manual_test.py` with the following.

```
import subprocess
from pathlib import Path

subprocess.run([
    "cargo", "run", "section", "polytwisters.json", "tetratwister", "--merged", "out.ply"
], check=True)

subprocess.run([
    "meshlab", "out.ply"
], check=True)
```