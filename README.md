# Polytwister Mesher

Command-line tool for creating 3D meshes of cross sections of [polytwisters](https://polytwisters.com/).

## Usage

### Building

Polytwister Mesher is written in Rust, and its few dependencies are all cross-platform. To build it, install the Rust toolchain and run:

```
cargo build --release
```

The executable is now at `./target/release/polytwister_mesher` (file name `polytwister_mesher.exe` on Windows). I will abbreviate this as just `polytwister_mesher` from now on.

The executable reads the database file at `./polytwisters.json` which has geometry and naming information on all polytwisters. It is assumed that your working directory contains `polytwisters.json` (if it doesn't, you can supply the `-d` option to the executable to set a custom location).

### Single cross section, merged mesh

For quick viewing of polytwisters, you can render a single cross section showing rings, strips, and twisters as a single PLY file:

```
polytwister_mesher tetter -w 0.1 --merged out.ply
```

MeshLab is a good tool for quick previewing of such files, because you can just run `meshlab out.ply`. All meshes have vertex normals, so make sure to configure MeshLab to shade using them.

In place of "tetter" you can use any Bowers acronym (`gaquapiditer`) or full name (`"cube twister"` or `cube-twister` or `cube_twister`) or index (`34`) or symbol `3.3`.

### Single cross section, split meshes

For production rendering, you usually want to produce individual PLY meshes for rings, strips, and twisters:

```
polytwister_mesher section sadtadoditer -w 0.1 --split sadtadoditer_meshes/
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

**NOTE:** Any of the PLY files produced by Polytwister Mesher might be empty meshes, such as if the W coordinate in question happens to be an empty cross section. The PLY files are not empty files, but they have zero vertices. Blender produces an error when importing such PLY files, but it is harmless. This may be the case with other 3D software as well. (You don't need to worry about this with the Blender scripts discussed below as they check for empty meshes.)

### Animation

Render evenly spaced W values to a directory of meshes:

```
polytwister_mesher animation sadtadoditer -n 24 out_dir
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

### Configuration

A configuration JSON file can be passed in with the `-c` option, for example:

```
polytwister_mesher -c mesher_configs/high_resolution.json tetter -w 0.1 --merged out.ply
```

This allows controlling the resolution of the meshing and the radius of the strip and ring cross sections. Two presets are provided in the `mesher_configs/` directory.

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

**NOTE:** these Python scripts are meant more as as examples than as production tools, so you might want to modify them for your own pipeline.

## Limitations

* Only uniform polytwisters are supported.
* Currently the discretization of twister cross sections is rather poor and especially has problems with twisters with spiky cross sections. This can be mitigated by cranking up the mesh resolution, but at the cost of a larger mesh and more compute time.
* The above problem is especially bad if you decide to create a polytwister section without any visualization of ring or strip cross sections, because the rings and strips hide the cracks where the twisters meet. So currently, this mesher is just meant for ball-and-tube polytwister visualizations.
* Overall, this is a research codebase, so it's a little janky.

## Development

To easily preview the results, install MeshLab and create a file at the root of this directory called `local_manual_test.py` with the following.

```
import subprocess
from pathlib import Path

subprocess.run([
    "cargo", "run", "section", "tetratwister", "--merged", "out.ply"
], check=True)

subprocess.run([
    "meshlab", "out.ply"
], check=True)
```

### Database file

`polytwisters.json` contains the 4D geometric and combinatorial information of all uniform polytwisters. It is checked in for convenience but you can generate it yourself from the polytwisters.com repo using `npm run export-geometry all polytwisters.json`.

## License

[MIT](./LICENSE)