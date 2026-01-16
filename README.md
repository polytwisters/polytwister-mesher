# Polytwister Mesher

Tool for creating 3D meshes of cross sections of [polytwisters](https://polytwisters.com/). Output is in Stanford PLY files as a triangle mesh with vertex normals.

This is a research codebase, so it might be a little janky.

## Usage

`polytwisters.json` contains the 4D geometric and combinatorial of all uniform polytwisters. It is checked in for convenience but you can generate it yourself from the polytwisters.com repo using `npm run export-geometry all polytwisters.json`.

```
cargo build --release
```

The executable is now at `./target/release/polytwister_mesher`.

### Merged mesh

The most convenient way to produce a mesh is to render a single cross section with rings, strips, and twisters merged together:

```
./target/release/polytwister_mesher polytwisters.json tetter -w 0.1 --merged out.ply
```

In place of "tetter" you can use any Bowers acronym (`gaquapiditer`) or full name (`"cube twister"` or `cube-twister` or `cube_twister`) or index (`34`) or symbol `3.3`.

The resulting mesh can be viewed in MeshLab.

### Individual meshes

For production rendering, you usually want to produce individual PLY meshes for rings, strips, and twisters:

```
./target/release/polytwister_mesher section polytwisters.json sadtadoditer \
    -w 0.1 \
    --rings rings.ply \
    --strips strips.ply \
    --twisters-1 twisters_1.ply \
    --twisters-2 twisters_2.ply
```

The two different twister files are for two different orbits of twisters. For regular twisters the second twister mesh will always be empty.

### Animation

Render evenly spaced W values to a directory of meshes:

```
./target/release/polytwister_mesher animation polytwisters.json sadtadoditer out_dir
```

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