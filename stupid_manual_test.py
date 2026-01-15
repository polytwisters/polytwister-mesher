import subprocess
from pathlib import Path

subprocess.run([
    Path.home() / ".cargo/bin/cargo.exe", "run", "section", "polytwisters.json", "tetratwister", "--merged", "out.ply"
], check=True)

subprocess.run([
    r"C:\Program Files\VCG\MeshLab\meshlab.exe", "out.ply"
], check=True)