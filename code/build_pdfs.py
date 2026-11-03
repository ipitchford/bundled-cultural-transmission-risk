"""Build the paper PDF from manuscript.md (requires pandoc and pdfLaTeX).

Usage: python3 -I -B code/build_pdfs.py [--out DIR]
The LaTeX intermediate is written to tmp/pdf/ (not distributed).
"""
import argparse, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def run(args, cwd):
    r = subprocess.run(args, cwd=cwd, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        sys.stdout.write(r.stdout[-3000:]); sys.stderr.write(r.stderr[-3000:])
        sys.exit(f"failed: {args[0]}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "output"))
    out = Path(ap.parse_args().out); out.mkdir(parents=True, exist_ok=True)
    for tool in ("pandoc", "pdflatex"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found on PATH")
    build = ROOT / "tmp" / "pdf"; build.mkdir(parents=True, exist_ok=True)
    header = build / "header.tex"
    header.write_text("\\usepackage{float}\\floatplacement{figure}{htbp}\n")
    run(["pandoc", str(ROOT / "manuscript.md"), "--standalone", "--from=markdown+tex_math_dollars",
         "-H", str(header), "-V", "fontsize=11pt", "-V", "papersize=a4", "-V", "geometry:margin=24mm",
         "-V", "colorlinks=true", "-V", "linkcolor=teal", "-V", "urlcolor=teal",
         "-V", "documentclass=article", "-o", str(build / "paper.tex")], ROOT)
    for _ in range(2):
        run(["pdflatex", "-interaction", "nonstopmode", "-halt-on-error",
             "-output-directory", str(build), str(build / "paper.tex")], ROOT)
    dest = out / "cultural_transmission_paper.pdf"
    shutil.copyfile(build / "paper.pdf", dest)
    print(dest)

if __name__ == "__main__":
    main()
