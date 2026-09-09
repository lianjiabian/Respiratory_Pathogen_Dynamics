"""One-command validation of S1-S11 and redraw of Figures 1-4 and S1-S5."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from build_public_figure_inputs import PUBLIC_COLUMNS


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_renderer(path: Path):
    spec = importlib.util.spec_from_file_location("public_figure_renderer", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_checked(command: list[str]) -> None:
    subprocess.run(command, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[1]
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--rscript", type=Path, default=None)
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else root / "reproduced_outputs"
    )
    table_output = output_dir / "tables"
    figure_output = output_dir / "figures"
    qa_output = output_dir / "qa"
    qa_output.mkdir(parents=True, exist_ok=True)

    run_checked(
        [
            sys.executable,
            str(root / "scripts" / "run_public_reproduction.py"),
            "--root",
            str(root),
            "--output-dir",
            str(table_output),
        ]
    )
    if args.rscript:
        run_checked(
            [
                sys.executable,
                str(root / "scripts" / "validate_rda_release.py"),
                "--root",
                str(root),
                "--manifest",
                str(root / "release_manifest.csv"),
                "--output",
                str(qa_output / "rda_python_and_r_validation.csv"),
                "--rscript",
                str(args.rscript.resolve()),
            ]
        )

    input_manifest = pd.read_csv(
        root / "qa" / "figure_inputs_manifest.csv", encoding="utf-8-sig"
    )
    if set(input_manifest["file"]) != set(PUBLIC_COLUMNS):
        raise ValueError("Figure-input manifest does not match the public input contract")
    input_checks: list[dict[str, object]] = []
    for item in input_manifest.itertuples(index=False):
        path = root / "figure_inputs" / item.file
        frame = pd.read_csv(path, encoding="utf-8-sig")
        passed = (
            sha256(path) == str(item.sha256)
            and frame.shape == (int(item.rows), int(item.columns))
            and list(frame.columns) == list(PUBLIC_COLUMNS[str(item.file)])
            and set(frame["branch"].astype(str)) == {"minimal_corrected"}
        )
        input_checks.append(
            {
                "file": item.file,
                "status": "PASS" if passed else "FAIL",
                "sha256": sha256(path),
                "rows": len(frame),
                "columns": len(frame.columns),
                "privacy_status": item.privacy_status,
            }
        )
        if not passed:
            raise ValueError(f"Figure-input validation failed: {item.file}")
    pd.DataFrame(input_checks).to_csv(
        qa_output / "figure_input_validation.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )

    renderer = load_renderer(root / "scripts" / "07_render_figures.py")
    tables = renderer.load_inputs(
        root / "figure_inputs", requested_branches=("minimal_corrected",)
    )
    input_qa = renderer.validate_inputs(
        tables, requested_branches=("minimal_corrected",)
    )
    render_qa = renderer.render_all(
        tables,
        figure_output,
        branches=("minimal_corrected",),
        figures=renderer.FIGURE_ORDER,
    )
    rendered_files = sorted(path for path in figure_output.rglob("*") if path.is_file())
    rendered_manifest = pd.DataFrame(
        [
            {
                "file": path.relative_to(output_dir).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in rendered_files
        ]
    )
    rendered_manifest.to_csv(
        qa_output / "rendered_figure_sha256_manifest.csv",
        index=False,
        encoding="utf-8-sig",
        lineterminator="\n",
    )
    result = {
        "status": "PASS",
        "scope": "11 aggregate tables and nine reported figures",
        "analysis_branch": "minimal_corrected",
        "table_validation": str(
            (table_output / "public_table_validation.csv").relative_to(output_dir)
        ).replace("\\", "/"),
        "figure_input_validation": input_qa,
        "figure_render_validation": render_qa,
        "rendered_files": len(rendered_manifest),
        "expected_rendered_files": 9 * 4,
        "formats": ["SVG editable", "EPS editable", "PDF editable", "TIFF 600 dpi RGB LZW"],
        "privacy_boundary": (
            "No event-level clinical records are used. Fig. 3 necessarily retains "
            "all observed descriptive monthly rates and the aggregate pair counts "
            "displayed in the author-designated figure; "
            "public upload requires author/institutional approval."
        ),
    }
    if result["rendered_files"] != result["expected_rendered_files"]:
        raise ValueError("Expected exactly 36 rendered figure files")
    (qa_output / "public_reproduction_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
