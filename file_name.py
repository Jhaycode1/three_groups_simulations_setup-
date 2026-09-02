from pathlib import Path
import re
import subprocess

project_dir = Path(__file__).resolve().parent
results_dir = project_dir / "test_run_results"


def next_save_name(results_path: Path, prefix: str = "sample_test_") -> str:
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+)(?:$|\.)")
    used_numbers = {
        int(match.group(1))
        for path in results_path.iterdir()
        if (match := pattern.match(path.name))
    }

    counter = 1
    while counter in used_numbers:
        counter += 1
    return f"{prefix}{counter}"


def run_simulation() -> None:
    save_name = next_save_name(results_dir)
    subprocess.run(
        ["Rscript", str(project_dir / "parent_script.R"), "1", save_name],
        cwd=project_dir,
        check=True,
    )


if __name__ == "__main__":
    run_simulation()


