from pathlib import Path
import runpy
if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "revision" / "make_si_figure_equilibration.py"), run_name="__main__")
