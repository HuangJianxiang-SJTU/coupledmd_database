"""Regenerate final Figure 2 from packaged tables."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'revision'))
import figure2_revision as plot
if __name__ == '__main__':
    plot.OUT = Path(__file__).resolve().parents[1] / 'figures' / 'regenerated'
    plot.style()
    plot.figure2()
