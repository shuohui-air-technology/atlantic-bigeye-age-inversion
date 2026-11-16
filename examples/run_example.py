"""Replay the frozen monthly fit and the three synthetic operations."""

from pathlib import Path

from bet_age_inversion.cli import main


if __name__ == "__main__":
    main(example_dir=Path(__file__).resolve().parent / "synthetic")
