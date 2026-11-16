"""Direct script and module entry point for the fragment tools."""

if __package__:
    from .cli import main
else:
    from pathlib import Path
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from fragment.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
