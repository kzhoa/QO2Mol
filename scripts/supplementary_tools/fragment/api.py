"""Fragmentation dispatch with explicitly selectable algorithms."""

DEFAULT_BACKEND = "f147"
BACKENDS = ("f147", "qq")


def process(
    record, families=None, include_parent=False, *, backend=DEFAULT_BACKEND, config=None
):
    if backend == "qq":
        from .qq import process as run

        return run(record, families, include_parent, config)
    if backend == "f147":
        from .modern.pipeline import process as run

        return run(record, families, include_parent, backend=backend, config=config)
    raise ValueError("Unknown fragment backend; choose f147 or qq.")
