"""Own batch files, SQLite aggregation and deterministic final exports."""

from contextlib import ExitStack
import json
import sqlite3
from .jsonl import dump


class BatchStore:
    """One output directory and its resources; no algorithm or batch policy."""

    def __init__(self, output):
        self.output = output

    def __enter__(self):
        self.output.mkdir(parents=True, exist_ok=False)
        self._resources = ExitStack()
        try:
            self.db = sqlite3.connect(self.output / "candidates.sqlite")
            self._resources.callback(self.db.close)
            self.streams = {
                name: self._resources.enter_context(
                    (self.output / f"{name}.jsonl").open("w")
                )
                for name in ("parents", "occurrences", "excluded", "errors")
            }
            self.db.execute(
                "CREATE TABLE candidates(identity_key TEXT PRIMARY KEY, data TEXT NOT NULL, occurrences INTEGER NOT NULL)"
            )
            self.db.execute(
                "CREATE TABLE sources(identity_key TEXT, input_line INTEGER, occurrence_index INTEGER, data TEXT, PRIMARY KEY(input_line, occurrence_index))"
            )
            self.db.execute("CREATE INDEX sources_identity ON sources(identity_key)")
            self.db.execute(
                "CREATE TABLE assignments(identity_key TEXT, family TEXT, classification TEXT, PRIMARY KEY(identity_key, family, classification))"
            )
        except BaseException:
            self._resources.close()
            raise
        return self

    def __exit__(self, exc_type, exc, traceback):
        return self._resources.__exit__(exc_type, exc, traceback)

    def _write(self, name, record):
        self.streams[name].write(dump(record) + "\n")

    def write_error(self, record):
        self._write("errors", record)

    def write_parent(self, result, input_line):
        """Persist every occurrence; identity deduplication never drops sources."""
        self._write(
            "parents",
            dict(
                {
                    k: v
                    for k, v in result.items()
                    if k not in ("occurrences", "candidates", "excluded")
                },
                input_line=input_line,
            ),
        )
        for failure in result["errors"]:
            self.write_error(
                dict(failure, input_line=input_line, input_id=result["id"])
            )
        for rejection in result["excluded"]:
            self._write(
                "excluded",
                dict(rejection, input_line=input_line, input_id=result["id"]),
            )
        for index, occurrence in enumerate(result["occurrences"]):
            item = dict(occurrence, input_line=input_line, occurrence_index=index)
            serialized = dump(item)
            self.streams["occurrences"].write(serialized + "\n")
            self.db.execute(
                "INSERT INTO sources VALUES(?,?,?,?)",
                (item["identity_key"], input_line, index, serialized),
            )
            self.db.execute(
                "INSERT OR IGNORE INTO assignments VALUES(?,?,?)",
                (
                    item["identity_key"],
                    item["frag_type"],
                    dump(item["classification_key"]),
                ),
            )
        for candidate in result["candidates"]:
            data = {
                k: v
                for k, v in candidate.items()
                if k not in ("occurrence_indices", "families")
            }
            self.db.execute(
                "INSERT INTO candidates VALUES(?,?,?) ON CONFLICT(identity_key) DO UPDATE SET occurrences=occurrences+excluded.occurrences",
                (
                    candidate["identity_key"],
                    dump(data),
                    len(candidate["occurrence_indices"]),
                ),
            )

    def commit(self):
        self.db.commit()

    def export(self):
        """Finish sorted candidate/group JSONL and return the unique count."""
        self.commit()
        unique = 0
        with (self.output / "candidates.jsonl").open("w") as stream:
            for key, data, count in self.db.execute(
                "SELECT * FROM candidates ORDER BY identity_key"
            ):
                candidate = json.loads(data)
                candidate["occurrence_count"] = count
                candidate["assignments"] = [
                    dict(family=f, classification_key=json.loads(g))
                    for f, g in self.db.execute(
                        "SELECT family,classification FROM assignments WHERE identity_key=? ORDER BY family,classification",
                        (key,),
                    )
                ]
                stream.write(dump(candidate) + "\n")
                unique += 1
        with (self.output / "groups.jsonl").open("w") as stream:
            for family, classification, count in self.db.execute(
                "SELECT family,classification,COUNT(DISTINCT identity_key) FROM assignments GROUP BY family,classification ORDER BY family,classification"
            ):
                stream.write(
                    dump(
                        dict(
                            family=family,
                            classification_key=json.loads(classification),
                            unique_identities=count,
                        )
                    )
                    + "\n"
                )
        return unique
