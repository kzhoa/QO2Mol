"""Disk-backed aggregation of elementary products with complete provenance."""
from collections import Counter
import json
import platform
import sqlite3
import time
from rdkit import rdBase
from .models import QQConfig, FAMILIES, VERSION, TOOL_VERSION, implementation_hash
from .pipeline import process, validate_request


def generate(args):
    from ..cli import dump, digest, records
    values = json.loads(args.qq_config.read_text()) if args.qq_config else {}
    if not isinstance(values, dict):
        raise ValueError('QQ config must be a JSON object')
    try:
        config = QQConfig(**values)
    except TypeError as exc:
        raise ValueError('Removed/unknown QQ options; see README.md: '+str(exc)) from exc
    families = None if args.families is None else tuple(args.families.split(','))
    validate_request(families, args.include_parent, config)
    input_sha256 = digest(args.input)  # validate readable input before output creation
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    counts = Counter()
    db = sqlite3.connect(args.output/'candidates.sqlite')
    try:
        db.execute('CREATE TABLE candidates (identity_key TEXT PRIMARY KEY, data TEXT NOT NULL, occurrences INTEGER NOT NULL)')
        db.execute('CREATE TABLE sources (identity_key TEXT, input_line INTEGER, occurrence_index INTEGER, data TEXT, PRIMARY KEY(input_line,occurrence_index))')
        db.execute('CREATE INDEX sources_identity ON sources(identity_key)')
        with (args.output/'parents.jsonl').open('w') as parents, (args.output/'occurrences.jsonl').open('w') as occurrences, (args.output/'excluded.jsonl').open('w') as excluded, (args.output/'errors.jsonl').open('w') as errors:
            for line, record, error in records(args.input):
                counts['input_records'] += 1
                if error:
                    counts['input_failures'] += 1
                    errors.write(dump(dict(line=line, stage='input', error=error))+'\n')
                    continue
                try:
                    result = process(record, families, False, config)
                except Exception as exc:
                    counts['input_failures'] += 1
                    errors.write(dump(dict(line=line, stage='parent', input_id=record.get('id'), error=f'{type(exc).__name__}: {exc}'))+'\n')
                    continue
                counts['processed_parents'] += 1
                counts['elementary_blocks'] += len(result['elementary_fragments'])
                parents.write(dump({k:v for k,v in result.items() if k not in ('candidates','occurrences','excluded')})+'\n')
                for failure in result['errors']:
                    counts['stage_errors'] += 1
                    errors.write(dump(dict(failure, input_line=line, input_id=result['id']))+'\n')
                for rejection in result['excluded']:
                    counts['excluded'] += 1
                    excluded.write(dump(dict(rejection, input_line=line, input_id=result['id']))+'\n')
                for candidate in result['candidates']:
                    data = {k:v for k,v in candidate.items() if k not in ('occurrence_indices','assembly_keys')}
                    db.execute('INSERT INTO candidates VALUES(?,?,?) ON CONFLICT(identity_key) DO UPDATE SET occurrences=occurrences+excluded.occurrences',
                        (candidate['identity_key'],dump(data),len(candidate['occurrence_indices'])))
                for i, occurrence in enumerate(result['occurrences']):
                    counts['occurrences'] += 1
                    item = dict(occurrence,input_line=line,occurrence_index=i)
                    occurrences.write(dump(item)+'\n')
                    db.execute('INSERT INTO sources VALUES(?,?,?,?)',(occurrence['candidate_key'],line,i,dump(item)))
                if counts['processed_parents'] % 100 == 0:
                    db.commit()
        db.commit()
        with (args.output/'candidates.jsonl').open('w') as stream:
            for key, data, count in db.execute('SELECT * FROM candidates ORDER BY identity_key'):
                candidate = json.loads(data)
                candidate['occurrence_count'] = count
                candidate['assembly_keys'] = sorted({json.loads(origin)['assembly']['module_key'] for (origin,) in db.execute('SELECT data FROM sources WHERE identity_key=?',(key,))})
                stream.write(dump(candidate)+'\n')
                counts['unique_identities'] += 1
        manifest = dict(backend='qq', profile=VERSION, tool_version=TOOL_VERSION,
            implementation_sha256=implementation_hash(),config=config.to_dict(),config_sha256=config.fingerprint,
            families=list(FAMILIES), include_parent=False, input_sha256=input_sha256,
            python=platform.python_version(),rdkit=rdBase.rdkitVersion,counts=dict(counts),
            selection='all_successful_strict_identities',elapsed_seconds=time.monotonic()-start,
            coordinates_generated=False,historical_reconstruction=False)
        (args.output/'run.json').write_text(dump(manifest)+'\n')
        print(json.dumps(dict(counts),indent=2))
        return 2 if (counts['input_failures'] or counts['stage_errors']) and not args.allow_errors else 0
    finally:
        db.close()
