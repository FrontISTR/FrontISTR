"""Compare two-rank shell integration-point output with an existing run.py result."""
import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from run import ROOT, read_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, default=ROOT/'build-shell-ep-mpi')
    parser.add_argument('--serial-run', type=Path, required=True)
    args = parser.parse_args()
    build = args.build.resolve()
    serial = args.serial_run.resolve()
    parent = build/'j2-validation'
    parent.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='mpi-', dir=parent))
    env = dict(os.environ, OMP_NUM_THREADS='1', OMPI_ALLOW_RUN_AS_ROOT='1',
               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
    summary = {'output': str(output), 'serial_run': str(serial), 'ranks': 2, 'checks': {}}

    def execute(command, directory, logname):
        completed = subprocess.run(command, cwd=directory, env=env,
                                   capture_output=True, text=True, timeout=120)
        text = completed.stdout+completed.stderr
        (directory/logname).write_text(text)
        assert completed.returncode == 0, text[-3000:]
        return text

    names = ('bend', 'mixed_layer_counts', 'separate_layer_counts',
             'mixed_elastic_output_on', 'mixed_elastic_layers',
             'mixed_elastic_strain_only', 'mixed_elastic_stress_only', 'elastic_output_on')
    for name in names:
        source = serial/name
        directory = output/name
        directory.mkdir()
        for filename in ('case.msh', 'case.cnt'):
            shutil.copy2(source/filename, directory/filename)
        (directory/'hecmw_ctrl.dat').write_text(
            '!MESH, NAME=fstrMSH, TYPE=HECMW-DIST\ncase.msh\n'
            '!CONTROL, NAME=fstrCNT\ncase.cnt\n!RESULT, NAME=fstrRES, IO=OUT\ncase.res\n'
            '!MESH, NAME=part_in, TYPE=HECMW-ENTIRE\ncase.msh\n'
            '!MESH, NAME=part_out, TYPE=HECMW-DIST\ncase.msh\n')
        (directory/'hecmw_part_ctrl.dat').write_text(
            '!PARTITION, TYPE=NODE-BASED, METHOD=RCB, DOMAIN=2, DEPTH=1\nx\n')
        execute(['mpiexec', '-n', '1', str(build/'hecmw1/tools/hecmw_part1')], directory, 'partition.log')
        log = execute(['mpiexec', '-n', '2', str(build/'fistr1/fistr1')], directory, 'solver.log')
        assert 'FrontISTR Completed !!' in log, log[-3000:]
        steps = sorted(int(path.name.rsplit('.', 1)[1]) for path in source.glob('case.res.0.*'))
        max_error = 0.0
        for step in steps:
            execute([str(build/'hecmw1/tools/rmerge'), '-n', '2', '-s', str(step),
                     '-e', str(step), 'case.res'], directory, f'merge-{step}.log')
            _, expected = read_result(source/f'case.res.0.{step}')
            _, actual = read_result(directory/f'case.res.{step}')
            assert set(expected) == set(actual)
            for ident, row in expected.items():
                assert set(row) == set(actual[ident])
                # Deliberately exclude nodal values: unequal-layer nodal averaging
                # is a separate known issue, not fixed by the output writer.
                for label, values in row.items():
                    assert len(values) == len(actual[ident][label])
                    error = max(abs(a-b) for a, b in zip(values, actual[ident][label]))
                    scale = max(map(abs, values))
                    assert error <= 1.e-9+1.e-8*scale, (name, step, ident, label, error, scale)
                    max_error = max(max_error, error/max(1.0, scale))
            rank_rows = [read_result(directory/f'case.res.{rank}.{step}')[1] for rank in range(2)]
            field_sets = [set(next(iter(rows.values()))) for rows in rank_rows]
            assert field_sets[0] == field_sets[1]
            if name == 'separate_layer_counts':
                assert all(len(rows) == 1 for rows in rank_rows), rank_rows
                assert set(rank_rows[0]).isdisjoint(rank_rows[1])
                assert 'GaussSTRESS16' in field_sets[0] and 'PLASTIC_GaussSTRAIN16' in field_sets[0]
        summary['checks'][name] = {'steps': steps, 'max_scaled_element_error': max_error}
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
