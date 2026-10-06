"""Compare a loaded-shell restart followed by unloading with uninterrupted runs."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from run import CASES, ROOT, integration_output, read_result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', default='build-shell-ep')
    parser.add_argument('--bend-unload-substeps', type=int, default=40)
    args = parser.parse_args()
    build = (ROOT / args.build).resolve()
    parent = build / 'j2-validation'
    parent.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='restart-', dir=parent))
    summary = {'output': str(output), 'checks': {}}

    def run(name, mesh, control, restart=None):
        directory = output / name
        directory.mkdir()
        shutil.copyfile(CASES / (mesh + '.msh'), directory / 'case.msh')
        (directory / 'case.cnt').write_text(integration_output(control))
        (directory / 'hecmw_ctrl.dat').write_text(
            '!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE\ncase.msh\n'
            '!CONTROL, NAME=fstrCNT\ncase.cnt\n!RESULT, NAME=fstrRES, IO=OUT\ncase.res\n'
            '!RESTART, NAME=restart_out, IO=INOUT\ncase.restart\n')
        if restart:
            files = list(restart.glob('case.restart*'))
            assert files, restart
            for file in files:
                shutil.copyfile(file, directory / file.name)
        proc = subprocess.run([str(build / 'fistr1/fistr1')], cwd=directory,
                              capture_output=True, text=True, timeout=120)
        log = proc.stdout + proc.stderr
        (directory / 'solver.log').write_text(log)
        assert proc.returncode == 0 and 'FrontISTR Completed !!' in log, log[-4000:]
        return directory

    for kind in ('MEMBRANE', 'BEND'):
        control = (CASES / f'K741J2{kind}.cnt').read_text()
        if kind == 'MEMBRANE':
            last_step = 14
            first_step = '!STEP, SUBSTEPS=12, MAXITER=30, CONVERG=1.0E-8\n  BOUNDARY, 1\n'
            marker = '!STEP, SUBSTEPS=2, MAXITER=30, CONVERG=1.0E-8\n  BOUNDARY, 2\n'
            loaded_control = control.replace(marker, '')
            assert loaded_control != control
        else:
            last_step = 12 + args.bend_unload_substeps
            first_step = ('!STEP, SUBSTEPS=12, MAXITER=40, CONVERG=1.0E-8\n'
                          '  BOUNDARY, 1\n  LOAD, 1\n')
            loaded_control = control
            control = control.replace('!MATERIAL, NAME=M1',
                '!CLOAD, GRPID=2\n  LOAD, 3, 0.04\n'
                f'!STEP, SUBSTEPS={args.bend_unload_substeps}, MAXITER=40, CONVERG=1.0E-8\n'
                '  BOUNDARY, 1\n  LOAD, 1\n  LOAD, 2\n!MATERIAL, NAME=M1')
        reference = run(kind.lower() + '_full', 'K741J2' + kind, control)
        loaded = run(kind.lower() + '_loaded', 'K741J2' + kind,
                     loaded_control.replace('!STATIC', '!RESTART, FREQUENCY=12, VERSION=6\n!STATIC'))
        assert first_step in control
        # A completed-step restart appends the steps in the new control file.
        resume_control = control.replace(first_step, '')
        resumed = run(kind.lower() + '_resumed', 'K741J2' + kind,
                      resume_control.replace('!STATIC', '!RESTART, FREQUENCY=-12, VERSION=6\n!STATIC'), loaded)
        errors = {}
        for old_block, new_block in zip(read_result(reference / f'case.res.0.{last_step}'),
                                       read_result(resumed / f'case.res.0.{last_step}')):
            assert old_block.keys() == new_block.keys()
            for ident, old_row in old_block.items():
                new_row = new_block[ident]
                assert old_row.keys() == new_row.keys()
                for label, values in old_row.items():
                    assert len(values) == len(new_row[label])
                    error = max(abs(a-b)/max(1.0, abs(a)) for a, b in zip(values, new_row[label]))
                    errors[label] = max(errors.get(label, 0.0), error)
                    assert error < 1.e-9, (kind, label, error)
        summary['checks'][kind] = {'max_scaled_error': max(errors.values()), 'fields': errors}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
