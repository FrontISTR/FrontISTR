"""Independent checks for a GNU/Linux serial build; never overwrite baselines."""
import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CASES = HERE.parent


def read_result(path):
    tokens = iter(path.read_text().split('*data\n', 1)[1].split())
    counts = [int(next(tokens)), int(next(tokens))]
    fields = [int(next(tokens)), int(next(tokens))]
    blocks = []
    for count, nfield in zip(counts, fields):
        sizes = [int(next(tokens)) for _ in range(nfield)]
        labels = [next(tokens) for _ in range(nfield)]
        rows = {}
        for _ in range(count):
            ident = int(next(tokens))
            rows[ident] = {label: [float(next(tokens)) for _ in range(size)]
                           for label, size in zip(labels, sizes)}
        blocks.append(rows)
    return blocks


def close(actual, expected, tol=1.e-8):
    assert abs(actual-expected) < tol, (actual, expected, tol)


def check_gauss_plastic_strain(row, nlayer, output_points=None):
    # Preserve the history-point mean; physical-volume averaging is separate.
    npoint = 4*nlayer*2
    labels = {name for name in row if name.startswith('PLASTIC_GaussSTRAIN')}
    assert labels == {f'PLASTIC_GaussSTRAIN{k}' for k in range(1, (output_points or npoint)+1)}
    average = 0.0
    for layer in range(nlayer):
        for thick, side in enumerate(('-', '+')):
            values = [row[f'PLASTIC_GaussSTRAIN{(ig*nlayer+layer)*2+thick+1}'][0]
                      for ig in range(4)]
            assert min(values) >= 0
            surface_average = sum(values)/4
            close(surface_average, row[f'ElementalPLSTRAIN_L{layer+1}{side}'][0], 1.e-11)
            average += surface_average/(2*nlayer)
    close(average, row['ElementalPLSTRAIN'][0], 1.e-11)
    for k in range(npoint+1, (output_points or npoint)+1):
        close(row[f'PLASTIC_GaussSTRAIN{k}'][0], 0.0, 1.e-15)


def integration_output(control):
    return control.replace('!OUTPUT_RES', '!OUTPUT_RES\n  ISTRAIN, ON\n  ISTRESS, ON')


def check_gauss_tensors(row, nlayer, output_points=None):
    # Existing thickness-rule average, without new layer or geometric weights.
    npoint = 8*nlayer
    for quantity in ('STRAIN', 'STRESS'):
        prefix = f'Gauss{quantity}'
        labels = {name for name in row if name.startswith(prefix)}
        assert labels == {f'{prefix}{k}' for k in range(1, (output_points or npoint)+1)}
        average = [0.0]*6
        for layer in range(nlayer):
            for thick, side in enumerate(('-', '+')):
                values = [row[f'{prefix}{(ig*nlayer+layer)*2+thick+1}'] for ig in range(4)]
                assert all(len(v) == 6 for v in values)
                for component in range(6):
                    surface_mean = sum(v[component] for v in values)/4
                    close(surface_mean, row[f'Elemental{quantity}_L{layer+1}{side}'][component])
                    average[component] += surface_mean/(2*nlayer)
        for actual, expected in zip(average, row[f'Elemental{quantity}']):
            close(actual, expected)
        for k in range(npoint+1, (output_points or npoint)+1):
            assert row[f'{prefix}{k}'] == [0.0]*6


def rotate_mesh(text):
    # Rotation about y: local normal becomes (0.6, 0, 0.8).
    output, nodes = [], False
    for line in text.splitlines():
        if line.startswith('!'):
            nodes = line.upper() == '!NODE'
        elif nodes and line.strip():
            ident, x, y, z = map(float, line.split(','))
            line = f'{int(ident)}, {0.8*x+0.6*z:.16g}, {y:.16g}, {-0.6*x+0.8*z:.16g}'
        output.append(line)
    return '\n'.join(output)+'\n'


def rotate_boundary(text):
    output, boundary = [], False
    for line in text.splitlines():
        if line.startswith('!'):
            boundary = line.upper().startswith('!BOUNDARY')
        elif boundary and line.strip():
            group, first, last, value = [v.strip() for v in line.split(',')]
            if group == 'ALLNODES':
                output.extend(['  ALLNODES, 2, 2, 0.0', '  ALLNODES, 4, 6, 0.0'])
                continue
            assert first == last == '1'
            output.extend([f'  {group}, 1, 1, {0.8*float(value):.16g}',
                           f'  {group}, 3, 3, {-0.6*float(value):.16g}'])
            continue
        output.append(line)
    return '\n'.join(output)+'\n'


def tensor_component(v, direction, engineering=False):
    x, y, z = direction
    factor = 1 if engineering else 2
    return (x*x*v[0]+y*y*v[1]+z*z*v[2]
            +factor*(x*y*v[3]+y*z*v[4]+z*x*v[5]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, default=ROOT/'build-shell-ep')
    args = parser.parse_args()
    build = args.build.resolve()
    parent = build/'j2-validation'
    parent.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='run-', dir=parent))
    summary = {'output': str(output), 'checks': {}}

    def run(name, case, cnt=None, mesh=None, expected_error=None):
        directory = output/name
        directory.mkdir()
        (directory/'case.cnt').write_text(cnt or (CASES/f'{case}.cnt').read_text())
        (directory/'case.msh').write_text(mesh or (CASES/f'{case}.msh').read_text())
        (directory/'hecmw_ctrl.dat').write_text(
            '!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE\ncase.msh\n'
            '!CONTROL, NAME=fstrCNT\ncase.cnt\n!RESULT, NAME=fstrRES, IO=OUT\ncase.res\n')
        completed = subprocess.run([str(build/'fistr1/fistr1')], cwd=directory,
                                   capture_output=True, text=True, timeout=120)
        log = completed.stdout+completed.stderr
        (directory/'solver.log').write_text(log)
        if expected_error:
            assert expected_error in log and 'FrontISTR Completed !!' not in log, log[-3000:]
            summary['checks'][name] = {'rejected': expected_error}
        else:
            assert completed.returncode == 0 and 'FrontISTR Completed !!' in log, log[-3000:]
            summary['checks'][name] = {
                'completed': True,
                'max_iterations': max(map(int, re.findall(r'iter:\s+(\d+)', log)), default=0)}
        return directory

    for name in ('check_j2_material', 'check_j2_output', 'check_j2_potential', 'check_shell_output'):
        executable = output/name
        subprocess.run(['gfortran', '-fcheck=all', '-I'+str(build/'fistr1'),
                        '-I'+str(build/'hecmw1'), str(HERE/f'{name}.f90'),
                        str(build/'fistr1/libfistr.a'), str(build/'hecmw1/libhecmw.a'),
                        '-llapack', '-lblas', '-lstdc++', '-o', str(executable)], check=True)
        proc = subprocess.run([str(executable)], capture_output=True, text=True, timeout=60)
        (output/f'{name}.log').write_text(proc.stdout+proc.stderr)
        # Fortran STOP with a character message can return status zero.
        assert proc.returncode == 0 and 'PASS ' in proc.stdout, proc.stdout+proc.stderr
        summary['checks'][name] = proc.stdout.strip()

    membrane = run('membrane', 'K741J2MEMBRANE',
                   integration_output((CASES/'K741J2MEMBRANE.cnt').read_text()))
    nodes, elements = read_result(membrane/'case.res.0.14')
    _, loaded = read_result(membrane/'case.res.0.12')
    e = elements[1]
    gap = max(abs(n['NodalSTRESS'][0]-e['ElementalSTRESS'][0]) for n in nodes.values())
    close(gap, 0)
    close(e['ElementalPLSTRAIN'][0], loaded[1]['ElementalPLSTRAIN'][0])
    close(e['ElementalSTRESS'][0], loaded[1]['ElementalSTRESS'][0]-210000/(1-.3**2)*.002)
    close(sum(e['ElementalSTRAIN'][:3]), (1-2*.3)/210000*sum(e['ElementalSTRESS'][:3]), 1.e-11)
    summary['checks']['membrane'].update(nodal_element_gap=gap, thickness_strain=e['ElementalSTRAIN'][2])
    for step in (0, 12, 14):
        _, rows = read_result(membrane/f'case.res.0.{step}')
        check_gauss_plastic_strain(rows[1], 1)
        check_gauss_tensors(rows[1], 1)
        for k in range(1, 9):
            for quantity in ('STRAIN', 'STRESS'):
                for actual, expected in zip(rows[1][f'Gauss{quantity}{k}'], rows[1][f'Elemental{quantity}']):
                    close(actual, expected)
        for k in range(1, 9):
            close(rows[1][f'PLASTIC_GaussSTRAIN{k}'][0], rows[1]['ElementalPLSTRAIN'][0], 1.e-11)
    assert e['PLASTIC_GaussSTRAIN1'][0] > 0
    summary['checks']['membrane']['plastic_gauss_points'] = 8

    cnt = integration_output((CASES/'K741J2MEMBRANE.cnt').read_text())
    cnt = cnt.replace('!OUTPUT_RES', '!OUTPUT_RES\n  NSTRAIN, ON')
    mesh = (CASES/'K741J2MEMBRANE.msh').read_text()
    for scale in (1.e-12, 1.e-6, 1.e6, 1.e12):
        elastic = f'{210000*scale:.16e}, 0.3'
        scaled_cnt = cnt.replace('210000.0, 0.3', elastic).replace('1000.0, 0.0', f'{1000*scale:.16e}, 0.0')
        scaled_mesh = mesh.replace('210000.0, 0.3', elastic)
        name = f'stress_units_{scale:g}'
        scaled = run(name, 'K741J2MEMBRANE', scaled_cnt, scaled_mesh)
        max_error = 0.0
        for step in (0, 12, 14):
            base_blocks = read_result(membrane/f'case.res.0.{step}')
            scaled_blocks = read_result(scaled/f'case.res.0.{step}')
            for base_rows, scaled_rows in zip(base_blocks, scaled_blocks):
                for ident, row in base_rows.items():
                    for label, values in row.items():
                        divisor = scale if 'STRESS' in label or 'MISES' in label else 1.0
                        normalized = [v/divisor for v in scaled_rows[ident][label]]
                        relative_error = max(abs(a-b) for a, b in zip(normalized, values))/max(1.0, *map(abs, values))
                        close(relative_error, 0.0, 1.e-9)
                        max_error = max(max_error, relative_error)
        summary['checks'][name]['max_normalized_error'] = max_error

    rotated = run('rotated_membrane', 'K741J2MEMBRANE', rotate_boundary(cnt), rotate_mesh(mesh))
    rnodes, relems = read_result(rotated/'case.res.0.14')
    normal = (0.6, 0, 0.8)
    close(tensor_component(relems[1]['ElementalSTRESS'], normal), 0)
    eps33 = tensor_component(relems[1]['ElementalSTRAIN'], normal, engineering=True)
    close(eps33, e['ElementalSTRAIN'][2], 1.e-11)
    for row in rnodes.values():
        close(tensor_component(row['NodalSTRESS'], normal), 0)
        close(tensor_component(row['NodalSTRAIN'], normal, engineering=True), eps33, 1.e-11)
    # All layer +/- results must use the same full tensor transformation.
    for name, values in relems[1].items():
        if name.startswith(('ElementalSTRESS', 'GaussSTRESS')):
            close(tensor_component(values, normal), 0)
        elif name.startswith(('ElementalSTRAIN', 'GaussSTRAIN')):
            close(tensor_component(values, normal, engineering=True), eps33, 1.e-11)
    summary['checks']['rotated_membrane']['local_thickness_strain'] = eps33

    bend_cnt = integration_output((CASES/'K741J2BEND.cnt').read_text())
    bend = run('bend', 'K741J2BEND', bend_cnt)
    bnodes, belems = read_result(bend/'case.res.0.12')
    close(sum(bnodes[i]['REACTION_FORCE'][2] for i in (1,4)), .17)
    assert max(r['ElementalPLSTRAIN'][0] for r in belems.values()) > 0
    assert summary['checks']['bend']['max_iterations'] > 1
    for row in belems.values():
        check_gauss_plastic_strain(row, 1)
        check_gauss_tensors(row, 1)

    auto_cnt = bend_cnt.replace(
        '!STEP, SUBSTEPS=12, MAXITER=40, CONVERG=1.0E-8',
        '!STEP, SUBSTEPS=1000, MAXITER=4, CONVERG=1.0E-8, INC_TYPE=AUTO\n'
        '  0.08333333333333333, 1.0, 0.00001, 0.08333333333333333')
    automatic = run('bend_cutback', 'K741J2BEND', auto_cnt)
    auto_log = (automatic/'solver.log').read_text()
    assert 'Fail to Converge' in auto_log, 'Cutback was not exercised'
    final = max(automatic.glob('case.res.0.*'), key=lambda p: int(p.name.rsplit('.', 1)[1]))
    anodes, aelems = read_result(final)
    close(sum(anodes[i]['REACTION_FORCE'][2] for i in (1, 4)), .17)
    assert max(row['ElementalPLSTRAIN'][0] for row in aelems.values()) > 0
    summary['checks']['bend_cutback'].update(
        failed_attempts=auto_log.count('Fail to Converge'), final_result=final.name,
        support_reaction_z=sum(anodes[i]['REACTION_FORCE'][2] for i in (1, 4)))

    bend_mesh = (CASES/'K741J2BEND.msh').read_text()
    quasi = run('bend_quasinewton', 'K741J2BEND', bend_cnt.replace(
        '!STATIC', '!NONLINEAR_SOLVER, METHOD=QUASINEWTON\n!STATIC'))
    qnodes, qelems = read_result(quasi/'case.res.0.12')
    differences = {}
    for expected, actual, labels in ((bnodes, qnodes, ('DISPLACEMENT', 'ROTATION')),
                                     (belems, qelems, tuple(belems[1]))):
        for label in labels:
            scale = max(abs(v) for row in expected.values() for v in row[label])
            error = max(abs(a-b) for ident, row in expected.items()
                        for a, b in zip(row[label], actual[ident][label]))
            assert error <= 1.e-9+5.e-4*scale, (label, error, scale)
            differences[label] = {'max_difference': error, 'reference_scale': scale}
    summary['checks']['bend_quasinewton']['newton_comparison'] = differences
    # The existing Quasi-Newton driver does not update the SPC reaction output.
    summary['checks']['bend_quasinewton']['reported_support_reaction_z'] = sum(
        qnodes[i]['REACTION_FORCE'][2] for i in (1, 4))
    layered_cnt = bend_cnt.replace('!CLOAD, GRPID=1\n  LOAD, 3, -0.085', '  LOAD, 3, 3, -0.05')
    layered_cnt = layered_cnt.replace('  LOAD, 1\n', '')
    layered_cnt = layered_cnt.replace('!SOLVER,',
        '!MATERIAL, NAME=M2\n!ELASTIC, INFINITESIMAL\n  40000.0, 0.3\n'
        '!PLASTIC, INFINITESIMAL\n  20.0, 0.0\n!SOLVER,')
    layered_mesh = bend_mesh.replace('!SECTION, TYPE=SHELL, EGRP=ALL, MATERIAL=M1\n  0.2, 2',
        '!EGROUP, EGRP=LAYERED\n  1\n!EGROUP, EGRP=SINGLE\n  2\n'
        '!SECTION, TYPE=SHELL, EGRP=LAYERED, MATERIAL=M2\n  0.2, 2\n'
        '!SECTION, TYPE=SHELL, EGRP=SINGLE, MATERIAL=M1\n  0.2, 2')
    layered_mesh = layered_mesh.replace('!NGROUP, NGRP=FIX',
        '!MATERIAL, NAME=M2, ITEM=2\n!ITEM=1, SUBITEM=7\n'
        '  0, 40000.0, 0.3, 0.25, 40000.0, 0.3, 0.75\n'
        '!ITEM=2, SUBITEM=1\n  8.01E-10\n!NGROUP, NGRP=FIX')
    layered = run('mixed_layer_counts', 'K741J2BEND', layered_cnt, layered_mesh)
    _, layer_rows = read_result(layered/'case.res.0.12')
    check_gauss_plastic_strain(layer_rows[1], 2, output_points=16)
    check_gauss_plastic_strain(layer_rows[2], 1, output_points=16)
    check_gauss_tensors(layer_rows[1], 2, output_points=16)
    check_gauss_tensors(layer_rows[2], 1, output_points=16)
    values = [layer_rows[1][f'PLASTIC_GaussSTRAIN{k}'][0] for k in range(1, 17)]
    assert max(values)-min(values) > 1.e-6
    summary['checks']['mixed_layer_counts'].update(plastic_gauss_points=16,
        min_plstrain=min(values), max_plstrain=max(values), unused_points_zero=True)

    for name, layer_data in (
            ('layer_young_mismatch', '40000.0, 0.3, 0.25, 80000.0, 0.3, 0.75'),
            ('layer_poisson_mismatch', '40000.0, 0.3, 0.25, 40000.0, 0.2, 0.75')):
        mismatched = layered_mesh.replace('40000.0, 0.3, 0.25, 40000.0, 0.3, 0.75', layer_data)
        assert mismatched != layered_mesh
        run(name, 'K741J2BEND', layered_cnt, mismatched,
            expected_error='J2 shell layers must share the !ELASTIC Young modulus and Poisson ratio')

    # Disconnected patches place different layer counts on separate MPI ranks.
    # Prescribe the same affine membrane strain on each patch.
    separate_mesh = (CASES/'K741J2MEMBRANE.msh').read_text().replace(
        '!ELEMENT, TYPE=741\n  1, 1, 2, 3, 4',
        '  5, 10.0, 0.0, 0.0\n  6, 11.0, 0.0, 0.0\n'
        '  7, 11.25, 1.0, 0.0\n  8, 10.25, 1.0, 0.0\n'
        '!ELEMENT, TYPE=741\n  1, 1, 2, 3, 4\n  2, 5, 6, 7, 8')
    separate_mesh = separate_mesh.replace(
        '!SECTION, TYPE=SHELL, EGRP=ALL, MATERIAL=M1\n  1.0, 2',
        '!EGROUP, EGRP=SINGLE\n1\n!EGROUP, EGRP=LAYERED\n2\n'
        '!SECTION, TYPE=SHELL, EGRP=SINGLE, MATERIAL=M1\n1.0, 2\n'
        '!SECTION, TYPE=SHELL, EGRP=LAYERED, MATERIAL=M2\n1.0, 2\n'
        '!MATERIAL, NAME=M2, ITEM=2\n!ITEM=1, SUBITEM=7\n'
        '0, 210000.0, 0.3, 0.25, 210000.0, 0.3, 0.75\n!ITEM=2, SUBITEM=1\n8.01E-10')
    separate_mesh = separate_mesh.replace('!NGROUP, NGRP=ALLNODES\n',
                                          '!NGROUP, NGRP=ALLNODES\n5,6,7,8\n')
    for node in range(1, 5):
        separate_mesh = separate_mesh.replace(f'!NGROUP, NGRP=N{node}\n  {node}',
                                              f'!NGROUP, NGRP=N{node}\n  {node}, {node+4}')
    separate_cnt = integration_output((CASES/'K741J2MEMBRANE.cnt').read_text()).replace(
        '!SOLVER,', '!MATERIAL, NAME=M2\n!ELASTIC, INFINITESIMAL\n210000.0, 0.3\n'
        '!PLASTIC, INFINITESIMAL\n1000.0, 0.0\n!SOLVER,')
    separate = run('separate_layer_counts', 'K741J2MEMBRANE', separate_cnt, separate_mesh)
    _, separate_rows = read_result(separate/'case.res.0.14')
    for ident, nlayer in ((1, 1), (2, 2)):
        check_gauss_tensors(separate_rows[ident], nlayer, output_points=16)
        check_gauss_plastic_strain(separate_rows[ident], nlayer, output_points=16)

    # The elastic element has no thickness history. Its tensors must be
    # evaluated, not confused with zero padding for nonexistent points.
    mixed_cnt = (CASES/'K741J2BEND.cnt').read_text().replace('!SOLVER,',
        '!MATERIAL, NAME=M2\n!ELASTIC, INFINITESIMAL\n  40000.0, 0.3\n!SOLVER,')
    mixed_mesh = bend_mesh.replace('!SECTION, TYPE=SHELL, EGRP=ALL, MATERIAL=M1\n  0.2, 2',
        '!EGROUP, EGRP=EP\n  1\n!EGROUP, EGRP=EL\n  2\n'
        '!SECTION, TYPE=SHELL, EGRP=EP, MATERIAL=M1\n  0.2, 2\n'
        '!SECTION, TYPE=SHELL, EGRP=EL, MATERIAL=M2\n  0.2, 2')
    mixed_mesh = mixed_mesh.replace('!END',
        '!MATERIAL, NAME=M2, ITEM=2\n!ITEM=1, SUBITEM=2\n  40000.0, 0.3\n'
        '!ITEM=2, SUBITEM=1\n  8.01E-10\n!END')
    mixed_off = run('mixed_elastic_output_off', 'K741J2BEND', mixed_cnt, mixed_mesh)
    mixed_on = run('mixed_elastic_output_on', 'K741J2BEND', integration_output(mixed_cnt), mixed_mesh)
    off_blocks = read_result(mixed_off/'case.res.0.12')
    on_blocks = read_result(mixed_on/'case.res.0.12')
    for before, after in zip(off_blocks, on_blocks):
        for ident, row in before.items():
            for label, values in row.items():
                assert values == after[ident][label], (ident, label)
    for side, indices in (('-', (1, 3, 5, 7)), ('+', (2, 4, 6, 8))):
        elastic_row = on_blocks[1][2]
        for quantity in ('STRAIN', 'STRESS'):
            for component in range(6):
                actual = sum(elastic_row[f'Gauss{quantity}{k}'][component] for k in indices)/4
                # Flat elastic bending: in-plane response is linear through
                # thickness; transverse shear is constant.
                expected = elastic_row[f'Elemental{quantity}_L1{side}'][component]
                if component < 4:
                    expected /= 3**0.5
                close(actual, expected)
    assert abs(elastic_row['GaussSTRESS1'][0]) > 1
    summary['checks']['mixed_elastic_output_on'].update(
        elastic_gauss_stress_xx=elastic_row['GaussSTRESS1'][0], existing_fields_unchanged=True)

    for quantity, keyword in (('STRAIN', 'ISTRAIN'), ('STRESS', 'ISTRESS')):
        control = mixed_cnt.replace('!OUTPUT_RES', f'!OUTPUT_RES\n  {keyword}, ON')
        single = run(f'mixed_elastic_{quantity.lower()}_only', 'K741J2BEND', control, mixed_mesh)
        _, rows = read_result(single/'case.res.0.12')
        for ident, row in rows.items():
            labels = {label for label in row if label.startswith('Gauss')}
            assert labels == {f'Gauss{quantity}{k}' for k in range(1, 9)}
            for label in labels:
                assert row[label] == on_blocks[1][ident][label]

    elastic_cnt = (CASES/'K741J2BEND.cnt').read_text().replace(
        '!PLASTIC, INFINITESIMAL\n  20.0, 0.0\n', '')
    assert '!PLASTIC' not in elastic_cnt
    elastic_off = run('elastic_output_off', 'K741J2BEND', elastic_cnt)
    elastic_on = run('elastic_output_on', 'K741J2BEND', integration_output(elastic_cnt))
    before_blocks = read_result(elastic_off/'case.res.0.12')
    after_blocks = read_result(elastic_on/'case.res.0.12')
    for before, after in zip(before_blocks, after_blocks):
        for ident, row in before.items():
            for label, values in row.items():
                assert values == after[ident][label], (ident, label)
    for row in after_blocks[1].values():
        assert abs(row['GaussSTRESS1'][0]) > 1
        assert {label for label in row if label.startswith('Gauss')} == {
            f'Gauss{quantity}{k}' for quantity in ('STRAIN', 'STRESS') for k in range(1, 9)}

    elastic_layers_mesh = mixed_mesh.replace(
        '!MATERIAL, NAME=M2, ITEM=2\n!ITEM=1, SUBITEM=2\n  40000.0, 0.3',
        '!MATERIAL, NAME=M2, ITEM=2\n!ITEM=1, SUBITEM=7\n'
        '  0, 40000.0, 0.3, 0.25, 40000.0, 0.3, 0.75')
    elastic_layers = run('mixed_elastic_layers', 'K741J2BEND', integration_output(mixed_cnt), elastic_layers_mesh)
    _, rows = read_result(elastic_layers/'case.res.0.12')
    check_gauss_tensors(rows[1], 1, output_points=16)
    for layer, (center, fraction) in enumerate(((-.75, .25), (.25, .75))):
        for thick, sign in enumerate((-1, 1)):
            mean = sum(rows[2][f'GaussSTRESS{(ig*2+layer)*2+thick+1}'][0] for ig in range(4))/4
            close(mean, 12.75*(center+fraction*sign/3**0.5))

    rotated_mixed_cnt = integration_output(mixed_cnt).replace(
        '  LOAD, 3, -0.085', '  LOAD, 1, -0.051\n  LOAD, 3, -0.068')
    rotated_mixed = run('mixed_elastic_rotated', 'K741J2BEND', rotated_mixed_cnt, rotate_mesh(mixed_mesh))
    _, rows = read_result(rotated_mixed/'case.res.0.12')
    for k in range(1, 9):
        for quantity in ('STRAIN', 'STRESS'):
            value = tensor_component(rows[2][f'Gauss{quantity}{k}'], (.8, 0, -.6), quantity == 'STRAIN')
            close(value, elastic_row[f'Gauss{quantity}{k}'][0])

    # Reject unavailable tensors during setup, before any result is written.
    beam_mesh = mixed_mesh.replace('!EGROUP, EGRP=EP',
        '!ELEMENT, TYPE=611\n  3, 1, 4\n!EGROUP, EGRP=BTEST\n  3\n'
        '!SECTION, TYPE=BEAM, EGRP=BTEST, MATERIAL=M2\n'
        '  0.0, 0.0, 1.0, 1.0, 0.0833333333333333, 0.0833333333333333, 0.1406\n'
        '!EGROUP, EGRP=EP')
    beam_output = run('mixed_unavailable_tensors', 'K741J2BEND', integration_output(mixed_cnt), beam_mesh,
                      expected_error='ISTRAIN/ISTRESS require shell thickness-point output support')
    assert not list(beam_output.glob('case.res.*'))

    # Preserve the established solid-only output, including the legacy
    # repeated values in slots beyond a tetrahedron's one Gauss point.
    solid_nodes = ((0,0,0), (1,0,0), (0,1,0), (0,0,1),
                   (2,0,0), (3,0,0), (3,1,0), (2,1,0),
                   (2,0,1), (3,0,1), (3,1,1), (2,1,1))
    solid_mesh = '!NODE\n'+''.join(f'{i}, {x}, {y}, {z}\n' for i,(x,y,z) in enumerate(solid_nodes, 1))
    solid_mesh += ('!ELEMENT, TYPE=341\n1, 1,2,3,4\n'
                   '!ELEMENT, TYPE=361\n2, 5,6,7,8,9,10,11,12\n'
                   '!SECTION, TYPE=SOLID, EGRP=ALL, MATERIAL=M1\n1.0\n'
                   '!MATERIAL, NAME=M1, ITEM=2\n!ITEM=1, SUBITEM=2\n210000.0, 0.3\n'
                   '!ITEM=2, SUBITEM=1\n1.0\n!NGROUP, NGRP=ALLNODES, GENERATE\n1,12,1\n')
    solid_cnt = '!VERSION\n3\n!SOLUTION, TYPE=STATIC, NONLINEAR\n!STATIC\n!BOUNDARY, GRPID=1\nALLNODES, 2,3,0.0\n'
    for i,(x,y,z) in enumerate(solid_nodes, 1):
        solid_mesh += f'!NGROUP, NGRP=N{i}\n{i}\n'
        solid_cnt += f'N{i}, 1,1,{.012*x}\n'
    solid_mesh += '!END\n'
    solid_cnt += ('!STEP, SUBSTEPS=12, MAXITER=30, CONVERG=1.0E-8\nBOUNDARY,1\n'
                  '!MATERIAL, NAME=M1\n!ELASTIC, INFINITESIMAL\n210000.0, 0.3\n'
                  '!PLASTIC, INFINITESIMAL\n1000.0, 0.0\n'
                  '!SOLVER, METHOD=CG, PRECOND=3, ITERLOG=NO, TIMELOG=NO\n10000,1\n1.0e-8,1.0,0.0\n'
                  '!WRITE, RESULT, FREQUENCY=9999\n!OUTPUT_RES\nISTRAIN,ON\nISTRESS,ON\nPL_ISTRAIN,ON\n!END\n')
    solid_output = run('solid_mixed_gauss_counts', 'K741J2BEND', solid_cnt, solid_mesh)
    _, rows = read_result(solid_output/'case.res.0.12')
    assert rows[1]['PLASTIC_GaussSTRAIN1'][0] > 0
    for prefix in ('GaussSTRAIN', 'GaussSTRESS', 'PLASTIC_GaussSTRAIN'):
        for k in range(2, 9):
            assert rows[1][f'{prefix}{k}'] == rows[1][f'{prefix}1'], (prefix, k)

    run('missing_nonlinear', 'K741J2MEMBRANE', cnt.replace('TYPE=STATIC, NONLINEAR', 'TYPE=STATIC'),
        expected_error='require !SOLUTION, TYPE=STATIC, NONLINEAR')
    run('temperature_table', 'K741J2MEMBRANE',
        cnt.replace('!ELASTIC, INFINITESIMAL\n  210000.0, 0.3',
                    '!ELASTIC, INFINITESIMAL, DEPENDENCIES=1\n  210000.0, 0.3, 0.0\n  200000.0, 0.3, 100.0'),
        expected_error='Temperature-dependent elasticity is not supported')
    run('temperature_load', 'K741J2MEMBRANE',
        cnt.replace('!END', '!TEMPERATURE, GRPID=1\n  ALLNODES, 100.0\n!END'),
        expected_error='Temperature loading is not supported')
    run('hardening', 'K741J2MEMBRANE', cnt.replace('1000.0, 0.0', '1000.0, 100.0'),
        expected_error='requires isotropic elasticity and perfect plasticity')
    # KIRCHHOFF selects TL; omitting both flags selects the default UL model.
    for name, control in (('tl', cnt.replace('INFINITESIMAL', 'KIRCHHOFF')),
                          ('ul', cnt.replace(', INFINITESIMAL', ''))):
        run(f'finite_strain_{name}', 'K741J2MEMBRANE', control,
            expected_error='require INFINITESIMAL material kinematics')
    for npoint in (1, 3, 5):
        run(f'thickness_points_{npoint}', 'K741J2BEND',
            mesh=bend_mesh.replace('  0.2, 2', f'  0.2, {npoint}'),
            expected_error='require 2 thickness integration points in !SECTION')
    orthotropic_mesh = mesh.replace('!ITEM=1, SUBITEM=2\n  210000.0, 0.3',
        '!ITEM=1, SUBITEM=9\n  1, 210000.0, 0.3, 100000.0, 80000.0, 60000.0, 70000.0, 0.0, 1.0')
    run('orthotropic_layer', 'K741J2MEMBRANE', mesh=orthotropic_mesh,
        expected_error='requires isotropic elasticity and perfect plasticity')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
