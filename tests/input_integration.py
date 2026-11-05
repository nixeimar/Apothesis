"""Run real CLI/simulation paths in temporary directories; leave the checkout clean."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
root = Path(sys.argv[2]).resolve()
original = json.loads((root / 'examples' / 'height-field.json').read_text())

with tempfile.TemporaryDirectory(prefix='apothesis-integration-') as temp:
    work = Path(temp)

    def run(args, success=True):
        result = subprocess.run([binary, *args], cwd=work, capture_output=True,
                                text=True, timeout=30)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result

    config = json.loads(json.dumps(original))
    config['simulation'] = {'start_time': 0, 'duration': 0.05, 'random_seed': 1234}
    path = work / 'input.json'
    path.write_text(json.dumps(config))
    run(['--validate-input', str(path)])
    assert list(work.iterdir()) == [path], 'validation wrote simulation output'
    run(['--input', str(path)])
    assert (work / 'Output.log').exists()
    assert 'co_adsorption' in (work / 'Output.log').read_text()
    run(['-i', str(path)])
    run([str(path)])

    # No arguments select input.json; legacy text is rejected.
    run([])
    legacy_path = work / 'old.kmc'
    legacy_path.write_text('lattice: SimpleCubic 10 10 10 X\n')
    run(['--input', str(legacy_path)], success=False)

    # Two independently configured channels with the same equation remain distinct.
    channel = json.loads(json.dumps(config['processes'][0]))
    channel['id'] = 'second_channel'
    config['processes'].append(channel)
    path.write_text(json.dumps(config))
    run([str(path)])
    log = (work / 'Output.log').read_text()
    assert 'co_adsorption' in log and 'second_channel' in log

    # Exercise the Arrhenius bridge, including the formerly incorrect diffusion indices.
    for kind in ['diffusion', 'desorption']:
        process = {
            'id': kind, 'type': kind,
            'reactants': [{'species': 'CO', 'phase': 'surface'}],
            'products': ([{'species': 'CO', 'phase': 'surface'}] if kind == 'diffusion'
                         else [{'vacancy': True}, {'species': 'CO', 'phase': 'gas'}]),
            'rate': {'model': 'arrhenius', 'parameters': {'prefactor': 1.0, 'activation_energy': 50000}}
        }
        if kind == 'diffusion':
            process['rate']['parameters']['migration_energy'] = 30000
        for all_neighbors in [False, True]:
            process['all_neighbors'] = all_neighbors
            config['processes'] = original['processes'] + [process]
            path.write_text(json.dumps(config))
            run([str(path)])

    # File paths are relative to input.json, matrices have y rows and x columns.
    config['lattice']['parameters'].update(dimensions=[4, 3], initial_height_file='heights.dat', initial_species_file='species.dat')
    del config['lattice']['parameters']['initial_height']
    config['processes'] = original['processes']
    (work / 'heights.dat').write_text('1 2 3 4\n5 6 7 8\n9 10 11 12\n')
    (work / 'species.dat').write_text('X X X X\nX X X X\nX X X X\n')
    path.write_text(json.dumps(config))
    run([str(path)])
    heights = (work / 'Height_0.dat').read_text().splitlines()[1:]
    assert [[int(x) for x in row.split()] for row in heights] == [[1,2,3,4],[5,6,7,8],[9,10,11,12]]

    # Atomistic geometry builds without invoking the height-field KMC engine.
    atomic = json.loads((root / 'examples' / 'atomistic-cu.json').read_text())
    atomic['lattice']['parameters']['vacancies'] = [0, 10]
    atomic_path = work / 'atomic.json'
    atomic_path.write_text(json.dumps(atomic))
    run(['--validate-input', str(atomic_path)])
    result = run(['--build-surface', str(atomic_path), '--output', 'copper.xyz'])
    assert '574 atoms on 576 crystal sites' in result.stdout
    xyz = (work / 'copper.xyz').read_text().splitlines()
    assert xyz[0] == '574' and len(xyz) == 576
    assert all(int(row.split()[4]) not in [0, 10] for row in xyz[2:])
    result = run(['--input', str(atomic_path)], success=False)
    assert 'geometry-only input' in result.stderr
    atomic['lattice']['parameters']['orientation'] = '111'
    atomic_path.write_text(json.dumps(atomic))
    result = run(['--validate-input', str(atomic_path)], success=False)
    assert 'orientation' in result.stderr

    # Full atomistic KMC: rate validation, reproducible output and atom balance.
    kinetic = json.loads((root / 'examples' / 'atomistic-kmc.json').read_text())
    kinetic_path = work / 'kinetic.json'
    kinetic_path.write_text(json.dumps(kinetic))
    run(['--validate-input', str(kinetic_path)])
    run(['--input', str(kinetic_path)])
    log = (work / 'atomistic.log').read_text()
    rows = [line.split() for line in log.splitlines() if line and not line.startswith('#')]
    last = rows[-1]
    counts = [int(n) for n in last[3:]]
    assert sum(counts) == int(last[2]) and int(last[1]) == 128 + counts[0] - counts[2]
    assert float(last[0]) <= kinetic['simulation']['duration']
    final = (work / 'atomistic_final.xyz').read_text()
    assert int(final.splitlines()[0]) == int(last[1])
    run([str(kinetic_path)])
    assert final == (work / 'atomistic_final.xyz').read_text()
    assert log == (work / 'atomistic.log').read_text()
    for mutate, field in [
        (lambda c: c['processes'][0].setdefault('parameters', {}).update(min_coordination=0), 'min_coordination'),
        (lambda c: c['processes'][0]['rate']['parameters'].update(value=-1), 'value'),
        (lambda c: c['processes'][0].update(type='reaction'), 'type'),
        (lambda c: c['lattice']['parameters'].update(periodic=[True,True,True],vacuum=0), 'periodic'),
        (lambda c: c['processes'][0]['rate'].update(model='simple'), 'rate.model'),
        (lambda c: c['processes'][0].update(all_neighbors=True), 'all_neighbors')
    ]:
        candidate = json.loads(json.dumps(kinetic))
        mutate(candidate)
        kinetic_path.write_text(json.dumps(candidate))
        result = run(['--validate-input', str(kinetic_path)], success=False)
        assert field in result.stderr, result.stderr

    # Invalid input produces an actionable diagnostic and a nonzero exit.
    config['system']['parameters']['temperature'] = -1
    path.write_text(json.dumps(config))
    result = run(['--validate-input', str(path)], success=False)
    assert 'system.parameters.temperature' in result.stderr
    run(['--input'], success=False)


    # Named height-field model: actual CLI dispatch, inventory, reproducibility.
    zgb = json.loads((root / 'examples/zgb-growth-2023.json').read_text())
    zgb['lattice']['parameters']['dimensions'] = [4, 4]
    zgb['simulation'].update(duration=100, random_seed=1234)
    zgb_path = work / 'zgb.json'
    zgb_path.write_text(json.dumps(zgb))
    run(['--validate-input', str(zgb_path)])
    run([str(zgb_path)])
    summary = json.loads((work / 'zgb_summary.json').read_text())
    assert summary['events'] > 0 and summary['reactions'] > 0
    assert abs(summary['mean_height'] - (10 + summary['reactions']/16)) < 1e-10
    assert abs(summary['coverage_ab'] - (summary['ab_adsorptions']-summary['reactions'])/16) < 1e-10
    run([str(zgb_path)])
    assert summary == json.loads((work / 'zgb_summary.json').read_text())
    for section, field, value in [('model', 'y_ab', 1.1), ('lattice', 'dimensions', [2,4])]:
        candidate = json.loads(json.dumps(zgb));candidate[section]['parameters'][field] = value
        zgb_path.write_text(json.dumps(candidate))
        run(['--validate-input',str(zgb_path)],success=False)

    # 3D chemistry must export real occupied sites and separate adsorbate roles.
    zgb3d = json.loads((root / 'examples/zgb-growth-3d.json').read_text())
    zgb3d['lattice']['parameters'].update(repetitions=[4,4,16], occupied_unit_cells=2)
    zgb3d_path = work / 'zgb3d.json'
    zgb3d_path.write_text(json.dumps(zgb3d))
    run(['--validate-input', str(zgb3d_path)])
    run([str(zgb3d_path)])
    summary = json.loads((work / 'zgb3d_summary.json').read_text())
    assert summary['representation'] == 'explicit_3d_occupancy'
    assert summary['solid_units'] == 32 + summary['reactions']
    xyz = (work / 'zgb3d_final.xyz').read_text().splitlines()
    assert 'chemical_state:S:1:role:S:1' in xyz[1]
    points = [row.split() for row in xyz[2:]]
    assert len(points) == int(xyz[0])
    assert sum(row[7] == 'solid' for row in points) == summary['solid_units']
    assert sum(row[6] == 'AB2_solid' for row in points) == summary['reactions']
    assert sum(row[6] == 'AB_adsorbate' for row in points) == summary['ab_adsorptions'] - summary['reactions']
    assert sum(row[6] == 'B_adsorbate' for row in points) == 2*summary['b2_adsorptions'] - summary['reactions']
    candidate = json.loads(json.dumps(zgb3d));candidate['lattice']['type']='FCC'
    zgb3d_path.write_text(json.dumps(candidate))
    run(['--validate-input',str(zgb3d_path)],success=False)

    # Exercise both ensemble modes and verify the runner really chooses 3D.
    for representation in ['height_field','3d']:
        sweep = work / ('sweep-' + representation)
        result = subprocess.run([sys.executable,str(root/'scripts/reproduce_zgb_2023.py'),
            '--binary',binary,'--output',str(sweep),'--representation',representation,
            '--size','4','--trials','2','--y','0.4','--variant','growth',
            '--growth-duration','100'],capture_output=True,text=True,timeout=30)
        assert result.returncode == 0, result.stdout+result.stderr
        summary_name = 'zgb3d_summary.json' if representation == '3d' else 'zgb_summary.json'
        assert len(list(sweep.glob('*/'+summary_name))) == 2
        if representation == '3d':
            assert not list(sweep.glob('*/zgb_summary.json'))
            assert not list(sweep.glob('*/*.xyz'))

print('All input integration tests passed.')
