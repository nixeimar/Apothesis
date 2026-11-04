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

print('All input integration tests passed.')
