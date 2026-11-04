"""Run real CLI/simulation paths in temporary directories; leave the checkout clean."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
root = Path(sys.argv[2]).resolve()
original = json.loads((root / 'input.json').read_text())

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

    # Invalid input produces an actionable diagnostic and a nonzero exit.
    config['system']['parameters']['temperature'] = -1
    path.write_text(json.dumps(config))
    result = run(['--validate-input', str(path)], success=False)
    assert 'system.parameters.temperature' in result.stderr
    run(['--input'], success=False)

print('All input integration tests passed.')
