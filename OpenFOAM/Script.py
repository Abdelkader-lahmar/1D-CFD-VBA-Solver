from shutil import copytree
from pathlib import Path
import subprocess


def main():
    base_case = Path('./Base/')
    for Scheme in ["linear", "upwind", "patankarLinearUpwind", "patankarQUICK", "patankarVanLeer", "QUICK", "vanLeer"]:
        if Scheme in ['linear', 'QUICK', "vanLeer"]:
            additional = ''
        else:
            additional = 'grad(T)'
        for Pe in [1, 2, 4, 10, 25, 50]:
            new_case = Path(f'./{Scheme}_{Pe}')
            copytree(base_case, new_case)
            change(new_case, Scheme, Pe, additional)
            subprocess.run(['scalarTransportFoam'],
                           cwd=new_case,
                           check=True,
                          )
            write_result(new_case, Scheme, Pe)


def write_result(case, Scheme, Pe):
    with open('Results.csv', 'a') as results_output:
        results_output.write(f"{Scheme}, {Pe}, ")
        results = Path(f'{get_latest(case)}/T').open('r')
        for line in results:
            if 'internalField' not in line:
                continue
            results_output.write(f'{line[line.index("(") +
                                         1:line.index(')')]}\n'.replace(' ', ','))


def get_latest(case):
    times = []
    for time in case.iterdir():
        if not time.is_dir():
            continue
        try:
            times.append(int(time.name))
        except ValueError:
            continue
    return Path(f'{case}/{max(times)}') 



def change(case, scheme, Pe, additional):
    files = {
        'initial_condition': f'{case}/0/U',
        'scheme': f'{case}/system/fvSchemes',
        'boundary_type': f'{case}/0/T',
    }
    changes = {
        'initial_condition': [
            'uniform (1 0 0)',
            f'uniform ({Pe} 0 0)'
        ],
        'scheme': [
            'div(phi,T)      Gauss',
            f'div(phi,T)      Gauss {scheme} {additional}'
        ],
        'boundary_type': [
            'patankarFixedValue',
            'fixedValue'
        ],
    }
    for key in files:
        if key == 'boundary_type' and scheme != 'linear':
            continue
        file = Path(files[key])
        content = file.read_text()
        content = content.replace(changes[key][0], changes[key][1])
        file.write_text(content)


if __name__ == "__main__":
    main()
