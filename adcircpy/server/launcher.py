
def bash(driver, server_config=None):
    launcher = []
    if server_config is not None:
        launcher.extend([
            '#!/bin/bash --login',
            f'#SBATCH -D {server_config.run_directory}',
            f'#SBATCH -J {server_config.run_name}',
            f'#SBATCH -A {server_config.account}'
        ])
    else:
        launcher.extend('#!/bin/bash')
