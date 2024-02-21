## Clone repo and cd into it

`git clone https://github.com/mitchellmcm27/tcg-ec.git && cd tcg-ec`

## Start a Docker container

Build an image using the provided Dockerfile, using `-t` to give the image a human-readable tag.

`docker build -t tcg-ec .`

Start an interactive container, binding the current directory to a volume inside the container

`docker run -it -v $PWD:/workspaces/tcg-ec tcg-ec`

## Build the reactions

Within the container, everything is in the `systems/ec` directory.

`cd /workspaces/tcg-ec/systems/ec`

Custom thermodynamic datbases are included as `.tar.gz` files.

Generate the `.rxml` files by running `scripts/generate-reactions`.

It would take too long to build all the reactions. To build a specific reaction, you can pass the path to the rxml, e.g., `scripts/build-reactions database/reactions/eclogitization_agu10_stx21_rx.rxml`.

## Calculations

Calculation scripts are in the `notebooks` directory.

`cd notebooks`

Run them as follows

- `python3 parallel_pd.py` generates a (p,T) pseudosection
- `python3 parallel_profile.py` generates a 1-d profile in (p,T)-space
- `python3 parallel_experiment2.py` is a simple geodynamic experiment, lithospheric thickening

In all cases, you can pass the name of any pre-defined composition that exists in the `notebooks/compositions` folder. 
For example, `python3 parallel_pd.py -c hacker_2015_md_xenolith`.
By default the `parallel_*` scripts use all available CPU cores.

Arguments can be passed as follows to customize the model runs:

| CLI argument    |  Purpose                           | Default |
|-----------------|------------------------------------|---------|
|   `-n [int]`    | Number of CPU processes to use     | `mp.cpu_count()` |
|   `-e [float]`  | Scaled ending time                 |  1               |
|   `-c [string]` | Bulk composition to use, by name   | an array of 4 compositions |
|   `-r [string]` | Reaction to use, by name           | `eclogitization_2024_stx21_rx` |
|   `-q`          | Run in "quick mode" (_Da_ ≤ 1e3)   | False |
|   `-f`          | Force model to recalculate results | False |
|   `-p [string]`  | String to prefix output directory  | None |

