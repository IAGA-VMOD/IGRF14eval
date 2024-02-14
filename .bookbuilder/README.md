Here is an extra conda environment, `bookbuilder`, that provides tools for executing the notebooks and assembling the html site using jupyter-book. They are kept separate from the notebook kernel for easier maintenance.

Example steps for local use:

1. Create the environment from the lockfile:
    ```
    conda create --name bookbuilder --file .bookbuilder/conda-linux-64.lock
    ```
2. Configure the jupyter kernel for the existing `igrf` environment if necessary:
    ```
    conda run --no-capture-output -n igrf python -m ipykernel install --user --name igrf
    ```
3. Execute the notebooks against the `igrf` kernel:
    ```
    conda run -n bookbuilder pytest --numprocesses auto --nbmake --overwrite --nbmake-kernel=igrf notebooks/*.ipynb
    ```
4. Build the jupyter book:
    ```
    conda run -n bookbuilder jupyter-book build .
    ```

Notes:

- jupyter-book is configured via `_config.yml` and `_toc.yml`
- We use pytest and the [nbmake plugin](https://github.com/treebeardtech/nbmake) to execute the notebooks first before running jupyter-book. This means that we can separate the steps for easier debugging (notebooks may fail to execute, or jupyter-book may fail to build). It also enables parallel execution for faster builds.
- `.github/workflows/publish.yml` configures the book building and push to github pages in CI
    - Github settings are at https://github.com/IAGA-VMOD/IGRF14eval/settings/pages
    - By default the site can only be published from the `main` branch. Alter the deployment protection rules in https://github.com/IAGA-VMOD/IGRF14eval/settings/environments to change this.