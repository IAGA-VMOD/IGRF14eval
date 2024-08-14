Generate with:

```
pipx run conda-lock --file environment-base.yml --kind env --platform linux-64
pipx run conda-lock --file environment-base.yml --kind explicit
```

The `env` format one is needed by binder (I think - `conda-linux-64.lock.yml`. Manually copy this into the `binder/environment.yml` so that binder will use it.

The other explicit-format ones are useful for local development. Install one with, e.g. `mamba create --name igrftest --file binder/lockfiles/conda-linux-64.lock`.
