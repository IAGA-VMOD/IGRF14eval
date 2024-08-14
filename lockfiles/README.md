Generate with:

```
pipx run conda-lock --file ../environment-base.yml --kind env --platform linux-64
pipx run conda-lock --file ../environment-base.yml --kind explicit
```

The `env` format one (`conda-linux-64.lock.yml`) is needed by repo2docker/binder. This is symlinked to `./environment.yml` so that repo2docker will see it.

The other explicit-format ones are useful for local development. Install one with, e.g. `mamba create --name igrftest --file lockfiles/conda-linux-64.lock`.
