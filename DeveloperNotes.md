# Developer Notes

This is only applicable if you want to develop the package.

## Coding Style

**Coding style:** Functions that are used externally are properly documented using Google docstring format (by default in VSCode): <https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings>

After changing the package entrypoints, you have to reinstall the package with
```{bash}
pip uninstall -y simreaduntil; pip install -e './[test,readfish,dev]'
```
This is also necessary when modifying the readfish dependency because it cannot easily be installed in editable mode with hatch: https://(github.com/pypa/hatch/issues/588).

These are some notes useful for developing:
```{bash}
# activate environment, then
pip install -e "./[test,readfish,dev]"

pip show simreaduntil

# copy pyproject.toml to tempdir, then replace root.uri with since tox does not recognize it, commands_pre/commands_post does not work because tox evaluates this before
sed -i'.backup' "s@{root:uri}@file://$(pwd)@g" pyproject.toml; tox run; mv pyproject.toml.backup pyproject.toml
sed -i'.backup' "s@{root:uri}@file://$(pwd)@g" pyproject.toml; tox run -e py38; mv pyproject.toml.backup pyproject.toml

# tox run
# To only run a specific test
# tox run -- -k test_covtracker_plotting

# test with pytest (instead of tox)
python -m pytest --cov .
python -m pytest --cov -n auto .
# for only a subset
python -m pytest ont_simulator/tests
python -m pytest --cov=. tests/simulator/gap_sampling/test_gap_sampling.py::test_run_replication
# test files should start with `test_*` (check online for more allowed patterns).

pydoctor "./src/simreaduntil"
# can also put one file to just compile it

# manually install ReadFish
# ReadFish imports its own files with `ru.*`, so it assumes that the ReadFish directory is in the `PYTHONPATH`. 
pip install -e ./external/ont_readfish
# test it with `python -c "import ru"`.

git submodule add <url> [<dir>]
git config --add oh-my-zsh.hide-dirty 1 # otherwise cd into NanoSim directory is slow

# Beware, they may destroy the git repository by accidentally changing a .gitfile (e.g. .gitmodules), so make a backup first
find . \( -type d -name .git -prune -o -type d -name external -prune -o -type d -name __pycache__ -prune -o -type d -name .tox -prune  -o -type d -name runs -prune -o -type d -name apidocs -prune \) -o -type f -print0 | pv | xargs -0 sed -i 's%simreaduntil%new_package_name%g'
# also rename the src/<packagename>, check name is no longer found (case insensitive search), reinstall the package under the new name, let vscode regenerate the caches
```

Test that the package can be imported:
```{python}
import simreaduntil.shared_utils.utils
# print(dir(simreaduntil.shared_utils.utils))
print(simreaduntil.shared_utils.utils.is_sorted([1, 2]))

# from simreaduntil.shared_utils.utils import is_sorted
# is_sorted([1, 2])
```

For profiling:

```{bash}
# first check no errors
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py
# then
scalene --html --profile-all ~/ont_project_all/ont_project/usecases/enrich_usecase.py

# Turn profiling on, only works on the main thread!
from scalene import scalene_profiler
scalene_profiler.start()
# code to profile
scalene_profiler.stop()
# scalene does not show a useful line profile when the code did not spent enough time, so try to call the function several times
```

## VSCode

### Testing in VSCode

VSCode nicely integrates with `pytest`.

Due to a bug with `conda run`, the `pytest` extension fails to discover tests in `VSCode` (as of June 2023). So use `venv` instead of `conda` (bug is related to `env/name/env/name` duplication, see <https://github.com/conda/conda/issues/11316>).

Use `Cmd+;+C` to run the current test under the cursor, `Cmd+;+F` for the current file. Breakpoints can be set!

Troubleshooting pytest Test Discovery:

- If it says `test result not found`, there is probably a syntax error. Check the `Output -> Python`. Or try to run the test file directly to detect any import errors.
Use `Clear all results` in VSCode to reset all test results, as pytest sometimes does not pick up the changes.
- If the error appears:
  ```{python}
  import file mismatch:
  imported module 'test_utils' has this __file__ attribute:
  ```
  then two identical files exist, so add `__init__.py` into every folder from the project root down until the test folder to make it a package.

If VSCode pytest stops to go into debugger, insert the following `breakpoint()`, run the test once, then you can remove this line again. Or start pytest in debug mode (it should usually enter by itself, but sometimes decides not to), or try clearing cache with `pytest --cache-clear .` or run the precise test from the command line with `pytest testfile.py::test_name` (afterwards, the debugging should work again in VSCode).

### VSCode Settings

In Settings, replace `Notebook File Root`: `${fileDirname} --> ${workspaceFolder}`
