# Instructions for updating the pip Package:

1. Modify the package
2. run ```python3
git tag v0.0.4` (use the actual version number)
```
3. update version in `setup.py` and `__init__.py`
4. delete `dist/` directory and run ```python3
python setup.py sdist bdist_wheel
```
5. add files, commit, and push to github
6. check distribution files using twine by running: ```python3
twine check dist/*
```
7. upload package to PyPI.org using: ```python3
twine upload dist/*
```
