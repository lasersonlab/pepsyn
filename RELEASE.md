To generate a release:

1. Generate a summary of all the commits since the last release

    ```bash
    git log $LAST_RELEASE_TAG..HEAD
    ```

2. Tag version number and summarize changes in the tag message

    ```bash
    git tag -a vX.Y.Z
    ```

3. Push the tag upstream

    ```bash
    git push upstream vX.Y.Z
    ```

    or

    ```bash
    git push upstream --tags
    ```

4. Create the distributions

    ```bash
    rm -rf dist/*
    python setup.py sdist bdist_wheel
    ```

5. Push to PyPI

    ```bash
    twine upload dist/*
    ```
