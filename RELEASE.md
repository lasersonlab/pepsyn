To generate a release:

1. Generate a summary of all the commits since the last release

    ```bash
    git log $LAST_RELEASE_TAG..HEAD
    ```

2. Set the release version in `setup.py` (remove the `.dev0` tag if applicable)
and commit the version number change.  Also set the new version number in the
readme if applicable.

3. Tag version number and summarize changes in the tag message

    ```bash
    git tag -a vX.Y.Z
    ```

4. Push the tag upstream

    ```bash
    git push upstream vX.Y.Z
    ```

    or

    ```bash
    git push upstream --tags
    ```

5. Create the distributions

    ```bash
    rm -rf dist/*
    python setup.py sdist bdist_wheel
    ```

6. Push to PyPI

    ```bash
    twine upload dist/*
    ```

7. If working on master, bump up to the next anticipated version with a `.dev0`
tag and commit
