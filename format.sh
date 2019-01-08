#! /usr/bin/env bash

if [ "$1" == "check" ]; then
    ISORT_CHECK_FLAG="--check-only"
    BLACK_CHECK_FLAG="--check"
elif [ "$1" == "format" ]; then
    ISORT_CHECK_FLAG=""
    BLACK_CHECK_FLAG=""
else
    echo "usage: format.sh [check | format]"
    exit 1
fi

FAILED=0
echo ""; echo "Running isort"; echo ""
isort -m 3 -i "    " -tc -rc $ISORT_CHECK_FLAG pepsyn || FAILED=1
echo ""; echo "Running black"; echo ""
black $BLACK_CHECK_FLAG . || FAILED=1

if [ "$FAILED" -ne 0 ]; then
    exit 1
fi
