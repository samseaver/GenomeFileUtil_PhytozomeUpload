#!/usr/bin/env bash
export PATH=/homes/chicago/seaver/kb_sdk/bin:$PATH
export PATH=/homes/chicago/seaver/jdk-15.0.1/bin:$PATH
source /homes/chicago/seaver/kb_sdk/src/sh/sdk-completion.sh
cat test.date
date
date > test.date
kb-sdk test -s > test.out 2>test.err
date >> test.date
date
