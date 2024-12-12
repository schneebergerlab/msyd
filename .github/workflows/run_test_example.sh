#!/bin/bash

# alias to a call launching the syri entrypoint from python
# necessary, as the hacky git install does not install the CLI entrypoints
alias syri='python <(echo "import syri.scripts.syri;syri.scripts.syri.main()")'

# run using source to preserve alias
source ./example/example_workflow.sh
