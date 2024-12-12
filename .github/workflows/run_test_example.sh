#!/bin/bash


# hacky way to hopefully alias the calls
# normal alias does not seem to work in GitHub CI
# necessary, as the hacky git install does not install the CLI entrypoints
echo "#!/bin/bash python <(echo 'import syri.scripts.syri;syri.scripts.syri.main()')" > syri
chmod +x ./syri
echo "minimap2.py" > ./minimap2
chmod +x ./minimap2
PATH=$PATH:./
syri --version
minimap2
# run using source to preserve alias
source ./example/example_workflow.sh
