import os

from random import choice
from sys import stdout
from uuid import uuid4

import yaml

GPU = {"01-beeler-reuter.jmd"}
LONG = {"02-workshop_solutions.jmd", "06-pendulum_bayesian_inference.jmd"}

cwd = os.getcwd()
os.chdir("tutorials")
folder = os.getenv("FOLDER")
if not folder:
    folder = choice([d for d in os.listdir() if os.path.isdir(d)])
file = os.getenv("FILE")
if not file:
    file = choice([f for f in os.listdir(folder) if f.endswith(".jmd") and f not in LONG])
os.chdir(cwd)

tags = ["nvidia"] if file in GPU else []

script = f"""
julia -e '
  using Pkg
  Pkg.instantiate()
  using DiffEqTutorials: weave_file
  weave_file("{folder}", "{file}")'

if [[ -z "$(git status -suno)" ]]; then
  echo "No changes"
  exit 0
fi

k="$(cat $SSH_KEY)"
echo "$k" > "$SSH_KEY"
chmod 400 "$SSH_KEY"
git config core.sshCommand "ssh -o StrictHostKeyChecking=no -i $SSH_KEY"
git config user.name "github-actions[bot]"
git config user.email "actions@github.com"
branch="rebuild/{str(uuid4())[:8]}"
git checkout -b "$branch"
git commit -am "Rebuild tutorials"
git remote add github "git@github.com:SciML/DiffEqTutorials.jl.git"
git push github "$branch"
"""

pipeline = {
    "include": "https://raw.githubusercontent.com/JuliaGPU/gitlab-ci/master/templates/v6.yml",
    "rebuild": {
        "extends": ".julia:1.4",
        "variables": {
            "CI_APT_INSTALL": "git python3-dev texlive-science texlive-xetex",
            "JULIA_NUM_THREADS": 4,
        },
        "tags": tags,
        "script": script,
    },
}

yaml.dump(pipeline, stdout)
