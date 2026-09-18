configfile: "snakemake.yml"

import os
from snakemake.logging import logger


EIC_SINGULARITY_CONTAINER = os.getenv("EIC_SINGULARITY_CONTAINER", "/cvmfs/eic.opensciencegrid.org/singularity/eicweb/eic_xl:nightly")


include: "benchmarks/pythia8/Snakefile"


rule metadata:
    output:
        "results/metadata.json"
    shell:
        """
cat > {output} <<EOF
{{
  "CI_COMMIT_REF_NAME": "${{CI_COMMIT_REF_NAME:-}}",
  "CI_COMMIT_SHA": "${{CI_COMMIT_SHA:-}}",
  "CI_PIPELINE_ID": "${{CI_PIPELINE_ID:-}}",
  "CI_PIPELINE_SOURCE": "${{CI_PIPELINE_SOURCE:-}}",
  "CI_PROJECT_ID": "${{CI_PROJECT_ID:-}}",
  "GITHUB_REPOSITORY": "${{GITHUB_REPOSITORY:-}}",
  "GITHUB_SHA": "${{GITHUB_SHA:-}}",
  "GITHUB_PR": "${{GITHUB_PR:-}}",
  "PIPELINE_NAME": $(echo "${{PIPELINE_NAME:-}}" | jq -Rs .)
}}
EOF
# validate JSON
jq '.' {output}
"""
