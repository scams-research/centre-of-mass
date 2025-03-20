rule sampling:
    input:
        "src/scripts/sampling.py"
    output:
        "src/data/sampling.txt"
    conda:
        "environment.yml"
    cache: True
    script:
        "src/scripts/sampling.py"