rule glswlsols:
    input:
        "src/code/sampling.py"
    output:
        "src/data/Fig3.npy"
    conda:
        "environment.yml"
    script:
        "src/code/sampling.py"