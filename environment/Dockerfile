FROM mambaorg/micromamba:2.1.0

# Added matplotlib and seaborn to package definitions
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    samtools=1.22 \
    python=3.11 \
    numpy=1.26.4 \
    pandas=2.2.2 \
    scipy=1.14.0 \
    scikit-learn=1.5.0 \
    statsmodels=0.14.2 \
    networkx=3.2.1 \
    pyarrow=16.1.0 \
    matplotlib=3.8.4 \
    pytest=8.2.2 \
    seaborn=0.13.2 && \
    micromamba clean --all --yes 

USER root
SHELL ["/bin/bash", "-lc"]

ENV WORKSPACE_DIR=/workspace
ENV DATA_DIR=/workspace/data
ENV OUTPUT_DIR=/workspace/output
ENV SOLUTION_DIR=/solution
ENV TESTS_DIR=/tests
ENV LOG_DIR=/logs/verifier

# Create the directories, including the target data mount point
RUN mkdir -p ${WORKSPACE_DIR} ${DATA_DIR} ${OUTPUT_DIR} ${SOLUTION_DIR} ${TESTS_DIR} ${LOG_DIR}

# --- REMOVED COPY LINE ---

WORKDIR ${SOLUTION_DIR}

COPY solution/cram_coverage.py .
RUN chmod +x cram_coverage.py

RUN echo 'import subprocess; subprocess.run(["python", "cram_coverage.py"], check=True)' > solve.py

CMD ["python", "solve.py"]

