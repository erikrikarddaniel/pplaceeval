FROM nfcore/base:1.12.1
LABEL authors="daniel.lundin@dbb.su.se" \
      description="Docker image containing requirements for R modules of the pplaceeval pipeline"

COPY environment.yml /
RUN conda install -c conda-forge mamba
RUN mamba env create -f /environment.yml && mamba clean -a
ENV PATH /opt/conda/envs/pplaceeval-r/bin:$PATH
