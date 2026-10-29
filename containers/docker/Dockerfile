FROM mambaorg/micromamba:2.6.2

ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN mkdir -p /opt/conda/envs /opt/celatlas-spatial-slam /workspace/{config,inputs,log,outputs,pipeline,reference,tmp} \
    && chown -R ${MAMBA_USER}:${MAMBA_USER} /opt/conda /opt/celatlas-spatial-slam /workspace
USER ${MAMBA_USER}

COPY --chown=${MAMBA_USER}:${MAMBA_USER} snakemake.yaml /opt/specs/snakemake.yaml
COPY --chown=${MAMBA_USER}:${MAMBA_USER} runtime.yaml /opt/specs/runtime.yaml
RUN micromamba create -y -p /opt/conda/envs/orchestrator -f /opt/specs/snakemake.yaml \
    && micromamba create -y -p /opt/conda/envs/runtime -f /opt/specs/runtime.yaml \
    && micromamba clean --all --yes

COPY --chown=${MAMBA_USER}:${MAMBA_USER} runtime_celatlas18/ /opt/celatlas18/
COPY --chown=${MAMBA_USER}:${MAMBA_USER} system_libs/ /opt/runtime/system_libs/
COPY --chown=${MAMBA_USER}:${MAMBA_USER} swin_tiny.pth /opt/celatlas-spatial-slam/celatlas/src/swin_tiny.pth
COPY --chown=${MAMBA_USER}:${MAMBA_USER} gedi-bin/ /opt/celatlas-spatial-slam/gedi-bin/
COPY --chown=${MAMBA_USER}:${MAMBA_USER} epi-logo.png /opt/celatlas-spatial-slam/assets/epi-logo.png
RUN chmod -R a+rX /opt/celatlas-spatial-slam /opt/celatlas18 /opt/conda \
    && sed -i 's#/opt/project/software/micromamba-root/envs/celatlas18-stardist-py311#/opt/celatlas18#g' /opt/celatlas18/lib/python3.11/site-packages/cv2/config.py /opt/celatlas18/lib/python3.11/site-packages/cv2/config-3.11.py 2>/dev/null || true
RUN ln -sf libGL.so.1.7.0 /opt/runtime/system_libs/libGL.so.1 \
    && ln -sf libGLdispatch.so.0.0.0 /opt/runtime/system_libs/libGLdispatch.so.0 \
    && ln -sf libGLX.so.0.0.0 /opt/runtime/system_libs/libGLX.so.0 \
    && ln -sf libX11.so.6.4.0 /opt/runtime/system_libs/libX11.so.6 \
    && ln -sf libxcb.so.1.1.0 /opt/runtime/system_libs/libxcb.so.1 \
    && ln -sf libXau.so.6.0.0 /opt/runtime/system_libs/libXau.so.6 \
    && ln -sf libXdmcp.so.6.0.0 /opt/runtime/system_libs/libXdmcp.so.6 \
    && ln -sf libbsd.so.0.11.5 /opt/runtime/system_libs/libbsd.so.0 \
    && ln -sf libmd.so.0.0.5 /opt/runtime/system_libs/libmd.so.0 \
    && printf '#!/bin/sh\nexec /opt/celatlas18/bin/python -c '\''from celatlas_spatial.runner.public_cli import main; raise SystemExit(main())'\'' "$@"\n' > /opt/conda/envs/runtime/bin/celatlas \
    && printf '#!/bin/sh\nexec /opt/celatlas18/bin/python -c '\''from celatlas_spatial.celatlas import main; raise SystemExit(main())'\'' "$@"\n' > /opt/conda/envs/runtime/bin/celatlas_spatial \
    && printf '#!/bin/sh\nexec /opt/celatlas18/bin/python -c '\''from celatlas_spatial.tools.spatial_report_generator import main; raise SystemExit(main())'\'' "$@"\n' > /opt/conda/envs/runtime/bin/celatlas_spatial_report \
    && chmod 755 /opt/conda/envs/runtime/bin/celatlas /opt/conda/envs/runtime/bin/celatlas_spatial /opt/conda/envs/runtime/bin/celatlas_spatial_report
RUN /opt/conda/envs/runtime/bin/python -m pip install --no-cache-dir \
      Cython 'leidenalg==0.10.2' \
    && /opt/conda/envs/runtime/bin/python -m pip check

COPY --chown=${MAMBA_USER}:${MAMBA_USER} --chmod=755 subread/bin/featureCounts /opt/subread/bin/featureCounts
COPY --chown=${MAMBA_USER}:${MAMBA_USER} subread/lib/libz.so.1 /opt/subread/lib/libz.so.1
RUN /opt/subread/bin/featureCounts -v

ENV PATH=/opt/conda/envs/runtime/bin:/opt/celatlas18/bin:/opt/subread/bin:/opt/conda/envs/orchestrator/bin:/opt/celatlas-spatial-slam/gedi-bin:${PATH}
ENV CELATLAS_MODEL=/opt/celatlas-spatial-slam/celatlas/src/swin_tiny.pth
ENV CELATLAS_PYTHON=/opt/celatlas18/bin/python
ENV PYTHONPATH=/opt/celatlas18/lib/python3.11/site-packages:/opt/celatlas-spatial-slam/celatlas
ENV LD_LIBRARY_PATH=/opt/runtime/system_libs:/opt/celatlas18/lib:/opt/conda/envs/runtime/lib
ENV SSL_CERT_FILE=/opt/celatlas18/ssl/cacert.pem
ENV CURL_CA_BUNDLE=/opt/celatlas18/ssl/cacert.pem
ENV REQUESTS_CA_BUNDLE=/opt/celatlas18/ssl/cacert.pem
ENV MPLBACKEND=Agg
WORKDIR /workspace
ENTRYPOINT ["/opt/conda/envs/orchestrator/bin/snakemake"]
CMD ["--help"]
