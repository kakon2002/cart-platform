# The design platform, as one container.
#
# What is NOT in here is the point: `.dockerignore` and `.gcloudignore` exclude
# the 8.3 GB single-cell matrix and the 2.6 GB archive it was expanded from. The
# loader returns a cached 5.7 MB summary whenever the fingerprint matches, so
# the serving path never opens either one — measured, not assumed: with the
# matrix renamed away a full screen completed in 7.3 seconds and reached the
# same end state. That takes the image from ~12 GB to ~1.1 GB.

FROM python:3.13-slim

# No compiler in the image: every dependency publishes manylinux wheels, and a
# build toolchain would add ~300 MB to carry nothing.
#
# PORT is set here rather than left to the runtime. The application falls back
# to 8000, so without this a plain `docker run -p 8080:8080` would refuse
# connections with nothing in the log to explain it. Cloud Run overrides it.
#
# CART_NO_MATRIX_FETCH turns a single-cell cache miss into a named error. The
# manifests for the excluded archives ship with the derived summaries, so
# without it a miss falls through to downloading 2.6 GB and expanding it to
# 8.3 GB onto a filesystem that is in memory here — an instance killed mid-job
# rather than anything a caller could read.
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8080 \
    CART_NO_MATRIX_FETCH=1

# Created before anything is copied. A `chown -R` afterwards would rewrite every
# one of the 648 MB of cache files into a second layer — layers store whole
# files, not metadata deltas — and roughly double the image.
RUN useradd --create-home --uid 10001 platform
WORKDIR /app

# Dependencies first, so editing the package does not re-resolve them.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=platform:platform car_pipeline/ ./car_pipeline/

# The caches. 648 MB, and the largest single item is the DepMap CSV at 413 MB.
# `data/stage5` comes with it deliberately: it is the retrieved binder set, and
# shipping it warm means the first screen skips one network call per pool member
# rather than spending five minutes on them. Delete `data/stage5` in a running
# container to force the real retrieval. (`data/stage4` rides along at 65 KB but
# nothing in the served path reads it — the pipeline recomputes pairing every
# run. It is here for parity with the local tree, not for speed.)
COPY --chown=platform:platform data/ ./data/

USER platform

# 0.0.0.0 because a container that binds loopback is unreachable from outside
# itself. The port comes from PORT above, which Cloud Run replaces with its own.
EXPOSE 8080
CMD ["python", "-m", "car_pipeline.api.server", "--host", "0.0.0.0"]
