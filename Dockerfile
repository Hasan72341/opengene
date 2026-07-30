FROM python:3.12-slim AS base
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1
WORKDIR /opt/opengene
COPY pyproject.toml README.md LICENSE ./
COPY apps ./apps
COPY packages ./packages
RUN pip install --no-cache-dir .


FROM base AS cli
RUN useradd --create-home --uid 10001 opengene
WORKDIR /work
USER opengene
ENTRYPOINT ["opengene"]
CMD ["--help"]


FROM base AS test
RUN apt-get update \
 && apt-get install --no-install-recommends -y nodejs \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/opengene
COPY tests ./tests
COPY scripts ./scripts
COPY index.html service-worker.js manifest.webmanifest package.json ./
COPY src ./src
COPY assets ./assets
CMD ["bash", "-lc", "python tests/run_all.py && node --test tests/web/*.test.mjs && python scripts/strip_comments.py --check"]


FROM nginx:1.27-alpine AS web
COPY deploy/nginx.conf /etc/nginx/conf.d/default.conf
WORKDIR /usr/share/nginx/html
COPY index.html service-worker.js manifest.webmanifest ./
COPY src ./src
COPY assets ./assets
EXPOSE 80
HEALTHCHECK --interval=10s --timeout=3s --start-period=3s --retries=5 \
  CMD wget --spider -q http://127.0.0.1/index.html || exit 1
