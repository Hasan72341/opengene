#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PORT="${PORT:-4174}"
step() { printf '\n=== %s ===\n' "$1"; }

step "build images"
docker compose --profile ci --profile tools build

step "test image: python + web suites, comment check"
docker compose --profile ci run --rm test

step "cli image: smoke workflows"
docker compose --profile tools run --rm cli info fixtures/cli-demo.fasta
docker compose --profile tools run --rm cli digest fixtures/cli-demo.fasta EcoRI BamHI
docker compose --profile tools run --rm cli align fixtures/overlap-left.fasta fixtures/overlap-right.fasta

step "web image: serve and probe"
docker compose up -d web
trap 'docker compose down --remove-orphans >/dev/null 2>&1 || true' EXIT

health=""
for _ in $(seq 1 40); do
  health=$(docker inspect --format '{{.State.Health.Status}}' "$(docker compose ps -q web)" 2>/dev/null || echo unknown)
  [ "$health" = "healthy" ] && break
  sleep 1
done
[ "$health" = "healthy" ] || { echo "FAILED: web container never became healthy (last status: ${health})"; exit 1; }
echo "container healthcheck: ${health}"

for path in / /index.html /src/app.js /src/styles.css /service-worker.js /manifest.webmanifest /assets/opengene-mark.svg; do
  code=$(curl -s -o /dev/null -w '%{http_code}' "http://localhost:${PORT}${path}")
  type=$(curl -s -o /dev/null -w '%{content_type}' "http://localhost:${PORT}${path}")
  printf '%-32s %s  %s\n' "$path" "$code" "$type"
  [ "$code" = "200" ] || { echo "FAILED: ${path} returned ${code}"; exit 1; }
done

curl -fsS "http://localhost:${PORT}/src/app.js" | grep -q "renderVerificationTool" \
  || { echo "FAILED: served app.js is missing expected symbols"; exit 1; }

step "all checks passed"
