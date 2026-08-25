# Correlate (V1) — the build served at correlate.cmm.se.
#
# Same shape as Greenlisted_v2's Dockerfile on purpose: static site, official
# Apache image, files copied into the web root. There is no backend and no build
# step here — this repo IS the built site (generated from the V2 repo by
# scripts/build_v1.py), so the image is just those files behind a web server.
#
#   podman build -t correlate:latest .
#   podman run --rm -p 8080:80 correlate:latest
FROM docker.io/library/httpd:alpine

COPY --chown=www-data:www-data . /usr/local/apache2/htdocs
