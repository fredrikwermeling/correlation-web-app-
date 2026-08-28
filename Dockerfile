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

# The image ships with mod_headers off and AllowOverride None, which would make
# the .htaccess beside this file dead weight. That file is what stops browsers
# serving a months-old index.html out of cache (see its own comments), so the
# two have to be enabled for it to mean anything. FileInfo is the narrowest
# override level that permits Header directives. Both seds are no-ops if the
# base image ever ships these already set.
RUN sed -i -e 's|^#LoadModule headers_module|LoadModule headers_module|' \
           -e 's|AllowOverride None|AllowOverride FileInfo|' \
    /usr/local/apache2/conf/httpd.conf

COPY --chown=www-data:www-data . /usr/local/apache2/htdocs
