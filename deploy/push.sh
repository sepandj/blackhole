#!/usr/bin/env bash
# copy the site into container 100 (run on the proxmox host)
set -euo pipefail
cd "$(dirname "$0")/.."
pct exec 100 -- mkdir -p /opt/blackhole/www
for f in index.html note.js render.js physics.js shaders.js robots.txt sitemap.xml; do
	pct push 100 "$f" "/opt/blackhole/www/$f" --perms 644
done
