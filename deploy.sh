#!/bin/bash
set -e

echo "Rendering site..."
quarto render

echo "Pushing to GitHub..."
git add .
git commit -m "update site $(date '+%Y-%m-%d')"
git push origin main

echo "✅ Done! Site will update in ~1 min at victoriasja18.github.io/eris-progress/"