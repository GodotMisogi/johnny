#!/bin/bash
cd public
(cd ..; hugo --theme=hugo-geo)
git add --all
echo "Commit message?"
read update
git commit -m "$update"
git push origin gh-pages