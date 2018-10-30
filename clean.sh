echo "cleaning GRIDOPT..."
find . -name \*~ -delete
find . -name \*.pyc -delete
find . -name __pycache__ -delete
rm -rf build
rm -rf dist
rm -rf GRIDOPT.egg-info
