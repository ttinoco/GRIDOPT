echo "cleaning..."
find . -name \*~ -delete
find . -name \*.pyc -delete
rm -rf build
rm -rf dist
rm -rf GRIDOPT.egg-info
