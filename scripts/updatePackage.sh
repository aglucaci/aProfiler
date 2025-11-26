

#pip install build
#pip install twine

python -m build

twine upload --repository testpypi dist/* --verbose


# pip install --index-url https://test.pypi.org/simple/ your_package_name


