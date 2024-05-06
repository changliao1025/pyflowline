
def test_installation():
    try:
        import pyflowline
        from pyflowline.classes.pycase import flowlinecase
        print('pyflowline is installed properly')
    except:
        print('pyflowline is not installed properly')

test_installation()

# This might be useful:
# import pkg_resources
# def test_pkg_installation():
#     try:
#         distribution = pkg_resources.get_distribution('pyflowline')
#         print("Pyflowline is installed at:", distribution.location)
#     except pkg_resources.DistributionNotFound:
#         print("Pyflowline is not installed.")
