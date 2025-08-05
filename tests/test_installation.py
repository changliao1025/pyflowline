def test_installation():
    try:
        import pyflowline
        from pyflowline.classes.pycase import flowlinecase
        print('pyflowline is installed properly')
    except:
        print('pyflowline is not installed properly')

test_installation()