import unittest


class InstallationTestCase(unittest.TestCase):
    def test_installation(self):
        try:
            import pyflowline
            from pyflowline.classes.pycase import flowlinecase

            print("pyflowline is installed properly")
            self.assertTrue(True)
        except Exception as e:
            print("pyflowline is not installed properly")
            self.fail(f"Installation check failed: {e}")


def test_installation():
    try:
        import pyflowline
        from pyflowline.classes.pycase import flowlinecase

        print("pyflowline is installed properly")
    except:
        print("pyflowline is not installed properly")


if __name__ == "__main__":
    unittest.main()

