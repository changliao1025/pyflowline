import unittest


class InstallationTestCase(unittest.TestCase):
    def test_installation(self):
        try:
            import pyflowline
            from pyflowline.classes.pycase import flowlinecase

            self.assertIsNotNone(pyflowline)
            self.assertIsNotNone(flowlinecase)
            print("pyflowline is installed properly")
        except Exception as e:
            print("pyflowline is not installed properly")
            self.fail(f"Installation check failed: {e}")


def test_installation():
    try:
        import pyflowline
        from pyflowline.classes.pycase import flowlinecase

        assert pyflowline is not None
        assert flowlinecase is not None
        print("pyflowline is installed properly")
    except Exception:
        print("pyflowline is not installed properly")


if __name__ == "__main__":
    unittest.main()
