import os
import shutil
import tempfile
import unittest
import subprocess

import rdd


class TestRdsParser(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        create_test_data_script_src = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "create_test_data.R"
        )
        shutil.copy(
            create_test_data_script_src,
            os.path.join(self.test_dir, "create_test_data.R"),
        )
        process = subprocess.run(
            ["Rscript", "create_test_data.R"],
            cwd=self.test_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if process.returncode != 0:
            raise Exception(
                f"Fail to setup tests. create_test_data.R returned code '{process.returncode}'\n"
                f"ERROR: {process.stderr}"
            )

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_atomic_vectors(self):
        vector_char = os.path.join(self.test_dir, "atomic_vector_char.rds")
        result = rdd.parse_rds(vector_char)
        expected = ["a", "b", "c", "d", "e"]
        self.assertEqual(result.values, expected)

        vector_num = os.path.join(self.test_dir, "atomic_vector_num.rds")
        result = rdd.parse_rds(vector_num)
        expected = [1, 2, 3, 4, 5]
        self.assertEqual(result.values, expected)

        vector_bool = os.path.join(self.test_dir, "atomic_vector_bool.rds")
        result = rdd.parse_rds(vector_bool)
        expcted = [True, True, False, False, True]
        self.assertEqual(result.values, expcted)

    def test_regular_sequence(self):
        vector_num = os.path.join(self.test_dir, "regular_sequence.rds")
        result = rdd.parse_rds(vector_num)
        expected = [1, 2, 3, 4, 5]
        self.assertEqual(result.values, expected)

    def test_data_frame(self):
        data_frame = os.path.join(self.test_dir, "data.frame.rds")
        result = rdd.parse_rds(data_frame)
        expected = {
            "CharVec": ["a", "b", "c"],
            "NumVec": [1, 2, 3],
            "BoolVec": [True, False, True],
        }
        self.assertEqual(result.values, {""})


if __name__ == "__main__":
    unittest.main
