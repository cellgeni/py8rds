import os
import shutil
import tempfile
import unittest
import subprocess

import pandas as pd
import py8rds

import logging
logging.disable(logging.CRITICAL) 

class TestRdsParser(unittest.TestCase):
    EXPECTED_RDS_FILES = {
        "atomic_vector_bool.rds",
        "atomic_vector_char.rds",
        "atomic_vector_num.rds",
        "data.frame_with_rownames.rds",
        "data.frame_without_rownames.rds",
        "environment.rds",
        "regular_sequence.rds",
        "seu_sketch.rds",
        "seu_sketch_no_cellnames.rds",
    }
    SEURAT_RDS_FILES = {"seu_sketch.rds", "seu_sketch_no_cellnames.rds"}

    @classmethod
    def setUpClass(cls):
        # Create a temporary directory and generate all test fixtures once.
        cls.test_dir = tempfile.mkdtemp()
        create_test_data_script_src = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "create_test_data.R"
        )
        shutil.copy(
            create_test_data_script_src,
            os.path.join(cls.test_dir, "create_test_data.R"),
        )
        process = subprocess.run(
            ["Rscript", "create_test_data.R"],
            cwd=cls.test_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if process.returncode != 0:
            raise Exception(
                f"Fail to setup tests. create_test_data.R returned code '{process.returncode}'\n"
                f"STDOUT: {process.stdout}\n"
                f"STDERR: {process.stderr}"
            )

    @classmethod
    def tearDownClass(cls):
        # Remove the directory after all tests in this class.
        shutil.rmtree(cls.test_dir)

    def _path(self, filename):
        return os.path.join(self.test_dir, filename)

    def test_expected_generated_files(self):
        generated = {f for f in os.listdir(self.test_dir) if f.endswith(".rds")}
        self.assertSetEqual(generated, self.EXPECTED_RDS_FILES)

    def test_all_generated_files_are_readable(self):
        for filename in sorted(self.EXPECTED_RDS_FILES):
            with self.subTest(filename=filename):
                result = py8rds.parse_rds(self._path(filename))
                self.assertIsNotNone(result)

    def test_atomic_vectors(self):
        vector_char = self._path("atomic_vector_char.rds")
        result = py8rds.parse_rds(vector_char)
        expected = ["a", "b", "c", "d", "e"]
        self.assertListEqual(list(result.value), expected)

        vector_num = self._path("atomic_vector_num.rds")
        result = py8rds.parse_rds(vector_num)
        expected = [1, 2, 3, 4, 5]
        self.assertListEqual(list(result.value), expected)

        vector_bool = self._path("atomic_vector_bool.rds")
        result = py8rds.parse_rds(vector_bool)
        expected = [True, True, False, False, True]
        self.assertListEqual(list(result.value), expected)

    def test_regular_sequence(self):
        vector_num = self._path("regular_sequence.rds")
        result = py8rds.parse_rds(vector_num)
        expected = [1, 2, 3, 4, 5]
        self.assertListEqual(list(result.value), expected)

    def test_data_frames(self):
        with_rownames = py8rds.as_data_frame(self._path("data.frame_with_rownames.rds"))
        self.assertEqual(with_rownames.shape, (4, 4))
        self.assertListEqual(
            with_rownames.columns.tolist(), ["CharVec", "NumVec", "BoolVec", "Factor"]
        )
        self.assertListEqual(with_rownames.index.tolist(), ["r1", "r2", "r3", "r4"])
        self.assertTrue(pd.isna(with_rownames.loc["r1", "CharVec"]))
        self.assertEqual(with_rownames.loc["r2", "CharVec"], "a")
        self.assertEqual(with_rownames.loc["r3", "CharVec"], "b")
        self.assertEqual(with_rownames.loc["r4", "CharVec"], "c")

        self.assertEqual(with_rownames.loc["r1", "NumVec"], 1)
        self.assertEqual(with_rownames.loc["r2", "NumVec"], 2)
        self.assertEqual(with_rownames.loc["r3", "NumVec"], 3)
        self.assertTrue(pd.isna(with_rownames.loc["r4", "NumVec"]))

        self.assertEqual(with_rownames.loc["r1", "BoolVec"], True)
        self.assertEqual(with_rownames.loc["r2", "BoolVec"], False)
        self.assertTrue(pd.isna(with_rownames.loc["r3", "BoolVec"]))
        self.assertEqual(with_rownames.loc["r4", "BoolVec"], True)

        self.assertEqual(with_rownames.loc["r1", "Factor"], "A")
        self.assertTrue(pd.isna(with_rownames.loc["r2", "Factor"]))
        self.assertEqual(with_rownames.loc["r3", "Factor"], "B")
        self.assertEqual(with_rownames.loc["r4", "Factor"], "C")

        without_rownames = py8rds.as_data_frame(
            self._path("data.frame_without_rownames.rds")
        )
        self.assertEqual(without_rownames.shape, (4, 4))
        self.assertListEqual(
            without_rownames.columns.tolist(),
            ["CharVec", "NumVec", "BoolVec", "Factor"],
        )
        self.assertListEqual(without_rownames.index.tolist(), [0, 1, 2, 3])
        self.assertTrue(pd.isna(without_rownames.loc[0, "CharVec"]))
        self.assertEqual(without_rownames.loc[1, "CharVec"], "a")
        self.assertEqual(without_rownames.loc[2, "CharVec"], "b")
        self.assertEqual(without_rownames.loc[3, "CharVec"], "c")

        self.assertEqual(without_rownames.loc[0, "NumVec"], 1)
        self.assertEqual(without_rownames.loc[1, "NumVec"], 2)
        self.assertEqual(without_rownames.loc[2, "NumVec"], 3)
        self.assertTrue(pd.isna(without_rownames.loc[3, "NumVec"]))

        self.assertEqual(without_rownames.loc[0, "BoolVec"], True)
        self.assertEqual(without_rownames.loc[1, "BoolVec"], False)
        self.assertTrue(pd.isna(without_rownames.loc[2, "BoolVec"]))
        self.assertEqual(without_rownames.loc[3, "BoolVec"], True)

        self.assertEqual(without_rownames.loc[0, "Factor"], "A")
        self.assertTrue(pd.isna(without_rownames.loc[1, "Factor"]))
        self.assertEqual(without_rownames.loc[2, "Factor"], "B")
        self.assertEqual(without_rownames.loc[3, "Factor"], "C")

    def test_environment_is_readable(self):
        environment = py8rds.parse_rds(self._path("environment.rds"))
        self.assertIsNotNone(environment)

    def test_seurat2adata_rna_assay_by_name_and_index(self):
        for filename in sorted(self.SEURAT_RDS_FILES):
            with self.subTest(filename=filename):
                seurat_path = self._path(filename)
                robj = py8rds.parse_rds(seurat_path)
                assay_names = robj.get(["assays", "names"]).value.tolist()
                rna_assay_index = assay_names.index("RNA")

                adata_by_name = py8rds.seurat2adata(seurat_path, assay="RNA")
                adata_by_index = py8rds.seurat2adata(seurat_path, assay=rna_assay_index)

                self.assertEqual(adata_by_name.shape, (2700, 13714))
                self.assertEqual(adata_by_index.shape, (2700, 13714))
                self.assertListEqual(
                    adata_by_name.obs_names.tolist(), adata_by_index.obs_names.tolist()
                )
                self.assertListEqual(
                    adata_by_name.var_names.tolist(), adata_by_index.var_names.tolist()
                )

    def test_seurat2adata_sketch_assay_by_name_and_index(self):
        for filename in sorted(self.SEURAT_RDS_FILES):
            with self.subTest(filename=filename):
                seurat_path = self._path(filename)
                robj = py8rds.parse_rds(seurat_path)
                assay_names = robj.get(["assays", "names"]).value.tolist()
                sketch_assay_index = assay_names.index("sketch")

                adata_by_name = py8rds.seurat2adata(seurat_path, assay="sketch")
                adata_by_index = py8rds.seurat2adata(
                    seurat_path, assay=sketch_assay_index
                )

                self.assertEqual(adata_by_name.shape, (500, 13714))
                self.assertEqual(adata_by_index.shape, (500, 13714))
                self.assertListEqual(
                    adata_by_name.obs_names.tolist(), adata_by_index.obs_names.tolist()
                )
                self.assertListEqual(
                    adata_by_name.var_names.tolist(), adata_by_index.var_names.tolist()
                )
                self.assertIn("pca", adata_by_name.obsm)
                self.assertEqual(adata_by_name.obsm["pca"].shape[0], 500)
                self.assertGreater(adata_by_name.obsm["pca"].shape[1], 0)


if __name__ == "__main__":
    unittest.main()
