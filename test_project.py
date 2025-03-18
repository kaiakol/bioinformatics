import os
import unittest
from project import reshape_data, postprocess_output

class TestPostprocessOutput(unittest.TestCase):
    def setUp(self):
        self.output_dir = 'meme_output'
        os.makedirs(self.output_dir, exist_ok=True)
        with open(os.path.join(self.output_dir, "meme.txt"), "w") as f:
            f.write("MOTIF TGTGANNDWKTTCAC\n")
            f.write("MOTIF ACGTACGTACGTACG\n")
            f.write("MOTIF GTCAGTCAGTCAGTC\n")
            f.write("MOTIF CAGTCAGTCAGTCAG\n")

    def tearDown(self):
        if os.path.exists(self.output_dir):
            for file in os.listdir(self.output_dir):
                os.remove(os.path.join(self.output_dir, file))
            try:
                os.rmdir(self.output_dir)
            except PermissionError:
                pass

    def test_reshape_data_single_feature(self):
        input_data = [1, 2, 3]
        reshaped_data = reshape_data(input_data)
        self.assertEqual(reshaped_data.shape, (3, 1))

    def test_reshape_data_single_sample(self):
        input_data = [[1, 2, 3]]
        reshaped_data = reshape_data(input_data)
        self.assertEqual(reshaped_data.shape, (1, 3))

    def test_output_directory_exists(self):
        self.assertTrue(os.path.exists(self.output_dir))
        
    def test_postprocess_output(self):
        postprocess_output(self.output_dir, (2, 10))

if __name__ == '__main__':
    unittest.main()