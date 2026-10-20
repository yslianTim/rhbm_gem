"""Regression checks for the retained fixed-seed study data generator."""
import hashlib
import unittest

import numpy as np

import fixed_seed_joint_noise_generator as generator


class FixedSeedJointNoiseGeneratorTest(unittest.TestCase):
    def test_retained_seeded_noise_streams(self):
        expected = {
            'iid': 'b4975f62fded3e787d44433e81f6b29234e3df69e5150d30dbf50952099f6375',
            'correlated': '86780b5eca52108530d0778ecdcd266f65c03860356b7c50a75f2e4bc6debf76',
        }
        self.assertEqual(generator.SEEDS[0], 20260921)
        self.assertEqual(generator.SEEDS[-1], 20260940)
        for kind, digest in expected.items():
            values = generator.noise_field(generator.SEEDS[0], kind)
            self.assertEqual(hashlib.sha256(values.astype('<f8').tobytes()).hexdigest(), digest)
            np.testing.assert_array_equal(values, generator.noise_field(generator.SEEDS[0], kind))

    def test_historical_generated_observation_hash(self):
        clean = generator.clean_map(1.2, 0.)
        scale = float(np.sqrt(np.mean(clean[generator.TARGET] ** 2)))
        values = clean + .01 * scale * generator.noise_field(generator.SEEDS[0], 'iid')
        digest = hashlib.sha256(values.astype('<f8').tobytes()).hexdigest()
        self.assertEqual(digest, 'c5147ebda1fe72d1012b49d91838743ce4e8ca00abf4b3f8d3dcfa233fba02e1')


if __name__ == '__main__':
    unittest.main()
