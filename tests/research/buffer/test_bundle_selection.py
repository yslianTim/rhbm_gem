import unittest

import numpy as np

from bundle_selection import choose, greedy, scopes
from geometry import equal_radius, gaussian


class BundleSelectionTest(unittest.TestCase):
    def test_complete_bundle_shared_gain_and_budget(self):
        # Three incidences from two output atoms share one blocking set.
        benefits = {(2, 4): 3., (1,): 1.}
        for bundles in ([(2, 4), (1,)], [(1,), (2, 4)]):
            self.assertEqual(greedy(bundles, benefits, 1, lambda _: True)[0], {1})
            self.assertEqual(greedy(bundles, benefits, 2, lambda _: True)[0], {2, 4})
            self.assertEqual(greedy(bundles, benefits, 3, lambda _: True)[0], {1, 2, 4})

    def test_shared_voxel_updates_both_output_targets(self):
        class Fixture:
            output=[1,2,3]
            retained=np.array([False,True,True,True])
            blockers={0:(10,11)}
            incidences={1:[(0,.09),(1,.52)],2:[(0,.24),(2,.81)],3:[(3,1.)],
                        10:[(0,.56),(1,.88),(2,1.)],11:[(0,.64),(1,1.),(3,2.)]}

            def atom_support(self, i):
                return np.array(self.incidences[i],dtype=[("row","u4"),("square","f8")])

            def rows_retained(self, native, selected):
                return np.array([all(i in selected for i in self.blockers.get(r,())) for r in native])

        bad=choose(Fixture(),.5,1,"G2b")
        good=choose(Fixture(),.5,2,"G2b")
        self.assertEqual(bad["buffer"],[])
        self.assertEqual(good["buffer"],[10,11])
        self.assertEqual(good["recovered_rows"],1)
        self.assertEqual(set(good["affected"]),{"1","2"})
        self.assertEqual(sum(s["core_recovered"] for s in good["affected"].values()),2)

    def test_rejects_unsupported_and_zero_benefit(self):
        found,steps,stop = greedy([(1,), (2,)], {(1,): 1., (2,): 0.}, 2, lambda _: False)
        self.assertEqual(found, set())
        self.assertEqual(steps, [])
        self.assertEqual(stop, "no-positive-feasible-gain")

    def test_equal_radius_keeps_machine_roundoff_ties(self):
        squares = np.array([.25, np.nextafter(.25, 1), .2500000001])
        np.testing.assert_array_equal(equal_radius(squares, .25), [True,True,False])

    def test_nearest_strategy_keeps_atoms_with_no_observed_core(self):
        class Fixture:
            output = [1]
            retained = np.array([False, True])

            def atom_support(self, _):
                return np.array([(0,.36),(1,.64)],dtype=[("row","u4"),("square","f8")])

        self.assertEqual(set(scopes(Fixture(), .5, None, "G2a")), {1})
        self.assertEqual(scopes(Fixture(), .5, None, "G2b"), {})

    def test_gaussian_proxy_references_the_same_original_rows(self):
        squares=np.array([.25,.25,.5,1.])
        mask=np.array([False,False,True,True])
        result=gaussian(squares,mask,1e-200)
        self.assertEqual(result["retention"],0.)
        self.assertEqual(result["original_effective"],2.)
        self.assertEqual(result["fitting_effective"],1.)

    def test_canonical_tie_for_reordered_bundles(self):
        for bundles in ([(8,), (2,)], [(2,), (8,)]):
            actual, _, _ = greedy(bundles, {(8,):1.,(2,):1.}, 1, lambda _: True)
            self.assertEqual(actual,{2})


if __name__ == "__main__":
    unittest.main()
