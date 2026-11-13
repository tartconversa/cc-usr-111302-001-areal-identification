import json
import tempfile
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

from src.areal_analysis import (
    Cell,
    InputError,
    Problem,
    Subtotal,
    analyse,
    build_model,
    check_witness,
    load_problem,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "fixtures"
FULL = [FIXTURES / "regions.json", FIXTURES / "crosswalk.json", FIXTURES / "cells.csv"]


def table(witness):
    return {(item["old"], item["new"]): item["persons"] for item in witness}


def make_problem(old_totals, new_regions, cells, **kwargs):
    return Problem(
        version="synthetic-test",
        unit="persons",
        old_totals=old_totals,
        new_regions=new_regions,
        cells=cells,
        **kwargs,
    )


def reported(problem):
    """Analyse and cross-check every witness against the model independently."""
    report = analyse(problem)
    model = build_model(problem)
    for query in report["queries"]:
        for side in ("min", "max"):
            witness = query[side]["witness"]
            if witness is not None:
                values = [item["persons"] for item in witness]
                assert check_witness(model, values) == [], (query["label"], side)
    return report, model


def find(report, label):
    return next(q for q in report["queries"] if q["label"] == label)


def assert_range(test, report, label, low, high):
    query = find(report, label)
    test.assertEqual((query["range"]["min"], query["range"]["max"]), (low, high))
    test.assertTrue(query["range"]["integer_certified"])
    test.assertEqual(query["range"]["kind"], "identified_feasible_range")
    for side in ("min", "max"):
        test.assertIsNotNone(query[side]["witness"])
        test.assertFalse(query[side]["unverified"])
    return query


def write_inputs(directory, regions, crosswalk, cells_text):
    directory = Path(directory)
    (directory / "regions.json").write_text(json.dumps(regions), encoding="utf-8")
    (directory / "crosswalk.json").write_text(json.dumps(crosswalk), encoding="utf-8")
    (directory / "cells.csv").write_text(cells_text, encoding="utf-8")
    return (
        directory / "regions.json",
        directory / "crosswalk.json",
        directory / "cells.csv",
    )


NO_CAP_CELLS = lambda: [  # noqa: E731 - small helper for scenario tables
    Cell("A", "X", None, None),
    Cell("A", "Y", None, None),
    Cell("B", "Y", None, None),
]


class FixtureExampleTests(unittest.TestCase):
    """The worked example stated in the task, using the shipped fixtures."""

    def setUp(self):
        self.problem = load_problem(*FULL)
        self.report, self.model = reported(self.problem)

    def test_ranges_match_worked_example(self):
        self.assertTrue(self.report["feasibility"]["feasible"])
        self.assertTrue(self.report["feasibility"]["certified"])
        self.assertEqual(len(self.report["feasibility"]["witness"]), 3)
        for label, low, high in (("X", 4, 6), ("Y", 4, 6), ("X - Y", -2, 2)):
            assert_range(self, self.report, label, low, high)

    def test_endpoint_witnesses_are_joint_assignments(self):
        diff = find(self.report, "X - Y")
        minimum = table(diff["min"]["witness"])
        maximum = table(diff["max"]["witness"])

        def difference(values):
            return (
                values[("A", "X")]
                - values[("A", "Y")]
                - values[("B", "Y")]
            )

        self.assertEqual(difference(minimum), -2)
        self.assertEqual(difference(maximum), 2)
        # X and Y are never at 5/5: the subtotal is not redistributed evenly.
        for values in (minimum, maximum):
            self.assertEqual(values[("A", "X")] + values[("A", "Y")] + values[("B", "Y")], 10)

    def test_witnesses_are_full_cross_cell_tables(self):
        keys = sorted(cell.key for cell in self.problem.cells)
        for query in self.report["queries"]:
            for side in ("min", "max"):
                witness = query[side]["witness"]
                self.assertEqual(sorted((item["old"], item["new"]) for item in witness), keys)

    def test_public_subtotal_is_a_separate_constraint(self):
        inventory = self.report["meta"]["constraint_inventory"]
        self.assertEqual(
            [item["kind"] for item in inventory].count("public_subtotal"), 1
        )
        self.assertEqual(
            sum(self.problem.old_totals.values()),
            sum(item["rhs"] for item in inventory if item["kind"] == "old_total"),
        )
        # 10 appears exactly once as a row; it is not folded into old totals.
        self.assertEqual(
            [item["rhs"] for item in inventory if item["kind"] == "public_subtotal"], [10]
        )

    def test_meta_is_aggregate_only_and_not_a_confidence_interval(self):
        meta = self.report["meta"]
        self.assertTrue(meta["no_individual_records"])
        self.assertEqual(meta["unit"], "persons")
        self.assertEqual(meta["region_version"], "synthetic-v1")
        self.assertIn("不是抽样置信区间", meta["range_semantics"])
        self.assertFalse(meta["area_weights_used_for_feasible_ranges"] if "area_weights_used_for_feasible_ranges" in meta else False)

    def test_constraints_cite_their_sources(self):
        inventory = self.report["meta"]["constraint_inventory"]
        self.assertEqual(
            {item["kind"] for item in inventory},
            {"old_total", "capacity", "public_subtotal"},
        )
        for item in inventory:
            self.assertTrue(item["source"])
        capacity = next(item for item in inventory if item["kind"] == "capacity")
        self.assertIn("<= 2", capacity["text"])
        self.assertTrue(capacity["source"].endswith("#row 3"))

    def test_area_baseline_is_noncompliant_and_not_projected(self):
        area = self.report["area_baseline"]
        self.assertEqual(area["status"], "noncompliant")
        self.assertFalse(area["compliant"])
        self.assertTrue(area["violations"])
        self.assertTrue(area["not_projected"])
        self.assertFalse(area["area_weights_used_for_feasible_ranges"])
        self.assertTrue(area["feasible_range_independent_of_area_weights"])
        self.assertTrue(area["not_a_confidence_interval"])
        declared = area["scenarios"][0]
        self.assertEqual(declared["name"], "declared_area_weights")
        # The non-compliant table is displayed verbatim: A->Y is 3, not snapped to 2.
        self.assertEqual(table(declared["allocation"])[("A", "Y")], 3)
        self.assertEqual(declared["violations"][0]["kind"], "capacity")


class RangeScenarioTests(unittest.TestCase):
    def test_no_capacity_no_subtotal(self):
        report, _ = reported(make_problem({"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS()))
        assert_range(self, report, "X", 0, 6)
        assert_range(self, report, "Y", 4, 10)

    def test_capacity_two(self):
        cells = [Cell("A", "X", None, None), Cell("A", "Y", 2, None), Cell("B", "Y", None, None)]
        report, _ = reported(make_problem(
            {"A": 6, "B": 4}, ["X", "Y"], cells,
            queries=[{"X": 1}, {"Y": 1}, {"X": 1, "Y": -1}],
        ))
        assert_range(self, report, "X", 4, 6)
        assert_range(self, report, "Y", 4, 6)
        assert_range(self, report, "X - Y", -2, 2)

    def test_public_subtotal_does_not_split_evenly(self):
        plain, _ = reported(make_problem({"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS()))
        subtotal, _ = reported(make_problem(
            {"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS(),
            public_subtotals=[Subtotal("all", ("X", "Y"), 10)],
        ))
        # With X+Y=10 fixed, X and Y are complementary; the range is still 0..6 and 4..10.
        self.assertEqual(plain["queries"][0]["range"], subtotal["queries"][0]["range"])
        self.assertEqual(plain["queries"][1]["range"], subtotal["queries"][1]["range"])
        self.assertNotEqual(subtotal["queries"][0]["range"]["min"], 5)
        self.assertNotEqual(subtotal["queries"][0]["range"]["max"], 5)

    def test_public_subtotal_narrows_range(self):
        report, _ = reported(make_problem(
            {"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS(),
            public_subtotals=[Subtotal("x", ("X",), 5)],
        ))
        assert_range(self, report, "X", 5, 5)
        assert_range(self, report, "Y", 5, 5)

    def test_group_comparison_is_joint_not_marginal(self):
        # Marginals are W in [0,1] and Y in [1,9]; the naive independent
        # combination would claim W-Y in [-9,0], but no single table reaches 0.
        cells = [
            Cell("A", "X", None, None), Cell("A", "Y", None, None), Cell("A", "Z", 2, None),
            Cell("B", "W", None, None), Cell("B", "Y", None, None), Cell("B", "Z", None, None),
            Cell("C", "X", 2, None), Cell("C", "Y", 1, None), Cell("C", "Z", 3, None),
        ]
        report, _ = reported(make_problem(
            {"A": 4, "B": 4, "C": 3}, ["W", "X", "Y", "Z"], cells,
            public_subtotals=[Subtotal("yz", ("Y", "Z"), 10)],
            queries=[{"W": 1}, {"Y": 1}, {"W": 1, "Y": -1}],
        ))
        query = assert_range(self, report, "W - Y", -9, -1)
        warning = query["marginal_combination_warning"]
        self.assertIsNotNone(warning)
        self.assertEqual(warning["naive_marginal_bounds"], [-9, 0])
        self.assertEqual(warning["joint_bounds"], [-9, -1])
        values = table(query["max"]["witness"])
        self.assertEqual(
            values[("B", "W")] - values[("A", "Y")] - values[("B", "Y")] - values[("C", "Y")],
            -1,
        )

    def test_zero_person_totals(self):
        report, _ = reported(make_problem({"A": 0, "B": 0}, ["X", "Y"], NO_CAP_CELLS()))
        self.assertTrue(report["feasibility"]["feasible"])
        self.assertEqual(
            table(report["feasibility"]["witness"]),
            {("A", "X"): 0, ("A", "Y"): 0, ("B", "Y"): 0},
        )
        assert_range(self, report, "X", 0, 0)
        assert_range(self, report, "Y", 0, 0)

    def test_zero_capacity_is_a_hard_constraint(self):
        cells = [Cell("A", "X", 0, None), Cell("A", "Y", None, None), Cell("B", "Y", None, None)]
        report, _ = reported(make_problem({"A": 6, "B": 4}, ["X", "Y"], cells))
        assert_range(self, report, "X", 0, 0)
        self.assertEqual(table(report["feasibility"]["witness"])[("A", "X")], 0)

    def test_relaxation_mode_reports_uncertified_endpoints(self):
        report, _ = reported(make_problem(
            {"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS(), solve_integer=False,
        ))
        query = report["queries"][0]
        self.assertEqual((query["range"]["min"], query["range"]["max"]), (0, 6))
        self.assertFalse(query["range"]["integer_certified"])
        self.assertEqual(query["range"]["kind"], "linear_relaxation_range")
        self.assertIsNone(query["min"]["witness"])
        self.assertIsNotNone(query["min"]["relaxation_bound"])

    def test_node_budget_keeps_endpoints_unverified(self):
        # Integrality is required but the root optimum is fractional, so a
        # one-node budget leaves the endpoint uncertified (relaxation bound kept).
        cells = [
            Cell("A", "X", None, None), Cell("A", "W", None, None),
            Cell("B", "Y", None, None), Cell("B", "W", None, None),
            Cell("C", "Z", None, None), Cell("C", "W", None, None),
        ]
        report, _ = reported(make_problem(
            {"A": 1, "B": 1, "C": 1}, ["W", "X", "Y", "Z"], cells,
            public_subtotals=[
                Subtotal("xy", ("X", "Y"), 1),
                Subtotal("yz", ("Y", "Z"), 1),
                Subtotal("xz", ("X", "Z"), 1),
            ],
            node_budget=1,
        ))
        self.assertIsNone(report["feasibility"]["feasible"])
        self.assertFalse(report["feasibility"]["certified"])
        self.assertIsNone(report["feasibility"]["conflict"])
        for query in report["queries"]:
            self.assertFalse(query["range"]["integer_certified"])
            for side in ("min", "max"):
                self.assertTrue(query[side]["unverified"])
                self.assertIsNotNone(query[side]["relaxation_bound"])


class InfeasibilityTests(unittest.TestCase):
    def test_capacity_conflict_core_keeps_old_total(self):
        report, _ = reported(make_problem({"A": 6}, ["X"], [Cell("A", "X", 5, None)]))
        feasibility = report["feasibility"]
        self.assertFalse(feasibility["feasible"])
        self.assertTrue(feasibility["certified"])
        conflict = feasibility["conflict"]
        kinds = [item["kind"] for item in conflict["core"]]
        self.assertIn("old_total", kinds)
        self.assertIn("capacity", kinds)
        self.assertEqual(
            next(item for item in conflict["core"] if item["kind"] == "old_total")["rhs"], 6
        )
        self.assertTrue(conflict["old_totals_retained"])
        self.assertFalse(conflict["integrality_required"])
        self.assertEqual(report["queries"], [])

    def test_integrality_only_conflict(self):
        # LP relaxation is feasible (every variable at 1/2) but no integer table exists.
        cells = [
            Cell("A", "X", None, None), Cell("A", "W", None, None),
            Cell("B", "Y", None, None), Cell("B", "W", None, None),
            Cell("C", "Z", None, None), Cell("C", "W", None, None),
        ]
        report, _ = reported(make_problem(
            {"A": 1, "B": 1, "C": 1}, ["W", "X", "Y", "Z"], cells,
            public_subtotals=[
                Subtotal("xy", ("X", "Y"), 1),
                Subtotal("yz", ("Y", "Z"), 1),
                Subtotal("xz", ("X", "Z"), 1),
            ],
        ))
        self.assertFalse(report["feasibility"]["feasible"])
        conflict = report["feasibility"]["conflict"]
        self.assertTrue(conflict["integrality_required"])
        self.assertTrue(conflict["old_totals_retained"])
        self.assertTrue(any(item["kind"] == "old_total" for item in conflict["core"]))
        self.assertGreater(len(conflict["log"]), 0)

    def test_infeasible_report_has_no_fake_endpoints(self):
        report, _ = reported(make_problem({"A": 6}, ["X"], [Cell("A", "X", 5, None)]))
        self.assertEqual(report["queries"], [])
        self.assertEqual(report["feasibility"]["witness"], None)


class AreaSensitivityTests(unittest.TestCase):
    def _model(self, weight):
        cells = [
            Cell("A", "X", None, Fraction(weight)),
            Cell("A", "Y", None, Fraction(1) - Fraction(weight)),
            Cell("B", "Y", None, Fraction(1)),
        ]
        return make_problem({"A": 6, "B": 4}, ["X", "Y"], cells)

    def test_area_weight_changes_point_estimate_not_feasible_range(self):
        heavy_x, _ = reported(self._model("3/4"))
        light_x, _ = reported(self._model("1/4"))
        self.assertEqual(heavy_x["queries"][0]["range"], light_x["queries"][0]["range"])
        heavy_points = heavy_x["area_baseline"]["query_point_values"]
        light_points = light_x["area_baseline"]["query_point_values"]
        self.assertEqual(heavy_points[0]["declared_area_weights"], 5)
        self.assertEqual(light_points[0]["declared_area_weights"], 2)
        self.assertNotEqual(heavy_points, light_points)
        self.assertFalse(heavy_x["area_baseline"]["area_weights_used_for_feasible_ranges"])
        self.assertTrue(heavy_x["area_baseline"]["not_a_confidence_interval"])

    def test_compliant_area_baseline_is_labelled_as_point_estimate(self):
        report, _ = reported(self._model("3/4"))
        area = report["area_baseline"]
        self.assertEqual(area["status"], "compliant_point_estimate")
        self.assertTrue(area["compliant"])
        self.assertEqual(area["violations"], [])
        self.assertEqual(area["scenarios"][1]["name"], "uniform_over_allowed_cells")
        self.assertNotEqual(
            area["scenarios"][0]["allocation"], area["scenarios"][1]["allocation"]
        )

    def test_missing_weights_leave_baseline_unavailable_without_projection(self):
        report, _ = reported(make_problem({"A": 6, "B": 4}, ["X", "Y"], NO_CAP_CELLS()))
        area = report["area_baseline"]
        self.assertEqual(area["status"], "unavailable")
        self.assertIsNone(area["compliant"])
        self.assertFalse(area["scenarios"][0]["available"])
        self.assertTrue(area["scenarios"][0]["not_projected"])

    def test_area_baseline_mode_none_is_explicit(self):
        cells = [Cell("A", "X", None, None), Cell("A", "Y", 2, None), Cell("B", "Y", None, None)]
        report, _ = reported(make_problem(
            {"A": 6, "B": 4}, ["X", "Y"], cells, area_baseline="none",
        ))
        self.assertTrue(report["area_baseline"]["skipped"])
        self.assertEqual(report["area_baseline"]["mode"], "none")


class InputContractTests(unittest.TestCase):
    REGIONS = {"version": "v1", "unit": "persons", "old_regions": {"A": 1}, "new_regions": ["X"]}
    CELLS = "old,new,capacity,area_weight\nA,X,,\n"

    def _load(self, regions=None, crosswalk=None, cells=None):
        regions = regions or self.REGIONS
        crosswalk = crosswalk if crosswalk is not None else {"public_subtotals": []}
        cells = cells or self.CELLS
        with tempfile.TemporaryDirectory() as directory:
            return load_problem(*write_inputs(directory, regions, crosswalk, cells))

    def test_unit_mismatch_rejected(self):
        with self.assertRaisesRegex(InputError, "单位不匹配"):
            self._load(crosswalk={"unit": "households", "public_subtotals": []})

    def test_zoning_version_mismatch_rejected(self):
        with self.assertRaisesRegex(InputError, "区划版本不匹配"):
            self._load(crosswalk={"region_version": "v2", "public_subtotals": []})

    def test_personal_level_column_rejected(self):
        cells = "old,new,capacity,area_weight,person_id\nA,X,,,P001\n"
        with self.assertRaisesRegex(InputError, "个人级/定位列"):
            self._load(cells=cells)

    def test_geolocation_column_rejected(self):
        cells = "old,new,capacity,area_weight,latitude\nA,X,,,40.7\n"
        with self.assertRaisesRegex(InputError, "个人级/定位列"):
            self._load(cells=cells)

    def test_unknown_column_rejected(self):
        cells = "old,new,capacity,area_weight,secret_note\nA,X,,,hello\n"
        with self.assertRaisesRegex(InputError, "未在契约中声明"):
            self._load(cells=cells)

    def test_non_integer_total_rejected(self):
        regions = dict(self.REGIONS, old_regions={"A": 1.5})
        with self.assertRaisesRegex(InputError, "非负整数"):
            self._load(regions=regions)

    def test_negative_capacity_rejected(self):
        cells = "old,new,capacity,area_weight\nA,X,-1,\n"
        with self.assertRaisesRegex(InputError, "不能小于 0"):
            self._load(cells=cells)

    def test_unknown_cell_endpoint_rejected(self):
        cells = "old,new,capacity,area_weight\nZ,X,,\n"
        with self.assertRaisesRegex(InputError, "未在 regions 中声明"):
            self._load(cells=cells)

    def test_unsupported_area_baseline_mode_rejected(self):
        crosswalk = {"public_subtotals": [], "area_baseline": "project_then_call_it_areal"}
        with self.assertRaisesRegex(InputError, "不支持的 area_baseline"):
            self._load(crosswalk=crosswalk)


class SolverAgreementTests(unittest.TestCase):
    """Deterministic sweep against exhaustive integer enumeration."""

    def _enumerate(self, olds, cells, subtotals):
        groups = []
        for old, total in olds.items():
            indices = [i for i, cell in enumerate(cells) if cell.old == old]
            options = []
            for combo in product(range(total + 1), repeat=len(indices)):
                if sum(combo) != total:
                    continue
                if any(
                    cells[index].capacity is not None and combo[j] > cells[index].capacity
                    for j, index in enumerate(indices)
                ):
                    continue
                options.append((indices, combo))
            groups.append(options)
        feasible = []
        for selection in product(*groups):
            vector = [0] * len(cells)
            for indices, combo in selection:
                for j, index in enumerate(indices):
                    vector[index] = combo[j]
            if all(
                sum(vector[i] for i, cell in enumerate(cells) if cell.new in subtotal.regions)
                == subtotal.count
                for subtotal in subtotals
            ):
                feasible.append(vector)
        return feasible

    CASES = [
        (
            {"A": 6, "B": 4}, ["X", "Y"],
            [Cell("A", "X", None, None), Cell("A", "Y", 2, None), Cell("B", "Y", None, None)],
            [],
        ),
        (
            {"A": 2, "B": 2, "C": 1}, ["X", "Y", "Z"],
            [Cell("A", "X", 1, None), Cell("A", "Y", None, None), Cell("A", "Z", None, None),
             Cell("B", "X", None, None), Cell("B", "Y", 1, None), Cell("B", "Z", None, None),
             Cell("C", "X", None, None), Cell("C", "Y", None, None), Cell("C", "Z", 0, None)],
            [Subtotal("xy", ("X", "Y"), 3)],
        ),
        (
            {"A": 1, "B": 1}, ["X", "Y"],
            [Cell("A", "X", None, None), Cell("A", "Y", None, None),
             Cell("B", "X", None, None), Cell("B", "Y", None, None)],
            [Subtotal("x", ("X",), 1), Subtotal("y", ("Y",), 1)],
        ),
        (
            {"A": 1, "B": 1, "C": 1}, ["W", "X", "Y", "Z"],
            [Cell("A", "X", None, None), Cell("A", "W", None, None),
             Cell("B", "Y", None, None), Cell("B", "W", None, None),
             Cell("C", "Z", None, None), Cell("C", "W", None, None)],
            [Subtotal("xy", ("X", "Y"), 1), Subtotal("yz", ("Y", "Z"), 1),
             Subtotal("xz", ("X", "Z"), 1)],
        ),
        (
            {"A": 4, "B": 4, "C": 3}, ["W", "X", "Y", "Z"],
            [Cell("A", "X", None, None), Cell("A", "Y", None, None), Cell("A", "Z", 2, None),
             Cell("B", "W", None, None), Cell("B", "Y", None, None), Cell("B", "Z", None, None),
             Cell("C", "X", 2, None), Cell("C", "Y", 1, None), Cell("C", "Z", 3, None)],
            [Subtotal("yz", ("Y", "Z"), 10)],
        ),
    ]

    def test_ranges_and_feasibility_agree_with_brute_force(self):
        for olds, news, cells, subtotals in self.CASES:
            with self.subTest(olds=olds, news=news):
                reference = self._enumerate(olds, cells, subtotals)
                problem = make_problem(
                    olds, news, cells,
                    public_subtotals=subtotals,
                    queries=[{region: 1} for region in news]
                    + [{news[0]: 1, news[-1]: -1}],
                )
                report, _ = reported(problem)
                self.assertEqual(report["feasibility"]["feasible"], bool(reference))
                if not reference:
                    self.assertIsNotNone(report["feasibility"]["conflict"])
                    continue
                for query in report["queries"]:
                    coefficients = query["expression"]
                    values = [
                        sum(
                            coefficient
                            * sum(row[i] for i, cell in enumerate(cells) if cell.new == region)
                            for region, coefficient in coefficients.items()
                        )
                        for row in reference
                    ]
                    self.assertEqual(query["range"]["min"], min(values))
                    self.assertEqual(query["range"]["max"], max(values))


if __name__ == "__main__":
    unittest.main()
