"""Areal (zone-system) conversion identification ranges.

Two notions are kept strictly apart:

* the **area baseline** - one deterministic apportionment of old-region
  totals along area weights.  It is an independent point scenario only and is
  never used to shrink or reinterpret the feasible set.  If it violates a
  capacity or a published subtotal it is reported as non-compliant; it is
  *not* projected onto the constraints and then still called an area baseline.
* the **feasible set** - every non-negative integer distribution of the old
  totals over the permitted old/new intersections that respects per-cell
  capacity upper bounds and any explicitly published public new-region
  subtotals.  No individual-level record is read or emitted.

Ranges are exact integer minima/maxima of a directly optimised linear
expression, each accompanied by a complete cross-cell witness.  A group
comparison (for example X-Y) is solved as one joint expression; independent
marginal ranges are never combined into a bound that no single assignment
attains.  If the integer search cannot certify an endpoint within the node
budget, the exact linear-relaxation bound is reported and the endpoint stays
unverified.  Results are identified feasible ranges, never sampling
confidence intervals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Sequence

MAXIMISE = "max"
MINIMISE = "min"

VERSION_KEYS = ("region_version", "zoning_version", "regions_version")
UNIT_KEYS = ("unit", "units", "region_unit")
ALLOWED_BASELINE_MODES = ("separate_without_projection", "independent", "none")
CELL_COLUMNS = {"old", "new", "capacity", "area_weight"}
# Columns that would indicate person-level or geolocated records.  Their
# presence is a contract violation: this analyser reads aggregates only.
FORBIDDEN_COLUMNS = {
    "person_id", "individual_id", "record_id", "household_id", "name",
    "address", "latitude", "longitude", "lat", "lon", "x_coord", "y_coord",
    "geom", "geometry", "wkt", "phone", "email",
}


# ---------------------------------------------------------------------------
# Errors and small value objects


class InputError(ValueError):
    """Input contradicts itself or cannot be interpreted safely."""


@dataclass(frozen=True)
class Cell:
    old: str
    new: str
    capacity: int | None
    area_weight: Fraction | None
    source: str = ""

    @property
    def key(self) -> tuple[str, str]:
        return (self.old, self.new)


@dataclass(frozen=True)
class Subtotal:
    key: str
    regions: tuple[str, ...]
    count: int
    source: str = ""


@dataclass
class Problem:
    version: str
    unit: str
    old_totals: dict[str, int]
    new_regions: list[str]
    cells: list[Cell]
    public_subtotals: list[Subtotal] = field(default_factory=list)
    queries: list[dict[str, int]] = field(default_factory=list)
    solve_integer: bool = True
    node_budget: int = 20000
    area_baseline: str = "separate_without_projection"
    source_labels: dict = field(default_factory=dict)

    def __post_init__(self):
        self.cells = sorted(self.cells, key=lambda cell: (cell.old, cell.new))
        self.new_regions = list(self.new_regions)
        self.old_totals = dict(self.old_totals)
        if not self.queries:
            self.queries = [{region: 1} for region in self.new_regions]


@dataclass(frozen=True)
class Constraint:
    """One mandatory aggregate constraint, with its source for auditing."""

    kind: str            # old_total | capacity | public_subtotal
    coeffs: tuple[int, ...]
    rhs: int
    text: str
    source: str
    detail: dict
    op: str = "="

    def row(self):
        return list(self.coeffs), self.op, self.rhs


@dataclass
class Model:
    problem: Problem
    cells: list[Cell]
    variables: dict[tuple[str, str], int]
    constraints: list[Constraint]

    @property
    def nvars(self):
        return len(self.cells)


# ---------------------------------------------------------------------------
# Exact two-phase simplex over the rationals


def _pivot(rows, basis, row_index, column):
    pivot = rows[row_index][column]
    rows[row_index] = [value / pivot for value in rows[row_index]]
    for i, row in enumerate(rows):
        if i == row_index:
            continue
        factor = row[column]
        if factor:
            target = rows[row_index]
            rows[i] = [value - factor * other for value, other in zip(row, target)]
    basis[row_index] = column


def _simplex(rows, basis, costs, forbidden, ncols, max_steps=200000):
    """Maximise costs.x from a feasible basis, Bland's rule for anti-cycling."""
    steps = 0
    while True:
        steps += 1
        if steps > max_steps:
            raise InputError("simplex 步数超限")
        basic = set(basis)
        entering = None
        for column in range(ncols):
            if column in basic or column in forbidden:
                continue
            reduced = costs[column] - sum(
                costs[basis[i]] * rows[i][column] for i in range(len(rows))
            )
            if reduced > 0:
                entering = column
                break
        if entering is None:
            value = sum(costs[basis[i]] * rows[i][-1] for i in range(len(rows)))
            return "optimal", value
        leaving = None
        best_ratio = None
        for i, row in enumerate(rows):
            if row[entering] > 0:
                ratio = row[-1] / row[entering]
                if (
                    best_ratio is None
                    or ratio < best_ratio
                    or (ratio == best_ratio and basis[i] < basis[leaving])
                ):
                    best_ratio, leaving = ratio, i
        if leaving is None:
            return "unbounded", None
        _pivot(rows, basis, leaving, entering)


def solve_lp(nvars, constraints, cost, sense):
    """Solve min/max cost.x subject to constraints, x >= 0, exact rationals.

    ``constraints`` is a sequence of ``(coeffs, op, rhs)`` with op in <=, >=, =.
    Returns ``(status, value, values)`` with status optimal/infeasible/unbounded.
    Verdicts are exact for the rational relaxation and therefore exact bounds
    for the integer problem.
    """
    normalised = []
    for coeffs, op, rhs in constraints:
        coeffs = list(coeffs)
        if rhs < 0:
            coeffs = [-value for value in coeffs]
            rhs = -rhs
            op = {"<=": ">=", ">=": "<=", "=": "="}[op]
        normalised.append((coeffs, op, rhs))
    m = len(normalised)
    artificial_base = nvars + m
    artificial_count = sum(1 for _, op, _ in normalised if op in (">=", "="))
    ncols = artificial_base + artificial_count
    rows = []
    basis = []
    artificials = set()
    artificial_cursor = artificial_base
    for i, (coeffs, op, rhs) in enumerate(normalised):
        row = [Fraction(0)] * (ncols + 1)
        for j, value in enumerate(coeffs):
            row[j] = Fraction(value)
        if op == "<=":
            slack = nvars + i
            row[slack] = Fraction(1)
            basis.append(slack)
        elif op == ">=":
            surplus = nvars + i
            row[surplus] = Fraction(-1)
            row[artificial_cursor] = Fraction(1)
            artificials.add(artificial_cursor)
            basis.append(artificial_cursor)
            artificial_cursor += 1
        elif op == "=":
            row[artificial_cursor] = Fraction(1)
            artificials.add(artificial_cursor)
            basis.append(artificial_cursor)
            artificial_cursor += 1
        else:
            raise InputError(f"未知约束方向: {op!r}")
        row[-1] = Fraction(rhs)
        rows.append(row)

    # Phase I: maximise -sum(artificials).
    phase1 = [Fraction(0)] * ncols
    for index in artificials:
        phase1[index] = Fraction(-1)
    status, objective = _simplex(rows, basis, phase1, set(), ncols)
    if status == "unbounded":
        raise InputError("phase I 无界；模型不合法")
    if objective < 0:
        return "infeasible", None, None

    # Drive artificials out; drop redundant zero rows.  Pivots rewrite every
    # row in place, so keep only indices until all pivots are done.
    drop = set()
    for i in list(range(len(rows))):
        if basis[i] not in artificials:
            continue
        if rows[i][-1] != 0:
            return "infeasible", None, None
        replacement = next(
            (column for column in range(ncols)
             if column not in artificials and rows[i][column] != 0),
            None,
        )
        if replacement is None:
            drop.add(i)
        else:
            _pivot(rows, basis, i, replacement)
    if drop:
        # Mutate in place: the caller holds references to both lists.
        rows[:] = [row for i, row in enumerate(rows) if i not in drop]
        basis[:] = [index for i, index in enumerate(basis) if i not in drop]

    # Phase II: original objective, artificial columns forbidden.
    if sense == MAXIMISE:
        phase2 = [Fraction(value) for value in cost] + [Fraction(0)] * (ncols - nvars)
    elif sense == MINIMISE:
        phase2 = [Fraction(-value) for value in cost] + [Fraction(0)] * (ncols - nvars)
    else:
        raise InputError(f"未知优化方向: {sense}")
    status, objective = _simplex(rows, basis, phase2, artificials, ncols)
    if status == "unbounded":
        return "unbounded", None, None
    if sense == MINIMISE:
        objective = -objective
    values = [Fraction(0)] * nvars
    for i, index in enumerate(basis):
        if index < nvars:
            values[index] = rows[i][-1]
    return "optimal", objective, values


# ---------------------------------------------------------------------------
# Input loading and validation


def _require_int(value, label, *, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int):
        raise InputError(f"{label} 必须是非负整数，实际为 {value!r}")
    if value < minimum:
        raise InputError(f"{label} 不能小于 {minimum}，实际为 {value}")
    return value


def _fraction(value, label):
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        parsed = Fraction(str(value).strip())
    except (ValueError, ZeroDivisionError) as exc:
        raise InputError(f"{label} 不是可解析的面积权重: {value!r}") from exc
    if parsed < 0:
        raise InputError(f"{label} 面积权重不能为负: {value!r}")
    return parsed


def _check_declarations(container, version, unit, source):
    """Reject explicit zoning-version or unit declarations that disagree."""
    if not isinstance(container, dict):
        return
    for key in VERSION_KEYS:
        if key in container and container[key] not in (None, ""):
            if str(container[key]) != str(version):
                raise InputError(
                    f"区划版本不匹配：{source} 的 {key}={container[key]!r}，"
                    f"regions 声明 {version!r}"
                )
    for key in UNIT_KEYS:
        if key in container and container[key] not in (None, ""):
            if str(container[key]) != str(unit):
                raise InputError(
                    f"单位不匹配：{source} 的 {key}={container[key]!r}，"
                    f"regions 声明 {unit!r}"
                )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_problem(regions_path, crosswalk_path, cells_path) -> Problem:
    """Read and validate aggregate inputs.  No personal data is read."""
    from .area_reader import load

    regions_path, crosswalk_path, cells_path = map(
        Path, (regions_path, crosswalk_path, cells_path)
    )
    regions = load(regions_path)
    crosswalk = load(crosswalk_path)
    rows = load(cells_path)
    if not isinstance(regions, dict):
        raise InputError("regions.json 必须是对象")
    if not isinstance(crosswalk, dict):
        raise InputError("crosswalk.json 必须是对象")
    if not isinstance(rows, list):
        raise InputError("cells.csv 必须是表格行列表")

    version = regions.get("version")
    unit = regions.get("unit")
    if not isinstance(version, str) or not version:
        raise InputError("regions.version 必须是非空字符串")
    if not isinstance(unit, str) or not unit:
        raise InputError("regions.unit 必须是非空字符串")
    _check_declarations(crosswalk, version, unit, str(crosswalk_path))

    old_totals_in = regions.get("old_regions")
    new_regions_in = regions.get("new_regions")
    if not isinstance(old_totals_in, dict) or not old_totals_in:
        raise InputError("regions.old_regions 必须是非空对象")
    if not isinstance(new_regions_in, list) or not new_regions_in:
        raise InputError("regions.new_regions 必须是非空数组")
    new_regions = []
    for region in new_regions_in:
        if not isinstance(region, str) or not region:
            raise InputError("new_regions 必须是非空字符串")
        if region in new_regions:
            raise InputError(f"new_regions 重复: {region}")
        new_regions.append(region)
    old_totals = {}
    for old, total in old_totals_in.items():
        if not isinstance(old, str) or not old:
            raise InputError("old_regions 的键必须是非空字符串")
        old_totals[old] = _require_int(total, f"regions.old_regions.{old}")

    cells = []
    seen = set()
    for index, row in enumerate(rows, start=2):  # header is line 1
        if not isinstance(row, dict) or not row.get("old") or not row.get("new"):
            raise InputError(f"cells.csv 第 {index} 行缺少 old/new")
        columns = {str(name).strip().lower() for name in row if name is not None}
        forbidden = columns & FORBIDDEN_COLUMNS
        if forbidden:
            raise InputError(
                f"cells.csv 包含疑似个人级/定位列 {sorted(forbidden)}；"
                "分析器只接受汇总交叉格资料"
            )
        _check_declarations(row, version, unit, f"cells.csv 第 {index} 行")
        unknown = columns - CELL_COLUMNS - set(UNIT_KEYS) - set(VERSION_KEYS)
        if unknown:
            raise InputError(
                f"cells.csv 含有未在契约中声明的列 {sorted(unknown)}；"
                "请显式扩充契约后再使用（避免把个人字段当汇总读入）"
            )
        old, new = row["old"].strip(), row["new"].strip()
        if old not in old_totals:
            raise InputError(f"cells.csv 第 {index} 行的旧区 {old!r} 未在 regions 中声明")
        if new not in new_regions:
            raise InputError(f"cells.csv 第 {index} 行的新区 {new!r} 未在 regions 中声明")
        if (old, new) in seen:
            raise InputError(f"cells.csv 重复的交叉格 ({old},{new})")
        seen.add((old, new))
        capacity_text = (row.get("capacity") or "").strip()
        capacity = None
        if capacity_text:
            try:
                capacity = int(capacity_text)
            except ValueError as exc:
                raise InputError(
                    f"cells.csv 第 {index} 行容量不是整数: {row.get('capacity')!r}"
                ) from exc
            capacity = _require_int(capacity, f"cells.csv 第 {index} 行容量")
        weight = _fraction(row.get("area_weight"), f"cells.csv 第 {index} 行面积权重")
        _check_declarations(row, version, unit, f"cells.csv 第 {index} 行")
        cells.append(Cell(old, new, capacity, weight, source=f"{cells_path}#row {index}"))

    for old, total in old_totals.items():
        if total and old not in {cell.old for cell in cells}:
            raise InputError(f"旧区 {old!r} 有 {total} 人但没有任何允许的交叉格")
    for new in new_regions:
        if new not in {cell.new for cell in cells}:
            raise InputError(f"新区 {new!r} 没有任何允许的交叉格")

    subtotals = []
    seen_keys = set()
    for index, entry in enumerate(crosswalk.get("public_subtotals", []) or []):
        if not isinstance(entry, dict):
            raise InputError(f"crosswalk.public_subtotals[{index}] 必须是对象")
        key = entry.get("key")
        if not isinstance(key, str) or not key:
            raise InputError(f"crosswalk.public_subtotals[{index}].key 必须是非空字符串")
        if key in seen_keys:
            raise InputError(f"public_subtotals 的 key 重复: {key}")
        seen_keys.add(key)
        regions_value = entry.get("regions")
        if not isinstance(regions_value, list) or not regions_value:
            raise InputError(f"public_subtotals[{index}].regions 必须是非空数组")
        members = []
        for region in regions_value:
            if region not in new_regions:
                raise InputError(f"public_subtotals[{index}] 引用未知新区 {region!r}")
            if region in members:
                raise InputError(f"public_subtotals[{index}] 重复新区 {region!r}")
            members.append(region)
        count = _require_int(entry.get("count"), f"public_subtotals[{index}].count")
        subtotals.append(
            Subtotal(key, tuple(members), count,
                     source=f"{crosswalk_path}#public_subtotals[{index}]")
        )

    queries_in = crosswalk.get("queries")
    queries = []
    if queries_in is not None:
        if not isinstance(queries_in, list) or not queries_in:
            raise InputError("crosswalk.queries 必须是非空数组（或省略以取默认单区查询）")
        for index, entry in enumerate(queries_in):
            if not isinstance(entry, dict) or not entry:
                raise InputError(f"crosswalk.queries[{index}] 必须是非空对象")
            expression = {}
            for region, coefficient in entry.items():
                if region not in new_regions:
                    raise InputError(f"crosswalk.queries[{index}] 引用未知新区 {region!r}")
                if isinstance(coefficient, bool) or not isinstance(coefficient, int):
                    raise InputError(f"crosswalk.queries[{index}][{region!r}] 系数必须是整数")
                expression[region] = coefficient
            if not any(expression.values()):
                raise InputError(f"crosswalk.queries[{index}] 至少需要一个非零系数")
            queries.append(expression)

    solve_integer = crosswalk.get("solve_integer", True)
    if not isinstance(solve_integer, bool):
        raise InputError("crosswalk.solve_integer 必须是布尔值")
    budget = _require_int(crosswalk.get("node_budget", 20000), "crosswalk.node_budget", minimum=1)
    baseline_mode = crosswalk.get("area_baseline", "separate_without_projection")
    if not isinstance(baseline_mode, str) or not baseline_mode:
        raise InputError("crosswalk.area_baseline 必须是非空字符串")
    if baseline_mode not in ALLOWED_BASELINE_MODES:
        raise InputError(
            f"不支持的 area_baseline={baseline_mode!r}；允许 {ALLOWED_BASELINE_MODES}。"
            "面积基准不会被投影到约束后再冒称面积摊分。"
        )

    return Problem(
        version=version,
        unit=unit,
        old_totals=old_totals,
        new_regions=new_regions,
        cells=cells,
        public_subtotals=subtotals,
        queries=queries,
        solve_integer=solve_integer,
        node_budget=budget,
        area_baseline=baseline_mode,
        source_labels={
            "regions": str(regions_path),
            "crosswalk": str(crosswalk_path),
            "cells": str(cells_path),
            "sha256": {
                "regions": _sha256(regions_path),
                "crosswalk": _sha256(crosswalk_path),
                "cells": _sha256(cells_path),
            },
        },
    )


# ---------------------------------------------------------------------------
# Model and witness checking


def build_model(problem: Problem) -> Model:
    cells = sorted(problem.cells, key=lambda cell: (cell.old, cell.new))
    variables = {cell.key: index for index, cell in enumerate(cells)}
    if len(variables) != len(cells):
        raise InputError("交叉格重复")
    constraints = []
    for old, total in problem.old_totals.items():
        coeffs = tuple(1 if cell.old == old else 0 for cell in cells)
        terms = [f"x[{cell.old}->{cell.new}]" for cell in cells if cell.old == old]
        constraints.append(Constraint(
            kind="old_total",
            coeffs=coeffs,
            rhs=total,
            text=f"{' + '.join(terms)} = {total}",
            source=problem.source_labels.get("regions", "regions.json") + f"#old_regions.{old}",
            detail={"old": old},
        ))
    for index, cell in enumerate(cells):
        if cell.capacity is None:
            continue
        coeffs = tuple(1 if j == index else 0 for j in range(len(cells)))
        constraints.append(Constraint(
            kind="capacity",
            coeffs=coeffs,
            rhs=cell.capacity,
            text=f"x[{cell.old}->{cell.new}] <= {cell.capacity}",
            source=cell.source,
            detail={"old": cell.old, "new": cell.new},
            op="<=",
        ))
    for subtotal in problem.public_subtotals:
        coeffs = tuple(1 if cell.new in subtotal.regions else 0 for cell in cells)
        terms = [f"x[{cell.old}->{cell.new}]" for cell in cells if cell.new in subtotal.regions]
        constraints.append(Constraint(
            kind="public_subtotal",
            coeffs=coeffs,
            rhs=subtotal.count,
            text=f"{' + '.join(terms)} = {subtotal.count}",
            source=subtotal.source or f"public_subtotals[{subtotal.key}]",
            detail={"key": subtotal.key, "regions": list(subtotal.regions)},
        ))
    return Model(problem, cells, variables, constraints)


def check_witness(model: Model, values: Sequence[int], constraints=None) -> list[str]:
    """Independently re-check an integer witness against declared constraints.

    ``constraints`` defaults to the full model; passing a subset is used by the
    conflict-core diagnostic so a candidate witness is judged only against the
    subset currently under test.
    """
    problems = []
    if len(values) != len(model.cells):
        return ["见证长度与交叉格数不一致"]
    active_capacities = None
    if constraints is not None:
        active_capacities = {
            (constraint.detail["old"], constraint.detail["new"])
            for constraint in constraints
            if constraint.kind == "capacity"
        }
    for index, value in enumerate(values):
        label = f"{model.cells[index].old}->{model.cells[index].new}"
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            problems.append(f"{label} 不是非负整数: {value!r}")
            continue
        capacity = model.cells[index].capacity
        checked = active_capacities is None or model.cells[index].key in active_capacities
        if checked and capacity is not None and value > capacity:
            problems.append(f"{label}={value} 超过容量上界 {capacity}")
    for constraint in (model.constraints if constraints is None else constraints):
        if constraint.kind == "capacity":
            continue
        actual = sum(coeff * value for coeff, value in zip(constraint.coeffs, values))
        if actual != constraint.rhs:
            problems.append(f"{constraint.text} 的见证合计为 {actual}")
    return problems


def _variable_upper_bound(model: Model, index: int) -> int:
    cell = model.cells[index]
    if cell.capacity is not None:
        return cell.capacity
    return model.problem.old_totals[cell.old]


def _subset_upper_bounds(model: Model, constraints):
    """Bounds that are valid for the *subset* under test.

    If a capacity row was removed for a diagnostic, its bound must not be
    silently reimposed; otherwise a conflict core could blame a constraint it
    never actually tested.
    """
    active = set()
    for constraint in constraints:
        if constraint.kind == "capacity":
            active.add((constraint.detail["old"], constraint.detail["new"]))
    global_bound = sum(model.problem.old_totals.values())
    bounds = []
    for index, cell in enumerate(model.cells):
        if cell.capacity is not None and cell.key in active:
            bounds.append(cell.capacity)
        elif any(
            constraint.kind == "old_total" and constraint.detail["old"] == cell.old
            for constraint in constraints
        ):
            bounds.append(model.problem.old_totals[cell.old])
        else:
            bounds.append(global_bound)
    return bounds


def _lp_solution(model: Model, cost, sense, extra_rows=()):
    rows = [constraint.row() for constraint in model.constraints]
    rows.extend(extra_rows)
    for index in range(model.nvars):
        coeffs = [1 if j == index else 0 for j in range(model.nvars)]
        rows.append((coeffs, "<=", _variable_upper_bound(model, index)))
    return solve_lp(model.nvars, rows, cost, sense)


def _constraints_feasible_lp(model: Model, constraints):
    bounds = _subset_upper_bounds(model, constraints)
    rows = [constraint.row() for constraint in constraints]
    for index in range(model.nvars):
        coeffs = [1 if j == index else 0 for j in range(model.nvars)]
        rows.append((coeffs, "<=", bounds[index]))
    status, _, _ = solve_lp(model.nvars, rows, [0] * model.nvars, MAXIMISE)
    return status != "infeasible"


# ---------------------------------------------------------------------------
# Branch and bound for exact integer endpoints


class _BudgetExceeded(Exception):
    pass


@dataclass
class SearchOutcome:
    status: str                 # optimal | optimal_relaxation | infeasible | budget
    value: Fraction | None
    vector: list[int] | None
    vertices: int
    lp_bound: Fraction | None
    certified: bool


def _integer_search(model: Model, cost, sense) -> SearchOutcome:
    budget = model.problem.node_budget
    if not model.problem.solve_integer:
        status, value, vector = _lp_solution(model, cost, sense)
        if status == "infeasible":
            return SearchOutcome("infeasible", None, None, 0, None, True)
        if status == "unbounded":
            raise InputError("可行域无界；输入缺少必要的容量上界")
        return SearchOutcome("optimal_relaxation", value, None, 0, value, False)

    vertices = [0]
    incumbent = {"value": None, "vector": None}

    def improves(value):
        if incumbent["value"] is None:
            return True
        return value > incumbent["value"] if sense == MAXIMISE else value < incumbent["value"]

    def search(extra_rows):
        if vertices[0] >= budget:
            raise _BudgetExceeded
        vertices[0] += 1
        status, value, vector = _lp_solution(model, cost, sense, extra_rows)
        if status == "unbounded":
            raise InputError("可行域无界；输入缺少必要的容量上界")
        if status == "infeasible":
            return
        if incumbent["value"] is not None:
            if sense == MAXIMISE and value <= incumbent["value"]:
                return
            if sense == MINIMISE and value >= incumbent["value"]:
                return
        fractional = next(
            (j for j, x in enumerate(vector) if x.denominator != 1), None
        )
        if fractional is None:
            integer_vector = [int(x) for x in vector]
            problems = check_witness(model, integer_vector)
            if problems:
                raise InputError("内部见证校验失败：" + "；".join(problems))
            if improves(value):
                incumbent["value"] = value
                incumbent["vector"] = integer_vector
            return
        floor_value = vector[fractional].numerator // vector[fractional].denominator
        coeffs = [1 if j == fractional else 0 for j in range(model.nvars)]
        rows_high = [(coeffs, ">=", floor_value + 1)]
        rows_low = [(coeffs, "<=", floor_value)]
        if sense == MAXIMISE:
            search(extra_rows + rows_high)
            search(extra_rows + rows_low)
        else:
            search(extra_rows + rows_low)
            search(extra_rows + rows_high)

    status, root_bound, _ = _lp_solution(model, cost, sense)
    if status == "infeasible":
        return SearchOutcome("infeasible", None, None, 1, None, True)
    if status == "unbounded":
        raise InputError("可行域无界；输入缺少必要的容量上界")
    try:
        search([])
    except _BudgetExceeded:
        # If the incumbent already attains the exact root relaxation bound, it
        # is provably optimal even though the tree was cut short.
        value = incumbent["value"]
        proven = value is not None and root_bound is not None and value == root_bound
        return SearchOutcome(
            "optimal" if proven else "budget",
            value,
            incumbent["vector"],
            vertices[0],
            root_bound,
            proven,
        )
    if incumbent["vector"] is None:
        return SearchOutcome("infeasible", None, None, vertices[0], root_bound, True)
    return SearchOutcome("optimal", incumbent["value"], incumbent["vector"],
                         vertices[0], root_bound, True)


def _integer_feasibility(model: Model, constraints):
    """Return (feasible, certified): True/False/None with exactness flag."""
    if not _constraints_feasible_lp(model, constraints):
        return False, True
    if not model.problem.solve_integer:
        return None, False
    vertices = [0]
    budget = model.problem.node_budget
    rows = [constraint.row() for constraint in constraints]
    bounds = _subset_upper_bounds(model, constraints)
    for index in range(model.nvars):
        coeffs = [1 if j == index else 0 for j in range(model.nvars)]
        rows.append((coeffs, "<=", bounds[index]))

    def search(extra_rows):
        if vertices[0] >= budget:
            raise _BudgetExceeded
        vertices[0] += 1
        status, _, vector = solve_lp(
            model.nvars, rows + list(extra_rows), [0] * model.nvars, MAXIMISE
        )
        if status == "infeasible":
            return False
        fractional = next((j for j, x in enumerate(vector) if x.denominator != 1), None)
        if fractional is None:
            integer_vector = [int(x) for x in vector]
            return not check_witness(model, integer_vector, constraints)
        floor_value = vector[fractional].numerator // vector[fractional].denominator
        coeffs = [1 if j == fractional else 0 for j in range(model.nvars)]
        return (
            search(extra_rows + [(coeffs, "<=", floor_value)])
            or search(extra_rows + [(coeffs, ">=", floor_value + 1)])
        )

    try:
        return bool(search([])), True
    except _BudgetExceeded:
        return None, False


# ---------------------------------------------------------------------------
# Conflict core extraction


def _describe(constraint: Constraint):
    return {
        "kind": constraint.kind,
        "text": constraint.text,
        "source": constraint.source,
        "detail": constraint.detail,
        "rhs": constraint.rhs,
    }


def _conflict_core(model: Model):
    """Deletion-filter core over capacities and public subtotals.

    Old-region totals are never candidates for removal: they appear in every
    reported core, and the corresponding rows are enforced in the full model
    for every range and witness.  The diagnostic only asks which of the other
    declared constraints still produce infeasibility on top of those totals.
    """
    feasible, certified = _integer_feasibility(model, model.constraints)
    if feasible is not False:
        return {
            "core": [],
            "certified": bool(certified),
            "integrality_checked": certified,
            "note": "未认证冲突：连续松弛可行但整数可行性未在预算内判定。"
            if feasible is None else "未发现不可行分配。",
            "log": [],
            "old_totals_retained": True,
        }

    core = list(model.constraints)
    log = []
    for candidate in [item for item in core if item.kind != "old_total"]:
        trial = [item for item in core if item is not candidate]
        value, check_certified = _integer_feasibility(model, trial)
        log.append({
            "removed": _describe(candidate),
            "still_infeasible": value,
            "integer_checked": check_certified,
            "diagnostic_only": True,
        })
        if value is False:
            core = trial
    for item in core:
        if item.kind == "old_total":
            log.append({
                "kept": _describe(item),
                "reason": "旧区总量从不作为删除候选；全模型与冲突检查均未放宽它。",
            })
    lp_relaxation_feasible = _constraints_feasible_lp(model, core)
    return {
        "core": [_describe(constraint) for constraint in core],
        "certified": True,
        "integrality_checked": True,
        "integrality_required": lp_relaxation_feasible,
        "old_totals_retained": True,
        "note": (
            "该约束子集（含全部旧区总量）本身可核对为整数不可行；"
            "未为求解而放宽任何旧区总量。"
            if lp_relaxation_feasible
            else "该约束子集（含全部旧区总量）在线性松弛下即不可行；未放宽任何旧区总量。"
        ),
        "log": log,
    }


# ---------------------------------------------------------------------------
# Area baseline: an independent point scenario only


def _largest_remainder(total: int, entries):
    weight_sum = sum(weight for _, weight in entries)
    if weight_sum == 0:
        return None
    quotas = [(cell, Fraction(total) * weight / weight_sum) for cell, weight in entries]
    floors = {cell.key: quota.numerator // quota.denominator for cell, quota in quotas}
    remainder = total - sum(floors.values())
    order = sorted(
        quotas,
        key=lambda item: (-(item[1] - floors[item[0].key]), item[0].key[0], item[0].key[1]),
    )
    allocation = dict(floors)
    for cell, _ in order[:remainder]:
        allocation[cell.key] += 1
    return allocation


def _baseline_allocation(model: Model, unified: bool):
    """Hamilton apportionment per old region; None unless every region is usable."""
    allocation = {cell.key: 0 for cell in model.cells}
    for old, total in model.problem.old_totals.items():
        entries = []
        for cell in model.cells:
            if cell.old != old:
                continue
            weight = Fraction(1) if unified else cell.area_weight
            if weight is None:
                return None, f"旧区 {old} 的交叉格缺少面积权重"
            entries.append((cell, weight))
        if not entries:
            return None, f"旧区 {old} 没有交叉格"
        shares = _largest_remainder(total, entries)
        if shares is None:
            if total == 0:
                shares = {cell.key: 0 for cell, _ in entries}
            else:
                return None, f"旧区 {old} 的面积权重全为 0，无法摊分 {total} 人"
        allocation.update(shares)
    return allocation, None


def _baseline_compliance(model: Model, allocation):
    violations = []
    for cell in model.cells:
        value = allocation[cell.key]
        if cell.capacity is not None and value > cell.capacity:
            violations.append({
                "kind": "capacity",
                "text": f"x[{cell.old}->{cell.new}]={value} 超过容量上界 {cell.capacity}",
                "source": cell.source,
            })
    for subtotal in model.problem.public_subtotals:
        actual = sum(allocation[cell.key] for cell in model.cells
                     if cell.new in subtotal.regions)
        if actual != subtotal.count:
            violations.append({
                "kind": "public_subtotal",
                "text": f"{subtotal.key}（{'/'.join(subtotal.regions)}）"
                        f"={actual}，公开小计为 {subtotal.count}",
                "source": subtotal.source,
            })
    return violations


def _table(model: Model, vector):
    return [
        {"old": cell.old, "new": cell.new, "persons": int(value),
         "capacity": cell.capacity, "area_weight": _free_bound(cell.area_weight)}
        for cell, value in zip(model.cells, vector)
    ]


def _baseline_scenarios(model: Model):
    scenarios = []
    for name, unified in (
        ("declared_area_weights", False),
        ("uniform_over_allowed_cells", True),
    ):
        allocation, reason = _baseline_allocation(model, unified)
        if allocation is None:
            scenarios.append({
                "name": name,
                "available": False,
                "compliant": None,
                "reason": reason,
                "not_projected": True,
            })
            continue
        violations = _baseline_compliance(model, allocation)
        scenarios.append({
            "name": name,
            "available": True,
            "allocation": _table(model, [allocation[cell.key] for cell in model.cells]),
            "compliant": not violations,
            "violations": violations,
            "not_projected": True,
        })
    return scenarios


def _area_block(model: Model, feasible: bool):
    if model.problem.area_baseline == "none":
        return {
            "skipped": True,
            "mode": "none",
            "role": "面积基准已停用；可行范围本来就不依赖面积权重。",
            "area_weights_used_for_feasible_ranges": False,
            "not_a_confidence_interval": True,
        }
    scenarios = _baseline_scenarios(model)
    declared = scenarios[0]
    block = {
        "role": "独立的面积点摊分基准；仅作假设展示，绝不收紧或解释可行范围。",
        "mode": model.problem.area_baseline,
        "scenarios": scenarios,
        "area_weights_used_for_feasible_ranges": False,
        "feasible_range_independent_of_area_weights": True,
        "not_a_confidence_interval": True,
    }
    if declared["available"]:
        block["compliant"] = declared["compliant"]
        block["violations"] = declared["violations"]
        block["not_projected"] = True
        if declared["compliant"]:
            block["status"] = "compliant_point_estimate"
            block["note"] = "面积基准满足约束，但它仍只是点估计，不是唯一分配。"
        else:
            block["status"] = "noncompliant"
            block["note"] = (
                "面积基准违反容量或公开小计，明确不合规；未将其投影到约束后再称为面积摊分。"
                "可行范围只由旧区总量、容量和公开小计决定。"
            )
    else:
        block["compliant"] = None
        block["status"] = "unavailable"
        block["note"] = "缺少完整面积权重，未形成基准摊分，也未做任何投影。"
    if feasible:
        comparisons = []
        for query in model.problem.queries:
            row = {"expression": query, "label": _expression_label(query)}
            for scenario in scenarios:
                if not scenario["available"]:
                    continue
                persons = {(item["old"], item["new"]): item["persons"]
                           for item in scenario["allocation"]}
                row[scenario["name"]] = sum(
                    coefficient * sum(
                        persons.get((cell.old, cell.new), 0)
                        for cell in model.cells if cell.new == region
                    )
                    for region, coefficient in query.items()
                )
            comparisons.append(row)
        block["query_point_values"] = comparisons
    return block


# ---------------------------------------------------------------------------
# Query evaluation


def _query_cost(model: Model, expression):
    cost = [0] * model.nvars
    for region, coefficient in expression.items():
        for index, cell in enumerate(model.cells):
            if cell.new == region:
                cost[index] += coefficient
    return cost


def _expression_label(expression):
    parts = []
    for region, coefficient in expression.items():
        magnitude = abs(coefficient)
        term = f"{magnitude}*{region}" if magnitude != 1 else region
        if not parts:
            parts.append(term if coefficient > 0 else f"-{term}")
        else:
            parts.append((" + " if coefficient > 0 else " - ") + term)
    return "".join(parts)


def _free_bound(value):
    if value is None:
        return None
    if isinstance(value, Fraction):
        return value.numerator if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return value


def _endpoint(model: Model, outcome: SearchOutcome):
    witness = None
    if outcome.vector is not None:
        problems = check_witness(model, outcome.vector)
        if problems:
            raise InputError("端点见证校验失败：" + "；".join(problems))
        witness = _table(model, outcome.vector)
    certified = outcome.status == "optimal" and outcome.certified
    if outcome.status == "optimal_relaxation":
        note = "连续松弛端点（未认证整数端点）；如需整数认证请启用 solve_integer。"
    elif outcome.status == "budget":
        note = "达到节点预算：给出当前最好可行见证与精确松弛界，该端点未认证。"
    elif outcome.status == "infeasible":
        note = "无可行分配。"
    else:
        note = "整数端点已认证，附完整交叉格见证。"
    return {
        "value": _free_bound(outcome.value),
        "witness": witness,
        "certified_integer": certified,
        "relaxation_bound": _free_bound(outcome.lp_bound),
        "unverified": outcome.status == "budget",
        "note": note,
    }


def _marginal_warning(model: Model, expression, minimum: SearchOutcome,
                      maximum: SearchOutcome):
    """Quantify where independent marginal combination would mislead.

    A single-region expression is its own marginal range, so no extra solves
    are needed there; only genuine multi-region comparisons are probed.
    """
    if minimum.value is None or maximum.value is None:
        return None
    if sum(1 for coefficient in expression.values() if coefficient) < 2:
        return None
    naive_low = Fraction(0)
    naive_high = Fraction(0)
    for region, coefficient in expression.items():
        if not coefficient:
            continue
        low = _integer_search(model, _query_cost(model, {region: 1}), MINIMISE)
        high = _integer_search(model, _query_cost(model, {region: 1}), MAXIMISE)
        if low.value is None or high.value is None:
            return None
        if coefficient > 0:
            naive_low += coefficient * low.value
            naive_high += coefficient * high.value
        else:
            naive_low += coefficient * high.value
            naive_high += coefficient * low.value
    joint_low, joint_high = minimum.value, maximum.value
    if naive_low != joint_low or naive_high != joint_high:
        return {
            "naive_marginal_bounds": [_free_bound(naive_low), _free_bound(naive_high)],
            "joint_bounds": [_free_bound(joint_low), _free_bound(joint_high)],
            "text": (
                "各变量边际范围的独立组合与联合可行范围不一致；"
                "该组合甚至可能不可能是任何单张交叉格表的结果。本报告直接优化联合表达式。"
            ),
        }
    return None


def _feasibility(model: Model):
    outcome = _integer_search(model, [0] * model.nvars, MAXIMISE)
    if outcome.status == "infeasible":
        return {
            "feasible": False,
            "certified": outcome.certified,
            "witness": None,
            "conflict": _conflict_core(model),
            "note": "不存在满足旧区总量、容量与公开小计的整数交叉格分配。",
        }
    if outcome.status == "optimal_relaxation":
        return {
            "feasible": True,
            "certified": False,
            "witness": None,
            "conflict": None,
            "note": "连续松弛可行；整数分配未求解，见各查询的松弛端点。",
        }
    if outcome.vector is None:
        return {
            "feasible": None,
            "certified": False,
            "witness": None,
            "conflict": None,
            "relaxation_feasible": True,
            "relaxation_bound": _free_bound(outcome.lp_bound),
            "note": "连续松弛可行但整数见证未在节点预算内找到；不报告冲突，保留未认证。",
        }
    # An integral witness is itself a proof of feasibility, even when the node
    # budget stopped before the objective was certified.
    return {
        "feasible": True,
        "certified": outcome.status == "optimal",
        "witness": _table(model, outcome.vector),
        "conflict": None,
        "note": "已找到满足全部声明约束的整数交叉格分配，作为可行性见证。",
    }


def analyse(problem: Problem) -> dict:
    model = build_model(problem)
    feasibility = _feasibility(model)
    report = {
        "meta": {
            "region_version": problem.version,
            "unit": problem.unit,
            "new_regions": problem.new_regions,
            "old_regions": problem.old_totals,
            "cell_count": model.nvars,
            "solve_integer": problem.solve_integer,
            "node_budget": problem.node_budget,
            "no_individual_records": True,
            "range_semantics": (
                "已识别的可行整数/松弛范围，不是抽样置信区间；面积权重不收缩该范围。"
            ),
            "area_baseline_mode": problem.area_baseline,
            "inputs": problem.source_labels,
            "constraint_inventory": [_describe(constraint) for constraint in model.constraints],
            "area_baseline_used_for_constraints": False,
            "results_are_sampling_confidence_intervals": False,
        },
        "feasibility": feasibility,
        "queries": [],
        "area_baseline": {},
    }
    if not feasibility["feasible"] and feasibility["feasible"] is not None:
        block = _area_block(model, False)
        block["feasibility_status"] = "infeasible"
        report["area_baseline"] = block
        return report

    for index, expression in enumerate(problem.queries):
        cost = _query_cost(model, expression)
        minimum = _integer_search(model, cost, MINIMISE)
        maximum = _integer_search(model, cost, MAXIMISE)
        min_endpoint = _endpoint(model, minimum)
        max_endpoint = _endpoint(model, maximum)
        certified = min_endpoint["certified_integer"] and max_endpoint["certified_integer"]
        report["queries"].append({
            "index": index,
            "expression": expression,
            "label": _expression_label(expression),
            "unit": f"{problem.unit}（表达式可有符号）",
            "min": min_endpoint,
            "max": max_endpoint,
            "range": {
                "min": min_endpoint["value"],
                "max": max_endpoint["value"],
                "integer_certified": certified,
                "kind": "identified_feasible_range" if certified else "linear_relaxation_range",
                "note": (
                    "两端分别由单张完整交叉格表联合实现；不是边际范围的独立组合。"
                    if certified
                    else "至少一端未获整数认证；请将 relaxation_bound 与端点值视为松弛界。"
                ),
            },
            "marginal_combination_warning": _marginal_warning(
                model, expression, minimum, maximum
            ),
        })

    report["area_baseline"] = _area_block(model, True)
    if feasibility.get("feasible") is None:
        report["area_baseline"]["feasibility_status"] = "unverified"
        report["area_baseline"]["note_unverified"] = (
            "整数可行性未在节点预算内认证；查询端点按松弛/未认证处理，"
            "面积基准仅作独立点估计展示。"
        )
    return report


# ---------------------------------------------------------------------------
# CLI


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="区划转换可行性、整数端点范围与完整交叉格见证（仅读取汇总合成资料）"
    )
    parser.add_argument("--regions", default="fixtures/regions.json")
    parser.add_argument("--crosswalk", default="fixtures/crosswalk.json")
    parser.add_argument("--cells", default="fixtures/cells.csv")
    parser.add_argument("--output", default=None, help="写入 JSON 的路径；省略则打印标准输出")
    parser.add_argument("--no-area-baseline", action="store_true", help="跳过面积基准情景")
    args = parser.parse_args(argv)
    problem = load_problem(args.regions, args.crosswalk, args.cells)
    report = analyse(problem)
    if args.no_area_baseline:
        report["area_baseline"] = {"skipped": True}
    text = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
