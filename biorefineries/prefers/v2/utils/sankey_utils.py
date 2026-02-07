# -*- coding: utf-8 -*-
"""
Sankey Utilities for PreFerS
============================

Step-wise Sankey data generation utilities for BioSTEAM systems.
"""

from __future__ import annotations

from collections import defaultdict
import re
from typing import Dict, Iterable, List, Tuple

import numpy as np


DEFAULT_STEP_MAP = {
    200: "Media Preparation",
    300: "Conversion",
    400: "Recovery",
    500: "Purification",
    600: "Formulation",
    900: "Facilities",
}


def _unit_area_from_id(unit_id: str | None) -> int | None:
    if not unit_id:
        return None
    match = re.search(r"(\d{2,4})", str(unit_id))
    if not match:
        return None
    try:
        value = int(match.group(1))
    except Exception:
        return None
    return int(value // 100) * 100


def get_unit_step(unit, step_map: Dict[int, str] | None = None) -> str:
    step_map = step_map or DEFAULT_STEP_MAP
    area = _unit_area_from_id(getattr(unit, "ID", None))
    if area in step_map:
        return step_map[area]
    return "Other"


def build_step_groups(system, step_map: Dict[int, str] | None = None) -> Dict[str, List]:
    step_map = step_map or DEFAULT_STEP_MAP
    groups: Dict[str, List] = {name: [] for name in step_map.values()}
    groups.setdefault("Other", [])
    for unit in system.units:
        step = get_unit_step(unit, step_map=step_map)
        groups.setdefault(step, []).append(unit)
    return groups


def _carbon_mass_flow(stream) -> float:
    try:
        carbon_kmol = stream.get_atomic_flow("C")
        return float(carbon_kmol) * 12.011
    except Exception:
        return 0.0


def _energy_flow(stream) -> float:
    try:
        lhv = getattr(stream, "LHV", None)
        if lhv is not None and np.isfinite(lhv):
            return float(lhv) * float(stream.F_mass)
    except Exception:
        pass
    try:
        enthalpy = float(stream.H)
        return max(enthalpy, 0.0)
    except Exception:
        return 0.0


def _stream_flow_value(stream, flow_property: str) -> float:
    if flow_property == "mass":
        return float(stream.F_mass)
    if flow_property == "carbon":
        return _carbon_mass_flow(stream)
    if flow_property == "energy":
        return _energy_flow(stream)
    if flow_property in stream.chemicals.IDs:
        return float(stream.imass[flow_property])
    try:
        return float(getattr(stream, flow_property))
    except Exception:
        return 0.0


def generate_step_sankey_data(
    system,
    flow_property: str = "mass",
    units: str = "kg/hr",
    step_map: Dict[int, str] | None = None,
    include_external: bool = True,
) -> Dict:
    """
    Generate step-wise sankey data for a BioSTEAM system.
    """
    step_map = step_map or DEFAULT_STEP_MAP
    unit_to_step = {unit: get_unit_step(unit, step_map=step_map) for unit in system.units}

    nodes: List[str] = []
    links: Dict[Tuple[int, int], float] = defaultdict(float)

    def node_index(name: str) -> int:
        if name not in nodes:
            nodes.append(name)
        return nodes.index(name)

    for unit in system.units:
        source_step = unit_to_step.get(unit, "Other")

        for stream in unit.outs:
            if not stream:
                continue
            if stream.sink:
                target_step = unit_to_step.get(stream.sink, "Other")
                if source_step == target_step:
                    continue
                value = _stream_flow_value(stream, flow_property)
                if value > 0:
                    src = node_index(source_step)
                    tgt = node_index(target_step)
                    links[(src, tgt)] += value
            elif include_external:
                value = _stream_flow_value(stream, flow_property)
                if value > 0:
                    src = node_index(source_step)
                    tgt = node_index("Outputs")
                    links[(src, tgt)] += value

        if include_external:
            for stream in unit.ins:
                if not stream:
                    continue
                if stream.source is None:
                    value = _stream_flow_value(stream, flow_property)
                    if value > 0:
                        src = node_index("Inputs")
                        tgt = node_index(source_step)
                        links[(src, tgt)] += value

    link_sources = []
    link_targets = []
    link_values = []
    link_labels = []

    for (src, tgt), value in links.items():
        link_sources.append(src)
        link_targets.append(tgt)
        link_values.append(value)
        link_labels.append(f"{value:.2f} {units}")

    return {
        "nodes": nodes,
        "links": {
            "source": link_sources,
            "target": link_targets,
            "value": link_values,
            "label": link_labels,
        },
        "flow_property": flow_property,
        "units": units,
    }


def compute_step_costs(system, step_map: Dict[int, str] | None = None) -> Dict[str, float]:
    step_map = step_map or DEFAULT_STEP_MAP
    costs = defaultdict(float)
    for unit in system.units:
        step = get_unit_step(unit, step_map=step_map)
        cost = 0.0
        try:
            cost = float(unit.installed_cost)
        except Exception:
            try:
                cost = float(unit.purchase_cost)
            except Exception:
                cost = 0.0
        if np.isfinite(cost):
            costs[step] += cost
    return dict(costs)


def compute_step_totals(step_costs: Dict[str, float]) -> Dict[str, float]:
    isbl = 0.0
    osbl = 0.0
    for step, cost in step_costs.items():
        if step.lower() == "facilities":
            osbl += cost
        else:
            isbl += cost
    return {
        "ISBL": isbl,
        "OSBL": osbl,
        "Total": isbl + osbl,
    }
