from __future__ import annotations

import json

import pytest
from cobra import Metabolite

from scripts.gem_annotate.metabolites import (
    _MNXM_CHEMISTRY_CONFLICT_NOTE,
    _apply_mnxm,
)


def _prop(*, formula: str, charge: object) -> dict:
    return {
        "MNXM_TEST": {
            "formula": formula,
            "charge": charge,
            "inchi": "",
            "inchikey": "",
        }
    }


def test_mnxm_complete_pair_is_applied_together() -> None:
    metabolite = Metabolite("m_c", formula=None, charge=0)

    _, status = _apply_mnxm(
        metabolite,
        "MNXM_TEST",
        {},
        _prop(formula="C3H2O6P", charge="-3"),
    )

    assert status == "applied"
    assert (metabolite.formula, metabolite.charge) == ("C3H2O6P", -3)


def test_mnxm_matching_formula_updates_charge_as_one_pair() -> None:
    metabolite = Metabolite("m_c", formula="C3H2O6P", charge=0)

    _, status = _apply_mnxm(
        metabolite,
        "MNXM_TEST",
        {},
        _prop(formula="C3H2O6P", charge="-3"),
    )

    assert status == "applied"
    assert (metabolite.formula, metabolite.charge) == ("C3H2O6P", -3)


def test_mnxm_formula_conflict_preserves_existing_pair_and_records_audit() -> None:
    metabolite = Metabolite("m_c", formula="C3H5O6P", charge=0)

    _, status = _apply_mnxm(
        metabolite,
        "MNXM_TEST",
        {},
        _prop(formula="C3H2O6P", charge="-3"),
    )

    assert status == "conflict"
    assert (metabolite.formula, metabolite.charge) == ("C3H5O6P", 0)
    conflict = json.loads(
        metabolite.notes[_MNXM_CHEMISTRY_CONFLICT_NOTE]
    )
    assert conflict == {
        "action": "preserved_existing_pair",
        "existing_charge": 0,
        "existing_formula": "C3H5O6P",
        "mnxm_id": "MNXM_TEST",
        "proposed_charge": -3,
        "proposed_formula": "C3H2O6P",
        "source": "MetaNetX",
    }


@pytest.mark.parametrize(
    ("formula", "charge"),
    [
        ("", "-3"),
        ("C3H2O6P", ""),
        ("C3H2O6P", "NA"),
        ("C3H2O6P", "not-a-charge"),
    ],
)
def test_mnxm_incomplete_pair_never_partially_updates(
    formula: str, charge: object
) -> None:
    metabolite = Metabolite("m_c", formula="CH4", charge=0)

    _, status = _apply_mnxm(
        metabolite,
        "MNXM_TEST",
        {},
        _prop(formula=formula, charge=charge),
    )

    assert status == "incomplete"
    assert (metabolite.formula, metabolite.charge) == ("CH4", 0)
    assert _MNXM_CHEMISTRY_CONFLICT_NOTE not in metabolite.notes
