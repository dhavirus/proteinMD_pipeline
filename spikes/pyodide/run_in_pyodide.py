"""Pyodide side of the ADR-0002 spike: rebuild the model from JSON and run the audit.

Only the pure core of simprep is importable here (no gemmi, PyYAML or jsonschema).
"""

import json
import time

from simprep.audit import run_audit
from simprep.config import default_config
from simprep.rules import Rule, RuleSet
from simprep.structure.model import (
    Atom,
    Link,
    LinkPartner,
    ModifiedResidue,
    Residue,
    ResidueClass,
    ResidueId,
    Structure,
    UnobservedResidue,
)


def rid(data):
    return ResidueId(**data)


def residue(data):
    atoms = tuple(Atom(**a | {"position": tuple(a["position"])}) for a in data["atoms"])
    return Residue(
        rid(data["id"]), data["name"], ResidueClass(data["residue_class"]), atoms, data["label_seq"]
    )


def partner(data):
    return LinkPartner(**data | {"residue": rid(data["residue"])})


def link(data):
    return Link(
        **data | {"partner1": partner(data["partner1"]), "partner2": partner(data["partner2"])}
    )


def structure(data):
    return Structure(
        name=data["name"],
        residues=tuple(residue(r) for r in data["residues"]),
        links=tuple(link(item) for item in data["links"]),
        modified_residues=tuple(
            ModifiedResidue(**m | {"residue": rid(m["residue"])}) for m in data["modified_residues"]
        ),
        unobserved_residues=tuple(
            UnobservedResidue(**u | {"residue": rid(u["residue"])})
            for u in data["unobserved_residues"]
        ),
        annotation_categories=frozenset(data["annotation_categories"]),
    )


def ruleset(data):
    rules = tuple(Rule.from_dict(r) for r in data["rules"])
    return RuleSet(data["version"], data["sha256"], rules, data["audit_defaults"])


def main(inputs_json):
    started = time.perf_counter()
    inputs = json.loads(inputs_json)
    model, rules = structure(inputs["structure"]), ruleset(inputs["ruleset"])
    rebuilt = time.perf_counter()
    findings = run_audit(model, rules, default_config(rules.audit_defaults))
    audited = time.perf_counter()
    return json.dumps(
        {
            "rebuild_ms": round((rebuilt - started) * 1000),
            "audit_ms": round((audited - rebuilt) * 1000),
            "atoms": sum(len(r.atoms) for r in model.residues),
            "findings": [f.to_dict() for f in findings],
        }
    )
