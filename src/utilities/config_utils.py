#!/usr/bin/env python3

from copy import deepcopy


SECTION_SCOPED_KEYS = (
    "OptimizeSoil",
    "PriorErrorOH",
    "OptimizeOH",
    "AdditionalDiagnostics",
    "SatelliteProduct",
    "UseWaterObs",
)

DIAGNOSTIC_KEY_MAP = {
    "OBSPACK": "DoObsPack",
    "PLANEFLIGHT": "PLANEFLIGHT",
    "GOSAT": "GOSAT",
    "TCCON": "TCCON",
    "AIRS": "AIRS",
}


def _active_species(config):
    species = config.get("Species")
    return str(species).upper() if species is not None else ""


def validate_hierarchical_config(config):
    """
    Validate section-scoped key policy for CH4/CO2 hierarchy.

    Returns a list of error strings.
    """
    errors = []

    if not isinstance(config, dict):
        return ["Top-level YAML must be a mapping (dict)."]

    for key in SECTION_SCOPED_KEYS:
        if key in config:
            species = _active_species(config)
            target_section = species if species in ("CH4", "CO2") else "CH4/CO2"
            errors.append(
                f"Top-level '{key}' is not allowed. Put it under the '{target_section}' section."
            )

    species = _active_species(config)
    if species in ("CH4", "CO2"):
        if species not in config or not isinstance(config.get(species), dict):
            errors.append(
                f"Missing required '{species}' section in config.yml for Species={species}."
            )

    return errors


def normalize_config(config):
    """
    Normalize CH4/CO2 hierarchical config to runtime flat keys.
    """
    if not isinstance(config, dict):
        return config

    cfg = deepcopy(config)
    species = _active_species(cfg)
    section = cfg.get(species, {}) if isinstance(cfg.get(species), dict) else {}

    if isinstance(section, dict):
        for key, value in section.items():
            cfg[key] = value

    # keep existing runtime expectations stable
    cfg.setdefault("OptimizeOH", False)
    cfg.setdefault("OptimizeSoil", False)
    cfg.setdefault("UseWaterObs", False)
    cfg.setdefault("AdditionalDiagnostics", [])

    diagnostics = cfg.get("AdditionalDiagnostics", [])
    if isinstance(diagnostics, str):
        diagnostics = [diagnostics]
    if not isinstance(diagnostics, list):
        diagnostics = []
    cfg["AdditionalDiagnostics"] = diagnostics

    for legacy_key in DIAGNOSTIC_KEY_MAP.values():
        cfg.setdefault(legacy_key, False)

    for item in diagnostics:
        if not isinstance(item, str):
            continue
        mapped = DIAGNOSTIC_KEY_MAP.get(item.strip().upper())
        if mapped is not None:
            cfg[mapped] = True

    return cfg
