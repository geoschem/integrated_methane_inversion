#!/usr/bin/env python3
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
try:
    from src.utilities.config_utils import load_config
except ModuleNotFoundError:
    from config_utils import load_config

"""
Utility script that validates the passed in IMI config file.
Arguments:
    config_path   [String]   : path to yaml config file
Usage:
    python sanitize_input_yaml.py <config.yml>

Validation rules:
 - Required keys come from `config_required` (dict: key -> rule)
 - Rules can be: a Python type (bool/int/float/str) or a list of allowed values
 - Booleans are flexible: true/false/yes/no/on/off/1/0 (any case) are accepted for bool keys
 - Conditional requirements are added ONLY when their controller is truthy
 - Optional keys are type-checked if present (via `optional_rules`)
 - Exits 1 if invalid
"""

# Sentinel meaning "required but no type/value validation"
ANY = object()
Rule = Union[type, List[str], object]  # object == ANY

# -------------------- REQUIRED KEYS (dict of rules) --------------------
config_required: Dict[str, Rule] = {
    "RunName": str,
    "Species": str,
    "SchedulerType": str,
    "SafeMode": bool,
    "S3Upload": bool,
    "UseGCHP": bool,
    "StartDate": int,
    "EndDate": int,
    "SpinupMonths": int,
    "SatelliteProduct": str,
    "UseWaterObs": bool,
    "isRegional": bool,
    "RegionID": str,
    "LonMin": float,
    "LonMax": float,
    "LatMin": float,
    "LatMax": float,
    "KalmanMode": bool,
    "CreateAutomaticRectilinearStateVectorFile": bool,
    "BufferRings": int,
    "BufferReductionFactor": float,
    "EmisThreshold": float,
    "ReducedDimensionStateVector": bool,
    "StateVectorFile": str,
    "ShapeFile": str,
    "PriorError": ANY,
    "ObsError": ANY,
    "Gamma": ANY,
    "PrecomputedJacobian": bool,
    "OffDiagonalPriorCov": bool,
    "LengthScalePriorCov": float,
    "Res": ["0.125x0.15625", "0.25x0.3125", "2.0x2.5", "4.0x5.0"],
    "Met": ["GEOSFP", "MERRA2"],
    "RunSetup": bool,
    "SetupTemplateRundir": bool,
    "SetupSpinupRun": bool,
    "SetupJacobianRuns": bool,
    "SetupInversion": bool,
    "SetupPosteriorRun": bool,
    "DoHemcoPriorEmis": bool,
    "DoSpinup": bool,
    "ReDoJacobian": bool,
    "DoJacobian": bool,
    "DoInversion": bool,
    "DoPosterior": bool,
    "DoPreview": bool,
    "DOFSThreshold": float,
    "RequestedMemory": ANY,
    "RequestedCPUs": int,
    "RequestedTime": str,
    "SchedulerPartition": str,
    "MaxSimultaneousRuns": int,
    "NumJacobianTracers": int,
    "PerturbValue": float,
    "HourlySpecies": bool,
    "PLANEFLIGHT": bool,
    "GOSAT": bool,
    "TCCON": bool,
    "AIRS": bool,
    "UseBCsForRestart": bool,
    "OutputPath": str,
    "DataPath": str,
    "PythonEnv": str,
    "RestartDownload": bool,
    "RestartFilePrefix": str,
    "BCpath": str,
    "BCversion": str,
    "HemcoPriorEmisDryRun": bool,
    "SpinupDryrun": bool,
    "ProductionDryRun": bool,
    "PosteriorDryRun": bool,
    "BCdryrun": bool,
    "LognormalErrors": bool,
    "MakePeriodsCSV": bool,
    "EnableOSSE": bool,
    "OptimizeSoil": bool,
}

# -------------------- OPTIONAL KEYS (validated if present) -------------
optional_rules: Dict[str, Rule] = {
    # KalmanMode-related
    "UpdateFreqDays": int,
    "NudgeFactor": float,
    "DynamicKFClustering": bool,
    "CustomPeriodsCSV": str,
    "FirstPeriod": int,
    "AutoAdvanceFirstPeriod": bool,
    # ReducedDimensionStateVector-related
    "ClusteringMethod": str,
    "NumberOfElements": int,
    "EmissionRateFilter": int,
    "PlumeCountFilter": int,
    "GroupByCountry": bool,
    "ForcedNativeResolutionElements": ANY,  # list of [lat, lon] pairs
    # S3 upload
    "S3UploadPath": str,
    "S3UploadFiles": ANY,  # list[str]
    # PrecomputedJacobian-related
    "ReferenceRunDir": str,
    # GCHP
    "CS_RES": int,
    "STRETCH_GRID": bool,
    "STRETCH_FACTOR": float,
    "TARGET_LAT": float,
    "TARGET_LON": float,
    # OSSE
    "DoOSSE": bool,
    "ObsErrorOSSE": int,
    "CreateAutomaticScaleFactorFileOSSE": bool,
    "EmisPerturbationOSSE": float,
    "ScaleFactorFileOSSE": str,
    # Other useful optional inputs
    "PointSourceDatasets": ANY,  # list[str]
    "InversionCPUs": int,
    "InversionMemory": ANY,
    "UseWaterObs": bool,  # already required above, here harmless if present
    "OptimizeBCs": bool,
    "OptimizeOH": bool,
    "PerturbValueOH": float,
    "PerturbValueBCs": float,
    "PriorErrorBCs": ANY,  # list[float]
    "PriorErrorBufferElements": ANY,  # list[float]
    "PriorErrorOH": ANY,  # list[float]
}

# -------------------- CONDITIONAL REQUIREMENTS -------------------------
conditional_dict: Dict[str, List[str]] = {
    "KalmanMode": [
        "UpdateFreqDays",
        "NudgeFactor",
        "DynamicKFClustering",
    ],
    "ReducedDimensionStateVector": [
        "ClusteringMethod",
        "NumberOfElements",
        "EmissionRateFilter",
        "PlumeCountFilter",
        "GroupByCountry",
    ],
    "PrecomputedJacobian": ["ReferenceRunDir"],
    "OffDiagonalPriorCov": ["LengthScalePriorCov"],
    "S3Upload": ["S3UploadPath", "S3UploadFiles"],
    "OptimizeBCs": ["PerturbValueBCs", "PriorErrorBCs"],
    "LognormalErrors": ["PriorErrorBufferElements"],
    "OptimizeOH": ["PerturbValueOH", "PriorErrorOH"],
    "EnableOSSE": ["DoOSSE", "ObsErrorOSSE", "CreateAutomaticScaleFactorFileOSSE"],
    "CreateAutomaticScaleFactorFileOSSE": ["EmisPerturbationOSSE"],
    "UseGCHP": ["CS_RES", "TOTAL_CORES", "NUM_NODES", "NUM_CORES_PER_NODE"],
    "STRETCH_GRID": ["STRETCH_FACTOR", "TARGET_LAT", "TARGET_LON"]
}

# -------------------- HELPERS -------------------------
def _is_boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return True
    if isinstance(value, (int, float)) and value in (0, 1):
        return True
    if isinstance(value, str) and value.strip().lower() in {
        "true",
        "false"
    }:
        return True
    return False


def _to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(int(value))
    if isinstance(value, str):
        s = value.strip().lower()
        if s in {"true"}:
            return True
        if s in {"false"}:
            return False
    return bool(value)


def _type_matches(value: Any, t: type) -> bool:
    if t is float:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if t is int:
        return isinstance(value, int) and not isinstance(value, bool)
    if t is str:
        return isinstance(value, str)
    if t is bool:
        return _is_boolish(value)
    return isinstance(value, t)


def _describe_rule(rule: Rule) -> str:
    if rule is ANY:
        return "any"
    if isinstance(rule, list):
        return " | ".join(repr(x) for x in rule)
    if isinstance(rule, type):
        return rule.__name__
    return repr(rule)


def _truthy_for_condition(value: Any) -> bool:
    return _to_bool(value) if _is_boolish(value) else bool(value)


def _missing_error_message(var: str) -> str:
    return (
        "Error: Missing input variable: "
        + var
        + ". Please add to config.yml file.\n"
        + "More information on config variables are available at: "
        + "https://imi.readthedocs.io/en/latest/getting-started/imi-config-file.html"
    )


# -------------------- VALIDATION -------------------------
def validate_config(cfg: Dict[str, Any]) -> Tuple[bool, List[str]]:
    errors: List[str] = []

    # 1) Base presence (required only)
    missing_base = [k for k in config_required.keys() if k not in cfg]
    if missing_base:
        errors.extend([_missing_error_message(k) for k in missing_base])

    # 2) Conditional expansions (only when controller is truthy)
    required_keys = set(config_required.keys())
    for controller, dependents in conditional_dict.items():
        if controller in cfg and _truthy_for_condition(cfg[controller]):
            required_keys.update(dependents)

    # 3) Presence after conditionals
    missing_after = [k for k in required_keys if k not in cfg]
    for k in missing_after:
        if k not in config_required:  # only add truly new misses
            errors.append(_missing_error_message(k))

    # 4) Value/type checks: validate required rules and optional rules for present keys
    def check_key(key: str, rule: Rule, val: Any) -> None:
        if rule is ANY:
            return
        if isinstance(rule, list):  # explicit allowed strings
            if not (isinstance(val, str) and val in rule):
                errors.append(
                    f"{key}: got {repr(val)} (type={type(val).__name__}), expected one of {rule}"
                )
            return
        if isinstance(rule, type):
            if not _type_matches(val, rule):
                errors.append(
                    f"{key}: got {repr(val)} (type={type(val).__name__}), expected {_describe_rule(rule)}"
                )
            # normalize bools if matched
            elif rule is bool:
                cfg[key] = _to_bool(val)
            return
        errors.append(f"Internal rule error for {key}: {_describe_rule(rule)}")

    # Required keys present -> check against their rules
    for key, rule in config_required.items():
        if key in cfg:
            check_key(key, rule, cfg[key])

    # Optional keys present -> check types too
    for key, rule in optional_rules.items():
        if key in cfg:
            check_key(key, rule, cfg[key])

    # 5) Sanity checks
    if (
        all(k in cfg for k in ("LonMin", "LonMax"))
        and isinstance(cfg["LonMin"], (int, float))
        and isinstance(cfg["LonMax"], (int, float))
    ):
        if not (-180 <= cfg["LonMin"] <= 180 and -180 <= cfg["LonMax"] <= 180):
            errors.append("LonMin/LonMax must be within [-180, 180].")
        if cfg["LonMin"] >= cfg["LonMax"]:
            errors.append("LonMin must be < LonMax.")

    if (
        all(k in cfg for k in ("LatMin", "LatMax"))
        and isinstance(cfg["LatMin"], (int, float))
        and isinstance(cfg["LatMax"], (int, float))
    ):
        if not (-90 <= cfg["LatMin"] <= 90 and -90 <= cfg["LatMax"] <= 90):
            errors.append("LatMin/LatMax must be within [-90, 90].")
        if cfg["LatMin"] >= cfg["LatMax"]:
            errors.append("LatMin must be < LatMax.")

    if (
        all(k in cfg for k in ("StartDate", "EndDate"))
        and isinstance(cfg["StartDate"], int)
        and isinstance(cfg["EndDate"], int)
    ):
        if cfg["StartDate"] > cfg["EndDate"]:
            errors.append("StartDate must be <= EndDate (yyyymmdd integers).")

    return (len(errors) == 0, errors)


# -------------------- CLI -------------------------
def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: validate_config.py <config.yml>", file=sys.stderr)
        sys.exit(2)

    config_path = Path(sys.argv[1])

    try:
        cfg = load_config(config_path)
    except ValueError as e:
        print("❌ Config is invalid.", file=sys.stderr)
        for err in str(e).split("\n"):
            print(err, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Failed to load YAML: {e}", file=sys.stderr)
        sys.exit(2)

    valid, errs = validate_config(cfg)

    if valid:
        print("✅ Config is valid.")
        sys.exit(0)
    else:
        print("❌ Config is invalid.")
        for e in errs:
            print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
