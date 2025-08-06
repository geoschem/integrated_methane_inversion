import sys
import yaml

"""
A simple utility script that tests whether all required IMI variables are present 
in the yaml config file given. Exits 1 if there are missing variables
Arguments
    config_path   [String]   : path to yaml config file
"""
# ************ Add required config variables to the corresponding list **************

# variables required on all systems
config_required = [
    "RunName",
    "UseSlurm",
    "SafeMode",
    "S3Upload",
    "StartDate",
    "EndDate",
    "SpinupMonths",
    "BlendedTROPOMI",
    "isRegional",
    "RegionID",
    "LonMin",
    "LonMax",
    "LatMin",
    "LatMax",
    "KalmanMode",
    "CreateAutomaticRectilinearStateVectorFile",
    "nBufferClusters",
    "BufferDeg",
    "LandThreshold",
    "OffshoreEmisThreshold",
    "ReducedDimensionStateVector",
    "StateVectorFile",
    "ShapeFile",
    "PriorError",
    "ObsError",
    "Gamma",
    "PrecomputedJacobian",
    "Res",
    "Met",
    "RunSetup",
    "SetupTemplateRundir",
    "SetupSpinupRun",
    "SetupJacobianRuns",
    "SetupInversion",
    "SetupPosteriorRun",
    "DoHemcoPriorEmis",
    "DoSpinup",
    "ReDoJacobian",
    "DoJacobian",
    "DoInversion",
    "DoPosterior",
    "DoPreview",
    "DOFSThreshold",
    "RequestedMemory",
    "RequestedCPUs",
    "RequestedTime",
    "SchedulerPartition",
    "MaxSimultaneousRuns",
    "NumJacobianTracers",
    "PerturbValue",
    "HourlyCH4",
    "PLANEFLIGHT",
    "GOSAT",
    "TCCON",
    "AIRS",
    "OutputPath",
    "DataPath",
    "CondaEnv",
    "CondaFile",
    "RestartDownload",
    "RestartFilePrefix",
    "BCpath",
    "BCversion",
    "HemcoPriorEmisDryRun",
    "SpinupDryrun",
    "ProductionDryRun",
    "PosteriorDryRun",
    "BCdryrun",
    "LognormalErrors",
    "MakePeriodsCSV",
    "UseWaterObs",
    "EnableOSSE",
]

# dict of variables that are required if another variable is set to true
# For example UpdateFreqDays is only required if KalmanMode is set to true
conditional_dict = {}
conditional_dict["KalmanMode"] = [
    "UpdateFreqDays",
    "NudgeFactor",
    "DynamicKFClustering",
]
conditional_dict["ReducedDimensionStateVector"] = [
    "ClusteringMethod",
    "NumberOfElements",
    "EmissionRateFilter",
    "PlumeCountFilter",
    "GroupByCountry",
]
conditional_dict["PrecomputedJacobian"] = ["ReferenceRunDir"]
conditional_dict["S3Upload"] = [
    "S3UploadPath",
    "S3UploadFiles",
]
conditional_dict["OptimizeBCs"] = ["PerturbValueBCs", "PriorErrorBCs"]
conditional_dict["LognormalErrors"] = ["PriorErrorBufferElements"]
conditional_dict["OptimizeOH"] = ["PerturbValueOH", "PriorErrorOH"]
conditional_dict["EnableOSSE"] = ["DoOSSE", "ObsErrorOSSE", "CreateAutomaticScaleFactorFileOSSE"]
conditional_dict["CreateAutomaticScaleFactorFileOSSE"] = ["EmisPerturbationOSSE"]


def raise_error_message(var):
    """
    Description: raise an error message about missing config variable
    """
    message = (
        "Error: Missing input variable: "
        + var
        + ". Please add to config.yml file."
        + "\n More information on config variables are available at:"
        + "https://imi.readthedocs.io/en/latest/getting-started/imi-config-file.html"
    )
    raise ValueError(message)


if __name__ == "__main__":
    config_path = sys.argv[1]
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    inputted_config = config.keys()

    # require additional variables if conditional dict key is set to true
    for key in conditional_dict.keys():
        if key not in inputted_config:
            raise_error_message(key)
        elif config[key]:
            config_required = config_required + conditional_dict[key]

    missing_input_vars = [x for x in config_required if x not in inputted_config]
    for var in missing_input_vars:
        raise_error_message(var)

    if len(missing_input_vars) > 0:
        sys.exit(1)
