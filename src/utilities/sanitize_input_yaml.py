import sys
import yaml

"""
A simple utility script that tests whether all required IMI variables are present 
in the yaml config file given. Exits 1 if there are missing variables
Arguments
    config_path   [String]   : path to yaml config file
"""
# ************ Add required config variables to the corresponding list **************

# variables only required by AWS
config_required_aws = [
    "CondaFile",
]

# variables only required by local cluster
config_required_local_cluster = [
    "DataPathTROPOMI",
    "GEOSChemEnv",
]

# variables required on all systems
config_required = [
    "RunName",
    "isAWS",
    "UseSlurm",
    "SafeMode",
    "StartDate",
    "EndDate",
    "SpinupMonths",
    "BlendedTROPOMI",
    "LonMin",
    "LonMax",
    "LatMin",
    "LatMax",
    "isRegional",
    "RegionID",
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
    "SetupTemplateRundir",
    "SetupSpinupRun",
    "SetupJacobianRuns",
    "SetupInversion",
    "SetupPosteriorRun",
    "RunSetup",
    "DoSpinup",
    "DoJacobian",
    "DoInversion",
    "DoPosterior",
    "DoPreview",
    "DOFSThreshold",
    "PerturbValue",
    "UseEmisSF",
    "UseOHSF",
    "HourlyCH4",
    "PLANEFLIGHT",
    "GOSAT",
    "TCCON",
    "AIRS",
    "OutputPath",
    "DataPath",
    "CondaEnv",
    "RestartDownload",
    "RestartFilePrefix",
    "RestartFilePreviewPrefix",
    "BCpath",
    "BCversion",
    "PreviewDryRun",
    "SpinupDryrun",
    "ProductionDryRun",
    "PosteriorDryRun",
    "BCdryrun",
    "SimulationMemory",
    "SimulationCPUs",
    "JacobianMemory",
    "JacobianCPUs",
    "RequestedTime",
    "SchedulerPartition",
    "KalmanMode",
    "S3Upload",
]

# dict of variables that are required if another variable is set to true 
# For example UpdateFreqDays is only required if KalmanMode is set to true
conditional_dict = {}
conditional_dict["KalmanMode"] = [
    "UpdateFreqDays",
    "NudgeFactor",
    "DynamicKFClustering"
]
conditional_dict["ReducedDimensionStateVector"] = [
    "ClusteringMethod",
    "NumberOfElements",
]
conditional_dict["PrecomputedJacobian"] = ["ReferenceRunDir"]
conditional_dict["S3Upload"] = [
    "S3UploadPath",
    "S3UploadFiles",
]
conditional_dict["OptimizeBCs"] = ["PerturbValueBCs", "PriorErrorBCs"]

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

    # update required vars based on system
    if config["isAWS"]:
        required_vars = config_required + config_required_aws
    else:
        required_vars = config_required + config_required_local_cluster

    missing_input_vars = [x for x in required_vars if x not in inputted_config]
    for var in missing_input_vars:
        raise_error_message(var)

    if len(missing_input_vars) > 0:
        sys.exit(1)
