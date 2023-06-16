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
config_required_local_cluster = []

# variables required on all systems
config_required = [
    "RunName",
    "isAWS",
    "UseSlurm",
    "SafeMode",
    "StartDate",
    "EndDate",
    "SpinupMonths",
    "LonMin",
    "LonMax",
    "LatMin",
    "LatMax",
    "NestedGrid",
    "NestedRegion",
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
]

clustering_vars = [
    "ClusteringMethod",
    "NumberOfElements",
]

if __name__ == "__main__":
    config_path = sys.argv[1]
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    inputted_config = config.keys()

    # only require clustering vars if reduced dimension state vector is true
    if config["ReducedDimensionStateVector"]:
        config_required = config_required + clustering_vars
        
    # update required vars based on system
    if config["isAWS"]:
        required_vars = config_required + config_required_aws
    else:
        required_vars = config_required + config_required_local_cluster

    missing_input_vars = [x for x in required_vars if x not in inputted_config]
    for var in missing_input_vars:
        message = (
            "Error: Missing input variable: "
            + var
            + ". Please add to config.yml file."
            + "\n More information on config variables are available at:"
            + "https://imi.readthedocs.io/en/latest/getting-started/imi-config-file.html"
        )
        raise ValueError(message)

    if len(missing_input_vars) > 0:
        sys.exit(1)
