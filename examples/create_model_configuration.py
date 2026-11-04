import argparse
import sys
from pathlib import Path

# %% set up workspace path
oPath_parent = Path(__file__).resolve().parents[1]
sys.path.append(str(oPath_parent))

from pyflowline.configuration import create_pyflowline_template_configuration_file

parser = argparse.ArgumentParser(
    description="Create a PyFlowline YAML model configuration template."
)
parser.add_argument(
    "--sMesh_type",
    help="Target mesh type (hexagon, square, latlon, mpas)",
    type=str,
    default="hexagon",
)
parser.add_argument("--iCase_index", help="Case simulation index", type=int, default=1)
parser.add_argument(
    "--dResolution_meter", help="Target metric resolution in meters", type=float, default=50000.0
)
parser.add_argument("--sDate", help="Simulation date (YYYYMMDD)", type=str, default="20230101")
parser.add_argument(
    "--sFilename_output", help="Output YAML configuration filepath", type=str, default=None
)

# Example usage:
# python create_model_configuration.py --sMesh_type hexagon --iCase_index 1 --dResolution_meter 50000 --sDate 20230101
pArgs = parser.parse_args()

sMesh_type = pArgs.sMesh_type
iCase_index = pArgs.iCase_index
dResolution_meter = pArgs.dResolution_meter
sDate = pArgs.sDate

oFolder_data = oPath_parent.joinpath("data", "susquehanna")
oFolder_input = oFolder_data.joinpath("input")
oFolder_output = oFolder_data.joinpath("output")

if pArgs.sFilename_output:
    oFilename_config = Path(pArgs.sFilename_output).resolve()
else:
    oFilename_config = oFolder_input.joinpath(f"template_{sMesh_type}.yaml")

oPyflowline = create_pyflowline_template_configuration_file(
    str(oFilename_config),
    sWorkspace_input=str(oFolder_input),
    sWorkspace_output=str(oFolder_output),
    iFlag_use_mesh_dem_in=1,
    sMesh_type_in=sMesh_type,
    iCase_index_in=iCase_index,
    dResolution_meter_in=dResolution_meter,
    sDate_in=sDate,
)

print(f"Generated YAML template configuration at: {oFilename_config}")
print("Finished")
