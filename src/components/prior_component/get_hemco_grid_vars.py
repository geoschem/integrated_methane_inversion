import sys
import xarray as xr
import numpy as np

def calc_edges(da):
    diff = (da[1].item()-da[0].item())/2
    edges = list(da.values - diff)
    edges.append(da[-1].item() + diff)
    return edges

def get_hemco_grid_vars(sv):
    YMID = list(sv.lat.values)
    YEDGE = calc_edges(sv.lat)
    XEDGE = calc_edges(sv.lon)
    
    if option == "YMID":
        print(" ".join([f"{num:.2f}" for num in YMID]))
    elif option == "YEDGE":
        print(" ".join([f"{num:.3f}" for num in YEDGE]))
    elif option == "XMIN":
        XMIN = XEDGE[0]
        print(XMIN)
    elif option == "XMAX":
        XMAX = XEDGE[-1]
        print(XMAX)
    elif option == "YMIN":
        YMIN = YEDGE[0]
        print(YMIN)
    elif option == "YMAX":
        YMAX = YEDGE[-1]
        print(YMAX)
    elif option == "NX":
        NX = len(list(sv.lon.values))
        print(NX)
    elif option == "NY":
        NY = len(YMID)
        print(NY)
    else:
        raise ValueError(f"Option '{option}' not recognized.")
    

if __name__ == "__main__":
    sv_path = sys.argv[1]
    option = sys.argv[2]
    sv = xr.load_dataset(sv_path)
    
    # print out relevant grid information
    get_hemco_grid_vars(sv)