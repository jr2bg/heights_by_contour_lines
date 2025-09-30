# -*- coding: utf-8 -*-
import os
import argparse
import pprint
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

parser = argparse.ArgumentParser(
    description="Generates the data from given points from contour lines",
    usage="%(prog)s contour_lines points"
)
parser.add_argument("contour_lines", help="csv output path from autocad that contains the lines of the contour map")
parser.add_argument("points", help="csv path with the desired points")
args = parser.parse_args()

def get_center_of_mass(points):
    """gets the center of mass"""
    x = 0
    y = 0
    for point in points:
        x += point[0]
        y += point[1]
    return x/len(points), y/len(points)

def get_directions(center_x, center_y, points):
    """gets the directions of the legs
    For this, a rest is done between the center and the points.
    Then the dot product is done between the points

    Parameters
    ----------
    center_x : float
    center_y : float
    points : list[(float,float)]
        points to get its height, in X,Y format

    Returns
    -------
    dict[key: tuple, value: List[float]
        directions and it corresponding points
    """
    # dictionary of the legs by each direction
    # the directions will be unitary and the legs will be in the 
    # non-modified origin
    legs_per_direction = {}

    for point in points:
    
        # taking as origin the center
        x,y = point
        nx,ny = (x - center_x,y - center_y)

        # iterate over the dictionary keys to determine direction
        dir_max_mag = None
        for dir in legs_per_direction.keys():
            magnitude_proj =nx*dir[0] + ny*dir[1]
            magnitude_vec = (nx**2 + ny**2)**(1/2)

            # the projection is at least 0.75 the magnitude of the vector
            # as they must be parallel
            if magnitude_proj > 0.75*magnitude_vec:
                dir_max_mag = dir

        # if it is the first element on the dict or the point doesnt follow
        # any of the previous directions, add the unitary direction
        if not legs_per_direction or dir_max_mag is None:
            len_np = (nx**2 + ny**2)**(1/2)
            unitary_dir = (nx/len_np, ny/len_np)
            legs_per_direction[unitary_dir]=[0 for _ in range(9)]

            # add the point in the position of its magnitude divided by 2, as 
            # they are separated by 2 meters
            legs_per_direction[unitary_dir][round(len_np/2)-1] = point
        else:
            # add the point in the position of its magnitude projection by 2, as 
            # they are separated by 2 meters
            legs_per_direction[dir_max_mag][round(magnitude_proj/2)-1] = point
    
    pprint.pprint(legs_per_direction)
    return legs_per_direction

def organize_points(st_x, st_y, end_x, end_y, center_x, center_y, points):
    """organizes the points with the given method
    
    the method consists in a line dividing the plane. The center is supposed 
    lie on it, but not in all cases. This line separates the plane in two. The
    first leg is the one on the rightest and below the line, the rest follow
    the clockwise direction.

    Parameters
    ----------
    st_x : float
        x beginning of the dividing line
    st_y : float
        y beginning of the dividing line
    end_x : float
        x end of the dividing line
    end_y : float
        y end of the dividing line
    center_x : float
        x position of the center
    center_y : float
        y position of the center
    points : list[(float,float)]
        points to get its height, in X,Y format
    
    Returns
    -------
    pandas.Dataframe
        data organized according to the selection rule
    """
    # get the 4 vectors of direction

    # from the center, get all 

    # get the 9 legs for each vector

    # 

    pass

def format_interpolated_result(output):
    """creates a new file with the appropriate format for the given points
    the points have to be in a specific order??

    File is saved in the folder where the script is run
    
    Parameter
    output: pandas.DataFrame
        dataframe containing the information of the interest points in X,Y,Z cols
    """
    return 

def main(contour_map_path, points_path):
    # read the passed file with the edges of contour map
    df = pd.read_csv(contour_map_path)

    # select only points
    ndf = df[["Inicial X", "Inicial Y", "Fin X", "Fin Y", "Fin Z"]]

    # renaming the columns to perform iteration over dataframe
    ndf.rename(
        columns={
            "Inicial X":"Inicial_X",
            "Inicial Y":"Inicial_Y",
            "Fin X":"Fin_X",
            "Fin Y":"Fin_Y",
            "Fin Z":"Fin_Z"
            },
        inplace=True)
    
    # read the passed file with the objective points
    df_points = pd.read_csv(points_path)
    
    # rename the columns
    df_points.rename(
        columns={"Posición X": "X", "Posición Y":"Y", "Posición Z": "Z"},
        inplace=True)
    
    points = []
    for row in df_points.itertuples():
        points.append((row.X,row.Y))

    cx, cy = get_center_of_mass(points)
    get_directions(cx,cy, points)

    # scatter points set
    points_contour_line = set()

    for row in ndf.itertuples():

        # as each row contains two points, then both have to be considered
        points_contour_line.add((row.Inicial_X, row.Inicial_Y, row.Fin_Z))
        points_contour_line.add((row.Fin_X, row.Fin_Y, row.Fin_Z))
    
    # create a numpy array
    np_scattered = np.array(list(points_contour_line))
    
    # get x,y grid
    grid_xy = np_scattered[:,:2]

    # get z data
    elevs = np_scattered[:,2]

    # use a numpy array for the points
    points = np.asarray(points)

    # interpolate with the provided data
    interpolated = griddata(grid_xy, elevs, points, method="linear")
    center_interpolated = griddata(grid_xy, elevs, (cx, cy), method="linear")

    # create dfs to store information
    df_grid = pd.DataFrame(points, columns = ["X", "Y"])
    df_elev = pd.DataFrame(interpolated, columns=["Z"])

    # save the information in a single dataframe for easier manipulation
    output = pd.concat([df_grid, df_elev], axis=1)

    # take the difference of the original and center's height
    output["Z"] = output["Z"] - center_interpolated
    print(output)
    print(f"center: x={cx}, y={cy}, z={center_interpolated}")
    return output

if __name__ == "__main__":
    # check that both files exist
    if (os.path.isfile(args.contour_lines) and os.path.isfile(args.points)):
        main(args.contour_lines,args.points)
    elif not os.path.isfile(args.contour_lines):
        print(f"[ERROR] Contour map file `{args.contour_lines}` was not found")
    else:
        print(f"[ERROR] Points file `{args.points}` was not found")