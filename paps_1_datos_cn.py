# -*- coding: utf-8 -*-
import os
import math
import argparse
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

parser = argparse.ArgumentParser(
    description="Generates the data from given points from contour lines",
    #usage="%(prog)s contour_lines points"
)
parser.add_argument("-cl","--contour_lines", help="csv output path from autocad that contains the lines of the contour map")
parser.add_argument("-p", "--points", help="csv path with the desired points")
parser.add_argument("-cx", "--center_x", type=float, help="x position of the center")
parser.add_argument("-cy", "--center_y", type=float, help="y position of the center")
parser.add_argument("-bx", "--beginning_x", type=float, help="x position of the beginning of the line")
parser.add_argument("-by", "--beginning_y", type=float, help="y position of the beginning of the line")
parser.add_argument("-ex", "--end_x", type=float, help="x position of the end of the line")
parser.add_argument("-ey", "--end_y", type=float, help="y position of the end of the line")
args = parser.parse_args()

def get_center_of_mass(points):
    """gets the center of mass"""
    x = 0
    y = 0
    for point in points:
        x += point[0]
        y += point[1]
    return x/len(points), y/len(points)

def get_polar(x,y):
    """given a vector, returns its degree"""
    theta = math.atan(y/x)

    if x < 0:
        theta += math.pi
    elif x > 0 and y < 0:
        theta += 2*math.pi
    
    return theta

def sort_section(mag_section):
    """mag_section is a list of tuple type (float, (float, float))
    
    """
    sorted_magsection =  sorted(mag_section)

    # return only the points
    return [p for _,p in sorted_magsection]

def get_directions_and_magnitudes(center_x, center_y, points):
    directions = []
    magnitudes = []

    unique_drs = []

    for point in points:
    
        # taking as origin the center
        x,y = point
        nx,ny = (x - center_x,y - center_y)

        # iterate over the dictionary keys to determine direction
        dir_max_mag = None
        for dir in unique_drs:
            magnitude_proj =nx*dir[0] + ny*dir[1]
            magnitude_vec = (nx**2 + ny**2)**(1/2)

            # the projection is at least 0.75 the magnitude of the vector
            # as they must be parallel
            if magnitude_proj > 0.75*magnitude_vec:
                dir_max_mag = dir

        # if it is the first element on the dict or the point doesnt follow
        # any of the previous directions, add the unitary direction
        if len(unique_drs) == 0 or dir_max_mag is None:
            len_np = (nx**2 + ny**2)**(1/2)
            unitary_dir = (nx/len_np, ny/len_np)

            unique_drs.append(unitary_dir)
            directions.append(unitary_dir)
            magnitudes.append(len_np)
        else:
            # append the point with the magnitude of its projection
            directions.append(dir_max_mag)
            magnitudes.append(magnitude_vec)
        
    return unique_drs, directions, magnitudes

def get_legs_with_index(legs_per_direction):
    """gets the legs and its associated index inside the direction"""
    legs = []
    for k in legs_per_direction.keys():
        mag_section = legs_per_direction[k]
        legs_per_direction[k] = sort_section(mag_section)

def get_sections_by_cardinal(center_x, center_y, points, line):
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
    unique_drs, directions, magnitudes = \
        get_directions_and_magnitudes(center_x, center_y, points)
    
    # find difference in radians from the directions and the line
    theta_l = get_polar(*line)

    # elements of the rotation matrix: sine and cos
    # negative to consider clockwise rotation
    R = {"s":math.sin(-theta_l), "c":math.cos(-theta_l)}

    # making a rotation for all the directions
    rotated_dirs = {
        dr:(dr[0]*R["c"] - dr[1]*R["s"], dr[0]*R["s"] + dr[1]*R["c"])  \
            for dr in unique_drs}
    
    cardinal_per_dr = {}
    for k,v in rotated_dirs.items():
        # nearest to the line to the right, corresponding to the 4th cuadrant
        if v[0] > 0 and v[1] < 0:
            cardinal_per_dr[k] = 1
        # 3th quadrant, following the clockwise convention for ennumeration
        elif v[0] < 0 and v[1] < 0:
            cardinal_per_dr[k] = 2
        # 2nd quadrant, following the clockwise convention for ennumeration
        elif v[0] < 0 and v[1] > 0:
            cardinal_per_dr[k] = 3
        # 1st quadrant, following the clockwise convention for ennumeration
        elif v[0] > 0 and v[1] > 0:
            cardinal_per_dr[k] = 4

    # converting directions to the appropriate cardinal
    cardinals = [cardinal_per_dr[d] for d in directions]
    return cardinals, magnitudes

def format_interpolated_result(output):
    """creates a new file with the appropriate format for the given points
    the points have to be in a specific order??

    File is saved in the folder where the script is run
    
    Parameter
    output: pandas.DataFrame
        dataframe containing the information of the interest points in X,Y,Z cols
    """
    base = """L.T. XXX XXXXXXXXXXXXXXXXXXXXXXXX        
Lev. TOP. ALGUIEN: MONTH/YY
TORRE XXXX cuerpo +XX            TORRE No.- DDD   KM. 
    19.200  3.70000     46.00      0.30
     14.00     14.00      14.00     14.00\n"""
    
    # create rows, in each row, the element considered belong to 1-4 sections
    # therefore, we require all the groupings at the same time
    s1 = output.loc[output["Section"]== 1].reset_index(drop=True)
    s2 = output.loc[output["Section"]== 2].reset_index(drop=True)
    s3 = output.loc[output["Section"]== 3].reset_index(drop=True)
    s4 = output.loc[output["Section"]== 4].reset_index(drop=True)

    # it is required that
    # 1. The index + 1 is present at the beginning, with format DD.
    #    If is only unity then ` D`
    # 2. 3 slots, one for the sign minus if needed and other for decimals
    # 3. digit, last digit is in index 9
    result = base
    for i in range(13):
        ind  = str(i+1).rjust(2, " ")
        if i < 9:
            s = ""
            # two decimal points always, and 7 spaces
            s += format(s1["Z"].iloc[i],".2f").rjust(7, " ")
            s += "   2.00"
            s += format(s2["Z"].iloc[i],".2f").rjust(7, " ")
            s += "   2.00"
            s += format(s3["Z"].iloc[i],".2f").rjust(7, " ")
            s += "   2.00"
            s += format(s4["Z"].iloc[i],".2f").rjust(7, " ")
            s += "   2.00\n"

        result += ind + s
    print(result)
    return 

def main(contour_map_path, points_path, line, cx, cy):
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

    #cx, cy = get_center_of_mass(points)
    sections, magnitudes = get_sections_by_cardinal(cx,cy, points, line)

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
    # round to two decimal places
    output["Z"] = output["Z"].round(2)
    output["Section"] = sections
    output["Magnitude"] = magnitudes
    output = output.sort_values(["Section", "Magnitude"])
    #print(output)
    #print(f"center: x={cx}, y={cy}, z={center_interpolated}")
    format_interpolated_result(output)
    return output

if __name__ == "__main__":
    args.contour_lines
    # start - end to have the appropriate vector
    line = (
        args.beginning_x-args.end_x,
        args.beginning_y-args.end_y
    )
    # check that both files exist
    if (os.path.isfile(args.contour_lines) and os.path.isfile(args.points)):
        main(
            args.contour_lines,
            args.points,
            line,
            args.center_x,
            args.center_y
            )
    elif not os.path.isfile(args.contour_lines):
        print(f"[ERROR] Contour map file `{args.contour_lines}` was not found")
    else:
        print(f"[ERROR] Points file `{args.points}` was not found")