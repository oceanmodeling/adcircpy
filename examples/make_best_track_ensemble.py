#! /usr/bin/env python3
"""
Script to extract an ATCF best track dataset and modify it 
e.g., Make it more intense or larger in size

Example for Sandy2012 and prescribed start and end dates 
that are valid. 
- central_pressure is decreased by 10%,
- max_sustained_wind_speed is increased by 10%
- radius_of_maximum_winds is increased by 10%

By William Pringle, Mar 2021 - 
"""

from copy import deepcopy
from datetime import datetime, timedelta

from adcircpy.forcing.winds.best_track import BestTrackForcing


def main():
    # set storm name
    storm_name = "Sandy2012"

    # set simulation dates
    start_date = datetime(2012, 10, 22)
    end_date = start_date + timedelta(days=5)

    # getting best track
    BT = BestTrackForcing(
            storm_name,
            start_date=start_date,
            end_date=end_date,
    )

    # write out original fort.22
    BT.write("original.22", overwrite=True)

    # extracting original dataframe   
    df_original = BT.df

    # modifying the neccessary variables and 
    # writing each to a new fort.22
    variable_list = ["central_pressure", "max_sustained_wind_speed",
                     "radius_of_maximum_winds"]
    alpha = [0.9, 1.1, 1.1]  # the multiplier for each variable
    for idx, var in enumerate(variable_list):
        print(var)
        # make a deepcopy to preserve the original dataframe
        df_modified = deepcopy(df_original)
        df_modified[var] = df_modified[var] * alpha[idx]
        # reset the dataframe
        BT._df = df_modified
        # write out the modified fort.22
        BT.write(var + ".22", overwrite=True)


if __name__ == '__main__':
    main()
