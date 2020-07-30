from pandas import read_csv


def atcf_id(storm_id: str):
    url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/storm.table'
    data_frame = read_csv(url, header=None)
    name = f'{storm_id[:-4].upper():>10}'
    year = f'{storm_id[-4:]:>5}'
    entry = data_frame[(data_frame[0].isin([name]) & data_frame[8].isin([year]))]
    if len(entry) == 0:
        return None
    else:
        return entry[20].tolist()[0].strip()
