#! /usr/bin/env python
import argparse

class GenerateHindcast(object):
  def __init__(self):
    self._parse_args()

  def _parse_args(self):
    # parser = argparse.
    print('gdhg')
    pass

def main():
  GenerateHindcast()

if __name__ == '__main__':
  main()
  



# import os
# import argparse
# from datetime import timedelta
# from AdcircPy.core.Winds import BestTrack

# def _parse_args():
#   global args
#   def __check_files():
#     _fort14 = args.project_dir+"/fort.14"
#     _fort13 = args.project_dir+"/fort.13"
#     if os.path.isfile(_fort14)==False or os.path.isfile(_fort13)==False:
#       raise FileNotFoundError("A fort.14 and a fort.13 are required to be present in the project directory.")
#   parser = argparse.ArgumentParser(description="Program generate ADCIRC input files for hindcast run.")
#   parser.add_argument("storm_id",  help="Storm id")
#   # parser.add_argument("project_dir",  help="Output path including the project's target directory.")
#   args =  parser.parse_args()
#   # __check_files()

# def _write_fort22():
#   fort22 = BestTrack(args.storm_id)
#   hotstart_time = fort22.end_time - fort22.start_time
#   coldstart_time = fort22.start_time - timedelta(days=15)
#   print(fort22.name)
#   print("Start time:     {}".format(coldstart_time.strftime("%-d %B %Y: %Hz")))
#   print("Wind Forcing:   {}".format(fort22.start_time.strftime("%-d %B %Y: %Hz")))
#   print("End time:       {}".format(fort22.end_time.strftime("%-d %B %Y: %Hz")))
#   _storm_duration = fort22.end_time - fort22.start_time
#   print("Storm Duration: {} days, {} hours".format(_storm_duration.days, int(_storm_duration.seconds/3600)))
#   _run_duration = fort22.end_time - coldstart_time
#   print("Run Duration:   {} days, {} hours".format(_run_duration.days, int(_run_duration.seconds/3600)))
#   # fort22.dump(args.project_dir+"/fort.22.best_track")

# def main():
#   _parse_args()

#   _write_fort22()

# if __name__ == "__main__":
#   main()