from datetime import timedelta

class PBS(object):
  def __init__(self, account, walltime, numprocs, queue='batch', description=None, email=None):
    self._account = account
    self._walltime = walltime
    self._numprocs = numprocs
    self._queue = queue
    self._description = description
    self._email = email
    self.__init_walltime()
    self.__init_PBS()

  def __init_walltime(self):
    if isinstance(self._walltime, timedelta)==False:
      raise Exception('walltime must be a datetime.timedelta instance.')
    total_seconds = int(self._walltime.total_seconds())
    hours, remainder = divmod(total_seconds,60*60)
    minutes, seconds = divmod(remainder,60)
    self._walltime = '{}:{}:{}'.format(hours, minutes, seconds)

  def __init_PBS(self):
    self._PBS = list()
    self._PBS.append('#!/bin/bash --login\n')
    self._PBS.append('#PBS -d .\n')
    if self._description is not None:
      self._PBS.append('#PBS -N {}\n'.format(self._description))
    self._PBS.append('#PBS -A {}\n'.format(self._account))
    if self._email is not None:
      self._PBS.append('#PBS -m be\n')
      self._PBS.append('#PBS -M {}\n'.format(self._email))
    self._PBS.append('#PBS -j oe\n')
    self._PBS.append('#PBS -o job.${PBS_JOBID}.oe\n')
    self._PBS.append('#PBS -l procs={}\n'.format(self._numprocs))
    self._PBS.append('#PBS -l walltime={}\n'.format(self._walltime))
    self._PBS.append('#PBS -q {}\n'.format(self._queue))
