from datetime import timedelta
import pathlib
from typing import Union
from adcircpy.server.confing import ServerConfig


class SlurmConfig(ServerConfig):

    def __init__(
        self,
        account: str,
        walltime: timedelta,
        modules: [str] = None,
        **kwargs
    ):
        super().__init__(**kwargs)
        self._account = account
        self._walltime = walltime
        self._modules = modules

    def write(self, path: Union(str, pathlib.Path)):
        raise NotImplementedError

    def _deploy_files_to_server(self):
        # avoids checking wdir on server.
        raise NotImplementedError
