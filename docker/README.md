
These folders contain the Dockerfiles and scripts needed to build Docker images.

The images are designed to run a particular type of STORM analysis. However each of them can also be used to run anything in the storm-analysis project, with the exception of Visualizer as the images do not include Qt5/PyQt5.

The images have the following structure:

* You need to source venv/bin/activate to start the Python3 virtual environment.
* The storm-analysis project is in the /storm-analysis directory.
* The /data directory is the volume that can be shared with the host os.

Linux example:

```sh
cd /directory/with/STORM/data
docker run -t -i -v $(pwd):/data -u `id -u`:`id -g` zhuanglab/base /bin/bash
source /venv/bin/activate
```

This part mounts the current directory as the /data directory in the docker image:

```sh
-v $(pwd):/data
```

This part sets user:group so that any files generated in the /data directory will be owned by the current user and not by root:

```sh
-u `id -u`:`id -g`
```

Note that any files that are generated inside the image will be lost when the image exits.
