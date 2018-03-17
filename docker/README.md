## Overview ##
These folders contain the Dockerfiles and scripts needed to build Docker images.

### Jupyter Image ###

This image allows you to use the storm-analysis package in Jupyter notebooks to analyze data on your computer.

For example:

`$ docker run -it --rm -p 8888:8888 -v /path/to/my/data/:/home/jovyan/work/share zhuanglab/jupyter-sa`

When Jupyter starts you will see:

```
--/sa_notebooks
  /share
  /image_date.txt 
  /sa_version.txt
```

* `sa_notebooks` - storm-analysis example notebooks.
* `share` - local folder on your computer
* `image_data.txt` - record of when this image was made.
* `sa_version.txt` - record of storm-analysis git version.

This container also supports all the options listed here:
[scipy-notebook](https://github.com/jupyter/docker-stacks/tree/master/scipy-notebook)


### Other Images ###

These images are designed to run a particular type of STORM analysis. However each of them can also be used to run anything in the storm-analysis project, with the exception of Visualizer as the images do not include Qt5/PyQt5.

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
