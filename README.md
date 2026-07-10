Steps you need to do, as per https://gitlab.com/openfoam/core/openfoam/-/wikis/precompiled/debian:

```
# Add the repository
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash

# Update the repository information
sudo apt-get update

# Install preferred package. Eg,
sudo apt-get install openfoam2512-default

# Use the openfoam shell session. Eg,
openfoam2512
```

create a python venv
```
sudo apt install python3.12-venv

python -m venv ~/your/venv/path

pip3 install foamlib numpy fluidfoam pillow scikit-image
```

I add a line to bashrc to autoload the environment variables, e.g. (in bashrc, your path may be different):

```source /usr/lib/openfoam/openfoam2512/etc/bashrc```

Then
```./Allwmake ```
from the top directory.

Now you should be able to use the package!
