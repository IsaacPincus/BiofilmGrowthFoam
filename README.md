Steps you need to do:

sudo apt install openfoam2512-default

create a python venv
sudo apt install python3.12-venv

python -m venv ~/your/venv/path

pip3 install foamlib numpy fluidfoam pillow scikit-image

I add a line to bashrc to autoload the environment variables, e.g. (in bashrc, your path may be different):
source /usr/lib/openfoam/openfoam2512/etc/bashrc

Then ./Allwmake 

Now you should be able to use the package!
