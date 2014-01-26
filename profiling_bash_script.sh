sudo apt-get install python-pip
sudo apt-get kcachegrind
sudo pip install pyprof2calltree
python -m cProfile -o outputfile_profile.pyprof ode3.py
pyprof2calltree -i outputfile_profile.pyprof -k
