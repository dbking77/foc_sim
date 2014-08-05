Field-Orientied-Controller (FOC) Simulator for Sinusoidal Brushless Motor.
Motor Simulator and FOC Controller is written in Python should be easy to modify
for different motors or differnt simulation scenarios.

Running ./motor_sim will run motor and FOC simulation, and generate a set of
plots using matplotlib.  The simulation settings and plots can be changed by 
editing the code.  The code has been tested on Ubuntu 14.04 (Trusty Tahr).

Files:
 motor_sim.py : Main motor simulator loop and graphical ploting functions
 foc.py : Python implementati of mostly standard FOC controller

To get matplotlib in Ubuntu:
sudo apt-get install python-matplotlib