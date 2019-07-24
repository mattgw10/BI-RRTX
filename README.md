# BI-RRTX
Batch informed real time planning algorithm.

BIRRTX.py and graph.py are the key components of the algorithm to look at the code. If you want to run the code everything you need is in the zip file. The algorithm runs with: run.sh [algorithm (-BIRRTX or -RRTX, -RRTX is default)] < map.txt. The program prints the time to start traversal and time to reach the goal to stdout and the rest of the data is sent to csv files. The .m file is a utility for viewing the robot plan from the csv files.
