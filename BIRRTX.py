import argparse
import importlib
import logging
import sys
import math
import copy
import random
import time
import numpy as np
from collections import namedtuple

import robot
import grid
import graph
import obstacles
import printing
import constants

# Matthew Westbrook

# Batch Informed Rapidly Exploring Randomizing Tree (BIRRT^x)
#  -Navigate obstacles with real time planner and 
#  time varying obstacles. Graph expands from random sampling
#  (batch informed or normal RRT random) and improves path
#  when possible to within epsilon consistency. Current path to
#  goal from each node makes a tree as in RRT. Planning is done
#  from goal to start and re-planning is done as neccessary.
#  Time to reach goal is affected by speed of connecting robot
#  to root of tree, amount that plan is able to be improved during 
#  traversal, and robot speed.

def get_batch(R, C, blocked, bound, p1, p2):
	QV = []
	for i in range(0,constants.m):
		[sx, sy] = [random.uniform(0, C-.1), random.uniform(0, R-.1)]
		if grid.coord2grid(C, sx, sy) not in blocked:
			d = grid.ed(p1[0],p1[1], sx, sy) + grid.ed(sx, sy, p2[0],p2[1])
			if d < bound:
				QV.append([d, [sx, sy]])
	QV.sort(key=lambda x: x[0])
	return QV
	
def search(alg, R, C, blocked, unknown, unknown_loc, x0, y0, gx, gy):
	# Tracks what edges and vertices are in each grid cell
	v_cell = [[] for i in range(0, R*C)]
	e_cell = [[] for i in range(0, R*C)]
	
	# Tracks vertices that could be connected without obstacles in the way
	e_blocked = [[] for i in range(0, R*C)]
	
	# Initialize: Create start conditions for planning
	
	# Start position
	x = x0
	y = y0
	v_bot = -1
	
	# Start time
	alg_start = time.time()
	# Is bot connected to tree
	v_connect = False
	
	# No movement until initial plan found
	moving = False
	
	# Put goal location in graph
	V = 1
	#  -Vertex locations as numpy array
	v_loc = np.array([[gx, gy]])
	#  -Vertex additional info structure is gval, lmcval, parent vertex, children, validity (incase later blocked)
	v_info = [[0.0, 0.0, -1, {}, True]]
	#  -Outgoing edge in structure is list of dictionaries for each vertex 
	#   where key is other vertex and value is edge length
	eoi = [{}]
	#  -Outgoing edge out
	eoo = [{}]
	#  -Running edge in
	eri = [{}]
	#  -Running edge out
	ero = [{}]
	
	v_cell[grid.coord2grid(R, x, y)].append(0)
	
	# Queue of vertices to rewire
	Q = []
	# List of disconnected vertices
	v_orphan = []
	
	# Batch of samples if using batch informed
	QV = []
	cost = float("inf")
	
	# Save first map after start connected
	first_saved = False
	block_remove = []
	block_add = []
				
	loops = 0
	t_print = time.time()
	
	# Loop until robot close to goal
	while (grid.ed(x, y, gx, gy) > constants.GOAL_RADIUS):
		# Neighborhood radius
		r = grid.shrinking_ball(V)
		# Update robot and obstacle if it has started moving
		if moving:
			cost = v_info[v_bot][0]
			v_bot_old = v_bot
			[v_bot, p, x, y, th, movement_start, movement_time] = robot.update_bot(movement_start, movement_time, v_bot, eri, p, v_loc, v_info, x, y, th)
			if p == -1:
				x = gx
				y = gy
				break
			[block_remove, block_add, unknown_loc, moving] = obstacles.check_obstacles(C, blocked, unknown_loc, unknown, Q, v_bot, v_cell, e_cell,\
			[x, y], e_blocked, v_loc, v_info, v_orphan, r, eoi, eoo, eri, ero)
			if len(block_remove) > 0 or len(block_add) > 0:
				printing.print_update(alg_start, block_remove, block_add)
			if not moving:
				cost = float("inf")
			else:
				if not v_bot == v_bot_old:
					printing.print_robot_path(v_bot, v_loc, v_info)
					t_print = time.time()
		# If not moving, see if robot is ready to move
		else:
			[success, v_loc] = robot.try_to_start(v_loc, x, y, C, blocked, e_blocked, v_info, r, V, eoi, eoo, eri, ero, v_cell, e_cell, Q, v_orphan)
			if success:
				if not first_saved:
					first_saved = True
					started_moving = time.time()
					print("STARTED MOVING IN: " + str(started_moving-alg_start))
					printing.print_map(v_loc, v_info, "first_path.csv")
				v_bot = V
				V += 1
				v_connect = True
				moving = True
				movement_start = time.time()
				p = v_info[v_bot][2]
				movement_time = grid.ed(v_loc[v_bot,0],v_loc[v_bot,1],v_loc[p,0],v_loc[p,1])/constants.s
				th = grid.angle(v_loc[v_bot,0], v_loc[v_bot,1], v_loc[p,0], v_loc[p,1])
		# Batch informed sampling or random
		if alg == 1:
			# New sample batch if empty or all remaining are useless
			if len(QV) == 0:
				QV = get_batch(R, C, blocked, cost, [x, y], [gx, gy])
			elif QV[0][0] > cost:
				QV = get_batch(R, C, blocked, cost, [x, y], [gx, gy])
			if len(QV) > 0:
				[d, [sx, sy]] = QV.pop(0)
			else:
				continue
		else:
			[sx, sy] = [random.uniform(0, C-.1), random.uniform(0, R-.1)]
			
		# Nearest vertex
		v_near = grid.nearest([sx, sy], v_loc)
		if grid.ed(v_loc[v_near,0], v_loc[v_near,1], sx, sy) > constants.delta:
			[sx, sy] = grid.saturate([sx, sy], v_loc[v_near])
		# Extend graph if possible
		if grid.coord2grid(C, sx, sy) not in blocked:
			[success, v_loc] = graph.extend(C, blocked, e_blocked, [sx, sy], v_loc, v_info, r, V, eoi, eoo, eri, ero, v_cell, e_cell)
			if success:
				graph.rewire_neighbors(Q, V, v_info, r, eoi[V], eri[V], ero[V])
				graph.reduce_inconsistency(Q, v_bot, v_info, eoi, eri, eoo, ero, r, v_orphan)
				V += 1
		block_remove = []
		block_add = []
		loops += 1
	
	print(str(V) + " NODES EXPANDED")
	print("REACHED GOAL IN: " + str(time.time()-alg_start))
	printing.print_map(v_loc, v_info, "final_path.csv")

def main():
	# Determine rrtx or birrtx
	alg = 0
	for line in sys.argv:
		if line == "-BIRRTX":
			alg = 1
	
	# Read in map
	#  -first two lines are rows and columns
	#  -third line is robot start location
	#  -forth line is goal location
	#  -rest of lines are initial guess of map and true map
	
	# Map size
	COLUMNS = int(sys.stdin.readline())
	ROWS = int(sys.stdin.readline())
	
	# Initial location
	arg = sys.stdin.readline().rstrip().split()
	[x0, y0] = [float(arg[0]), float(arg[1])]
	
	# Goal location
	arg = sys.stdin.readline().rstrip().split()
	[gx, gy] = [float(arg[0]), float(arg[1])]

	# Blocked cells list (current knowledge of map)
	blocked = {}
	blocked_coord = []
	for i in range(0,ROWS):
		cells = sys.stdin.readline().rstrip()
		for j in range(0, COLUMNS):
			if (cells[j] == "#"):
				# space is blocked
				blocked[(ROWS-i-1)*COLUMNS+j] = 1
				blocked_coord.append([j, (ROWS-i-1)])
			# do nothing if space is empty
	printing.print_grid(blocked_coord, "first_grid.csv")
	
	# Cells that are different from initial knowledge
	sys.stdin.readline()
	unknown = {}
	unknown_loc = np.array([[2*ROWS, 2*COLUMNS]])
	blocked_coord = []
	for i in range(0,ROWS):
		cells = sys.stdin.readline().rstrip()
		for j in range(0, COLUMNS):
			if (cells[j] == "#"):
				blocked_coord.append([j, (ROWS-i-1)])
			# Check if different from initial knowledge
			if (cells[j] == "#") and (((ROWS-i-1)*COLUMNS+j) not in blocked):
				unknown[(ROWS-i-1)*COLUMNS+j] = 1
				unknown_loc = np.append(unknown_loc, [[j, ROWS-i-1]], axis = 0)
			elif (cells[j] == "_") and (((ROWS-i-1)*COLUMNS+j) in blocked):
				unknown[(ROWS-i-1)*COLUMNS+j] = 1
				unknown_loc = np.append(unknown_loc, [[j, ROWS-i-1]], axis = 0)
				
	search(alg, ROWS, COLUMNS, blocked, unknown, unknown_loc, x0, y0, gx, gy)
	printing.print_grid(blocked_coord, "final_grid.csv")

if __name__ == '__main__':
    main()