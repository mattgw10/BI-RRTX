import math
import numpy as np

import grid
import constants

# Matthew Westbrook

# Functions for connecting/rewiring/propagating vertices
	
def connect(C, blocked, x0, y0, x1, y1):
	# Create edge from 0->1 and 1->0
	cells = grid.intersected_cells(C, x0, y0, x1, y1)
	
	connected = True
	interfered = []
	for i in range(0, len(cells)):
		if cells[i] in blocked:
			connected = False
			interfered.append(cells[i])
	
	return [connected, interfered, cells]
	
def reconnect(v, v_info, eo_out, er_out):
	lmc = float("inf")
	p = -1
	for u in eo_out:
		if lmc > eo_out[u]+v_info[u][1]:
			p = u
			lmc = eo_out[u]+v_info[u][1]
	for u in er_out:
		if lmc > er_out[u]+v_info[u][1]:
			p = u
			lmc = er_out[u]+v_info[u][1]
	if p == -1:
		return False
	else:
		make_parent_of(v, p, v_info)
		return True
	
def find_parent(C, blocked, v, v_loc, v_info, v_near, r):
	# Find parent of v
	lmc = float("inf")
	p = -1
	for u in v_near[0]:
		d = grid.ed(v[0], v[1], v_loc[u,0], v_loc[u,1])
		[connected, interfered, cells] = connect(C, blocked, v[0], v[1],  v_loc[u,0], v_loc[u,1])
		if (d <= r) and (d > 0) and (lmc > d + v_info[u][1]) and connected:
			p = u
			lmc = d + v_info[u][1]
	return [p, lmc]
	
def make_parent_of(u, v, v_info):
	# Change u's parent to v, add u to v's children, and delete
	# u as child of old parent
	if u in v_info[v_info[u][2]][3]:
		del v_info[v_info[u][2]][3][u]
	v_info[u][2] = v
	v_info[v][3][u] = 1

def updateLMC(v, v_info, r, er_out, eo_out, er_in, v_orphan):
	cull_neighbors(v, v_info, r, er_out, er_in)
	for u in er_out:
		if (v_info[v][1] > er_out[u]+v_info[u][1]) and not (v_info[u][2] == v) \
			and (u not in v_orphan):
			v_info[v][1] = er_out[u]+v_info[u][1]
			make_parent_of(v, u, v_info)
	for u in eo_out:
		if (v_info[v][1] > eo_out[u]+v_info[u][1]) and not (v_info[u][2] == v) \
			and (u not in v_orphan):
			v_info[v][1] = eo_out[u]+v_info[u][1]
			make_parent_of(v, u, v_info)
	
def verify_queue(v, g, Q):
	# Verify Queue
	index = [i for i, x in enumerate(Q) if x[0] == v]
	if len(index) > 0:
		Q[index[0]][1] = g
	else:
		Q.append([v, g])
	Q.sort(key=lambda x: x[1])
	
def verify_orphan(Q, v, v_orphan):
	# Verify Orphan List
	index = [i for i, x in enumerate(Q) if x[0] == v]
	for i in index:
		Q.pop(i)
	if not v in v_orphan:
		v_orphan.append(v)

def extend(C, blocked, e_blocked, v, v_loc, v_info, r, V, eoi, eoo, eri, ero, v_cell, e_cell):
	# Extend vertex to neighbors and get parent
	v_near = grid.near(v, v_loc, r)
	if len(v_near) > 0:
		[p, lmc] = find_parent(C, blocked, v, v_loc, v_info, v_near, r)
	else:
		return [False, v_loc]
	
	if p == -1:
		return [False, v_loc]
	v_loc = np.append(v_loc, [[v[0], v[1]]], axis = 0)
	v_info.append([lmc, lmc, p, {}, True])
	v_info[p][3][V] = 1
	v_cell[grid.coord2grid(C, v[0], v[1])].append(V)
	
	eoi.append({})
	eoo.append({})
	eri.append({})
	ero.append({})

	# Make connections to all neighbors
	for i in range(0, len(v_near)):
		u = v_near[i][0]
		# Connect
		[connected, interfered, cells] = connect(C, blocked, v[0], v[1], v_loc[u,0], v_loc[u,1])
		if connected:
			d = grid.ed(v[0], v[1], v_loc[u,0], v_loc[u,1])
			eoo[V][u] = d
			eri[u][V] = d
			ero[u][V] = d
			eoi[V][u] = d
			for cell in cells:
				e_cell[cell].append([V, u])
		else:
			for cell in interfered:
				e_blocked[cell].append([V, u])
	return [True, v_loc]

def cull_neighbors(v, v_info, r, e_out, e_in):
	# Remove edges from v->u if distance > radius
	remove = []
	for u in e_out:
		# If edge long and not parent edge
		if (r < e_out[u]) and not (v_info[v][2] == u) and not (v_info[u][2] == v):
			remove.append(u)
	for u in remove:
		del e_out[u]
		del e_in[u]

def rewire_neighbors(Q, v, v_info, r, eo_in, er_in, er_out):
	# Rewire neighbors to make better path
	if v_info[v][0]-v_info[v][1] > constants.eps:
		cull_neighbors(v, v_info, r, er_out, er_in)
	for u in eo_in:
		if (v_info[u][1] > eo_in[u]+v_info[v][1]) and not (v_info[u][2] == v):
			v_info[u][1] = eo_in[u]+v_info[v][1]
			make_parent_of(u, v, v_info)
			if v_info[u][0]-v_info[u][1] > constants.eps:
				verify_queue(u, v_info[u][0], Q)
	for u in er_in:
		if (v_info[u][1] > er_in[u]+v_info[v][1]) and not (v_info[u][2] == v):
			v_info[u][1] = er_in[u]+v_info[v][1]
			make_parent_of(u, v, v_info)
			if v_info[u][0]-v_info[u][1] > constants.eps:
				verify_queue(u, v_info[u][0], Q)
					
def reduce_inconsistency(Q, v_bot, v_info, eoi, eri, eoo, ero, r, v_orphan):
	# Cascading repairs
	cascade = True
	while (len(Q) > 0) and cascade:
		if not v_bot == -1:
			if not ((Q[0][1] > v_info[v_bot][1]) or not (v_info[v_bot][0] == v_info[v_bot][1]) \
				or (v_info[v_bot][0] == float("inf"))):
				cascade = False
				break
		[v, d] = Q.pop(0)
		
		if (v_info[v][0]-v_info[v][1] > constants.eps):
			updateLMC(v, v_info, r, ero[v], eoo[v], eri[v], v_orphan)
			rewire_neighbors(Q, v, v_info, r, eoi[v], eri[v], ero[v])
		v_info[v][0] = v_info[v][1]
		
def propagate_descendents(Q, v_orphan, v_info, eo_out, er_out):
	# If vertex cut off, all children are cut off
	orphan_new = []
	for v in v_orphan:
		for c in v_info[v][3]:
			if c not in v_orphan:
				v_orphan.append(c)
	
	for v in v_orphan:
		for u in eo_out[v]:
			if u not in v_orphan:
				v_info[u][0] = float("inf")
				verify_queue(u, float("inf"), Q)
		for u in er_out[v]:
			if u not in v_orphan:
				v_info[u][0] = float("inf")
				verify_queue(u, float("inf"), Q)
				
	for v in v_orphan:
		v_orphan.remove(v)
		v_info[v][0] = float("inf")
		v_info[v][1] = float("inf")
		if not v_info[v][2] == -1:
			del v_info[v_info[v][2]][3][v]
			v_info[v][2] = -1