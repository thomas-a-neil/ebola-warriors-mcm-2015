from __future__ import division

import sys

import numpy as np
import time
import os


#simple_rho_vec = [0.0, 0.3, 0.6, 0.9]
rho_vec = np.append(np.append(np.linspace(0, 0.3, 3), np.linspace(0.35, 0.65, 10)), np.linspace(0.7, 0.95, 3))
p_vec = np.logspace(-1, 0, 10)

def runTests(N, K, p_vec, rho_vec, test_N):
	curDir = os.getcwd()
	for p in p_vec:
		for rho in rho_vec:
			death_rates = np.empty(test_N)
			start = time.clock()
			for i in range(test_N):
				s = SmallWorld(N, K, p)
				if rho != 0:
					# s.vaccinateRandomly(rho)
					s.vaccinatePrioritize(rho)
				s.seed_infection()
				s.run_percolation_all_infect()
				# path = os.path.normpath("%s/model_results/rho=%g_p=%g_N=%g_K=%g_test%g.csv" % (curDir, rho, p, N, K, i))
				#np.savetxt(path, np.transpose([s.t, s.I_t]), fmt = '%g', delimiter=',', comments = '')
				death_rates[i] = s.numRemov/(N*N - s.numVacc)
			surv_path = os.path.normpath("%s/model_results/death_rates_rho=%g_p=%g_N=%g_K=%g.csv" % (curDir, rho, p, N, K))
			np.savetxt(surv_path, np.transpose(death_rates), fmt = '%g', delimiter=',', comments = '')
			print("%s trials ran in %g seconds") % (test_N, (time.clock() - start))
	



class SmallWorld:
	"""
	A class, each instance of which is a randomly generated world of Nodes and Edges
	with some connectivity
	"""
	def __init__(self, N, K, p):
		self.N = N
		self.K = K
		self.p = p

		self.curT = 0

		self.numVacc = 0
		self.numInf = 0
		self.numSusc = N*N
		self.numRemov = 0

		self.t = []
		self.I_t = []

		self.infectedNodes = set()
		# if (N, K) not in self.gen_2d_lattices:
		# 	self.gen_2d_lattices[(N, K)] = self.make_2d_graph(N, K)
		# self.nodes, self.edges = self.gen_2d_lattices[(N, K)]
		# 
		# if i can figure out how to memoize and copy, I can come back to this. Copying is no better than creating! Same process.

		self.lattice, self.edges = make_2d_graph(N, K)
		self.nodes = self.lattice.values()
		self.reconnect_graph(self.edges, p)
		self.setNeighbors()

		while not areConnected(self.nodes):
			#remake graph if disconnected
			# print("disconnected") #TESTING
			self.lattice, self.edges = make_2d_graph(N, K)
			self.nodes = self.lattice.values()
			self.reconnect_graph(self.edges, p)
			self.setNeighbors()

	def reconnect_graph(self, edges, p):
		"""
		Reconnects the small-world network with list of Edges with randomness
		parameter p. Keeps total number of connections the same.
		"""
		# print("reconnecting") #TESTING
		if p == 0:
			return
		for e in edges:
			curNode = e.fromNode
			if np.random.random() < p:
 			    i, j = np.random.randint(self.N, size = 2)
 			    while self.lattice[(i, j)] in curNode.getNearestNeighbors():
    				i, j = np.random.randint(self.N, size = 2)
    			    e.reconnect(self.lattice[(i, j)])
				

	def setNeighbors(self):
		"""called when graph is no longer mutating"""
		for node in self.nodes:
			node.setNearestNeighbors()

	def vaccinateRandomly(self, rho):
		"""Called before the percolation, vaccinates randomly within the graph"""
		indices = set()

		while len(indices) <= self.N * self.N * rho:
			indices.add(tuple(np.random.randint(self.N, size = 2)))

		for i, j in indices:
			self.lattice[(i, j)].remove()
		self.numVacc = len(indices) + 1
		self.numSusc -= len(indices) + 1

	def vaccinatePrioritize(self, rho):
		l = list(self.nodes)
		sorted(l, key =lambda x: len(x.getNearestNeighbors()))
		num = int(np.floor(self.N*self.N*rho))
		for i in range(num):
			l[-i].remove()
		self.numVacc = num
		self.numSusc -= num

	def vaccinateIndividual(self, i, j):
		self.lattice[(i, j)].remove()
		self.numVacc += 1
		self.numSusc -= 1

	def seed_infection(self):
		"""
		Chooses a random first susceptible node to be infected. Assumes curT = 0
		"""
		i, j = np.random.randint(self.N, size = 2)
		while not self.try_infect(self.lattice[(i, j)]):
			i, j = np.random.randint(self.N, size = 2)

	def try_infect(self, node):
		"""
		Try to infect a Node at the current time, if successful, add to to list of infected nodes.
		Returns if successful
		"""
		worked = node.infect(self.curT)
		if worked:
			self.numInf += 1
			self.numSusc -= 1
			self.infectedNodes.add(node)
		return worked

	def stepTime(self):
		"""Steps time according to t = t + A/N_I(t), A = 1"""
		self.curT = self.curT + 1/self.numInf


	def run_percolation_all_infect(self, T = 3):
		"""After vaccinating and seeding infection, called to run percolation until no infected nodes remain"""
		self.curT = 0
		self.t = [0, ]
		self.I_t = [self.numInf, ]
		while len(self.infectedNodes) > 0:
			self.curT += 1
			# print("time is now %s" % self.curT) #TEST
			# print("there are %s infected (%s)" % (self.numInf, len(self.infectedNodes))) #TEST
			cur_infected = list(self.infectedNodes)
			for curNode in cur_infected:
				if (self.curT - curNode.t_inf) >= T:
					# print("removing node %s" % curNode) #TEST
					curNode.remove()
					self.numInf -= 1
					self.numRemov += 1
					self.infectedNodes.remove(curNode)
				else:
					n_n = curNode.getNearestNeighbors()
					p = 0.8
					#infects each neighbor with probability p - 
					for chosenNode in n_n:
						if np.random.rand() <= p:
							self.try_infect(chosenNode)
			self.t += [self.curT]
			self.I_t += [self.numInf]

	def run_percolation_rand_choice(self, T = 3):
		"""After vaccinating and seeding infection, called to run percolation until no infected nodes remain"""
		self.curT = 0
		self.t = [0, ]
		self.I_t = [self.numInf, ]
		while len(self.infectedNodes) > 0:
			self.stepTime()
			# print("time is now %s" % self.curT) #TEST
			# print("there are %s infected (%s)" % (self.numInf, len(self.infectedNodes))) #TEST
			curNode = np.random.choice(tuple(self.infectedNodes))
			if (self.curT - curNode.t_inf) >= T:
				# print("removing node %s" % curNode) #TEST
				curNode.remove()
				self.numInf -= 1
				self.numRemov += 1
				self.infectedNodes.remove(curNode)
			else:
				n_n = curNode.getNearestNeighbors()
				chosenNode = n_n[np.random.randint(len(n_n))]
				self.try_infect(chosenNode)
				# if self.try_infect(chosenNode):
				# 	print("infecting node %s" % chosenNode) #TEST
			self.t += [self.curT]
			self.I_t += [self.numInf]


def make_2d_graph(N, K):
	"""creates a network of size N**2 in 2d, where nodes are connected
	to their K nearest neighbors, N>>K. Returns a tuple of the list
	of Nodes and the list of Edges"""
	mat = {}
	edges = []
	for i in range(0, N):
		for j in range(0, N):
		#create the nodes in this row
			n = Node((i, j))
			mat[(i, j)] = n
	for i in range(0, N):
		for j in range(0, N):
			for k in range(1, K+1):
				for c in range(k+1):
					e = mat[(i, j)].addNewEdge(mat[ ((i + k - c) % N, (j + c) % N) ])
					edges += [e]
				for c in range(1, k):
					e = mat[(i, j)].addNewEdge(mat[ ((i- c) % N, (j + k - c) % N) ])
					edges += [e]
	return (mat, edges)

def areConnected(nodes):
	"""
	Takes a list of nodes and checks if they are connected
	"""
	stack = []
	stack.append(nodes[0])
	seen = []
	while len(stack) > 0:
		node = stack.pop()
		seen.append(node)
		for next in node.getNearestNeighbors():
			if next not in seen:
				stack.append(next)
	missing = False
	for node in nodes:
		if node not in seen:
			missing = True
	return not missing

class Node:
	"""
	Nodes can be infected, removed, etc.
	These are the building blocks of the simulation, they store their time of
	infection, they can find their nearest neighbors
	"""
	def __init__(self, coord):
		self.susc = True
		self.inf = False
		self.rem = False
		self.t_inf = False
		self.edges = set()
		self.nearestNeighbors = []
		self.inStatic = False
		self.coord = coord

	def __str__(self):
		return str(self.coord)

	def addNewEdge(self, toNode):
		"""Returns the added edge"""
		e = Edge(self, toNode)
		return e

	def addEdge(self, edge):
		self.edges.add(edge)

	def deleteEdge(self, edge):
		self.edges.remove(edge)

	def setNearestNeighbors(self):
		"""
		Called after a graph of which Node is a part has been modified to its final form, 
		and the list of nearest neighbors will no longer mutate
		"""
		self.nearestNeighbors = self.getNearestNeighbors()
		self.inStatic = True
		
	def unsetNearestNeighbors(self):
		self.inStatic = False

	def getNearestNeighbors(self):
		if self.inStatic:
			return self.nearestNeighbors
		else:
			neighbors = []
			for e in self.edges:
				neighbors += [e.getNeighbor(self)]
			return neighbors

	def infect(self, curTime):
		"""Returns the success of an infection, and stores the time of infection"""
		if self.susc:
			self.susc = False
			self.inf = True
			self.t_inf = curTime
			return True
		else:
			return False

	def remove(self):
		"""For removal of a node by death or vaccination"""
		self.rem = True
		self.susc = False
		self.inf = False

	def isInfected(self):
		return self.inf

	def isSusceptible(self):
		return self.susc

	def isRemoved(self):
		return self.rem

class Edge:
	def __init__(self, fromNode, toNode):
		self.fromNode = fromNode
		self.toNode = toNode
		fromNode.addEdge(self)
		toNode.addEdge(self)

	def reconnect(self, new_b):
		self.toNode.deleteEdge(self)
		self.toNode = new_b
		new_b.addEdge(self)

	def getNeighbor(self, curNode):
		if self.fromNode == curNode:
			return self.toNode
		else: 
			return self.fromNode


### ----- uncomment if running from command line ------
#if __name__ == '__main__':
#	N = int(sys.argv[1])
#	K = int(sys.argv[2])
#	p = float(sys.argv[3])
#	rho = float(sys.argv[4])
#	test_N = int(sys.argv[5])
#	runTests(N, K, [p], [rho], test_N)