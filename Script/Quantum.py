import numpy as np

class Quantum:
	def __init__(self,n_fermi,TBME):
		"""
		n_fermi - Number of electrons (number of states below fermi level)
		TBME - function which returns the two-body matrix elements of the Hamiltonian.
		"""
		self.n_fermi = n_fermi
		self.TBME = TBME
	def set_Z(self,Z):
		"""
		Z - Atomic number / coloumb interaction
		"""
		self.Z = Z
		self.u = self.TBME(Z)
	def set_states(self,states):
		"""
		states - list containing lists [n,m_s] (basis states), where n refers to quantum number n of the state
				 and m_s is the secondary spin number of the state.
		"""
		self.n_basis = len(states)
		self.states = states
	def ihj(self,i,j):
		if self.states[i,1] == self.states[j,1] and self.states[i,0] == self.states[j,0]:
			return(-self.Z**2/(2*self.states[i,0]**2))
		else:
			return(0)
	def ijvkl(self,i,j,k,l):
		if self.states[i,1] == self.states[k,1] and self.states[j,1] == self.states[l,1]:
			if self.states[i,0] == self.states[j,0] and self.states[i,1] == self.states[j,1]:
				return(0)
			elif self.states[k,0] == self.states[l,0] and self.states[k,1] == self.states[l,1]:
				return(0)
			else:
				return(self.u[self.states[i,0]-1,self.states[j,0]-1,\
						self.states[k,0]-1,self.states[l,0]-1])
		else:
			return(0)
	def cHc(self):
		expval = 0
		for i in range(self.n_fermi):
			expval += self.ihj(i,i)
			for j in range(self.n_fermi):
				expval += 0.5*(self.ijvkl(i,j,i,j) - self.ijvkl(i,j,j,i))
		return(expval)
	def cHia(self,i,a):
		expval = self.ihj(i,a)
		for j in range(self.n_fermi):
			expval += self.ijvkl(i,j,a,j) - self.ijvkl(i,j,j,a)
		return(expval)
	def iaHjb(self,i,a,j,b):
		expval = self.ijvkl(j,a,b,i) - self.ijvkl(j,a,i,b)
		if i == j:
			expval += self.ihj(a,b)
			for k in range(self.n_fermi):
				expval += self.ijvkl(a,k,b,k) - self.ijvkl(a,k,k,b)
			if a == b:
				for k in range(self.n_fermi):
					expval += self.ihj(k,k)
					for l in range(self.n_fermi):
						expval += 0.5*(self.ijvkl(k,l,k,l) - self.ijvkl(k,l,l,k))
		if a == b:
			expval -= self.ihj(j,i)
			for k in range(self.n_fermi):
				expval += self.ijvkl(j,k,k,i) - self.ijvkl(j,k,i,k)
		return(expval)

class CI(Quantum):
	"""
	Find CI (singles) approximation of eigenstates of the Hamiltonian
	"""
	def solve(self):
		n_fermi = self.n_fermi
		n_basis = self.n_basis
		n = n_fermi*(n_basis - n_fermi) + 1
		H = np.zeros((n,n))
		H[0,0] = self.cHc()
		for i in range(n_fermi):
			for j in range(n_fermi,n_basis):
				cHa = self.cHia(i,j)
				H[0,1+i*(n_basis-n_fermi)+j-n_fermi] = cHa
				H[1+i*(n_basis-n_fermi)+j-n_fermi,0] = cHa
				for k in range(n_fermi):
					for l in range(n_fermi,n_basis):
						H[j+1-n_fermi+i*(n_basis-n_fermi),k*(n_basis-n_fermi)+l+1-n_fermi] = self.iaHjb(i,j,k,l)
		self.H_CI = H
		self.e_CI, self.u_CI = np.linalg.eigh(H)
		return(self.e_CI,self.u_CI,self.H_CI)

class HF(Quantum):
	"""
	Find the Hartree-Fock approximation of the ground state energy
	"""
	def iter(self):
		n_basis = self.n_basis
		self.H_hf = np.zeros((n_basis,n_basis))
		for alpha in range(n_basis):
			for beta in range(n_basis):
				self.H_hf[alpha,beta] += self.ihj(alpha,beta)
				for j in range(self.n_fermi):
					for gamma in range(n_basis):
						for delta in range(n_basis):
							self.H_hf[alpha,beta] += self.C[j,gamma]*self.C[j,delta]*\
												(self.ijvkl(alpha,gamma,beta,delta) - self.ijvkl(alpha,gamma,delta,beta))
		self.e_1,C = np.linalg.eigh(self.H_hf)
		self.C = C.T
	def calculate_energies(self):
		n_basis = self.n_basis
		n_fermi = self.n_fermi
		e = 0
		for i in range(n_fermi):
			for alpha in range(n_basis):
				for beta in range(n_basis):
					e += self.C[i,alpha]*self.C[i,beta]*self.ihj(alpha,beta)
					for j in range(n_fermi):
						for gamma in range(n_basis):
							for delta in range(n_basis):
								e += 0.5*self.C[i,alpha]*self.C[j,beta]*self.C[i,gamma]*self.C[j,delta]\
									*(self.ijvkl(alpha,beta,gamma,delta) - self.ijvkl(alpha,beta,delta,gamma))
		return(e)

	def solve(self,tol = 1e-8,max_iter=1000):
		n_basis = self.n_basis
		n_fermi = self.n_fermi
		self.C = np.identity(n_basis)
		e_0 = np.zeros(n_basis)
		err=1
		converge=False
		for i in range(max_iter):
			self.iter()
			err=np.sum(np.abs(self.e_1 - e_0))/n_basis
			e_0 = self.e_1
			if err < tol:
				converge = True
				break
		if converge != True:
			print('HF did not converge')
		E = self.calculate_energies()
		return(self.e_1,self.C,self.H_hf, E)

