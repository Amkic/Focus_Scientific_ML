import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as linalg
import copy as cp
import scipy.linalg as spla
import seaborn as sns
import scipy.linalg as splin

def deim(Phi):
    """
    Discrete Empirical Interpolation Method (DEIM).
    Computes interpolation indices to accelerate non-linear term evaluation.
    """
    n, m = Phi.shape
    i = np.argmax(abs(Phi[:, 0]))

    U = Phi[:, [0]]
    P = np.zeros((n, 1))
    P[i] = 1.0

    i_list = [i]
    r_list = [Phi[:, [0]]]
    uc_list = [None]

    for l in range(1, m):
        u_l = Phi[:, l]

        c = splin.solve(P.T @ U, P.T @ u_l)
        r = u_l - U @ c
        i = np.argmax(abs(r))

        i_list.append(i)
        r_list.append(r)
        uc_list.append(U @ c)

        U = Phi[:, : l + 1]
        P = np.concatenate([P, np.zeros((n, 1))], 1)
        P[i, -1] = 1.0

    return P, i_list, r_list, uc_list

def greedy_basis(snapshots,K):
    """
    Constructs a reduced basis using a Greedy algorithm approach.
    """
    basis = np.zeros((snapshots.shape[0],K))
    basis[:,0] =snapshots[:,0]/np.linalg.norm(snapshots[:,0])
    error= np.zeros(snapshots.shape[1])
    projector = np.zeros_like(snapshots)
    for k  in range(1,K):
        error= np.zeros(snapshots.shape[1])
        projector = ((basis.T@snapshots).T)@(basis.T) 
        error = np.linalg.norm(snapshots.T - projector, axis=1)
        imax=np.argmax(error)
       
        current_phi =snapshots[:,imax]
        
        proj = (basis)@ (basis.T@ current_phi[:])
        basis[:,k]=(current_phi[:] - proj)/np.linalg.norm(current_phi[:] - proj)
    return basis

def apply_basis(basis,vec):
    return (basis)@ (basis.T@ vec)

class Eq_full:
    """
    Full Order Model (FOM) class for PDE solving.
    """
    def __init__(self,Nx,Re,f,df):
        self.Nx=Nx  # number of cells
        self.a=0.0  # left boundary
        self.b=1.0  # right boundary
        self.Re=Re     # Reynolds number
        self.dt=0.0
        self.f=f
        self.df=df
                
        # Mesh setup with Nx cells and 2 ghost cells on each side for boundary conditions
        self.h=(self.b-self.a)/self.Nx  # mesh step
        self.begin = 2  # index of the first non-ghost cell
        self.end=  self.Nx+2  # index of the last non-ghost cell  
       
        self.cells=np.linspace(self.a-1.5*self.h,self.b+1.5*self.h,self.Nx+4)
        # list of cell centers: n centers and 2 ghost centers
        self.rho= np.zeros((self.Nx,1)) ## rho on mesh + ghost cells
        self.flux= None
        self.rho_next= np.zeros((self.Nx,1)) 
    
    def init(self):
        b= self.begin
        e= self.end
        x0=(self.a+self.b)/2.0-0.3
        self.rho[:,0]=(1.0/10.0)*np.exp(-(self.cells[b:e]-x0)*(self.cells[b:e]-x0)/0.004)+0.5
    
        # store initial rho
        self.rho_0=cp.copy(self.rho)
        
        
    def slope(self,p):
        """Slope limiter function."""
        return 2.0*(p+np.abs(p))/(p+3)#max(0.0,min(1.0,p))
    

    def flot(self,rho_c):
        """Calculates numerical fluxes."""
        eps=0.0000000001
        #### Bc
        rho=np.zeros((self.Nx+4,(rho_c.shape)[1]))
        rho[0,:]=rho_c[0,:]
        rho[1,:]=rho_c[0,:]
        rho[self.begin:self.end,:]=rho_c[:,:]
        rho[-1,:]=rho_c[-1,:]
        rho[-2,:]=rho_c[-1,:]
        
        diff_im = rho[self.begin-1:self.end-1]-rho[self.begin-2:self.end-2]
        diff_i = rho[self.begin:self.end]-rho[self.begin-1:self.end-1]
        diff_ip = rho[self.begin+1:self.end+1]-rho[self.begin:self.end]
        diff_ipp = rho[self.begin+2:self.end+2]-rho[self.begin+1:self.end+1]
                    
        pl=  diff_im/ (diff_i + eps)
        p=  diff_i/(diff_ip + eps)
        pr=  diff_ip/(diff_ipp + eps)
        
        ##################################
        rl = (rho[self.begin-1:self.end-1]+0.5*self.slope(pl) *diff_i)
        rr = (rho[self.begin:self.end]-0.5*self.slope(p) * diff_ip)
        
        center_L = 0.5*(self.f(rl)+self.f(rr))
        vel_L = np.maximum(np.abs(self.df(rl)),np.abs(self.df(rr)))
        vis_L= (vel_L*0.5)*(rl-rr) -(1.0/(self.Re*self.h))* diff_i
        
        ##################################
        rl = rho[self.begin:self.end]+0.5*self.slope(p)* diff_ip
        rr = rho[self.begin+1:self.end+1]-0.5*self.slope(pr)*diff_ipp
    
        center_R = 0.5*(self.f(rl)+self.f(rr))
        vel_R = np.maximum(np.abs(self.df(rl)),np.abs(self.df(rr)))
        vis_R= (vel_R*0.5)*(rl-rr) -(1.0/(self.Re*self.h))*(diff_ip)
                        
        self.flux= -(center_R+vis_R - (center_L+vis_L))/self.h


    # plot rho_0 and rho
    def plot(self):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        axes[0].plot(self.cells[self.begin:self.end],self.rho_0,color="red",linewidth=3)
        axes[0].title.set_text('rho(0)')
        axes[1].plot(self.cells[self.begin:self.end],self.rho,color="blue",linewidth=3)
        axes[1].title.set_text('rho(t)')
        plt.show()
        
    def solve(self,T):
        """Standard FOM solver."""
        sol=[]
        nt = 0
        time = 0.0
        self.dt = 0.4*min(self.h,self.Re*self.h*self.h)
        sol.append(cp.copy(self.rho))
        while time < T:
            if time+self.dt>T:
                self.dt=T-time
      
            # explicit flux calculation
            self.flot(self.rho)
            self.rho_next=self.rho +0.5*self.dt*self.flux
            
            self.flot(self.rho_next)
            self.rho_next=self.rho +self.dt*self.flux

            self.rho=cp.copy(self.rho_next)
            sol.append(self.rho)
            time=time+self.dt
            nt=nt+1
        sol= np.concatenate(sol,axis=1)
        return sol
    
def f_adv(u):
    return u

def df_adv(u):
    return np.ones_like(u)

def f_burgers(u):
    return 0.5*u*u

def df_burgers(u):
    return u

class Reduced:
    """
    Reduced Order Model (ROM) class using POD and DEIM.
    """
    def __init__(self,snapshots,model,n_modes,red,hyp,n_modes_DEIM=10,affine=False):
        self.snap= snapshots
        self.snap_0= cp.copy(snapshots)
        self.model =model 
        self.hyp = hyp
        self.n_modes=n_modes
        self.rho_ref= np.zeros_like(self.snap[:,0])
        self.red = red
        
        if self.red=="POD":
            self.s_ref = np.zeros_like(self.snap[:,0])
            if affine==True:
                self.s_ref = np.zeros_like(self.snap[:,0])
                self.s_ref=np.sum(self.snap,axis=1)/self.snap.shape[1]
                for i in range(0,self.snap.shape[1]):
                    self.snap[:,i]= self.snap_0[:,i]-self.s_ref[:]

            podmodes, svals, _ = spla.svd(self.snap, full_matrices=False)              
            self.phi= podmodes[:,:self.n_modes]
        
            if self.hyp=="DEIM":
                self.n_modes_DEIM=n_modes_DEIM
                self.fsnap = np.zeros_like(self.snap)
                self.model.flot(self.snap_0)
                self.fsnap=self.model.flux

                podmodes, svals, _ = spla.svd(self.fsnap, full_matrices=False) 
                self.phif= podmodes[:,:self.n_modes_DEIM]
                self.Z,_ ,_ ,_ = deim(self.phif)
                self.Pi_ob = self.phif @ np.linalg.inv(self.Z.T @ self.phif) @ self.Z.T
        else:
            self.phi=greedy_basis(self.snap,self.n_modes)
            self.s_ref = np.zeros_like(self.snap[:,0])
            
        self.rho_hat= np.zeros(n_modes)
        self.rho_hat_next= np.zeros(n_modes)
      
    def init(self):
        self.rho_hat= self.phi.T @ (self.snap_0[:,0] - self.s_ref[:])
        
    def full_sol(self,sol_n):
        res= self.phi @ sol_n + self.s_ref
        return res
    
    def solve(self,Tf):
        """ROM solver with optional Hyper-reduction (DEIM)."""
        sol=[]
        nt = 0
        time = 0.0
        rho_hat_temp = np.zeros(self.n_modes)
        self.dt = 0.4*min(self.model.h,self.model.Re*self.model.h*self.model.h)
        sol.append(cp.copy(self.rho_hat))
        if self.hyp=='No' or self.red=="Greedy":
            while time < Tf:
                if time+self.dt> Tf:
                    self.dt=Tf-time            
            
                self.model.flot((self.phi @ self.rho_hat+ self.s_ref)[:,np.newaxis] )
                self.rho_hat_next=self.rho_hat +0.5*(self.dt
                    * self.phi.T @ self.model.flux[:,0]) 
                
                self.model.flot((self.phi @ self.rho_hat_next+ self.s_ref)[:,np.newaxis] )
                self.rho_hat_next=self.rho_hat +(self.dt
                    * self.phi.T @ self.model.flux[:,0]) 
                  
                self.rho_hat=cp.copy(self.rho_hat_next)
                sol.append(self.rho_hat)
                time=time+self.dt
                nt=nt+1
        else:
            while time < Tf:
                if time+self.dt> Tf:
                    self.dt=Tf-time            
            
                self.model.flot((self.phi @ self.rho_hat+ self.s_ref)[:,np.newaxis] )
 
                self.rho_hat_next=self.rho_hat +(self.dt
                    * self.phi.T @ self.Pi_ob @ self.model.flux[:,0]) 
                  
                self.rho_hat=cp.copy(self.rho_hat_next)
                sol.append(self.rho_hat)
                time=time+self.dt
                t=nt+1
            
        sol= np.transpose(np.array(sol))
        return sol    
        
    def plot(self,sol):
        b=self.model.begin
        e=self.model.end
        print("Error >>",np.linalg.norm(self.full_sol(sol[:,-1])-self.snap_0[:,-1]))
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes[0].plot(self.model.cells[b:e],self.full_sol(sol[:,-1]),c="blue",linewidth=3)
        axes[0].plot(self.model.cells[b:e],self.snap_0[:,-1],c="red",linewidth=3)
        axes[0].title.set_text('rho(t)')
        for k in range(0,self.n_modes):
            axes[1].plot(self.model.cells[b:e],self.phi[:,k],linewidth=2)
        axes[1].title.set_text('modes')
        sns.set_palette("Spectral",10)
        plt.show() 
        
        
def plot_reduced(model,Re,n_modes,red,hyp,affine):
    Tf=0.8
    if model=="Burgers":
        sim = Eq_full(1000,Re,f_burgers,df_burgers)
    else:
        sim = Eq_full(1000,Re,f_adv,df_adv)
    sim.init()
    snapshots=sim.solve(Tf)

    model_red = Reduced(snapshots,sim,n_modes,red,hyp,n_modes_DEIM=n_modes,affine=affine)
    model_red.init()
    res=model_red.solve(Tf)
    model_red.plot(res)
