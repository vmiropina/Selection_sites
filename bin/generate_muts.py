import math
import mpmath as mp
import numpy as np
import scipy.special as sp

def f(n, r, s_obs, theta, model):
    L = 1
    #r = w_s*m/l_s
    if model == 1.0:
        a,b = theta
        def numerator(s, thr):
            return np.exp(n*np.log(r) + (n + s)*np.log(b) + (-n - s - a)*np.log(1 + b + b*r) + sp.gammaln(n + s + a) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) ) 
        def pofs(s, L, thr):
            return np.exp( s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) 
    elif model == 2.0:
        a,b = theta
        def numerator(s, thr):
            if thr:
                return  mp.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + a))*math.log(b) -(0.5*(n+s-a)*math.log(1+r)) + mp.log(mp.besselk(n + s - a, 2. * math.sqrt(b+b*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) )
            else:
                return np.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + a))*math.log(b) -(0.5*(n+s-a)*math.log(1+r)) + np.log(sp.kv(n + s - a, 2. * math.sqrt(b+b*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) )
        def pofs(s, L, thr):
            if thr:
                return 2. *  mp.exp( ((s + a)/2.)*np.log(L*b) + mp.log(mp.besselk(-s + a, 2*mp.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
            else:
                return  2. * np.exp( ((s + a)/2.)*np.log(L*b) + np.log(sp.kv(-s + a, 2*np.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
    elif model==3.0:
        a,b,t,w = theta
        def numerator(s, thr):
            return w * t* np.exp((-1. -n -s)* np.log(1 + r + t) + n* np.log(r)  + sp.gammaln(1 + n + s ) - sp.gammaln(s+1) - sp.gammaln(n+1)) + (1-w) * np.exp(n*np.log(r) + (n + s)*np.log(b) + (-n - s - a)*np.log(1 + b + b*r) + sp.gammaln(n + s + a) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) ) 
        def pofs(s, L, thr):
            return np.exp( np.log(w) + s*np.log(L) + np.log(t) + (-1 - s)*np.log(L + t) ) + np.exp( np.log(1.-w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) )
    elif model==4.0:
        a,b,t,w = theta
        def numerator(s, thr):
            if thr:
                return w * t* np.exp((-1. -n -s)* np.log(1 + r + t) + n* np.log(r)  + sp.gammaln(1 + n + s ) - sp.gammaln(s+1) - sp.gammaln(n+1)) + (1-w) * mp.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + a))*math.log(b) -(0.5*(n+s-a)*math.log(1+r)) + mp.log(mp.besselk(n + s - a, 2. * math.sqrt(b+b*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) )
            else:
                return w * t* np.exp((-1. -n -s)* np.log(1 + r + t) + n* np.log(r)  + sp.gammaln(1 + n + s ) - sp.gammaln(s+1) - sp.gammaln(n+1)) + (1-w) * np.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + a))*math.log(b) -(0.5*(n+s-a)*math.log(1+r)) + np.log(sp.kv(n + s - a, 2. * math.sqrt(b+b*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) )
        def pofs(s, L, thr):
            if thr:
                return (w * t * mp.exp( s*np.log(L) + (-1 - s)*np.log(L + t) )) + mp.exp( np.log(1.-w) + np.log(2.) + ((s + a)/2.)*np.log(L*b) + mp.log(mp.besselk(-s + a, 2*math.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
            else:
                return (w * L**s * t * (L + t)**(-1 - s)) + np.exp( np.log(1.-w) + np.log(2.) + ((s + a)/2.)*np.log(L*b) + np.log(sp.kv(-s + a, 2*math.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )

    elif model == 5.0:
        a,b,g,d,w = theta
        def numerator(s, thr):
            return w * np.exp(n*np.log(r) + (n + s)*np.log(b) + (-n - s - a)*np.log(1 + b + b*r) + sp.gammaln(n + s + a) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) ) + (1-w) * np.exp(n*np.log(r) + (n + s)*np.log(d) + (-n - s - g)*np.log(1 + d + d*r) + sp.gammaln(n + s + g) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(g) ) 
        def pofs(s, L, thr):
            return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1-w) + s*np.log(L*d) + (-s-g)*np.log(1 + L*d) + sp.gammaln(s + g) - sp.gammaln(s+1) - sp.gammaln(g) )

    elif model == 6.0:
        a,b,g,d,w = theta
        def numerator(s, thr):
            if thr:
                return w * np.exp(n*np.log(r) + (n + s)*np.log(b) + (-n - s - a)*np.log(1 + b + b*r) + sp.gammaln(n + s + a) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) ) + (1-w) * mp.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + g))*math.log(d) -(0.5*(n+s-g)*math.log(1+r)) + mp.log(mp.besselk(n + s - g, 2. * math.sqrt(d+d*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(g) )
            else:
                return w * np.exp(n*np.log(r) + (n + s)*np.log(b) + (-n - s - a)*np.log(1 + b + b*r) + sp.gammaln(n + s + a) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(a) ) + (1-w) * np.exp( n*np.log(r) + math.log(2.) + (0.5 * (n + s + g))*math.log(d) -(0.5*(n+s-g)*math.log(1+r)) + np.log(sp.kv(n + s - g, 2. * math.sqrt(d+d*r))) - sp.gammaln(s+1) - sp.gammaln(n+1) - sp.gammaln(g) )
        def pofs(s, L, thr):
            if thr:
                return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + mp.exp( np.log(1.-w) + np.log(2.) + ((s + g)/2.)*np.log(L*d) + mp.log(mp.besselk(-s + g, 2*mp.sqrt(L*d))) - sp.gammaln(s+1) - sp.gammaln(g) )
            else:
                return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1.-w) + np.log(2.) + ((s + g)/2.)*np.log(L*d) + np.log(sp.kv(-s + g, 2*np.sqrt(L*d))) - sp.gammaln(s+1) - sp.gammaln(g) )
    
    if s_obs>25:
        thr = 1
    else:
        thr = 0

    return numerator(s_obs, thr)/pofs(s_obs,L,thr)


def generate_n(r, s_obs, theta, model, reps, max):
	possible_num_muts = list(range(max))
	probs_num_muts = np.array([float(f(n,r, s_obs, theta, model)) for n in possible_num_muts])
	#print(np.nan_to_num(probs_num_muts, nan = np.nan, neginf = 0, posinf = 0))    
	multinomial = np.random.multinomial(1, np.nan_to_num(probs_num_muts, nan = np.nan, neginf = 0, posinf = 0), size=reps)
	n = np.where(multinomial==1)[1]
	return n

def generate_s(lambda_s, mi, ls, reps):
    return np.random.poisson(lambda_s*mi/ls, reps)
