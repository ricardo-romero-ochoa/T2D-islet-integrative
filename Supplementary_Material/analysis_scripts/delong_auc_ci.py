import numpy as np, pandas as pd
from scipy.stats import norm

def compute_midrank(x):
    J=np.argsort(x); Z=x[J]; N=len(x); T=np.zeros(N); i=0
    while i<N:
        j=i
        while j<N and Z[j]==Z[i]: j+=1
        T[i:j]=0.5*(i+j-1)+1
        i=j
    out=np.empty(N); out[J]=T; return out

def fast_delong(preds, label_1_count):
    m=label_1_count; n=preds.shape[1]-m; k=preds.shape[0]
    tx=np.empty((k,m)); ty=np.empty((k,n)); tz=np.empty((k,m+n))
    for r in range(k):
        tx[r]=compute_midrank(preds[r,:m]); ty[r]=compute_midrank(preds[r,m:]); tz[r]=compute_midrank(preds[r])
    aucs=tz[:,:m].sum(axis=1)/(m*n)- (m+1)/(2*n)
    v01=(tz[:,:m]-tx)/n; v10=1-(tz[:,m:]-ty)/m
    sx=np.cov(v01); sy=np.cov(v10)
    cov=sx/m+sy/n
    return aucs, np.atleast_2d(cov)

def ci(y,score,alpha=.95):
    y=np.asarray(y).astype(int); score=np.asarray(score,float)
    order=np.argsort(-y); y=y[order]; score=score[order]
    m=int(y.sum()); aucs,cov=fast_delong(score[None,:],m); auc=float(aucs[0]); se=float(np.sqrt(cov[0,0])); z=norm.ppf(1-(1-alpha)/2)
    return auc,max(0,auc-z*se),min(1,auc+z*se),se
p='/mnt/data/_audit_reanalysis/final_repos_run/t2d-repo-run/results/external_validation/gse50244_fixed_threshold/external_reduced_scores.csv'
d=pd.read_csv(p); y=(d.group=='T2D').astype(int).values; s=d.reduced_score.values
print(ci(y,s))
