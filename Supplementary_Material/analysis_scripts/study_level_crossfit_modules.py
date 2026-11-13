from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
from rds_miniparser import read_matrix

BASE=Path('/mnt/data/_audit_reanalysis/final_repos_run/t2d-islet')
OUT=Path('/mnt/data/_audit_reanalysis/new_results/crossfit_modules'); OUT.mkdir(parents=True,exist_ok=True)
STUDIES=['GSE25724','GSE20966','GSE38642','GSE164416']
rng=np.random.default_rng(20261021)

def clean(x): return str(x).strip().upper()

def dl_meta3(tabs):
    # intersect genes among three training cohorts
    maps=[]
    for d in tabs:
        dd=d.copy(); dd['gene2']=dd.gene.map(clean); dd=dd.drop_duplicates('gene2').set_index('gene2')
        maps.append(dd)
    genes=sorted(set(maps[0].index)&set(maps[1].index)&set(maps[2].index))
    Y=np.column_stack([m.loc[genes,'logFC'].to_numpy(float) for m in maps])
    S=np.column_stack([m.loc[genes,'SE'].to_numpy(float) for m in maps])
    ok=np.isfinite(Y).all(1)&np.isfinite(S).all(1)&(S>0).all(1)
    genes=np.array(genes)[ok]; Y=Y[ok]; S=S[ok]
    w=1/S**2; sw=w.sum(1); muf=(w*Y).sum(1)/sw
    Q=(w*(Y-muf[:,None])**2).sum(1); C=sw-(w**2).sum(1)/sw; tau=np.maximum(0,(Q-2)/C)
    wr=1/(S**2+tau[:,None]); swr=wr.sum(1); mu=(wr*Y).sum(1)/swr; se=np.sqrt(1/swr); z=mu/se; p=2*norm.sf(np.abs(z))
    same=(Y>0).all(1)|(Y<0).all(1)
    return pd.DataFrame({'gene':genes,'meta_effect':mu,'meta_se':se,'z':z,'p':p,'tau2':tau,'direction_consistent':same})

def cohens_d(x,g):
    a=np.asarray(x)[np.asarray(g)=='ND']; b=np.asarray(x)[np.asarray(g)=='T2D'];
    a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return (b.mean()-a.mean())/sp if np.isfinite(sp) and sp>0 else np.nan

def load_expr_meta(study):
    m,genes,samples=read_matrix(BASE/'data/processed/expression'/f'{study}_expression_genelevel.rds')
    genes=np.array([clean(x) for x in genes]); samples=np.array(samples)
    meta=pd.read_csv(BASE/'data/processed/metadata'/f'{study}_metadata_harmonized.csv').set_index('sample_id')
    keep=np.array([s in meta.index and bool(meta.loc[s,'include_primary_analysis']) and meta.loc[s,'disease_group'] in ['ND','T2D'] for s in samples])
    m=m[:,keep]; samples=samples[keep]; groups=np.array([meta.loc[s,'disease_group'] for s in samples])
    # unique gene rows
    _,idx=np.unique(genes,return_index=True); idx=np.sort(idx); genes=genes[idx]; m=m[idx]
    return m,genes,samples,groups

def zrows(m):
    mu=np.nanmean(m,1,keepdims=True); sd=np.nanstd(m,1,ddof=1,keepdims=True); sd[~np.isfinite(sd)|(sd==0)]=1
    z=(m-mu)/sd; z[~np.isfinite(z)]=0; return z

degs={s:pd.read_csv(BASE/'data/processed/deg'/f'{s}_DE_full.csv') for s in STUDIES}
rows=[]; members=[]; nullrows=[]
for held in STUDIES:
    train=[s for s in STUDIES if s!=held]
    m,genes,samples,g=load_expr_meta(held); Z=zrows(m); gi={x:i for i,x in enumerate(genes)}
    mt=dl_meta3([degs[s] for s in train]); mt=mt[(mt.direction_consistent)&(mt.p<0.05)&(mt.gene.isin(set(genes)))]
    up=mt[mt.meta_effect>0].assign(absz=lambda x:abs(x.z)).sort_values('absz',ascending=False).head(14)
    dn=mt[mt.meta_effect<0].assign(absz=lambda x:abs(x.z)).sort_values('absz',ascending=False).head(14)
    if len(up)<14 or len(dn)<14: raise RuntimeError((held,len(up),len(dn)))
    for module,dd in [('CrossFitImmuneStress',up),('CrossFitIdentitySecretion',dn)]:
        for _,r in dd.iterrows(): members.append({'held_out':held,'training_cohorts':';'.join(train),'module':module,**r.drop(labels='absz').to_dict()})
    upg=[x for x in up.gene if x in gi]; dng=[x for x in dn.gene if x in gi]
    su=Z[[gi[x] for x in upg]].mean(0); sd=Z[[gi[x] for x in dng]].mean(0); sc=su-sd
    obs={'ImmuneStress':cohens_d(su,g),'IdentitySecretion':cohens_d(sd,g),'IsletDysfunctionScore':cohens_d(sc,g)}
    for mod,dval in obs.items(): rows.append({'held_out':held,'training_cohorts':';'.join(train),'module':mod,'n_ND':int((g=='ND').sum()),'n_T2D':int((g=='T2D').sum()),'n_genes_up_available':len(upg),'n_genes_down_available':len(dng),'cohens_d':dval})
    # null candidates matched on mean expression and nonzero/finite detection rate
    mean=np.nanmean(m,1); det=np.mean(np.isfinite(m)&(m!=0),1); valid=np.isfinite(mean)&np.isfinite(det)
    pool_idx=np.where(valid & ~np.isin(genes,upg+dng))[0]
    pm=mean[pool_idx]; pdet=det[pool_idx]; sm=np.std(pm); sd_det=np.std(pdet); sm=sm if sm>0 else 1; sd_det=sd_det if sd_det>0 else 1
    zpm=(pm-pm.mean())/sm; zpd=(pdet-pdet.mean())/sd_det
    def candidate_lists(targets):
        out=[]
        for tg in targets:
            ix=gi[tg]; tm=(mean[ix]-pm.mean())/sm; td=(det[ix]-pdet.mean())/sd_det
            dist=(zpm-tm)**2+(zpd-td)**2; near=pool_idx[np.argsort(dist)[:100]]; out.append(near)
        return out
    cup=candidate_lists(upg); cdn=candidate_lists(dng)
    B=2000; nd_up=np.empty(B); nd_dn=np.empty(B); nd_comp=np.empty(B)
    def draw(cands):
        chosen=[]
        for c in cands:
            avail=np.setdiff1d(c,np.array(chosen,dtype=int),assume_unique=False)
            if len(avail)==0: avail=c
            chosen.append(int(rng.choice(avail)))
        return chosen
    for b in range(B):
        ru=draw(cup); rd=draw(cdn); a=Z[ru].mean(0); bsc=Z[rd].mean(0)
        nd_up[b]=cohens_d(a,g); nd_dn[b]=cohens_d(bsc,g); nd_comp[b]=cohens_d(a-bsc,g)
    for mod,nv in [('ImmuneStress',nd_up),('IdentitySecretion',nd_dn),('IsletDysfunctionScore',nd_comp)]:
        o=obs[mod]; emp=(1+np.sum(np.abs(nv)>=abs(o)))/(B+1)
        nullrows.append({'held_out':held,'module':mod,'observed_d':o,'null_mean':np.nanmean(nv),'null_sd':np.nanstd(nv,ddof=1),'null_q025':np.nanquantile(nv,.025),'null_q975':np.nanquantile(nv,.975),'empirical_p_two_sided':emp,'n_perm':B})
    print(held,obs,'available',len(upg),len(dng))

pd.DataFrame(rows).to_csv(OUT/'crossfit_heldout_effect_sizes.csv',index=False)
pd.DataFrame(members).to_csv(OUT/'crossfit_module_membership.csv',index=False)
pd.DataFrame(nullrows).to_csv(OUT/'matched_gene_set_permutation_null.csv',index=False)
print('\nEFFECTS'); print(pd.DataFrame(rows).to_string(index=False)); print('\nNULL'); print(pd.DataFrame(nullrows).to_string(index=False))
