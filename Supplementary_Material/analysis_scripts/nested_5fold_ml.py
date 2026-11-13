import os, warnings, json, math
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import LinearSVC, SVC
from sklearn.feature_selection import RFECV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, f1_score, matthews_corrcoef, accuracy_score
warnings.filterwarnings('ignore')

BASE='/mnt/data/_audit_reanalysis/final_repos_run/t2d-repo-run'
OUT='/mnt/data/_audit_reanalysis/new_results/nested_ml'
os.makedirs(OUT, exist_ok=True)
RANDOM_STATE=42
DEG_FDR=0.01
DEG_LOGFC=1.5
LASSO_FOLDS=10
SVM_FOLDS=5
SVM_STEP=0.1
RF_TREES=500
TOP_N=50
FINAL_N=10

expr=pd.read_csv(BASE+'/data/processed/GSE164416_expr_normalized.csv',index_col=0)
ldf=pd.read_csv(BASE+'/data/processed/GSE164416_labels.csv')
labels=pd.Series(ldf.label.values,index=ldf.sample_id.values)
common=[c for c in expr.columns if c in labels.index and labels[c] in [0,1]]
expr=expr[common]
y=labels[common].astype(int).values
samples=np.array(common)
G=np.array(expr.index)
Xall=expr.values.T.astype(float)
print('data',Xall.shape,'T2D',y.sum(),'ND',(y==0).sum())

def bh(p):
    p=np.asarray(p); n=len(p); order=np.argsort(p); out=np.empty(n,float)
    ranked=p[order]*n/np.arange(1,n+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]
    out[order]=np.minimum(ranked,1.0); return out

def deg_genes(Xtr,ytr):
    a=Xtr[ytr==1]; b=Xtr[ytr==0]
    mean1=a.mean(0); mean0=b.mean(0); fc=mean1-mean0
    var1=a.var(0,ddof=1); var0=b.var(0,ddof=1)
    df=len(a)+len(b)-2
    sp2=((len(a)-1)*var1+(len(b)-1)*var0)/df
    prior_var=np.exp(np.median(np.log(sp2+1e-8))); prior_df=max(3,len(sp2)/10)
    shr=(prior_df*prior_var+df*sp2)/(df+prior_df)
    se=np.sqrt(shr*(1/len(a)+1/len(b)))
    t=fc/(se+1e-10); p=2*stats.t.sf(np.abs(t),df=df+prior_df); q=bh(p)
    idx=np.where((q<DEG_FDR)&(np.abs(fc)>=DEG_LOGFC))[0]
    return idx, fc, q

def select_panel(X,y,genes, seed):
    # LASSO -- training-only scaling and internal CV
    sc=StandardScaler(); Xs=sc.fit_transform(X)
    min_class=min(np.bincount(y))
    lfold=min(LASSO_FOLDS,max(3,min_class))
    cv=StratifiedKFold(n_splits=lfold,shuffle=True,random_state=seed)
    lasso=LogisticRegressionCV(Cs=np.logspace(-4,2,30),cv=cv,penalty='l1',solver='liblinear',max_iter=2000,
                               random_state=seed,scoring='roc_auc',n_jobs=-1)
    lasso.fit(Xs,y); lsel=np.abs(lasso.coef_[0])>0; lscore=np.abs(lasso.coef_[0])
    # SVM-RFE training only
    sfold=min(SVM_FOLDS,max(3,min_class))
    cv2=StratifiedKFold(n_splits=sfold,shuffle=True,random_state=seed)
    step=max(1,int(X.shape[1]*SVM_STEP))
    rfe=RFECV(LinearSVC(C=1.0,max_iter=5000,random_state=seed),step=step,cv=cv2,scoring='roc_auc',min_features_to_select=min(5,X.shape[1]),n_jobs=-1)
    rfe.fit(Xs,y); ssel=rfe.support_; srank=rfe.ranking_.astype(float)
    # RF importance averaged across internal folds
    cv3=StratifiedKFold(n_splits=min(5,max(3,min_class)),shuffle=True,random_state=seed)
    imps=[]
    for k,(tr,_) in enumerate(cv3.split(Xs,y)):
        rf=RandomForestClassifier(n_estimators=RF_TREES,max_features='sqrt',n_jobs=-1,random_state=seed+k,class_weight='balanced')
        rf.fit(Xs[tr],y[tr]); imps.append(rf.feature_importances_)
    imp=np.mean(imps,0); order=np.argsort(-imp); rsel=np.zeros(len(genes),bool); rsel[order[:min(TOP_N,len(genes))]]=True
    votes=lsel.astype(int)+ssel.astype(int)+rsel.astype(int)
    lnorm=lscore/(lscore.max()+1e-10)
    snorm=1-(srank-1)/(srank.max()+1e-10)
    rnorm=imp/(imp.max()+1e-10)
    comp=(lnorm+snorm+rnorm)/3
    cand=np.where(votes>=2)[0]
    if len(cand)<FINAL_N:
        cand=np.argsort(-comp)[:min(FINAL_N,len(comp))]
    else:
        cand=cand[np.argsort(-comp[cand])[:FINAL_N]]
    return cand, votes, comp

def build_ensemble(seed):
    return Pipeline([('scaler',StandardScaler()),('clf',VotingClassifier(estimators=[
        ('svm',SVC(kernel='rbf',probability=True,class_weight='balanced',random_state=seed)),
        ('rf',RandomForestClassifier(n_estimators=300,class_weight='balanced',random_state=seed,n_jobs=-1)),
        ('lr',LogisticRegression(C=0.1,class_weight='balanced',max_iter=2000,random_state=seed)),
        ('gb',GradientBoostingClassifier(n_estimators=100,max_depth=3,random_state=seed))],voting='soft'))])

probs=np.zeros(len(y)); pred=np.zeros(len(y),int); fold_rows=[]; sel_rows=[]; fold_metrics=[]
outer=StratifiedKFold(n_splits=5,shuffle=True,random_state=20261021)
for fold,(tr,te) in enumerate(outer.split(Xall,y),1):
    Xtr=Xall[tr]; ytr=y[tr]; Xte=Xall[te]; yte=y[te]
    didx,fc,q=deg_genes(Xtr,ytr)
    if len(didx)<10:
        didx=np.argsort(q)[:max(10,min(184,len(q)))]
    Xd=Xtr[:,didx]; gd=G[didx]
    pidx,votes,comp=select_panel(Xd,ytr,gd,RANDOM_STATE+fold)
    selected=gd[pidx]
    model=build_ensemble(RANDOM_STATE+fold)
    model.fit(Xtr[:,didx][:,pidx],ytr)
    pr=model.predict_proba(Xte[:,didx][:,pidx])[:,1]
    probs[te]=pr; pred[te]=(pr>=0.5).astype(int)
    f_auc=roc_auc_score(yte,pr) if len(np.unique(yte))==2 else np.nan
    fold_metrics.append({'fold':fold,'n_test':len(te),'n_deg_train':len(didx),'auc':f_auc,'panel_genes':';'.join(selected)})
    for j,ix in enumerate(te):
        fold_rows.append({'fold':fold,'test_sample':samples[ix],'true_label':int(y[ix]),'prob_t2d':float(pr[j]),'n_deg_train':len(didx),'n_panel':len(pidx),'panel_genes':';'.join(selected)})
    for g in selected: sel_rows.append({'fold':fold,'gene':g})
    print('fold',fold,'nDEG',len(didx),'AUC',round(f_auc,4),'panel',','.join(selected[:3]),'...')

auc=roc_auc_score(y,probs)
cm=confusion_matrix(y,pred,labels=[0,1]); tn,fp,fn,tp=cm.ravel()
metrics={'nested_loocv_auc':auc,'sensitivity':tp/(tp+fn),'specificity':tn/(tn+fp),'f1':f1_score(y,pred),'mcc':matthews_corrcoef(y,pred),'accuracy':accuracy_score(y,pred),'n':len(y),'n_t2d':int(y.sum()),'n_nd':int((y==0).sum())}
# stratified bootstrap distribution of OOF AUC
rng=np.random.default_rng(20261021); aucs=[]
pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
for b in range(5000):
    idx=np.r_[rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]
    aucs.append(roc_auc_score(y[idx],probs[idx]))
aucs=np.array(aucs)
metrics.update({'bootstrap_auc_median':float(np.median(aucs)),'auc_ci_2_5':float(np.quantile(aucs,.025)),'auc_ci_97_5':float(np.quantile(aucs,.975))})

fold_df=pd.DataFrame(fold_rows); fold_df.to_csv(OUT+'/nested_5fold_predictions_and_panels.csv',index=False); pd.DataFrame(fold_metrics).to_csv(OUT+'/nested_5fold_fold_metrics.csv',index=False)
sel=pd.DataFrame(sel_rows); stab=sel.groupby('gene').agg(n_folds=('fold','nunique')).reset_index(); stab['selection_frequency']=stab.n_folds/5; stab=stab.sort_values(['selection_frequency','gene'],ascending=[False,True]); stab.to_csv(OUT+'/nested_feature_selection_stability.csv',index=False)
pd.DataFrame({'auc':aucs}).to_csv(OUT+'/nested_auc_bootstrap_distribution.csv',index=False)
metrics['fold_auc_mean']=float(np.nanmean([x['auc'] for x in fold_metrics])); metrics['fold_auc_sd']=float(np.nanstd([x['auc'] for x in fold_metrics],ddof=1)); metrics['fold_auc_min']=float(np.nanmin([x['auc'] for x in fold_metrics])); metrics['fold_auc_max']=float(np.nanmax([x['auc'] for x in fold_metrics])); pd.DataFrame([metrics]).to_csv(OUT+'/nested_5fold_metrics.csv',index=False)
with open(OUT+'/nested_5fold_summary.json','w') as f: json.dump(metrics,f,indent=2)
print('\nMETRICS',json.dumps(metrics,indent=2))
print('\nTOP STABILITY\n',stab.head(25).to_string(index=False))
