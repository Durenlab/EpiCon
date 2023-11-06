import pandas as pd
import numpy as np
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
def list2mat(df,TF_motif_set,RE_set,i_n,j_n,x_n):
    TFs = TF_motif_set
    REs = RE_set
#Initialize matrix as numpy array 
#Map row and col indices for lookup
    row_map = {r:i for i,r in enumerate(REs)}
    col_map = {c:i for i,c in enumerate(TFs)}
    row_indices = np.array([row_map[row] for row in df[i_n]])
    col_indices = np.array([col_map[col] for col in df[j_n]])
    from scipy.sparse import coo_matrix
    matrix = coo_matrix((df[x_n], (row_indices, col_indices)), shape=(len(REs), len(TFs)))
    mat=matrix.toarray()
    return mat
def TF_RE_correlation(TF_psedubulk,RE_psedubulk):
    import pandas as pd
    import numpy as np
    #TF_RE_cor = TF_psedubulk.T.corrwith(RE_psedubulk.T)
    X=TF_psedubulk.values
    Y=RE_psedubulk.values
    X_variance = np.std(X, axis=1)*np.sqrt(X.shape[1])
    Y_variance = np.std(Y, axis=1)*np.sqrt(X.shape[1])
    import warnings
    row_mean = X.mean(axis=1)
    # Subtract the row-wise mean from each element in X
    X = X - row_mean[:, np.newaxis]
    row_mean = Y.mean(axis=1)
    Y=Y-row_mean[:, np.newaxis]
    warnings.filterwarnings("ignore")
    TF_RE_cor=Y.dot(X.T)/(Y_variance[:, np.newaxis].dot(X_variance[:,np.newaxis].T))
    warnings.resetwarnings()
    TF_RE_cor=pd.DataFrame(TF_RE_cor,index=RE_psedubulk.index,columns=TF_psedubulk.index)
    sample_size=TF_psedubulk.shape[1]
    idxRE=(Y_variance>0)
    idxTF=(X_variance>0)
    return TF_RE_cor,idxRE,idxTF
def rm_outlayer(X):
    X[X<0]=-np.log(1-X[X<0])
    X[X>0]=np.log(1+X[X>0])
    X[X>4]=4
    X[X<-4]=-4
    return X

def impute(TFexp,ATAC):
    K=np.floor(np.sqrt(TFexp.shape[1]))+1
    from sklearn.impute import KNNImputer
    imputer = KNNImputer(n_neighbors=int(K))
    import scanpy as sc
    TF_psedubulk = imputer.fit_transform(sc.pp.log1p(TFexp.values.T, base=2))
    TF_psedubulk=pd.DataFrame(TF_psedubulk.T,index=TFexp.index,columns=TFexp.columns)
    RE_psedubulk = imputer.fit_transform(sc.pp.log1p(ATAC.values.T, base=2))
    RE_psedubulk=pd.DataFrame(RE_psedubulk.T,index=ATAC.index,columns=ATAC.columns)
    TF_RE_cor,idxRE,idxTF=TF_RE_correlation(TF_psedubulk,RE_psedubulk)
    return TF_RE_cor,idxRE,idxTF
def extend_region(outdir,N):
    region=pd.read_csv(outdir+'/Region.bed',sep='\t',header=None)
    region.columns=['chr','str','end','Name']
    region['str']=region['str']-N
    region['end']=region['end']+N
    region.loc[region['str'] < 0, 'str'] = 0
    region.to_csv(outdir+'/Region_200k.bed',sep='\t',header=None,index=None)

## region.bed extend 200k
def find_distance(outdir,N,code_dir,provide_dir,genome):
    import subprocess
    extend_region(outdir,N)
    subprocess.run(["sh", code_dir+"distance.sh",outdir,provide_dir, genome])
    data0=pd.read_csv(outdir+'/Distance0.txt',sep='\t',header=None)
    data0.columns=['RE','gene','TSS']
    data0['Start'] = data0['RE'].str.split(':|-').str[1].astype(int)
    data0['End'] = data0['RE'].str.split(':|-').str[2].astype(int)
    data0['dis']=np.abs(np.floor((data0['Start']+data0['End'])/2)-data0['TSS'])
    data0[['RE','gene','dis']].to_csv(outdir+'/distance.txt',sep='\t',header=None,index=None)

def RE_dis(distance_file,ATAC,RNA):
    distance=pd.read_csv(distance_file,sep='\t',header=None);
    distance.columns=['RE','TG','dis']
    TG_set=distance['TG'].unique()
    RE_set=distance['RE'].unique()
    distance=list2mat(distance,TG_set,RE_set,'RE','TG','dis')
    RE_set1= RE_set
    distance=pd.DataFrame(distance,index=RE_set1,columns=TG_set)
    TG_set=list(set(RNA.index)&(set(TG_set)))
    distance=distance[TG_set]
    distance1=np.zeros((ATAC.shape[0],len(TG_set)))
    index_RE=pd.DataFrame(range(ATAC.shape[0]),index=ATAC.index)
    distance1[index_RE.loc[RE_set1][0].values,:]=distance.values
    distance1=pd.DataFrame(distance1,index=ATAC.index,columns=TG_set)
    ATAC1=(np.exp(-distance1.values/25000-0.5)*(distance1>0)).dot(RNA.loc[TG_set])
    ATAC1=pd.DataFrame(ATAC1,index=ATAC.index,columns=ATAC.columns)
    return ATAC1


def ECS(outdir,ATAC,TF_motif_match,ATAC1,TF_RE_cor):
    MotifTarget=pd.read_csv(outdir+'/MotifTarget_u.txt',header=None,sep='\t')
    MotifTarget.columns=['RE','motif']
    MotifTarget['score']=1
    motifset=list(set(MotifTarget['motif']))
    RE_set=ATAC.index
    MotifTarget_m=list2mat(MotifTarget,motifset,RE_set,'RE','motif','score')
    TF_motif_match=TF_motif_match.loc[TF_motif_match['motif'].isin(motifset)]
    ATAC_E=ATAC.values.sum(axis=0).reshape(-1, 1).dot(ATAC.values.sum(axis=1).reshape( 1,-1))/ATAC.values.sum(axis=0).sum()
    ATAC_E1=ATAC1.values.sum(axis=0).reshape(-1, 1).dot(ATAC1.values.sum(axis=1).reshape( 1,-1))/ATAC1.values.sum(axis=0).sum()
    MotifTarget_m[MotifTarget_m>0]=1
    MotifTarget_m=pd.DataFrame(MotifTarget_m,index=RE_set,columns=motifset)
    TF_motif_match=TF_motif_match[TF_motif_match['TF'].isin(TF_RE_cor.columns)]
    MotifTarget_m1=MotifTarget_m.copy()
    MotifTarget_m1=MotifTarget_m1[TF_motif_match['motif'].values]
    #MotifTarget_m1.replace(0, np.nan, inplace=True)
    TF_RE_cor1=TF_RE_cor[TF_motif_match['TF'].values]
    TF_motif_RE_cor=TF_RE_cor1.values*MotifTarget_m1.values
    TF_motif_RE_cor_0=TF_motif_RE_cor.copy()
    TF_motif_RE_cor_0[TF_motif_RE_cor_0<0]=0
    ATAC_Y_TF_motif_n=(TF_motif_RE_cor_0.T.dot(ATAC.values.astype('float32'))-TF_motif_RE_cor_0.T.dot(ATAC_E.T))/TF_motif_RE_cor_0.T.dot(ATAC_E.T)
    ATAC_Y_TF_motif_n1=(TF_motif_RE_cor_0.T.dot(ATAC1.values.astype('float32'))-TF_motif_RE_cor_0.T.dot(ATAC_E1.T))/TF_motif_RE_cor_0.T.dot(ATAC_E1.T)
    ATAC_Y_TF_motif_n=rm_outlayer(ATAC_Y_TF_motif_n)
    ATAC_Y_TF_motif_n1=rm_outlayer(ATAC_Y_TF_motif_n1)
    ATAC_Y_TF_motif_n=pd.DataFrame(ATAC_Y_TF_motif_n,columns=ATAC.columns)
    ATAC_Y_TF_motif_n1=pd.DataFrame(ATAC_Y_TF_motif_n1,columns=ATAC.columns)
    #ATAC_Y_TF_motif_n.to_csv('ATAC_Y_TF_motif_n_0_.txt',sep='\t')
    #ATAC_Y_TF_motif_n1.to_csv('ATAC_Y_TF_motif_n_0_1.txt',sep='\t')
    (ATAC_Y_TF_motif_n+ATAC_Y_TF_motif_n1).to_csv(outdir+'/EpiCon_TF_motif.txt',sep='\t')
    ATAC_Y_TF_motif_n.index=TF_motif_match['TF'].values
    ATAC_Y_TF_motif_n1.index=TF_motif_match['TF'].values
    ATAC_Y_TF=ATAC_Y_TF_motif_n.groupby(ATAC_Y_TF_motif_n.index).mean()
    ATAC_Y_TF1=ATAC_Y_TF_motif_n1.groupby(ATAC_Y_TF_motif_n.index).mean()
    (ATAC_Y_TF+ATAC_Y_TF1).to_csv(outdir+'/EpiCon.txt',sep='\t')

def EpiCon(Input_dir,RNA_file,ATAC_file,label_file,provide_dir,code_dir,genome,outdir,N):
    import subprocess
    subprocess.run(["sh", code_dir+"motif_Target.sh", Input_dir, genome,provide_dir,outdir])
    find_distance(outdir,N,code_dir,provide_dir,genome)
    RNA=pd.read_csv(Input_dir+RNA_file,sep='\t',index_col=0)
    ATAC=pd.read_csv(Input_dir+ATAC_file,sep='\t',index_col=0)
    TFlist=pd.read_csv(provide_dir+'TFName715.txt',sep='\t',header=None)
    TF_motif_match=pd.read_csv(provide_dir+'Match2.txt',sep='\t',header=0)
    TF_motif_match.columns=['motif','TF']
    TF=list(set(TF_motif_match['TF'].values)&(set(RNA.index)))
    TFexp=RNA.loc[TF]
#TF_psedubulk,RE_psedubulk,indexfinal,clusternew=pseudo_bulk(TFexp,ATAC,label_file,Input_dir)
#TF_RE_cor,idxRE,idxTF=TF_RE_correlation(TF_psedubulk,RE_psedubulk)
    import scanpy as sc
    import numpy as np
    from sklearn.impute import KNNImputer
    TF_RE_cor,idxRE,idxTF=impute(TFexp,ATAC)
    ATAC=ATAC.loc[idxRE]
    TFexp=TFexp.loc[idxTF]
    ATAC1=RE_dis(outdir+'/distance.txt',ATAC,RNA)
    ECS(outdir,ATAC,TF_motif_match,ATAC1,TF_RE_cor)