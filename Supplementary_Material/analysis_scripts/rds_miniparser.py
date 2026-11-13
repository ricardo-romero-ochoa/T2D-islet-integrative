import gzip, struct, io, numpy as np
HAS_ATTR=1<<9; HAS_TAG=1<<10
NIL=254; REF=255
class RDSReader:
    def __init__(self,path):
        raw=gzip.open(path,'rb').read() if open(path,'rb').read(2)==b'\x1f\x8b' else open(path,'rb').read()
        self.f=io.BytesIO(raw); self.refs=[]
        magic=self.f.read(2)
        if magic not in (b'X\n',b'B\n'): raise ValueError(magic)
        self.version=self.i(); self.writer=self.i(); self.minreader=self.i()
        # version 3 has encoding string: int length then bytes
        if self.version>=3:
            n=self.i(); self.encoding=self.f.read(n).decode('ascii','replace')
    def i(self): return struct.unpack('>i',self.f.read(4))[0]
    def u(self): return struct.unpack('>I',self.f.read(4))[0]
    def d(self,n): return np.frombuffer(self.f.read(8*n),dtype='>f8').astype(float)
    def read(self): return self.sexp()
    def sexp(self):
        flags=self.u(); typ=flags&0xff
        if typ==NIL: return None
        if typ==REF:
            idx=flags>>8
            if idx==0: idx=self.i()
            return self.refs[idx-1]
        tag=None
        if flags & HAS_TAG:
            tag=self.sexp()
        if typ==1: # SYMSXP
            val=self.sexp(); self.refs.append(val); return val
        if typ==9: # CHARSXP
            n=self.i();
            if n<0: val=None
            else: val=self.f.read(n).decode('utf-8','replace')
            self.refs.append(val); return val
        if typ in (10,13):
            n=self.i(); val=np.frombuffer(self.f.read(4*n),dtype='>i4').astype(int)
        elif typ==14:
            n=self.i(); val=self.d(n)
        elif typ==16:
            n=self.i(); val=[self.sexp() for _ in range(n)]
        elif typ in (19,20):
            n=self.i(); val=[self.sexp() for _ in range(n)]
        elif typ in (2,6,17):
            car=self.sexp(); cdr=self.sexp(); val=('pair',tag,car,cdr)
            return val
        else:
            raise NotImplementedError(f'type {typ} flags {hex(flags)} pos {self.f.tell()}')
        attrs=None
        if flags & HAS_ATTR:
            attrs=self.sexp()
        return {'type':typ,'value':val,'attrs':attrs,'tag':tag}

def pairlist_to_dict(x):
    out={}
    while isinstance(x,tuple) and x and x[0]=='pair':
        _,tag,car,cdr=x
        out[tag]=car
        x=cdr
    return out

def read_matrix(path):
    r=RDSReader(path); obj=r.read(); attrs=pairlist_to_dict(obj['attrs'])
    dim=np.array(attrs['dim']['value'],int)
    dn=attrs.get('dimnames'); rownames=colnames=None
    if dn:
        vals=dn['value']; rownames=vals[0]['value'] if vals[0] else None; colnames=vals[1]['value'] if vals[1] else None
    mat=np.array(obj['value']).reshape(tuple(dim),order='F')
    return mat, list(rownames) if rownames else None, list(colnames) if colnames else None
if __name__=='__main__':
 import sys
 m,r,c=read_matrix(sys.argv[1]); print(m.shape,r[:3],r[-3:],c[:3],c[-3:]); print(m[:2,:3])
