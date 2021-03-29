const int inf=0x20202020;
namespace mincost {
  const int V=101000,E=10010000,_inf=0x20;
  int dis[V],q[V*30],inq[V],fst[V],pre[V],nxt[E],y[E],f[E],c[E],S,T,flow,cost,tot,tn;
  void init(int s,int t,int Tn) {
    tot=1,tn=Tn,S=s,T=t;
  }
  void add(int u,int v,int ff,int cc) {
    tot++,y[tot]=v,nxt[tot]=fst[u],f[tot]=ff,c[tot]=cc,fst[u]=tot,
    tot++,y[tot]=u,nxt[tot]=fst[v],f[tot]=0,c[tot]=-cc,fst[v]=tot;
  }
  bool spfa() {
    rep(i,0,tn) dis[i]=inf;
    dis[S]=0,q[0]=S,inq[S]=1;
    int t=1;
    rep(i,0,t) {
      int u=q[i];
      for (int j=fst[u];j;j=nxt[j]) {
        int v=y[j];
        if (f[j]&&dis[v]>dis[u]+c[j]) {
          dis[v]=dis[u]+c[j];
          pre[v]=j;
          if (!inq[v]) inq[v]=1,q[t++]=v;
        }
      }
      inq[u]=0;
    }
    return dis[T]!=inf;
  }
  void augment() {
    int p=T,_f=inf;
    while (pre[p]) _f=min(_f,f[pre[p]]),p=y[pre[p]^1];
    flow+=_f,cost+=_f*dis[T],p=T;
    while (pre[p]) f[pre[p]]-=_f,f[pre[p]^1]+=_f,p=y[pre[p]^1];
  }
  void solve() {
    while (spfa()) augment();
  }
}