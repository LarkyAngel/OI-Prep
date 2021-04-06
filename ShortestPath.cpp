ll dis[N];
set<pair<ll,int>> hs;
void dijkstra(int S) {
  rep(i,1,n+1) dis[i]=1ll<<60;
  dis[S]=0;
  rep(i,1,n+1) hs.insert(mp(dis[i],i));
  rep(i,0,n) {
    int u=hs.begin()->se; hs.erase(hs.begin());
    rep(j,0,SZ(e[u])) {
      int v=e[u][j].fi;
      if (dis[v]>dis[u]+e[u][j].se) {
        hs.erase(mp(dis[v],v));
        dis[v]=dis[u]+e[u][j].se;
        hs.insert(mp(dis[v],v));
      }
    }
  }
}
vector<PII> e[N];
ll dis[N];
int q[N*N],cnt[N];
bool inq[N];
bool spfa(int S) {
  rep(i,1,n+1) dis[i]=1ll<<60;
  dis[S]=0,q[0]=S,inq[S]=1;
  int t=1;
  rep(i,0,t) {
    int u=q[i];
    for (auto p:e[u]) {
      int v=p.fi;
      if (dis[v]>dis[u]+p.se) {
        dis[v]=dis[u]+p.se;
        if (!inq[v]) {
          inq[v]=1,q[t++]=v,cnt[v]++;
          if (cnt[v]>n) return 1;
        }
      }
    }
    inq[u]=0;
  }
  return 0;
}
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
struct edge {
  int u,v,w;
};
vector<edge> E;
VI cyc;
ll dis[N];
int p[N];
bool bellman(int S) {
  rep(i,1,n+1) dis[i]=1ll<<60;
  dis[S]=0;
  int x;
  rep(i,0,n) {
    x=0;
    for (auto p:E) {
      int u=p.u,v=p.v;
      if (dis[v]>dis[u]+p.w) {
        dis[v]=dis[u]+p.w;
        p[v]=u;
        x=v;
      }
    }
  } 
  if (x) {
    rep(i,0,n) x=p[x];
    int y=x;
    while (x!=y||SZ(cyc)<1) {
      cyc.pb(y);
      y=p[y];
    }
    return 1;
  }
  return 0;
}
VI path;
ll dis[N][N];
int src[N][N];
void get_path(int i,int j) {
  int k=src[i][j];
  if (k!=-1) {
    get_path(i,k);
    get_path(k,j);
  } else path.pb(j);
}
ll run() {
  ll res=1e9+7;
  rep(i,0,n) rep(j,0,n) src[i][j]=-1,dis[i][j]=e[i][j];
  rep(k,0,n) {
    rep(i,0,k) rep(j,i+1,k) {
      int tmp=e[k][i]+e[j][k];
      if (dis[i][j]+tmp>res) continue;
      path.clear();
      get_path(i,j);
      path.pb(k);
      path.pb(i);
      res=dis[i][j]+tmp;
    }
    rep(i,0,n) rep(j,0,n) {
      ll tmp=dis[i][k]+dis[k][j];
      if (tmp<dis[i][j]) {
        dis[i][j]=tmp;
        src[i][j]=k;
      }
    }
  }
  return res;
}