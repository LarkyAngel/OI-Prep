//computational-geometry
const db eps=1e-8,pi=acos(-1.0);
inline int sgn(db x,db e=eps) { return x<-e?-1:x>e;}
struct Point {
  db x,y;
  Point(db x=0,db y=0):x(x),y(y) {}
  bool operator < (const Point &o) const { return x<o.x||(x==o.x&&y<o.y);}
  bool operator == (const Point &o) const { return x==o.x&&y==o.y;}
  Point operator + (const Point &o) const { return Point(x+o.x,y+o.y);}
  Point operator - (const Point &o) const { return Point(x-o.x,y-o.y);}
  Point operator * (db k) const { return Point(x*k,y*k);}
  Point operator / (db k) const { return Point(x/k,y/k);}
  db dot(const Point &o) const { return x*o.x+y*o.y;}
  db det(const Point &o) const { return x*o.y-y*o.x;}
  db norm() const { return hypot(x,y);}
  db ang() const { return atan2(y,x);}
};
bool on_seg(const Point &a,const Point &b,const Point &o) {
  return sgn((a-o).det(b-o))==0&&sgn((a-o).dot(b-o))<=0;
}
int in_polygon(const vector<Point> &p,const Point &o) {
  int cnt=0,n=SZ(p);
  rep(i,0,n) {
    const auto &a=p[i],&b=p[(i+1)%n];
    if (on_seg(a,b,o)) return 2;
    int k=sgn((b-a).det(o-a)),d1=sgn(a.y-o.y),d2=sgn(b.y-o.y);
    cnt+=(k>0&&d1<=0&&d2>0);
    cnt-=(k<0&&d2<=0&&d1>0);
  }
  return cnt!=0;
}
db area(const vector<Point> &p,Point &o) {
  db sum=0;
  o={0,0};
  int n=SZ(p);
  rep(i,0,n) {
    const auto &a=p[i],&b=p[(i+1)%n];
    sum+=b.det(a);
    o=o+(a+b)*b.det(a);
  }
  sum=abs(sum);
  o=o/(3.0*sum);
  return sum*0.5;
}
struct ConvexHull {
  vector<Point> ps;
  void build(vector<Point> &u) {
    sort(all(u));
    u.erase(unique(all(u),u.end()));
    if (SZ(u)<3) {
      ps=u;
      return;
    }
    for (int i=0,o=1,m=1;~i;i+=o) {
      while (SZ(ps)>m) {
        Point A=ps.back()-ps[SZ(ps)-2],B=ps.back()-u[i];
        if (sgn(A.det(B))<0) break;
        ps.pop_back();
      }
      ps.pb(u[i]);
      if (i+1==SZ(u)) m=SZ(ps),o=-1;
    }
    ps.pop_back();
  }
  db diameter(int &first,int &second) const {
    int si=0,sj=0,n=SZ(ps);
    rep(i,0,n) {
      if (ps[si].x>ps[i].x) si=i;
      if (ps[sj].x<ps[i].x) sj=i;
    }
    first=si,second=sj;
    db mx=0;
    for (int i=si,j=sj;i!=sj||j!=si;) {
      db tmp=(ps[i]-ps[j]).norm();
      if (tmp>mx) mx=tmp,first=i,second=j;
      if ((ps[(i+1)%n]-ps[i]).det(ps[(j+1)%n]-ps[j])<0) i=(i+1)%n;
      else j=(j+1)%n;
    }
    return mx;
  }
  db area(Point &center) const {
    return ::area(ps,center);
  }
  int contain(const Point &P) const {
    return in_polygon(ps,P);
  }
};
//data-structure
struct node {
  bool rev;
  ll fg0,fg1;
  int cnt0,cnt1;
  ll s;
}nd[4*N];
void upd(int p) {
  nd[p].s=nd[p+p].s+nd[p+p+1].s;
  nd[p].cnt0=nd[p+p].cnt0+nd[p+p+1].cnt0;
  nd[p].cnt1=nd[p+p].cnt1+nd[p+p+1].cnt1;
}
void setr(int p) {
  nd[p].rev^=1; swap(nd[p].fg0,nd[p].fg1); swap(nd[p].cnt0,nd[p].cnt1);
}
void setf(int p,ll f0,ll f1) {
  nd[p].fg0+=f0; nd[p].fg1+=f1;
  nd[p].s+=f0*nd[p].cnt0+f1*nd[p].cnt1;
}
void build(int p,int l,int r) {
  if (l==r) {
    if (s[l]==0) nd[p].cnt0++; else nd[p].cnt1++;
    nd[p].s=a[l];
  } else {
    int md=(l+r)>>1;
    build(p+p,l,md);
    build(p+p+1,md+1,r);
    upd(p);
  }
}
void push(int p) {
  if (nd[p].rev) {
    setr(p+p);
    setr(p+p+1);
    nd[p].rev=0;
  }
  if (nd[p].fg0||nd[p].fg1) {
    setf(p+p,nd[p].fg0,nd[p].fg1);
    setf(p+p+1,nd[p].fg0,nd[p].fg1);
    nd[p].fg0=nd[p].fg1=0;
  }
}
ll query(int p,int l,int r,int tl,int tr) {
  if (tl==l&&tr==r) return nd[p].s;
  push(p);
  int md=(l+r)>>1;
  if (tr<=md) return query(p+p,l,md,tl,tr);
  if (tl>md) return query(p+p+1,md+1,r,tl,tr);
  return query(p+p,l,md,tl,md)+query(p+p+1,md+1,r,md+1,tr);
}
void modifyr(int p,int l,int r,int tl,int tr) {
  if (tl>tr) return;
  if (tl==l&&tr==r) return setr(p);
  push(p);
  int md=(l+r)>>1;
  if (tr<=md) modifyr(p+p,l,md,tl,tr);
  else if (tl>md) modifyr(p+p+1,md+1,r,tl,tr);
  else modifyr(p+p,l,md,tl,md),modifyr(p+p+1,md+1,r,md+1,tr);
  upd(p);
}
void modifyo(int p,int l,int r,int tl,int tr,int v) {
  if (tl>tr) return;
  if (tl==l&&tr==r) return setf(p,0,v);
  push(p);
  int md=(l+r)>>1;
  if (tr<=md) modifyo(p+p,l,md,tl,tr,v);
  else if (tl>md) modifyo(p+p+1,md+1,r,tl,tr,v);
  else modifyo(p+p,l,md,tl,md,v),modifyo(p+p+1,md+1,r,md+1,tr,v);
  upd(p);
}
int tot;
struct node {
  int mx,fg,l,r;
}nd[100*N];
void upd(int p) {
  nd[p].mx=max(nd[nd[p].l].mx,nd[nd[p].r].mx);
}
void setf(int p) {
  nd[p].mx++;
  nd[p].fg++;
}
void build(int &p,int l,int r) {
  nd[tot]=nd[p],p=tot++;
  if (l==r) {
    if (l==1) nd[p].mx=1;
  } else {
    int md=(l+r)>>1;
    build(nd[p].l,l,md);
    build(nd[p].r,md+1,r);
    upd(p);
  }
}
void modify(int &p,int l,int r,int tl,int tr) {
  nd[tot]=nd[p],p=tot++;
  if (l==tl&&tr==r) return setf(p);
  int md=(l+r)>>1;
  if (tr<=md) modify(nd[p].l,l,md,tl,tr);
  else if (tl>md) modify(nd[p].r,md+1,r,tl,tr);
  else modify(nd[p].l,l,md,tl,md),modify(nd[p].r,md+1,r,md+1,tr);
  upd(p);
}
int query(int p,int l,int r,int tl,int tr,int v) {
  if (l==tl&&tr==r) return nd[p].mx+v;
  v+=nd[p].fg;
  int md=(l+r)>>1;
  if (tr<=md) return query(nd[p].l,l,md,tl,tr,v);
  if (tl>md) return query(nd[p].r,md+1,r,tl,tr,v);
  return max(query(nd[p].l,l,md,tl,md,v),query(nd[p].r,md+1,r,md+1,tr,v));
}
//GraphUtility
namespace biconnect {
int dep[N],dsu[N],p[N],bcp[N],m;
bool vis[N],v2[N];
VI e[N];
void dfs(int u) {
  vis[u]=1,dep[u]=dep[p[u]]+1;
  for (auto v:e[u]) if (!vis[v]) p[v]=u,dfs(v);
}
int get(int x) { return dsu[x]==x?x:dsu[x]=get(dsu[x]);}
void solve(int n,vector<PII> E) {
  rep(i,1,n+1) dsu[i]=i;
  for (auto p:E) e[p.fi].pb(p.se),e[p.se].pb(p.fi);
  rep(i,1,n+1) if (!vis[i]) dfs(i);
  for (auto r:E) {
    int u=r.fi,v=r.se;
    if (p[u]==v) swap(u,v);
    if (p[v]==u&&!v2[v]) {
      v2[v]=1; 
      continue;
    }
    u=get(u); v=get(v);
    while (u!=v) {
      if (dep[u]>dep[v]) swap(u,v);
      v=get(v); dsu[v]=get(p[v]);
      v=get(v);
    }
  }
  rep(i,1,n+1) if (get(i)==i) bcp[i]=++m;
  rep(i,1,n+1) bcp[i]=bcp[get(i)];
}
void build() {
  rep(i,1,n+1) if (get(i)==i&&p[i]!=0) {
    int u=bcp[i],v=bcp[get(p[i])];
    vE[u].pb(v);
    vE[v].pb(u);
  }
}
};
int cnt,dfn[N],low[N],st[N],bcp[N],top;
VI keyV,keyE;
void dfs(int u,int dep) {
  int src=0,out=1<dep;
  dfn[u]=low[u]=dep;
  for (auto p:e[u]) {
    int v=p.fi,id=p.se;
    if (!dfn[v]) {
      st[++top]=id;
      dfs(v,dep+1);
      low[u]=min(low[u],low[v]);
      if (low[v]>dfn[u]) keyE.pb(id);
      if (low[v]>=dfn[u]) {
        if (++out==2) keyV.pb(u);
        cnt++;
        while (1) {
          bcp[st[top]]=cnt;
          if (st[top--]==id) break;
        }
      }
    } else if (dfn[v]!=dfn[u]-1||src++) low[u]=min(low[u],dfn[v]);
  }
}
int cnt,dfn[N],low[N],st[N],top;
VI keyV,keyE,bcp[N];
void dfs(int u,int dep) {
  int src=0,out=1<dep;
  dfn[u]=low[u]=dep;
  st[++top]=u;
  for (auto p:e[u]) {
    int v=p.fi;
    if (!dfn[v]) {
      dfs(v,dep+1);
      low[u]=min(low[u],low[v]);
      if (low[v]>dfn[u]) keyE.pb(p.se);
      if (low[v]>=dfn[u]) {
        if (++out==2) keyV.pb(u);
        cnt++;
        while (st[top]!=u) bcp[st[top--]].pb(cnt);
        bcp[u].pb(cnt);
      }
    } else if (dfn[v]!=dfn[u]-1||src++) low[u]=min(low[u],dfn[v]);
  }
}
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
      if (dis[i][j]>res-tmp) continue;
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
VI toposort(int n) {
  VI deg(n,0);
  rep(i,0,n) for (auto v:e[i]) deg[v]++;
  VI topo;
  rep(i,0,n) if (deg[i]==0) topo.pb(i);
  for (auto u:topo) for (auto v:e[u]) if (--deg[v]==0) topo.pb(v);
  if (SZ(topo)!=n) topo.clear();
  return topo;
}
int p[N][22],dep[N];
void dfs(int u,int f) {
  dep[u]=dep[f]+1;
  p[u][0]=f;
  for (auto v:e[u]) {
    if (v==f) continue;
    dfs(v,u);
  }
}
#define LOGN 19
int lca(int u,int v) {
  if (dep[u]>dep[v]) swap(u,v);
  per(i,0,LOGN) if (dep[p[v][i]]>=dep[u]) v=p[v][i];
  if (u==v) return u;
  per(i,0,LOGN) if (p[v][i]!=p[u][i]) u=p[u][i],v=p[v][i];
  return p[u][0];
}
int dis(int u,int v) {
  int w=lca(u,v);
  return dep[u]+dep[v]-2*dep[w];
}
void run() {
  dfs(1,0);
  rep(j,1,LOGN) rep(i,1,n+1) p[i][j]=p[p[i][j-1]][j-1];
}
VI e[N];
int cnt,dfn[N],low[N],st[N],bel[N],top,ind;
bool ins[N];
void tarjan(int u) {
  dfn[u]=low[u]=++ind;
  ins[u]=1;
  st[++top]=u;
  rep(i,0,SZ(e[u])) {
    int v=e[u][i];
    if (!dfn[v]) tarjan(v),low[u]=min(low[u],low[v]);
    else if (ins[v]) low[u]=min(low[u],dfn[v]);
  } 
  if (dfn[u]==low[u]) {
    while (1) {
      bel[st[top]]=cnt;
      ins[st[top]]=0;
      if (st[top--]==u) break;
    }
    cnt++;
  }
}
namespace kosaraju {
int st[N],bel[N],cnt,top;
bool vis[N];
void dfs1(int u,int f) {
  for (auto v:te[u]) if (v!=f) dfs(v);
  st[++top]=u;
}
void dfs2(int u) {
  vis[u]=1,bel[u]=cnt;
  for (auto v:e[u]) {
    if (vis[v]) {
      if (bel[u]!=bel[v]) vE[bel[u]].pb(bel[v]);
      continue;
    }
    dfs2(v);
  }
}
void build() {
  rep(i,0,n) if (!vis[i]) dfs(i);
  fill(vis,vis+n,0);
  while (top) {
    int u=st[top--];
    if (!vis[u]) {
      dfs2(v);
      cnt++;
    }
  }
}
}
struct TwoSAT {
  VI e[N],scc[N];
  int low[N],dfn[N],bel[N],st[N],sol[N],top,cnt,ind;
  void init(int n) {
    this->n=n;
    rep(i,0,n*2) e[i].clear();
  }
  void tarjan(int u) {
    dfn[u]=low[u]=++ind;
    ins[u]=1;
    st[++top]=u;
    rep(i,0,SZ(e[u])) { 
      int v=e[u][i];
      if (!dfn[v]) tarjan(v),low[u]=min(low[u],low[v]);
      else if (ins[v]) low[u]=min(low[u],dfn[v]);
    }
    if (dfn[u]==low[u]) {
      while (1) {
        scc[cnt].pb(st[top]);
        bel[st[top]]=cnt;
        ins[st[top]]=0;
        if (st[top--]==u) break;
      }
      cnt++;
    }
  }
  bool solve() {
    rep(i,0,n*2) if (!dfn[i]) dfs(i);
    for (int i=0;i<n*2;i+=2) {
      if (bel[i]==bel[i^1]) return false;
      sol[i]=-1;
    }
    rep(i,0,cnt) {
      int val=1;
      for (auto x:scc[i]) {
        if (sol[x]==0) val=0;
        for (auto y:e[x]) if (sol[y]==0) val=0;
        if (val==0) break;
      }
      for (auto x:scc[i]) {
        sol[x]=val;
        sol[x^1]=!val;
      }
    }
    return true;
  }
  void add_clause(int x,int xv,int y,int yv) {//x=xv or y=yv
    x=x<<1|xv,y=y<<1|yv;
    e[x^1].pb(y);
    e[y^1].pb(x);
  }
  void add_var(int x,int xv) {//x=xv
    x=x<<1|xv;
    e[x^1].pb(x);
  }
};
const int inf=0x20202020;
typedef int flowt;
namespace flow {
  const int M=100000,N=1000;
  int y[M],nxt[M],gap[N],fst[N],c[N],pre[N],q[N],dis[N];
  flowt f[M];
  int S,T,tot,Tn;
  void init(int s,int t,int tn) {
    tot=1; assert(tn<N);
    rep(i,0,tn) fst[i]=0;
    S=s;T=t;Tn=tn;
  }
  void add(int u,int v,flowt c1,flowt c2=0) {
    tot++;y[tot]=v;f[tot]=c1;nxt[tot]=fst[u];fst[u]=tot;
    tot++;y[tot]=u;f[tot]=c2;nxt[tot]=fst[v];fst[v]=tot;
  }
  flowt sap() {
    int u=S,t=1;flowt flow=0;
    rep(i,0,Tn) c[i]=fst[i],dis[i]=Tn,gap[i]=0;
    q[0]=T;dis[T]=0;pre[S]=0;
    rep(i,0,t) {
      int u=q[i];
      for (int j=fst[u];j;j=nxt[j]) if (dis[y[j]]>dis[u]+1&&f[j^1]) 
        q[t++]=y[j],dis[y[j]]=dis[u]+1;
    }
    rep(i,0,Tn) gap[dis[i]]++;
    while (dis[S]<=Tn) {
      while (c[u]&&(!f[c[u]]||dis[y[c[u]]]+1!=dis[u])) c[u]=nxt[c[u]];
      if (c[u]) {
        pre[y[c[u]]]=c[u]^1;
        u=y[c[u]];
        if (u==T) {
          flowt minf=inf;
          for (int p=pre[T];p;p=pre[y[p]]) minf=min(minf,f[p^1]);
          for (int p=pre[T];p;p=pre[y[p]]) f[p^1]-=minf,f[p]+=minf;
          flow+=minf;u=S;
        }
      } else {
        if (!(--gap[dis[u]])) break;
        int mind=Tn;
        c[u]=fst[u];
        for (int j=fst[u];j;j=nxt[j]) if (f[j]&&dis[y[j]]<mind) 
          mind=dis[y[j]],c[u]=j;
        dis[u]=mind+1;
        gap[dis[u]]++;
        if (u!=S) u=y[pre[u]];
      }
    }
    return flow;
  }
};
const int inf=0x20202020;
namespace mincost {
  const int V=101000,E=10010000,_inf=0x20;
  int dis[V],q[V*30],vis[V],fst[V],pre[V],nxt[E],y[E],f[E],c[E],S,T,flow,cost,tot,tn;
  void init(int s,int t,int Tn) {
    tot=1; tn=Tn;
    rep(i,0,tn) fst[i]=0;
    S=s;T=t;
  }
  void add(int u,int v,int ff,int cc) {
    tot++;y[tot]=v;nxt[tot]=fst[u];f[tot]=ff;c[tot]=cc;fst[u]=tot;
    tot++;y[tot]=u;nxt[tot]=fst[v];f[tot]=0;c[tot]=-cc;fst[v]=tot;
  }
  bool spfa() {
    rep(i,0,tn) dis[i]=inf,vis[i]=0,pre[i]=0;
    dis[S]=0,q[0]=S,vis[S]=1;
    int t=1;
    rep(i,0,t) {
      int u=q[i];
      for (int j=fst[u];j;j=nxt[j]) {
        int v=y[j];
        if (f[j]&&dis[v]>dis[u]+c[j]) {
          dis[v]=dis[u]+c[j];
          pre[v]=j;
          if (!vis[v]) vis[v]=1,q[t++]=v;
        }
      }
      vis[u]=0;
    }
    return dis[T]!=inf;
  }
  void augment() {
    int p=T,_f=inf;
    while (pre[p]) _f=min(_f,f[pre[p]]),p=y[pre[p]^1];
    flow+=_f;cost+=_f*dis[T];
    p=T;
    while (pre[p]) f[pre[p]]-=_f,f[pre[p]^1]+=_f,p=y[pre[p]^1];
  }
  void solve() {
    flow=0,cost=0;
    while (spfa()) augment();
  }
}
//mathematics
ll pow_mod(ll a,ll n,ll m) {
  ll res=1;
  for (a%=m;n;n>>=1) {
    if (n&1) res=res*a%m;
    a=a*a%m;
  }
  return res;
}
int gcd(int a,int b) { return b?gcd(b,a%b):a;}
void exgcd(int a,int b,int &g,int &x,int &y) {
  if (!b) x=1,y=0,g=a;
  else {
    exgcd(b,a%b,g,y,x);
    y-=x*(a/b);
  }
}
int mod_inv(int a,int mod) {
  if (gcd(a,mod)!=1) return -1;
  int b=mod,s=1,t=0;
  while (b) {
    int q=a/b;
    swap(a-=q*b,b);
    swap(s-=q*t,t);
  }
  return s<0?s+mod:s;
}
int crt2(int r1,int mod1,int r2,int mod2) {
  int inv=mod_inv(mod1,mod2);
  return (r2+mod2-r1)*inv%mod2*mod1+r1;
}
int crt(int n,int *c,int *m) {
  int M=1,ans=0;
  rep(i,0,n) M*=m[i];
  rep(i,0,n) {
    int x,y,g,tm=M/m[i];
    exgcd(tm,m[i],g,x,y);
    ans=(ans+tm*x*c[i]%M)%M;
  }
  return (ans+M)%M;
}
int solve(int n,int *c,int *m) {
  auto mod=[](int x,int y) {
    x%=y;
    if (x<0) x+=y;
    return x;
  };
  int ans=c[0],LCM=m[0];
  rep(i,1,n) {
    int g,x,y;
    exgcd(LCM,m[i],g,x,y);
    if ((c[i]-ans)%g) return -1;
    int tmp=mod((c[i]-ans)/g*x,m[i]/g);
    ans=mod(ans+LCM*tmp,LCM/g*m[i]);
    LCM=LCM/g*m[i];
  }
  return mod(ans,LCM);
}
namespace BinomialCoefficient {
  VI fac,inv,ifac;
  void lucas_init(int n,int mod) {
    fac.assign(n+1,1);
    inv.assign(n+1,1);
    ifac.assign(n+1,1);
    rep(i,2,n+1) {
      fac[i]=(ll)i*fac[i-1]%mod;
      inv[i]=ll(mod-mod/i)*inv[mod%i]%mod;
      ifac[i]=(ll)ifac[i-1]*inv[i]%mod;
    }
  }
  inline int lucas_binom(int n,int m,int p) {
    return (ll)fac[n]*inv[fac[m]]*inv[fac[n-m]]%p;
  }
}
VI pl,phi,spf;
void fast_sieve(int n) {
  pl.clear();
  phi.assign(n,1);
  spf.assign(n,0);
  rep(i,2,n) {
    if (!spf[i]) {
      pl.pb(i);
      spf[i]=i;
      phi[i]=i-1;
    }
    for (int j=0;j<SZ(pl)&&i*pl[j]<n;j++) {
      int p=pl[j];
      spf[i*p]=p;
      phi[i*p]=phi[i]*(p-!!(i%p));
      if (i%p==0) break;
    }
  }
}
vector<ll> factorize(ll n) {
  vector<ll> u;
  for (int i=0,t=sqrt(n+1);pl[i]<=t;i++) if (n%pl[i]==0) {
    while (n%pl[i]==0) {
      n/=pl[i];
      u.pb(pl[i]);
    }
    t=sqrt(n+1);
  }
  if (n>1) u.pb(n);
  return u;
}
bool is_prime(int n) {
  if (n<SZ(spf)) return spf[n]==n;
  for (int i=0,sq=sqrt(n+1);pl[i]<=sq;i++) {
    if (n%pl[i]==0) return false;
  }
  return true;
}
//RangeMinimumQuery
struct SparseTable {
  int u[N][20];
  void build(int *a,int n) {
    rep(i,0,n) u[i][0]=a[i];
    for (int j=1;(1<<j)<=n;j++) for (int i=0;i+(1<<j)<=n;i++) 
      u[i][j]=min(u[i][j-1],u[i+(1<<(j-1))][j-1]);
  }
  int ask(int l,int r) {
    if (l>r) swap(l,r);
    int k=31-__builtin_clz(r-l+1);
    return min(u[l][k],u[r-(1<<k)+1][k]);
  }
};
//others
#include <bits/stdc++.h>
using namespace std;
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
#define pb push_back
#define mp make_pair
#define all(x) (x).begin(),(x).end()
#define fi first
#define se second
#define SZ(x) ((int)(x).size())
typedef vector<int> VI;
typedef long long ll;
typedef pair<int,int> PII;
typedef double db;
template<class T> inline void read(T &x) {
    x=0; int c=getchar(),f=1;
    for (;!isdigit(c);c=getchar()) if (c==45) f=-1;
    for (;isdigit(c);c=getchar()) (x*=10)+=f*(c-'0');
}
VI lis(int n,int *a) {
  VI g,res(n);
  rep(i,0,n) {
    res[i]=lower_bound(all(g),a[i])-g.begin();
    if (res[i]==SZ(g)) g.pb(a[i]);
    else g[res[i]]=a[i];
  }
  return res;
}
typedef long long matrix[N][N];
matrix ret,tmp,base;
void mul(matrix &a,matrix &b) {
  rep(i,0,101) rep(j,0,101) tmp[i][j]=0;
  rep(i,0,101) rep(j,0,101) rep(k,0,101) tmp[i][j]=(tmp[i][j]+a[i][k]*b[k][j])%mod;
  rep(i,0,101) rep(j,0,101) a[i][j]=tmp[i][j];
}
void powmod(ll b) {
  rep(i,0,101) ret[i][i]=1;
  for(;b;b>>=1) {
    if (b&1) mul(ret,base);
    mul(base,base);
  }
}
//StringUtility
VI kmp(char *s,int n) {
  VI fail(n,-1);
  for (int i=1,j=-1;i<n;i++) {
    while (j>=0&&s[j+1]!=s[i]) j=fail[j];
    fail[i]=(s[j+1]==s[i])?++j:j;
  }
  return fail;
}
VI manacher(char *s,int n) {
  VI u(n<<=1,0);
  for (int i=0,j=0,k;i<n;i+=k,j=max(j-k,0)) {
    while (i>=j&&i+j+1<n&&s[(i-j)>>1]==s[(i+j+1)>>1]) j++;
    for (u[i]=j,k=1;i>=k&&u[i]>=k&&u[i-k]!=u[i]-k;k++) u[i+k]=min(u[i-k],u[i]-k);
  }
  return u;
}
struct Hashing {
  typedef unsigned long long ull;
  ull seed;
  vector<ull> hs,pw;
  Hashing(ull _seed=9875321):seed(_seed) {}
  void build(char *s,int n) {
    hs.assign(n+1,0);
    pw.assign(n+1,1);
    rep(i,1,n+1) pw[i]=seed*pw[i-1];
    per(i,0,n) hs[i]=seed*hs[i+1]+s[i];
  }
  ull get_hash(int l,int r) {// [l, r)
    return hs[l]-hs[r]*pw[r-l];
  }
};
namespace SA {
  int cnt[N],tr[2][N],ts[N],sa[N],rk[N],ht[N];
  SparseTable rmq;
  void build(char *s,int n,int m=256) {
    int *x=tr[0],*y=tr[1];
    memset(cnt,0,sizeof(cnt[0])*m);
    rep(i,0,n) cnt[s[i]-'a']++;
    partial_sum(cnt,cnt+m,cnt);
    rep(i,0,n) rk[i]=cnt[s[i]-'a']-1;
    for (int k=1;k<=n;k<<=1) {
      rep(i,0,n) {
        x[i]=rk[i];
        y[i]=i+k<n?rk[i+k]+1:0;
      }
      fill(cnt,cnt+n+1,0);
      rep(i,0,n) cnt[y[i]]++;
      partial_sum(cnt,cnt+n+1,cnt);
      per(i,0,n) ts[--cnt[y[i]]]=i;
      fill(cnt,cnt+n+1,0);
      rep(i,0,n) cnt[x[i]]++;
      partial_sum(cnt,cnt+n+1,cnt);
      per(i,0,n) sa[--cnt[x[ts[i]]]]=ts[i];
      rep(i,rk[sa[0]]=0,n-1) rk[sa[i+1]]=rk[sa[i]]+
          (x[sa[i]]!=x[sa[i+1]]||y[sa[i]]!=y[sa[i+1]]);
    }
  }
  void get_height(char *s,int n) {
    for (int i=0,l=0,j;i<n;i++) if (rk[i]) {
      j=sa[rk[i]-1];
      while (i+l<n&&j+l<n&&s[i+l]==s[j+l]) l++;
      ht[rk[i]]=l; if (l) l--;
    }
  }
  void init_rmq(int n) {
    rmq.build(ht,n);
  }
  inline int lcp(int a,int b) {
    a=rk[a],b=rk[b];
    if (a>b) swap(a,b);
    if (a==b) return 1e9;
    return rmq.ask(a+1,b);
  }
};
