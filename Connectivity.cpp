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
VI toposort(int n) {
  VI deg(n,0);
  rep(i,0,n) for (auto v:e[i]) deg[v]++;
  VI topo;
  rep(i,0,n) if (deg[i]==0) topo.pb(i);
  for (auto u:topo) for (auto v:e[u]) if (--deg[v]==0) topo.pb(v);
  if (SZ(topo)!=n) topo.clear();
  return topo;
}
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
  int low[N],dfn[N],bel[N],st[N],sol[N],top,cnt,ind,n;
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