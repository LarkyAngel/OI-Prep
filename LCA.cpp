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