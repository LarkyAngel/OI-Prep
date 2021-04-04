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