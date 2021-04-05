VI kmp(int n,char *s) {
  VI fail(n,-1);
  for (int i=1,j=-1;i<n;i++) {
    while (j>=0&&s[j+1]!=s[i]) j=fail[j];
    fail[i]=(s[j+1]==s[i])?++j:j;
  }
  return fail;
}
VI manacher(int n,char *s) {
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
  void build(int n,char *s) {
    hs.assign(n+1,0);
    pw.assign(n+1,1);
    rep(i,1,n+1) pw[i]=seed*pw[i-1];
    per(i,0,n) hs[i]=seed*hs[i+1]+s[i];
  }
  ull get_hash(int l,int r) {
    return hs[l]-hs[r]*pw[r-l];
  }
};
namespace SA {
  int cnt[N],tr[2][N],ts[N],sa[N],rk[N],ht[N];
  SparseTable rmq;
  void build(int n,char *s,int m=256) {
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
  void get_height(int n,char *s) {
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