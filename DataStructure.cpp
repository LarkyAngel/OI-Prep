namespace centroid {
  VI sz,mark;
  int total,mins,rt;
  void init(int n) {
    sz.resize(n);
    mark.assign(n,0);
  }
  void get_center(int u,int f) {
    int mx=0; sz[u]=1;
    for (auto v:e[u]) if (v!=f&&!mark[v]) {
      get_center(v,u);
      sz[u]+=sz[v];
      mx=max(mx,sz[v]);
    }
    mx=max(mx,total-sz[u]);
    if (mx<mins) mins=mx,rt=u;
  }
  void run(int u,int tot) {
    total=tot,mins=tot*2,rt=u;
    get_center(u,-1);// find centroid
    mark[u=rt]=true;
    get_center(u,-1);// recalculate subtree size
    for (auto v:e[u]) if (!mark[v]) run(v,sz[v]);
  }
}
namespace Normal {
  const int N = 100000 + 10, P = 333;
  struct Q {
    int l, r, id;
    bool operator < (const Q &rhs) const {
      return l / P == rhs.l / P ? r < rhs.r : l < rhs.l;
    }
  } seq[N];
  int A[N], ans[N], n, m;
  namespace mo {
    int u[N], now;
    void init() {}
    void add(int x) {}
    void del(int x) {}
  }
  void run() {
    scanf("%d%d", &n, &m);
    for (int i = 0; i < n; ++i) scanf("%d", A + i);
    for (int i = 0; i < m; ++i) {
      scanf("%d%d", &seq[i].l, &seq[i].r);
      seq[i].id = i;
    }
    std::sort(seq, seq + m);
    mo::init();
    for (int i = 0, l = 0, r = -1; i < m; ++i) {
      int L = seq[i].l, R = seq[i].r;
      while (r < R) mo::add(A[++r]);
      while (r > R) mo::del(A[r--]);
      while (l < L) mo::del(A[l++]);
      while (l > L) mo::add(A[--l]);
      ans[seq[i].id] = mo::now;
    }
  }
}

namespace Tree {
  using int64 = long long;
  const int N = 100000 + 10, K = 17;
  std::vector<int> G[N];
  int dep[N], f[N][K + 1];
  int st[N], ed[N], loc[N << 1];
  int val[N], n, m, P, dfn;
  void dfs(int x, int par = -1) {
    loc[st[x] = dfn++] = x;
    for (int i = 1; i <= K; ++i) f[x][i] = f[f[x][i - 1]][i - 1];
    for (auto &y: G[x]) if (y != par) {
      dep[y] = dep[f[y][0] = x] + 1;
      dfs(y, x);
    }
    loc[ed[x] = dfn++] = x;
  }
  int lca(int x, int y) {
    if (x == y) return x;
    if (dep[x] < dep[y]) std::swap(x, y);
    for (int i = K; ~i; --i) {
      if (dep[f[x][i]] >= dep[y]) x = f[x][i];
    }
    if (x == y) return x;
    for (int i = K; ~i; --i) {
      if (f[x][i] != f[y][i]) x = f[x][i], y = f[y][i];
    }
    return f[x][0];
  }
  namespace mo {
    int cnt[N], vis[N], sum;
    void deal(int x) {//这部分维护序列
      int c = val[x]; vis[x] ^= 1;
      sum -= !!cnt[c];
      if (!vis[x]) cnt[c]--;
      else cnt[c]++;
      sum += !!cnt[c];
    }
  }
  struct Node {
    int l, r, z, id;
    bool operator < (const Node &rhs) {
      return l / P == rhs.l / P ? r < rhs.r : l / P < rhs.l / P;
    }
  } Q[N];
  int ans[N];
  void run() {
    scanf("%d", &n);
    for (int i = 1; i <= n; ++i) scanf("%d", val + i);
    for (int i = 1; i < n; ++i) {
      int u, v; scanf("%d%d", &u, &v);
      G[u].push_back(v);
      G[v].push_back(u);
    }
    dfs(dep[1] = 1);
    P = sqrt(n * 2);
    scanf("%d", &m);
    for (int i = 0; i < m; ++i) {
      int x, y; scanf("%d%d", &x, &y);
      if (st[x] > st[y]) std::swap(x, y);
      int z = lca(x, y); Q[i].id = i;
      if (z == x) Q[i].l = st[x], Q[i].r = st[y];
      else Q[i].l = ed[x], Q[i].r = st[y], Q[i].z = z;
    }
    std::sort(Q, Q + m);
    for (int i = 0, l = 0, r = -1; i < m; ++i) {
      while (r < Q[i].r) mo::deal(loc[++r]);
      while (r > Q[i].r) mo::deal(loc[r--]);
      while (l < Q[i].l) mo::deal(loc[l++]);
      while (l > Q[i].l) mo::deal(loc[--l]);
      if (Q[i].z) mo::deal(Q[i].z);
      ans[Q[i].id] = mo::sum;
      if (Q[i].z) mo::deal(Q[i].z);
    }
    for (int i = 0; i < m; ++ i) printf("%d\n", ans[i]);
  }
}

namespace Modified {
  using int64 = long long;
  using pii = std::pair<int, int>;

  constexpr int N = 1e5 + 10;
  constexpr int M = 1e5 + 10;

  struct query_t {
    int l, r, t;
  } ask[N];

  std::vector<int> query[60][60];
  int a[N], value[N];
  int angry[N];
  int n, q, B;

  namespace mo {
    int cnt[M], mark[N];
    int64 ret;
    void init() {
      memset(cnt, 0, sizeof(cnt));
      memset(mark, 0, sizeof(mark));
      ret = 0;
    }
    void add(int x) {
      int v = value[x];
      if (mark[x]) ret -= (int64)v * angry[cnt[v]--];
      else ret += (int64)v * angry[++cnt[v]];
      mark[x] ^= 1;
    }
  }

  void run() {
    scanf("%d", &n);
    for (int i = 0; i < n; ++i) {
      scanf("%d", &a[i]);
    }
    for (int i = 1; i <= n; ++i) {
      scanf("%d", &angry[i]);
    }
    for (B = 1; B * B * B < n; ++B);
    B = B * B;
    std::vector<pii> modify;
    scanf("%d", &q);
    for (int i = 0; i < q; ++i) {
      int op, x, y;
      scanf("%d%d%d", &op, &x, &y);
      if (op == 2) {
        ask[i].t = -2;
        modify.emplace_back(x - 1, y);
      } else {
        --x, --y;
        ask[i].t = modify.size() - 1;
        ask[i].l = x;
        ask[i].r = y;
        query[x / B][y / B].emplace_back(i);
      }
    }
    std::vector<int64> ret(q);
    int m = (n - 1) / B + 1;
    for (int i = 0; i < m; ++i) {
      for (int j = i; j < m; ++j) {
        if (query[i][j].empty()) continue;
        memcpy(value, a, sizeof(*a) * n);
        mo::init();
        int t = 0, l = 0, r = -1;
        for (auto &&v: query[i][j]) {
          auto qr = ask[v];
          while (r < qr.r) mo::add(++r);
          while (r > qr.r) mo::add(r--);
          while (l < qr.l) mo::add(l++);
          while (l > qr.l) mo::add(--l);
          for (int x, y; t <= qr.t; ++t) {
            std::tie(x, y) = modify[t];
            bool flag = qr.l <= x && x <= qr.r;
            if (flag) mo::add(x);
            value[x] = y;
            if (flag) mo::add(x);
          }
          ret[v] = mo::ret;
        }
      }
    }
    for (int i = 0; i < q; ++i) {
      if (ask[i].t == -2) continue;
      printf("%lld\n", ret[i]);
    }
  }
}
struct BitVector {
  using uint=uint32_t;
  vector<uint> blocks,ranktable;
  int len,blocks_len,count;
  BitVector(int n=0) {init(n);}
  void set(int i) { blocks[i>>5]|=1<<(i&31);}
  bool get(int i) const { return blocks[i>>5]>>(i&31)&1;}
  void init(int n) {
    len=n,blocks_len=(len+31)/32+1;
    blocks.assign(blocks_len,0);
  }
  void build() {
    count=0;
    if (len==0) return;
    ranktable.assign(blocks_len+1,0);
    rep(i,0,blocks_len) {
      ranktable[i]=count;
      count+=popcount(blocks[i]);
    }
    ranktable[blocks_len]=count;
  }
  int rank(int p) const {// number of 1s in [0, p)
    int idx=p>>5;
    return ranktable[idx]+popcount(blocks[idx]&(1u<<(p&31))-1);
  }
  int rank(bool b,int p) const { return b?rank(p):p-rank(p);}
  int rank(bool b,int l,int r) const { return rank(b,r)-rank(b,l);}
 private:
  static inline int popcount(uint x) {
    x=x-((x>>1)&0x55555555);
    x=(x&0x33333333)+((x>>2)&0x33333333);
    return ((x+(x>>4)&0xF0F0F0F)*0x1010101)>>24;
  }
};
struct WaveletMatrix {
  using value_t=unsigned int;
  int height,len;
  value_t maxval;
  vector<BitVector> B;
  VI mids;
  WaveletMatrix() { init({});}
  void init(const vector<value_t> &data) {
    len=data.size();
    maxval=len?*max_element(all(data)):0;
    for (height=1;maxval>>1>>(height-1);height++);
    B.assign(height,BitVector(len));
    mids.resize(height);
    vector<value_t> now(data),next(len);
    rep(i,0,height) {
      rep(j,0,len) mids[i]+=(now[j]>>(height-1-i)&1)==0;
      BitVector &bv=B[i];
      int zero=0,one=mids[i];
      rep(j,0,len) {
        bool b=now[j]>>(height-1-i)&1;
        if (b) next[one++]=now[j],bv.set(j);
        else next[zero++]=now[j];
      }
      bv.build();
      next.swap(now);
    }
  }
  value_t get(int p) const {
    value_t ret=0;
    rep(i,0,height) {
      const BitVector &bv=B[i];
      bool dir=bv.get(p);
      ret=ret<<1|dir;
      p=bv.rank(dir,p);
      if (dir) p+=mids[i];
    }
    return ret;
  }
  // k-th element in position [left, right)
  value_t quantile(int left,int right,int k) const {
    if(k<0||right-left<=k) { return -1;}
    value_t val=0;
    rep(i,0,height) {
      const BitVector &bv=B[i];
      int count=bv.rank(true,left,right);
      bool dir=k<count;
      val=val<<1|(dir?1:0);
      if(!dir) k-=count;
      left=bv.rank(dir,left),right=bv.rank(dir,right);
      if(dir) left+=mids[i],right+=mids[i];
    }
    return val;
  }
  // number of element less/greater than val in position [left, right), return the rank?
  int rank_all(value_t val,int left,int right,int &res_lt,int &res_gt) const {
    if(val>maxval) {
      res_lt=right-left;
      res_gt=0;
      return 0;
    }
    res_lt=res_gt=0;
    rep(i,0,height) {
      const BitVector &bv=B[i];
      bool dir=val>>(height-i-1)&1;
      int leftcount=bv.rank(dir,left),rightcount=bv.rank(dir,right);
      (dir?res_lt:res_gt)+=(right-left)-(rightcount-leftcount);
      left=leftcount,right=rightcount;
      if(dir) left+=mids[i],right+=mids[i];
    }
    return right-left;
  }
};