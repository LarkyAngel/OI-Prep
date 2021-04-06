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
  void work(int u,int tot) {
    total=tot,mins=tot*2,rt=u;
    get_center(u,-1);
    mark[u=rt]=true;
    get_center(u,-1);
    for (auto v:e[u]) if (!mark[v]) work(v,sz[v]);
  }
}
struct wavelet_node {
  int lo,hi;
  wavelet_node *l,*r;
  VI b;
  wavelet_node(int *from,int *to,int x,int y) {
    lo=x,hi=y;
    if (lo==hi||from>=to) return;
    int md=(lo+hi)/2;
    auto f=[md](int x) {
      return x<=md;
    };
    b.reserve(to-from+1);
    b.pb(0);
    for (auto it=from;it!=to;it++) b.pb(b.back()+f(*it)); 
    auto pivot=stable_partition(from,to,f);
    l=new wavelet_node(from,pivot,lo,md);
    r=new wavelet_node(pivot,to,md+1,hi);
  }
  int kth(int l,int r,int k) {
    if (l>r) return 0;
    if (lo==hi) return lo;
    int inLeft=b[r]-b[l-1],lb=b[l-1],rb=b[r];
    if (k<=inLeft) return this->l->kth(lb+1,rb,k);
    return this->r->kth(l-lb,r-rb,k-inLeft);
  }
  int LTE(int l,int r,int k) {
    if (l>r||k<lo) return 0;
    if (hi<=k) return r-l+1;
    int lb=b[l-1],rb=b[r];
    return this->l->LTE(lb+1,rb,k)+this->r->LTE(l-lb,r-rb,k);
  }
  int count(int l,int r,int k) {
    if (l>r||k<lo||k>hi) return 0;
    if (lo==hi) return r-l+1;
    int lb=b[l-1],rb=b[r],md=(lo+hi)/2;
    if (k<=md) return this->l->count(lb+1,rb,k);
    return this->r->count(l-lb,r-rb,k);
  }
  ~wavelet_node() {
    delete l;
    delete r;
  }
};