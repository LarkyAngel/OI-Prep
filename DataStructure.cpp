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
namespace Normal {
  const int N=100010,P=333;
  struct Q {
    int l,r,id;
    bool operator < (const Q &rhs) const {
      return l/P==rhs.l/P?r<rhs.r:l<rhs.l;
    }
  } seq[N];
  int A[N],ans[N],n,m;
  namespace mo {
    int u[N],now;
    void init() {}
    void add(int x) {}
    void del(int x) {}
  }
  void run() {
    scanf("%d%d",&n,&m);
    for (int i=0;i<n;i++) scanf("%d",A+i);
    for (int i=0;i<m;i++) {
      scanf("%d%d",&seq[i].l,&seq[i].r);
      seq[i].id =i;
    }
    sort(seq,seq+m);
    mo::init();
    for (int i=0,l=0,r=-1;i<m;i++) {
      int L=seq[i].l,R=seq[i].r;
      while (r<R) mo::add(A[++r]);
      while (r>R) mo::del(A[r--]);
      while (l<L) mo::del(A[l++]);
      while (l>L) mo::add(A[--l]);
      ans[seq[i].id]=mo::now;
    }
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