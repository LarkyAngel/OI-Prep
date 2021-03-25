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
        Point A=ps.back()-ps[SZ(ps)-2];
        Point B=ps.back()-u[i];
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
  flt area(Point &center) const {
    return ::area(ps,center);
  }
  int contain(const Point &P) const {
    return in_polygon(ps,P);
  }
};
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
