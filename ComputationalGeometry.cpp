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
  return (a-o).det(b-o)==0&&(a-o).dot(b-o)<=0;
}
int in_polygon(const vector<Point> &p,const Point &o) {
  int cnt=0,n=SZ(p);
  rep(i,0,n) {
    const auto &a=p[i],&b=p[(i+1)%n];
    if (on_seg(a,b,o)) return 2;
    int k=(b-a).det(o-a),d1=a.y-o.y,d2=b.y-o.y;
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
        if (A.det(B)<0) break;
        ps.pop_back();
      }
      ps.pb(u[i]);
      if (i==SZ(u)-1) m=SZ(ps),o=-1;
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