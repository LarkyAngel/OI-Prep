ll pow_mod(ll a,ll n,ll m) {
  ll res=1;
  for (a%=m;n;n>>=1) {
    if (n&1) res=res*a%m;
    a=a*a%m;
  }
  return res;
}
ll gcd(ll a,ll b) { return b?gcd(b,a%b):a;}
void exgcd(ll a,ll b,ll &g,ll &x,ll &y) {
  if (!b) x=1,y=0,g=a;
  else {
    exgcd(b,a%b,g,y,x);
    y-=x*(a/b);
  }
}
ll mod_inv(ll a,ll mod) {
  if (gcd(a,mod)!=1) return -1;
  ll b=mod,s=1,t=0;
  while (b) {
    ll q=a/b;
    swap(a-=q*b,b);
    swap(s-=q*t,t);
  }
  return s<0?s+mod:s;
}
bool linear_equation(ll a,ll b,ll c,ll &x,ll &y) {
  ll g;
  exgcd(a,b,g,x,y);
  if (c%g) return 0;
  a/=g,b/=g,c/=g; 
  x=(x%b*(c%b)%b+b)%b;
  y=(c-a*x)/b;
  return 1;
}
vector<ll> congruence_equation(ll a,ll b,ll m) {
  vector<ll> ret;
  ll g=gcd(a,m),x;
  if (b%g!=0) return ret;
  a/=g,b/=g;
  x=mod_inv(a,m/g);
  rep(k,0,g) ret.pb((x*b+m/g*k)%m);
  return ret;
}
typedef unsigned long long ull;
ull crt2(ull r1,ull mod1,ull r2,ull mod2) {
  ull inv=mod_inv(mod1,mod2);
  return (r2+mod2-r1)*inv%mod2*mod1+r1;
}
ll crt(int n,ll *c,ll *m) {
  ll M=1,ans=0;
  rep(i,0,n) M*=m[i];
  rep(i,0,n) {
    ll g,x,y,tm=M/m[i];
    exgcd(tm,m[i],g,x,y);
    ans=(ans+tm*x*c[i]%M)%M;
  }
  return (ans+M)%M;
}
ll solve(int n,ll *c,ll *m) {
  auto mod=[](ll x,ll y) {
    x%=y;
    if (x<0) x+=y;
    return x;
  };
  ll ans=c[0],LCM=m[0];
  rep(i,1,n) {
    ll g,x,y;
    exgcd(LCM,m[i],g,x,y);
    if ((c[i]-ans)%g) return -1;
    ll tmp=mod((c[i]-ans)/g*x,m[i]/g);
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
      ifac[i]=(ll)inv[i]*ifac[i-1]%mod;
    }
  }
  inline int lucas_binom(int n,int m,int p) {
    return (ll)fac[n]*inv[fac[m]]*inv[fac[n-m]]%p;
  }
}
VI pl,phi,spf;
void fast_sieve(int n) {
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