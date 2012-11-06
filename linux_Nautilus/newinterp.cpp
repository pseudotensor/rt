struct Base_interp
{int n, mm, jsav, cor, dj;
const doub *xx, *yy;
Base_interp(Vecdoub_I &x, const doub *y, int m)
: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
dj = MIN(1,(int)pow((doub)n,0.25));
}
doub interp(doub x) {
int jlo = cor ? hunt(x) : locate(x);
return rawinterp(jlo,x);
}
int locate(const doub x); See definitions below.
int hunt(const doub x);
doub virtual rawinterp(int jlo, doub x) = 0;
};