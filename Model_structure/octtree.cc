#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <assert.h>
#include <set>
using namespace std;
typedef vector< vector <bool> > bits;

//------------------------------------------------------------------------------------
struct index{
        bits indices;
        index(){indices.resize(3);}
        void add(bool i,bool j,bool k){
            indices[0].push_back(i);
            indices[1].push_back(j);
            indices[2].push_back(k);
        }
        index up(int w){
            index result;
            bits B=indices;
            unsigned p=B[w].size();
            while (p>0){
                if (B[w][p-1]==1){B[w][p-1]=0;p--;} else {B[w][p-1]=1;p=0;}
            }
            result.indices=B;
            return result;
        }
        index down(int w){
            index result;
            bits B=indices;
            unsigned p=B[w].size();
            while (p>0){
                if (B[w][p-1]==0){B[w][p-1]=1;p--;} else {B[w][p-1]=0;p=0;}
            }
            result.indices=B;
            return result;
        }
        index neibr(int w,int d){
            if (d>0)return up(w);
            if (d<0)return down(w);
            index r;
            r.indices=indices;
            return r;
        }
};
//------------------------------------------------------------------------------------
struct Cohort{
    double x,y,z;
    index idx;
    Cohort(double a,double b,double c):x(a),y(b),z(c){;}
};
///////////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------------
class octTree{
    Cohort* c;
    bool leaf;
    index idx;
    vector< vector <vector <octTree*> > > child;
    octTree* parent;
    double splitu,splitv,splitw;
    double sx,sy,sz;
    void makeChildren();
    static int maxDepth;
public:
    octTree(octTree*,int,int,int,double,double,double,double,double,double);
    ~octTree();
    void add(Cohort*);
    void searchTree(int , index& , set<Cohort*>& );
    set<Cohort*> findNeighbours(Cohort& );
};
//------------------------------------------------------------------------------------
int octTree:: maxDepth=32;
octTree::octTree(octTree* _parent, int i,int j, int k, double x, double y, double z,double u, double v, double w)
: c(0), leaf(true), parent(_parent) {
    splitu = u;
    splitv = v;
    splitw = w;
    sx=x/2;
    sy=y/2;
    sz=z/2;
    if (_parent!=0){idx=_parent->idx;idx.add(i,j,k);}else{idx.add(0,0,0);}

}
//------------------------------------------------------------------------------------
octTree::~octTree() {
    for (int i = 0;i < 2;i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                if (child.size() != 0)delete child[i][j][k];
            }
        }
    }
}
//------------------------------------------------------------------------------------
void octTree::add(Cohort* a) {
    if (c == 0 && leaf) {
        c = a;
        a->idx=idx;
    } else {
        if (leaf) {
            leaf = false;
            makeChildren();
            if (idx.indices[0].size()<maxDepth)child[c->x > splitu][c->y > splitv][c->z > splitw]->add(c);
            c = 0;
            if (idx.indices[0].size()<maxDepth)child[a->x > splitu][a->y > splitv][a->z > splitw]->add(a);
        } else {
            if (idx.indices[0].size()<maxDepth)child[a->x > splitu][a->y > splitv][a->z > splitw]->add(a);
        }
    }
}
//------------------------------------------------------------------------------------
void octTree::makeChildren() {
    child.resize(2);
    for (int i = 0; i < 2; i++) {
        child[i].resize(2);
        for (int j = 0; j < 2; j++) {
            child[i][j].resize(2);
            for (int k = 0; k < 2; k++) {
                child[i][j][k] = new octTree(this, i,j,k,sx,sy,sz, splitu+(2*i-1)*sx/2, splitv+(2*j-1)*sy/2, splitw+(2*k-1)*sz/2);
            }
        }
    }
}
//------------------------------------------------------------------------------------
void octTree::searchTree(int level, index& sdx, set<Cohort*>& N) {
    if (sdx.indices[0][0] == 0 && sdx.indices[1][0] == 0 && sdx.indices[2][0] == 0) {

        if (!leaf){
            child[sdx.indices[0][level]]
                 [sdx.indices[1][level]]
                 [sdx.indices[2][level]]->searchTree(level + 1, sdx, N);
        }else{
            if (c != 0)N.insert(c);}
    }
}
//------------------------------------------------------------------------------------
set <Cohort* > octTree::findNeighbours(Cohort& c) {
    set<Cohort*> Nbrs;
    index ix;
    for (int k = -1; k < 2; k++) {
        for (int j = -1; j < 2; j++) {
            for (int i = -1; i < 2; i++) {
                if (!(k == 0 && j == 0 && i == 0)) {
                    ix = c.idx.neibr(2, k).neibr(1, j).neibr(0, i);
                    searchTree(1, ix, Nbrs);
                }
            }
        }
    }
    return Nbrs;
}
//------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////
void test1(){
        octTree o(0,0,0,0,1.,1.,1.,0.5,0.5,0.5);
        Cohort c(0.91,0.91,0.91);
        Cohort d(0.92,0.92,0.91);
        o.add(&c);
        o.add(&d);
        cout<<c.idx.indices[0].size()<<endl;
        cout<<c.idx.indices[0][7]<<c.idx.indices[1][7]<<c.idx.indices[2][7]<<endl;
        cout<<d.idx.indices[0][7]<<d.idx.indices[1][7]<<d.idx.indices[2][7]<<endl;
        set<Cohort*> Nbrs;
        index ix=c.idx.up(0).up(1).up(2).down(2);
        o.searchTree(1,ix,Nbrs);
        for (auto n:Nbrs)cout<<n->x<<endl;
        ix=c.idx.neibr(0,1).neibr(1,1);
        o.searchTree(1,ix,Nbrs);
        for (auto n:Nbrs)cout<<n->x<<endl;
}
void test2(){
    vector<Cohort*> Cs;
    int sz=1000;
    for (int i=0;i<sz;i++){
        Cohort *c=new Cohort(float(rand())/RAND_MAX,float(rand())/RAND_MAX,float(rand())/RAND_MAX);
        Cs.push_back(c);
    }
    double mxx=-1,mnx=2;
    double mxy=-1,mny=2;
    double mxz=-1,mnz=2;
    for (int i=0;i<sz;i++){mxx=max(Cs[i]->x,mxx);mnx=min(Cs[i]->x,mnx);}
    for (int i=0;i<sz;i++){mxy=max(Cs[i]->y,mxy);mny=min(Cs[i]->y,mny);}
    for (int i=0;i<sz;i++){mxz=max(Cs[i]->z,mxz);mnz=min(Cs[i]->z,mnz);}
    
    for (int kk=0;kk<10;kk++){


    octTree o(0,0,0,0,mxx-mnx,mxy-mny,mxz-mnz,(mxx+mnx)/2,(mxy+mny)/2,(mxz+mnz)/2);
    for (int i=0;i<sz;i++){o.add(Cs[i]);}
    }
}
void test3(){
    index idx;
    idx.add(0,0,0);
    idx.add(1,0,0);
    idx.add(1,0,1);
    idx.add(1,1,1);
    for (int i=0;i<idx.indices[0].size();i++)cout<<idx.indices[0][i];
    cout<<endl;
    index ix;
    ix=idx.up(0);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[0][i];
    cout<<endl;
    ix=ix.up(0).up(0);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[0][i];
    cout<<endl;
    ix=ix.down(0);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[0][i];
    cout<<endl;
    ix=ix.down(0);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[0][i];
    cout<<endl;
    ix =ix.down(1);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[1][i];
    cout<<endl;
    ix=ix.up(2);
    for (int i=0;i<ix.indices[0].size();i++)cout<<ix.indices[2][i];
    cout<<endl;
}
double dist(Cohort& c,Cohort&d){
    return pow((c.x-d.x),2)+pow((c.y-d.y),2)+pow((c.z-d.z),2);
}
void test4() {
    octTree o(0, 0, 0, 0, 1., 1., 1., 0.5, 0.5, 0.5);
    Cohort c(0.91, 0.91, 0.91);
    Cohort d(0.92, 0.92, 0.91);
    Cohort e(0.93, 0.92, 0.91);
    Cohort f(0.912, 0.92, 0.91);
    Cohort g(0.97, 0.92, 0.91);
    o.add(&c);
    o.add(&d);
    o.add(&e);
    o.add(&f);
    o.add(&g);
    set<Cohort*>   Nbrs=o.findNeighbours(c);
    for (auto n : Nbrs)cout << "c "<<dist(c,*n) << endl;
    Nbrs.clear();  Nbrs=o.findNeighbours(d);
    for (auto n : Nbrs)cout << "d "<<dist(d,*n) << endl;
    Nbrs.clear();   Nbrs=o.findNeighbours(e);
    for (auto n : Nbrs)cout << "e "<<dist(e,*n) << endl;
    Nbrs.clear();   Nbrs=o.findNeighbours(f);
    for (auto n : Nbrs)cout << "f "<<dist(f,*n) << endl;
    Nbrs.clear();   Nbrs=o.findNeighbours(g);
    for (auto n : Nbrs)cout << "g "<<dist(g,*n) << endl;
    
    cout<<"foo"<<endl;
    cout << dist(c,d) << endl;
    cout << dist(c,e) << endl;
    cout << dist(c,f) << endl;
    cout << dist(c,g) << endl;
    cout << dist(d,e) << endl;
    cout << dist(d,f) << endl;
    cout << dist(d,g) << endl;
    cout << dist(e,f) << endl;
    cout << dist(e,g) << endl;
    cout << dist(f,g) << endl;
}
//------------------------------------------------------------------------------------
int main(){

    //test1();
    //test2();
    //test3();
    test4();
    return 0;
}
