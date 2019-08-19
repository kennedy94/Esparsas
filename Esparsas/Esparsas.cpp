#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>

using namespace std;

//Variáveis globais
int N = 1000;
vector<float> w_full(N, 0);


struct VETOR{
	vector<int> ind;
	vector<float> vals;

	VETOR(vector<int> id, vector<float> va){
		ind = id;
		vals = va;
	}

	int size() { return ind.size(); }

};

float produto_interno(vector<float> v_full, VETOR v_e);

float produto_interno_ord(VETOR v, VETOR u);

void soma_empac(VETOR &x, float alpha, VETOR y);

void soma_empac_seq(VETOR &x, vector<float> alpha, vector<VETOR> y);

struct MATRIZ{
	vector<int> irn, jcn;
	vector<float> val;

	int size() { return irn.size(); }

	MATRIZ(vector<int> i, vector<int> j, vector<float> v) {
		irn = i;
		jcn = j;
		val = v;
	}
};

vector<float> prod_matriz_vetor(MATRIZ A, vector<float> x);




int main() {

	vector<int> uind = {2, 4}, vind = {0, 2};
	vector<float> uval = {1.0, 2.0}, vval = { 4.0, 3.0 };

	VETOR u(uind, uval), v(vind, vval);

	//cout << produto_interno_ord(u,v) << endl;
	
	//soma_empac(u, 0.5, v);

	vector<int>
		Ai = { 0, 1 },
		Aj = { 0, 1 };
	vector<float> Av = { 1.0, 1.0 };

	vector<float> x = { 2, 0 };

	MATRIZ A(Ai, Aj, Av);

	vector<float> y = prod_matriz_vetor(A, x);

	return 0;
}

float produto_interno (vector<float> v_full, VETOR v_e){

	float prod = 0.0;

	int n = v_full.size();

	for(int i = 0; i < n; ++i)
		prod += v_e.vals[i] * v_full[v_e.ind[i]];

	return prod;
}

float produto_interno_ord(VETOR v, VETOR u){

	float prod = 0.0;
	int i = 0, j = 0, nzv = v.size(),nzu = u.size();

	while(true){
		if(v.ind[j] < u.ind[i]){
			j++;
		}
		else {
			if(v.ind[j] > u.ind[i])
				i++;
			else{
				prod += u.vals[i]* v.vals[j];
				i++; j++;
			}
		}
		if(i == nzu || j == nzv) break;
	}

	return prod;
}

void soma_empac(VETOR & x, float alpha, VETOR y){
	int nzy = y.size(),
		nzx = x.size();

	for (int i = 0; i < nzy; i++)
		w_full[y.ind[i]] = y.vals[i];

	for (int i = 0; i < nzx; i++){
		if (w_full[x.ind[i]] != 0) {
			x.vals[i] += alpha * w_full[x.ind[i]];
			w_full[x.ind[i]] = 0;
		}
	}

	for (int i = 0; i < nzy; i++){
		if (w_full[y.ind[i]] != 0) {
			x.ind.push_back(y.ind[i]);
			x.vals.push_back(alpha * w_full[y.ind[i]]);
			w_full[y.ind[i]] = 0;
		}
	}
}

void soma_empac_seq(VETOR & x, vector<float> alpha, vector<VETOR> y){


}

vector<float> prod_matriz_vetor(MATRIZ A, vector<float> x) {
	int annz = A.size();
	int n = x.size();
	vector<float> y(n, 0);

	int i, j;
	for (int k = 0; k < annz; k++)
		y[A.irn[k]] += A.val[k] * x[A.jcn[k]];

	return y;
}
