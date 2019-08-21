#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>

using namespace std;

//Variáveis globais
int N = 1000;
vector<float> w_full(N, 0);
vector<int> p_full(N, -1);


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

void soma_empac_alternativa(VETOR &x, float alpha, VETOR y);

void soma_empac_seq(VETOR &x, vector<float> alpha, vector<VETOR> y);

void ord_vetor(VETOR &x);

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

struct MATRIZ_LINKED_COL {
	vector<int> col_start, row_index, link;
	vector<float> value;
};


/*Matriz como coceção de vetores esparsos coluna*/
struct MATRIZ_VET_LIN{
	vector<int> len_col, col_start, row_index;
	vector<float> value;
};

struct MATRIZ_VET_COL{
	vector<int> lenrow, irowst, jcn;
	vector<float> value;
};

vector<float> prod_matriz_vetor(MATRIZ A, vector<float> x);

MATRIZ_LINKED_COL converter_matriz(MATRIZ A);



int main() {

	vector<int> uind = {2, 4}, vind = {0, 2};
	vector<float> uval = {1.0, 2.0}, vval = { 4.0, 3.0 };

	VETOR u(uind, uval), v(vind, vval), z(uind,vval);

	vector<float> alpha_vet = {0.25, 0.5};

	//cout << produto_interno_ord(u,v) << endl;
	
	//soma_empac(u, 0.5, v);
	//soma_empac_alternativa(u, 0.5, v);

	soma_empac_seq(u, alpha_vet, vector<VETOR> {v, z});

	/*vector<int>
		Ai = { 0, 1 },
		Aj = { 0, 1 };
	vector<float> Av = { 1.0, 1.0 };

	vector<float> x = { 2, 0 };

	MATRIZ A(Ai, Aj, Av);

	vector<float> y = prod_matriz_vetor(A, x);*/

	//ord_vetor(u);

	//TESTES PARA MATRIZ------------------------------------------------------------

	vector<int> irow = { 4,5,1,1,5,2,4,3,3,2,1 },
		jcol = { 1,2,2,1,5,3,4,5,2,4,5 };

	for (int i = 0; i < irow.size(); i++){
		irow[i]--;
		jcol[i]--;
	}

	vector<float> val = { -1,3,2,1,6,-3,-4,-5,-2,4,5 };

	MATRIZ A(irow, jcol, val);

	converter_matriz(A);

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

void soma_empac_alternativa(VETOR & x, float alpha, VETOR y) {
	//guarda indíces usados em x em um vetor denso
	for (int i = 0; i < x.size(); i++)
		p_full[x.ind[i]] = i;


	//para todo yi soma alpha*yi correspondente em x
	//se nao houver posicao em x, adiciona e atualiza vetor denso
	for (int i = 0; i < y.size(); i++) {
		if (p_full[y.ind[i]] != -1) {
			x.vals[p_full[y.ind[i]]] += alpha * y.vals[i];
		}
		else {
			p_full[y.ind[i]] = x.size();
			x.ind.push_back(y.ind[i]);
			x.vals.push_back(alpha * y.vals[i]);
		}
	}
	
	//recupera valores originais do vetor auxiliar denso de indices
	for (int i = 0; i < x.size(); i++) 
		p_full[x.ind[i]] = -1;
}

void soma_empac_seq(VETOR & x, vector<float> alpha, vector<VETOR> y){
	for (int i = 0; i < x.size(); i++)
		p_full[x.ind[i]] = i;

	for (int iy = 0; iy < y.size(); iy++)
	{
		for (int i = 0; i < y.size(); i++) {
			if (p_full[y[iy].ind[i]] != -1) {
				x.vals[p_full[y[iy].ind[i]]] += alpha[iy] * y[iy].vals[i];
			}
			else {
				p_full[y[iy].ind[i]] = x.size();
				x.ind.push_back(y[iy].ind[i]);
				x.vals.push_back(alpha[iy] * y[iy].vals[i]);
			}

		}
	}

	for (int i = 0; i < x.size(); i++)
		p_full[x.ind[i]] = -1;
}

void ord_vetor(VETOR &x) {
	vector<int> work(N, 0);

	for (int i = 0; i < x.size(); i++) {
		int j = x.ind[i];
		work[j]++;
	}

	work[0]++;
	for (int i = 1; i < N; i++)
		work[i] += work[i - 1];
	
	vector<float> value_new(x.size(), 0);
	vector<int> map(x.size(), 0);

	for (int k = 0; k < x.size(); k++){
		int j = x.ind[k];
		int k_new = work[j] - 1;
		work[j] = k_new;
		value_new[k_new-1] = x.vals[k];
		map[k_new-1] = j;
	}
	x.ind = map;
	x.vals = value_new;
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

MATRIZ_LINKED_COL converter_matriz(MATRIZ A)
{
	MATRIZ_LINKED_COL B;

	vector<int> col_start(N, -1), row_index(A.size()), link(A.size(),-1);

	for (int k = 0; k < A.size(); k++){
		//if (col_start[A.jcn[k]] == -1) {
		//	col_start[A.jcn[k]] = k;
		//}

		row_index[k] = A.irn[k];
		link[k] = col_start[A.jcn[k]];
		col_start[A.jcn[k]] = k;
	}




	return B;
}
