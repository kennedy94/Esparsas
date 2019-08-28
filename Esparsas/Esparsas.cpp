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

struct VETOR_link {
	vector<float> value;
	vector<int> link;
	int header;
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
struct MATRIZ_VET_COL{
	vector<int> len_col, col_start, row_index;
	vector<float> value;
};

struct MATRIZ_VET_LIN{
	vector<int> lenrow, irowst, jcn;
	vector<float> value;
};

struct VETOR_double_link {
	vector<int> f_link, b_link;
	vector<float> values;
	int f_header, b_header;
};

vector<float> prod_matriz_vetor(MATRIZ A, vector<float> x);

MATRIZ_LINKED_COL converter_matriz(MATRIZ A);

VETOR questao2_1(VETOR x, float alpha, VETOR y);

void questao2_2(VETOR &x, float alpha, VETOR y);

void questao2_3(VETOR &x, float alpha, VETOR y);

void questao2_6(VETOR_link &x, int i);

void questao2_7(VETOR_double_link &x, int i);

int main() {

	/*vector<int> uind = {2, 4, 0}, vind = {0, 2};
	vector<float> uval = {1.0, 2.0, 9.0}, vval = { 4.0, 3.0 };

	VETOR u(uind, uval), v(vind, vval), z(uind,vval);

	vector<float> alpha_vet = {0.25, 0.5};*/

	//cout << produto_interno_ord(u,v) << endl;
	
	//soma_empac(u, 0.5, v);
	/*VETOR Z = questao2_1(u, 0.5, v);

	questao2_3(u, 0.5, v);*/

	//soma_empac_seq(u, alpha_vet, vector<VETOR> {v, z});

	/*vector<int>
		Ai = { 0, 1 },
		Aj = { 0, 1 };
	vector<float> Av = { 1.0, 1.0 };

	vector<float> x = { 2, 0 };

	MATRIZ A(Ai, Aj, Av);

	vector<float> y = prod_matriz_vetor(A, x);*/

	//ord_vetor(u);

	//TESTES PARA MATRIZ------------------------------------------------------------

	/*vector<int> irow = { 4,5,1,1,5,2,4,3,3,2,1 },
		jcol = { 1,2,2,1,5,3,4,5,2,4,5 };

	for (int i = 0; i < irow.size(); i++){
		irow[i]--;
		jcol[i]--;
	}

	vector<float> val = { -1,3,2,1,6,-3,-4,-5,-2,4,5 };

	MATRIZ A(irow, jcol, val);

	converter_matriz(A);*/

	//Questao 2.6------------------------------------------------------------
	vector<float> values = {10, 3, 5, 2};
	/*vector<int> links = {1, 2, 3, -1};
	int header = 0;

	VETOR_link vet;
	vet.header = header;
	vet.link = links;
	vet.value = values;
	
	questao2_6(vet, 1);*/

	//Questao 2.7------------------------------------------------------------
	vector<int> f_links = { 3, 0, -1, 2 },
		b_links = { 1, -1, 3, 0 };
	values = { 3, 10, 2, 5 };

	VETOR_double_link vet_d;
	vet_d.b_header = 2;
	vet_d.f_header = 1;
	vet_d.f_link = f_links;
	vet_d.b_link = b_links;
	vet_d.values = values;

	questao2_7(vet_d, 3);


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
	vector<int> col_start(N, -1), row_index(A.size()), link(A.size(),-1);

	for (int k = 0; k < A.size(); k++){
		row_index[k] = A.irn[k];
		link[k] = col_start[A.jcn[k]];
		col_start[A.jcn[k]] = k;
	}

	MATRIZ_LINKED_COL B;
	B.col_start = col_start;
	B.link = link;
	B.row_index = row_index;
	B.value = A.val;

	return B;
}

VETOR questao2_1(VETOR x, float alpha, VETOR y)
{
	vector<int> zind;
	vector<float> zval;

	int i = 0, j = 0;

	while (true){
		if (i == x.size() && j == y.size()) break;

		if (i < x.size() && j < y.size()){
			if (x.ind[i] == y.ind[j] && i < x.size() && j < y.size()) {
				zval.push_back(x.vals[i] + alpha * y.vals[j]);
				zind.push_back(x.ind[i]);
				i++; j++;
			}
			else {
				if (x.ind[i] < y.ind[j]) {
					zind.push_back(x.ind[i]);
					zval.push_back(x.vals[i]);
					i++;
				}
				else
				{
					zind.push_back(y.ind[j]);
					zval.push_back(alpha * y.vals[j]);
					j++;
				}
			}
		}
		else{
			if (i < x.size()) {
				zind.push_back(x.ind[i]);
				zval.push_back(x.vals[i]);
				i++;
			}
			else{
				zind.push_back(y.ind[j]);
				zval.push_back(alpha * y.vals[j]);
				j++;
			}
		}
	}

	VETOR z(zind, zval);
	return z;
}

void questao2_2(VETOR & x, float alpha, VETOR y){

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

void questao2_3(VETOR & x, float alpha, VETOR y){

	//guarda indíces usados em x em um vetor denso
	for (int i = 0; i < x.size(); i++)
		p_full[x.ind[i]] = i;


	//para todo yi soma alpha*yi correspondente em x
	//se nao houver posicao em x, adiciona e atualiza vetor denso
	for (int i = 0; i < y.size(); i++) {
		x.vals[p_full[y.ind[i]]] += alpha * y.vals[i];
	}

	//recupera valores originais do vetor auxiliar denso de indices
	for (int i = 0; i < x.size(); i++)
		p_full[x.ind[i]] = -1;
}

void questao2_6(VETOR_link &x, int i) {
	if (i == 0) {
		x.header = x.link[x.header];
		x.link[i] = -1;
		x.value[i] = NULL;

		return;
	}

	int cont = x.header,
		anterior;

	for (int j = 0; j < i; j++){
		anterior = cont;
		cont = x.link[cont];
	}

	x.link[anterior] = x.link[cont];

	x.link[cont] = -1;
	x.value[cont] = NULL;
}



void questao2_7(VETOR_double_link &x, int i){
	if (i == 0) {
		int antigo_header = x.f_header;
		x.f_header = x.f_link[antigo_header];

		x.b_link[x.f_header] = -1;
		
		x.f_link[antigo_header] = -1;
		x.values[antigo_header] = NULL;
		x.b_link[antigo_header] = -1;

		return;
	}

	int cont = x.f_header,
		anterior;

	for (int j = 0; j < i; j++) {
		anterior = cont;
		cont = x.f_link[anterior];
	}

	if (cont != x.b_header)	{
		x.f_link[anterior] = x.f_link[cont];

		x.b_link[x.f_link[cont]] = anterior;

		x.f_link[cont] = -1;
		x.b_link[cont] = -1;
		x.values[cont] = 0;
	}
	else {
		x.f_link[anterior] = -1;

		x.f_link[cont] = -1;
		x.b_link[cont] = -1;
		x.values[cont] = 0;

		x.b_header = anterior;

	}

	



}