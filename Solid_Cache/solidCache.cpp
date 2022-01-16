#include <iostream>


#include <fstream>

#define Tags(i,j,k)  Tags[(i*n_ways*2)+(j*2)+k]

#define my_optimize1 false

#include <sstream>

//#include <pthread.h>
#include<thread>		//To compile programs with std::thread support use:
//g++ -std=c++11 -pthread

using namespace std;

typedef unsigned long ulong;
//typedef unsigned long long int ulong;			// if using windows OS uncomment this.
typedef unsigned int uint;

string filename00 = "_result";
string filename = "_result.csv";



class Way {
public:
	int ID = 0;
	int rank = 0;
	ulong Tag = 0;
	ulong BlockSize = 0;
	ulong Miss = 0;
	ulong Access = 0;

	Way(int _ID, int _rank, ulong _blocksize) {
		ID = _ID;
		rank = _rank;
		BlockSize = _blocksize;
		Miss = 0;
		Access = 0;
	}

	bool check(ulong _Tag) {
		Access++;
		if (_Tag == Tag) {
			return true;
		}
		Miss++;
		return false;
	}
	void update(ulong _Tag, int _rank = 0) {
		Tag = _Tag;
		rank = _rank;
	}
};

class CSet {
public:
	uint Index = 0;
	uint Number_of_ways = 0;
	Way** Ways;

	ulong Miss = 0, Access = 0;

	CSet(uint _Index, uint _Number_of_ways, ulong block_size) {
		Index = _Index;
		Number_of_ways = _Number_of_ways;
		Miss = 0;
		Access = 0;

		Ways = new Way*[_Number_of_ways];
		for (int i = 0; i < _Number_of_ways; i++) {
			Ways[i] = new Way(i, i, block_size);
		}
	}

	bool check(ulong _Tag) {
		Access++;
		for (int i = 0; i < Number_of_ways; i++) {
			if (Ways[i]->check(_Tag) == true) {
				return true;
			}
		}
		Miss++;
		return false;
	}
	void Put(ulong _Tag, int _rank = 0) {
		//LRU:
		for (int i = 0; i < Number_of_ways; i++) {
			if (Ways[i]->rank >= Number_of_ways) {
				cout << "\nErrorrrrr!!!!\n";
			}

			if (Ways[i]->rank == Number_of_ways - 1) {
				Ways[i]->update(_Tag, _rank);
			}
			else {
				Ways[i]->rank++;
			}
		}
	}
};

class Cache
{
private:
	string Cache_hitmiss_penalty()
	{
		string txt = "";
		txt += to_string(Miss) + ",";
		txt += to_string(Access) ;
		return txt;
	}

	static ulong rwmiss(Cache* a, Cache* b, Cache* c) {
		return ((a->Miss + b->Miss + c->Miss)/3);
	}

public:
	ulong Miss = 0;
	ulong Access = 0;
	uint n_sets = 0;
	uint block_size = 0;
	uint n_ways = 0;
	CSet** Sets;


	Cache(uint setSize, uint block_size_in_byte, uint num_ways)
	{
		n_sets = setSize;
		block_size = block_size_in_byte;
		n_ways = num_ways;
		Miss = 0;
		Access = 0;


		Sets = new CSet*[setSize];
		for (int i = 0; i < setSize; i++) {
			Sets[i] = new CSet(i, num_ways, block_size_in_byte);
		}
	}

	bool check_and_put_Data(ulong address)
	{
		ulong index_bits = (address / block_size) % n_sets;
		ulong tag = (address / block_size) / n_sets;

		if (Sets[index_bits]->Index != index_bits) {
			cout << "\nError!!!!!!!!!!!!!\n";
			return false;
		}

		Access++;
		if (Sets[index_bits]->check(tag) == true) {
			return true;
		}
		else {
			Sets[index_bits]->Put(tag, 0);
			Miss++;
			return false;
		}
	}

	static void info(Cache* A0, Cache* A1, Cache* A2, Cache* B0, Cache* B1, Cache* B2, Cache* C0, Cache* C1, Cache* C2, Cache* A01, Cache* A11, Cache* A21, Cache* B01, Cache* B11, Cache* B21, Cache* C01, Cache* C11, Cache* C21, uint dm, uint dk, uint dn, uint _bsize, uint _asso)
	{
		char buffer[200];
		snprintf(buffer, sizeof(buffer), "%d,%d,%d,%d,%d,", dm, dk, dn, _bsize, _asso);

		string I = "";

		I = string(buffer);


		I += A0->Cache_hitmiss_penalty() + "," + B0->Cache_hitmiss_penalty() + "," + C0->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A0, B0, C0)) + ",";
		I += A1->Cache_hitmiss_penalty() + "," + B1->Cache_hitmiss_penalty() + "," + C1->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A1, B1, C1)) + ",";
		I += A2->Cache_hitmiss_penalty() + "," + B2->Cache_hitmiss_penalty() + "," + C2->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A2, B2, C2)) + ",0,";

		I += A01->Cache_hitmiss_penalty() + "," + B01->Cache_hitmiss_penalty() + "," + C01->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A01, B01, C01)) + ",";
		I += A11->Cache_hitmiss_penalty() + "," + B11->Cache_hitmiss_penalty() + "," + C11->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A11, B11, C11)) + ",";
		I += A21->Cache_hitmiss_penalty() + "," + B21->Cache_hitmiss_penalty() + "," + C21->Cache_hitmiss_penalty() + "," + to_string(rwmiss(A21, B21, C21)) + "\n";
		cout << filename << ":\t" << I << endl;

		ofstream ofs;
		ofs.open(filename, std::ofstream::out | std::ofstream::app);
		ofs << I.c_str();
		ofs.close();
	}
	static void Trash_info(Cache* ccc) {
		string txt = "";
		for (uint i = 0; i < ccc->n_sets; i++) {
			txt += to_string(i) + ",>>>>,";
			for (uint j = 0; j < ccc->n_ways; j++) {
				txt += to_string(ccc->Sets[i]->Ways[j]->Miss) + "," + to_string(ccc->Sets[i]->Ways[j]->Access) + ",|,";
			}
			txt += "\n";
		}
		txt += "***\n\n";
		ofstream ofs;
		ofs.open(filename + "__Trash.csv", std::ofstream::out | std::ofstream::app);
		ofs << txt.c_str();
		ofs.close();
	}

};

inline ulong getElementAddress(ulong a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
{
	ulong k0 = a0 + (ulong)((i*column_in_row) + j)*floatsize;
	return k0;
}

void inner_product0(Cache* A0, Cache* B0, Cache* C0, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_k = 0;
	int step_k = 1;

	int init_n = 0;
	int step_n = 1;

	//inner product
	for (int m = 0; m < M; m++)
	{

		for (int n = init_n; n < N && n >=0; n+= step_n)
		{

			C0->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
																		 //float REGISTER = C[m, n];
			for (int k = init_k; k < K && k >= 0; k+= step_k)
			{
				A0->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
				B0->check_and_put_Data(getElementAddress(B, n, k, K));       //read B[k,n] =====> B.T[n,k]

																			 //REGISTER += A[m, k] * B[k, n]; ===> REGISTER += A[m, k] * B.T[n,k];

			}
			if (step_k == 1 && my_optimize1) {
				step_k = -1;
				init_k = K - 1;
			}
			else 
			{
				step_k = 1;
				init_k = 0;
			}


			C0->check_and_put_Data(getElementAddress(C, m, n, N));       //write C[m,n]
		}
		if (step_n == 1 && my_optimize1) {
			step_n = -1;
			init_n = N - 1;
		}
		else
		{
			step_n = 1;
			init_n = 0;
		}
	}
}
void outer_product0(Cache* A1, Cache* B1, Cache* C1, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_m = 0;
	int step_m = 1;

	int init_n = 0;
	int step_n = 1;

	//outer product
	for (int k = 0; k < K; k++)
	{
		for (int m = init_m; m < M &&  m >= 0 ; m+=step_m)
		{
			A1->check_and_put_Data(getElementAddress(A, k, m, M));       //read A[m,k] =====> A.T[k,m]
																		 //float REGISTER = A.T[k,m];
			for (int n = init_n; n < N && n >= 0; n += step_n)
			{

				C1->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]      
				B1->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																			 //C[m, n] += REGISTER * B[k, n];

				C1->check_and_put_Data(getElementAddress(C, m, n, N));       //write C[m,n]

			}
			if (step_n == 1 && my_optimize1) {
				step_n = -1;
				init_n = N - 1;
			}
			else
			{
				step_n = 1;
				init_n = 0;
			}
			//Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
		}
		if (step_m == 1 && my_optimize1) {
			step_m = -1;
			init_m = M - 1;
		}
		else
		{
			step_m = 1;
			init_m = 0;
		}
	}
}
void gustavson0(Cache* A2, Cache* B2, Cache* C2, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_k = 0;
	int step_k = 1;

	int init_n = 0;
	int step_n = 1;
	//Gustavson product
	for (int m = 0; m < M; m++)
	{
		for (int k = init_k; k < K && k >= 0; k += step_k)
		{
			A2->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																		 //float REGISTER = A[m, k];
			for (int n = init_n; n < N && n >= 0; n += step_n)
			{
				C2->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
				B2->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																			 //C[m, n] += REGISTER * B[k, n];

				C2->check_and_put_Data(getElementAddress(C, m, n, N));       //write C[m,n]

			}
			if (step_n == 1 && my_optimize1) {
				step_n = -1;
				init_n = N - 1;
			}
			else
			{
				step_n = 1;
				init_n = 0;
			}
			//Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

		}
		if (step_k == 1 && my_optimize1) {
			step_k = -1;
			init_k = K - 1;
		}
		else
		{
			step_k = 1;
			init_k = 0;
		}
	}
}

void inner_product1(Cache* A0, Cache* B0, Cache* C0, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_k = 0;
	int step_k = 1;

	int init_m = 0;
	int step_m = 1;

	//inner product


	for (int n = 0; n < N; n++)
	{
		for (int m = init_m; m < M && m >= 0; m += step_m)
		{

			C0->check_and_put_Data(getElementAddress(C, n, m, M));       //read C[m,n] =======> C.T[n,m]
																		 //float REGISTER = C[m, n];
			for (int k = init_k; k < K && k >= 0; k += step_k)
			{
				A0->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
				B0->check_and_put_Data(getElementAddress(B, n, k, K));       //read B[k,n] =====> B.T[n,k]

																			 //REGISTER += A[m, k] * B[k, n]; ===> REGISTER += A[m, k] * B.T[n,k];

			}
			if (step_k == 1 && my_optimize1) {
				step_k = -1;
				init_k = K - 1;
			}
			else
			{
				step_k = 1;
				init_k = 0;
			}


			C0->check_and_put_Data(getElementAddress(C, n, m, M));       //write C[m,n] =======> C.T[n,m]
		}
		if (step_m == 1 && my_optimize1) {
			step_m = -1;
			init_m = M - 1;
		}
		else
		{
			step_m = 1;
			init_m = 0;
		}
	}
}
void outer_product1(Cache* A1, Cache* B1, Cache* C1, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_m = 0;
	int step_m = 1;

	int init_n = 0;
	int step_n = 1;

	//outer product
	for (int k = 0; k < K; k++)
	{
		for (int n = init_n; n < N && n >= 0; n += step_n)
		{
			B1->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]
																		// REGISTER = B[k,n]

			for (int m = init_m; m < M && m >= 0; m += step_m)
			{

				C1->check_and_put_Data(getElementAddress(C, n, m, M));       //read C[m,n] =======> C.T[n,m]      
				A1->check_and_put_Data(getElementAddress(A, k, m, M));       //read A[m,k] =====> A.T[k,m]

																			 // C.T[n,m] += REGISTER *  A.T[k,m];

				C1->check_and_put_Data(getElementAddress(C, n, m, M));       //write C.T[n,m]

			}
			if (step_m == 1 && my_optimize1) {
				step_m = -1;
				init_m = M - 1;
			}
			else
			{
				step_m = 1;
				init_m = 0;
			}

			//Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
		}
		if (step_n == 1 && my_optimize1) {
			step_n = -1;
			init_n = N - 1;
		}
		else
		{
			step_n = 1;
			init_n = 0;
		}
	}
}
void gustavson1(Cache* A2, Cache* B2, Cache* C2, int M, int K, int N, ulong A, ulong B, ulong C)
{
	int init_k = 0;
	int step_k = 1;

	int init_m = 0;
	int step_m = 1;
	//Gustavson product
	for (int n = 0; n < N; n++)
	{
		for (int k = init_k; k < K && k >= 0; k += step_k)
		{
			B2->check_and_put_Data(getElementAddress(B, n, k, K));       //read B[k,n] =====> B.T[n,k]		
																		 //float REGISTER = B.T[n,k];
			for (int m = init_m; m < M && m >= 0; m += step_m)
			{
				C2->check_and_put_Data(getElementAddress(C, n, m, M));       //read C.T[n,m]

				A2->check_and_put_Data(getElementAddress(A, k, m, M));       //read A.T[k,m]

																			 //C.T[n,m] += REGISTER * A.T[k,m];

				C2->check_and_put_Data(getElementAddress(C, n, m, M));       //write C[m,n]

			}
			if (step_m == 1 && my_optimize1) {
				step_m = -1;
				init_m = M - 1;
			}
			else
			{
				step_m = 1;
				init_m = 0;
			}
			//Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

		}
		if (step_k == 1 && my_optimize1) {
			step_k = -1;
			init_k = K - 1;
		}
		else
		{
			step_k = 1;
			init_k = 0;
		}
	}
}



static void delta(uint dm, uint dk, uint dn, uint _bsize, uint _asso, int* Sizes)
{
	int M = dm;
	int K = dk;
	int N = dn;

	ulong k0, k1, k2;
	k0 = 0xF0000000;
	k1 = k0 + (ulong)(M * K * sizeof(float));		//k1 = k0 + (ulong)(M * K * sizeof(float) + (rand() % 16) * 256);
	k2 = k1 + (ulong)(K * N * sizeof(float));		//k2 = k1 + (ulong)(K * N * sizeof(float) + (rand() % 16) * 256);

	ulong A = k0;
	ulong B = k1;
	ulong C = k2;


	int sizeAi = Sizes[0];
	int sizeBi = Sizes[1];
	int sizeCi = Sizes[2];

	Cache A0 = Cache(sizeAi, _bsize, _asso);
	//Cache B0 = Cache(sizeBi, _bsize, _asso);
	//Cache C0 = Cache(sizeCi, _bsize, _asso);
	Cache A01 = Cache(sizeAi, _bsize, _asso);
	//Cache B01 = Cache(sizeBi, _bsize, _asso);
	//Cache C01 = Cache(sizeCi, _bsize, _asso);
	///////////////////////////////////////

	int sizeAo = Sizes[3];
	int sizeBo = Sizes[4];
	int sizeCo = Sizes[5];


	Cache A1 = Cache(sizeAo, _bsize, _asso);
	//Cache B1 = Cache(sizeBo, _bsize, _asso);
	//Cache C1 = Cache(sizeCo, _bsize, _asso);
	Cache A11 = Cache(sizeAo, _bsize, _asso);
	//Cache B11 = Cache(sizeBo, _bsize, _asso);
	//Cache C11 = Cache(sizeCo, _bsize, _asso);
	///////////////////////////////////////

	int sizeAg = Sizes[6];
	int sizeBg = Sizes[7];
	int sizeCg = Sizes[8];


	Cache A2 = Cache(sizeAg, _bsize, _asso);
	//Cache B2 = Cache(sizeBg, _bsize, _asso);
	//Cache C2 = Cache(sizeCg, _bsize, _asso);
	Cache A21 = Cache(sizeAg, _bsize, _asso);
	//Cache B21 = Cache(sizeBg, _bsize, _asso);
	//Cache C21 = Cache(sizeCg, _bsize, _asso);
	///////////////////////////////////////

	thread t0(inner_product0, &A0, &A0, &A0, M, K, N, A, B, C);
	thread t1(outer_product0, &A1, &A1, &A1, M, K, N, A, B, C);
	thread t2(gustavson0, &A2, &A2, &A2, M, K, N, A, B, C);
	
	thread t3(inner_product1, &A01, &A01, &A01, M, K, N, A, B, C);
	thread t4(outer_product1, &A11, &A11, &A11, M, K, N, A, B, C);
	thread t5(gustavson1, &A21, &A21, &A21, M, K, N, A, B, C);

	t0.join();
	t1.join();
	t2.join();
	t3.join();
	t4.join();
	t5.join();

	////////////Save results:///////////////////////////////////////////////////////////////////////////////////////////
	Cache::info(&A0, &A1, &A2, &A0, &A1, &A2, &A0, &A1, &A2, &A01, &A11, &A21, &A01, &A11, &A21, &A01, &A11, &A21, dm, dk, dn, _bsize, _asso);

	//Cache::Trash_info(&A0);
	//Cache::Trash_info(&B0);
	//Cache::Trash_info(&C0);

	//Cache::Trash_info(&A1);
	//Cache::Trash_info(&B1);
	//Cache::Trash_info(&C1);

	//Cache::Trash_info(&A2);
	//Cache::Trash_info(&B2);
	//Cache::Trash_info(&C2);
}

int main(int argc, char *argv[])
{
	if (argc != 6) {
		cout << "Invalid input arguments error!\nYou should input 5 integer arguments as:\nM K N  Block_size(B) #_of_Ways\nExit -1" << endl;
		//cin >> argc;
		return -1;
	}

	uint Data[5];
	for (int i = 0; i < 5; i++) {
		std::stringstream ss(argv[i + 1]);
		if (ss >> Data[i])
			std::cout << "Arg" << i << " is: " << Data[i] << endl;
		else
			std::cout << "error";
	}

	/////////////////
	string Sizes = "";
	int k[9];
	int ar_c = 0, li_c = 0;
	ifstream ifs;
	ifs.open("Sizes.txt", std::ofstream::in);
	while (getline(ifs, Sizes)) {
		li_c++;
		cout << "line " << li_c << endl;
		stringstream stream(Sizes);
		ar_c = 0;
		while (stream >> k[ar_c]) {
			//cout << k[ar_c] << endl;
			ar_c++;
		}
		filename = filename00 + "_line_" + to_string(li_c) + ".csv";

		std::ifstream file(filename.c_str());
		if (!file) {
			string I = "M,K,N,BlockSize(B),Ways,Inner_A_miss,Inner_A_access,Inner_B_miss,Inner_B_access,Inner_C_miss,Inner_C_access,Inner_total_miss,";
			I += "Outer_A_miss,Outer_A_access,Outer_B_miss,Outer_B_access,Outer_C_miss,Outer_C_access,Outer_total_miss,";
			I += "Gus_A_miss,Gus_A_access,Gus_B_miss,Gus_B_access,Gus_C_miss,Gus_C_access,Gus_total_miss,";

			I += "XXX,iInner_A_miss,iInner_A_access,iInner_B_miss,iInner_B_access,iInner_Ci_miss,iInner_C_access,iInner_total_miss,";
			I += "iOuter_A_miss,iOuter_A_access,iOuter_B_miss,iOuter_B_access,iOuter_C_miss,iOuter_C_access,iOuter_total_miss,";
			I += "iGus_A_miss,iGus_A_access,iGus_B_miss,iGus_B_access,iGus_C_miss,iGus_C_access,iGus_total_miss\n";
			ofstream ofs;
			ofs.open(filename, std::ofstream::out);
			ofs << I.c_str();
			ofs.close();
		}

		delta(Data[0], Data[1], Data[2], Data[3], Data[4], k);
	}
	ifs.close();
	//////////////////////////

	cout << "finished successfully!" << endl << endl;
	cin >> Data[0];
	return 0;
}

