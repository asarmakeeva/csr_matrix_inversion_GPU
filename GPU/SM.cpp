#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <math.h>
#include <sys/time.h>

using namespace std;

void read1(double **AV1, int **IA1, int **JA1, double **buf1, int *dim,
		int *NumEl) {
	int n, m, nc, i, ik, jk, k, l; //,
	double *AV, *buf;
	int *IA, *JA;
	ifstream cout1;
	cout1.open("new.mtx");
	char s[255];
	cout1.getline(s, 255);
	cout1 >> n >> nc >> m;
	AV = new double[m];
	IA = new int[m];
	JA = new int[n + 1];
	buf = new double[n];
	l = -1;
	for (i = 0; i < m; i++) {
		cout1 >> ik >> jk >> k;
		AV[i] = k;
		IA[i] = jk - 1;
		if (l != (ik - 1)) {
			JA[++l] = i;
		}
	}
	for (i = 0; i < n; i++) {
		buf[i] = 0;
		for (int j = JA[i]; j < JA[i + 1]; j++)
			buf[i] += AV[j];
	}
	cout << "check2" << endl;
	JA[n] = m;
	*dim = n;
	*NumEl = m;
	*AV1 = AV;
	*IA1 = IA;
	*JA1 = JA;
	*buf1 = buf;
	cout1.close();

}
void read3(double **AV1, int **IA1, int **JA1, double **buf1, int *dim,
		int *NumEl) {
	double *AV, *V, *buf;
	int *IA, *NC, *JA, *NR;

	ifstream data, data1;
	int i, j, ik, jk, n, nc, m, m1;
	double k;
	data.open("new.mtx");
	char s[255], c[255];
	data.getline(s, 255);
	short flag = 1;
	while ((flag) && (!data.eof())) {
		data.getline(c, 255);
		if (c[0] != '%')
			flag = 0;
	}

	cout << "check2" << endl;
	int len = strlen(c);
	cout << len << "  " << c << endl;
	int size[3];
	j = 0;
	i = 0;
	while (i < len) {
		n = 0;
		while ((c[i] != ' ') && (i < len)) {
			n = n * 10 + c[i] - 0x30;
			i++;
		}
		size[j++] = n;
		i++;
	}

	n = size[0];
	nc = size[1];
	m1 = size[2];
	m = 2 * m1 - n;

	cout << "!!" << n << "   " << nc << "   " << m1 << "   " << m << endl;

	V = new double[m1];
	NC = new int[m1];
	NR = new int[m1];
	AV = new double[m];
	IA = new int[m];
	JA = new int[n + 1];
	buf = new double[n];
	for (i = 0; i < n; i++) {
		JA[i] = 0;
	}
	JA[n] = 0;
	i = 0;
	while (!data.eof()) {
		data >> ik >> jk >> k;
		V[i] = k;
		NC[i] = jk - 1;
		NR[i] = ik - 1;
		JA[ik]++;
		if (ik != jk)
			JA[jk]++;
		i++;
	}

	cout << "check2" << endl;

	ik = 0;
	jk = 0;
	for (i = 1; i <= n; i++) {
		ik = JA[i];
		JA[i] = jk;
		jk += ik;
	}

	cout << "check3" << endl;

	for (i = 0; i < m1; i++) {
		int col = NC[i];
		int row = JA[NR[i] + 1];
		double val = V[i];
		AV[row] = val;
		IA[row] = col;
		JA[NR[i] + 1]++;
		if (NR[i] != col) {
			AV[JA[col + 1]] = val;
			IA[JA[col + 1]] = NR[i];
			JA[col + 1]++;
		}
	}

	cout << "check4" << endl;

/*	for (i = 0; i < n; i++) {
		buf[i] = 0;
		for (j = JA[i]; j < JA[i + 1]; j++)
			buf[i] += AV[j];
	}*/
//        delete [] V;
//        delete [] NC;
//        delete [] NR;

	*dim = n;
	*NumEl = m;
	*AV1 = AV;
	*IA1 = IA;
	*JA1 = JA;
	*buf1 = buf;
	data.close();
}

void read2(double **AV1, int **IA1, int **JA1, double **buf1, int *dim,
		int *NumEl) {
	double *AV, *V, *buf;
	int *IA, *NC, *JA, *NR;

	ifstream data, data1;
	int i, j, ik, jk, n, nc, m, m1;
	double k;
	data.open("new.mtx");
	char s[255], c[255];
	data.getline(s, 255);
	short flag = 1;
	while ((flag) && (!data.eof())) {
		data.getline(c, 255);
		if (c[0] != '%')
			flag = 0;
	}

	cout << "check2" << endl;
	int len = strlen(c);
	cout << len << "  " << c << endl;
	int size[3];
	j = 0;
	i = 0;
	while (i < len) {
		n = 0;
		while ((c[i] != ' ') && (i < len)) {
			n = n * 10 + c[i] - 0x30;
			i++;
		}
		size[j++] = n;
		i++;
	}

	n = size[0];
	nc = size[1];
	m = size[2];
	m1 = (m+n)/2;

	cout << "!!" << n << "   " << nc << "   " << m1 << "   " << m << endl;

	V = new double[m1];
	NC = new int[m1];
	NR = new int[m1];
	AV = new double[m];
	IA = new int[m];
	JA = new int[n + 1];
	buf = new double[n];
	for (i = 0; i < n; i++) {
		JA[i] = 0;
	}
	JA[n] = 0;
	i = 0;
	while (!data.eof()) {
		data >> ik >> jk >> k;
		V[i] = k;
		NC[i] = ik - 1;
		NR[i] = jk - 1;
		JA[ik]++;
		if (ik != jk)
			JA[jk]++;
		i++;
	}

	cout << "check2" << endl;

	ik = 0;
	jk = 0;
	for (i = 1; i <= n; i++) {
		ik = JA[i];
		JA[i] = jk;
		jk += ik;
	}

	cout << "check3" << endl;

	for (i = 0; i < m1; i++) {
		int col = NC[i];
		int row = JA[NR[i] + 1];
		double val = V[i];
		AV[row] = val;
		IA[row] = col;
		JA[NR[i] + 1]++;
		if (NR[i] != col) {
			AV[JA[col + 1]] = val;
			IA[JA[col + 1]] = NR[i];
			JA[col + 1]++;
		}
	}

	cout << "check4" << endl;

/*	for (i = 0; i < n; i++) {
		buf[i] = 0;
		for (j = JA[i]; j < JA[i + 1]; j++)
			buf[i] += AV[j];
	}*/
//        delete [] V;
//        delete [] NC;
//        delete [] NR;

	*dim = n;
	*NumEl = m;
	*AV1 = AV;
	*IA1 = IA;
	*JA1 = JA;
	*buf1 = buf;
	data.close();
}
void SM(double *AV, int *IA, int *JA, int *dim, double ***AV_11, int *NumEl) //, int **IA_11, int **JA_11,int str, double **mass, double **&LU,double **&A_1 )
		{
	int n, m;
	n = *dim;
	double ttime;
	double **SM,**A_1, **check;
	m = *NumEl;
	double *v = new double[n];
	double *buf2 = new double[n];
	double *buf = new double[n];

	SM = new double*[n];
	for (int i = 0; i < n; i++)
		SM[i] = new double[n];

	A_1 = new double*[n];
	for (int i = 0; i < n; i++)
		A_1[i] = new double[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				A_1[i][j] = 1;
			else
				A_1[i][j] = 0;
		}
	}
	cout << "check in JG " << n << " " << m << endl;


		for(int i = 0; i<n; i++ )
	        for (int j = 0;j<n;j++)
		    SM[i][j]=0;
	for(int i = 0; i<n; i++ )
	{    for (int j = JA[i];j<JA[i+1];j++)
		   SM[i][IA[j]]=AV[j];}


	timeval tim;
	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
for(int k=0;k<n;k++)
{
	for(int i=0; i<n;i++)	    v[i]=SM[k][i];
	    v[k]-=1;
	for(int j=0;j<n;j++)
	{	
	    double bufer=0;
	    for(int i=0;i<n;i++)
	    {
		bufer+=v[i]*A_1[i][j];
	    }
	buf[j]=bufer;
	}
        for(int i=0;i<n;i++)
		buf2[i]=A_1[i][k];

    for(int j=0;j<n;j++)
    {
	for(int i=0;i<n;i++)
	{
		A_1[i][j]-=(buf[j]*buf2[i])/(buf[k]+1);
//cout<<A_1[i][j]<<" ";
	}
//cout<< endl;
    }
//cout<< endl;
}

//
//END METHOD_POPLNENUA
//
	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	ttime = t2 - t1;
	cout << "time = " << ttime << endl;

	check = new double*[n];
	for (int i = 0; i < n; i++)
		check[i] = new double[n];

	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++){
		check[i][k]=0;
}
}

	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++) {
			check[i][k] = 0;
			for (int j = JA[k]; j < JA[k + 1]; j++)
				check[i][k] += A_1[i][IA[j]] * AV[j];
		}
	}

/*
double  sum_buff,buff=0;
//vichitaem 1
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++) {
		if(i==k) check[i][k]-=1;
		}
	}

	for (int i = 0; i < n; i++) {
	sum_buff=0;
		for (int k = 0; k < n; k++) {
		sum_buff+=fabs(check[i][k]);
	    }
	if (buff<sum_buff) buff=sum_buff;
	}


cout<<"norma!!!= "<<buff<<endl;

*/
	ofstream write;
	write.open("check.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			write << check[i][j] << " ";
		}
		write << endl;
	}
	delete[] v;
	delete[] buf;
	delete[] buf2;
	write.close();
	delete[] SM;
	*AV_11 = A_1;
}

void write(double **AV_1, int dim) {
	ofstream wr;
	cout << "check9" << endl;
	wr.open("result.txt");
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			wr << AV_1[i][j] << " "; //������ ���������� � ����
		}
		wr << endl;
	}
	//	*AV_11 = *AV_1;
	cout << "check8" << endl;
	wr.close();
	cout << "check20" << endl;
}
int main() {
	double *AV, *buf, **AV_1;
	int *IA, *JA, dim, NumEl;
//	read1(&AV, &IA, &JA, &buf, &dim, &NumEl);
	read3(&AV, &IA, &JA, &buf, &dim, &NumEl);
//	read2(&AV, &IA, &JA, &buf, &dim, &NumEl);
	cout << "read2" << endl;
	SM(AV, IA, JA, &dim, &AV_1, &NumEl);
	cout << "read3" << endl;
	write(AV_1, dim);
	cout << "read4" << endl;
	delete[] IA;
	delete[] AV;
	delete[] JA;
	delete[] buf;
	return 0;
}
