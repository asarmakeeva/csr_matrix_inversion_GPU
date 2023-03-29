#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <math.h>
#include <sys/time.h>

using namespace std;

void read(double **AV1,int **ANC1,int **ANL1,double **res1,int *dim,int *NumEl)
{
    double *AV, *b, *res;
    int *ANC, *ANL, *ANR;
cout<<"Point1"<<endl;
    ifstream data;//, data1;
    int i,j,ik,jk,l,n,nc,m;
    double k;
    data.open("matrix.mtx");
//	data1.open("rhs");
    char s[255];
	data.getline(s,255);//ñ÷èòûâàåì áëàíê ìàòðèöû
        short flag=1;
    while ((flag)&&(!data.eof()))
    {
        data.getline(s,255);
        if (s[0]!='%')
            flag=0;
    }

        int len = strlen(s);
        int size[3];
        j=0;
        i=0;
        while (i<len)
        {
            n=0;
            while( (s[i]!=' ') && (i<len) )
            {
                n=n*10 + s[i] - 0x30;
                i++;
            }
            size[j++]=n;
            i++;
        }

        n=size[0];
        nc=size[1];
        m=size[2];

//	data>>n>>nc>>m; //êîë-âî ñòðîê, ñòîëáöîâ è íå íóëåâûõ ýëëåìåíòîâ ìàòðèöû
    AV=new double [m];
    ANC=new int [m];
    ANR=new int [m];
    ANL=new int [n+1];
    res=new double [n];
    b=new double [n];
    for(i=0;i<n;i++)
    {
	res[i]=0;		//èíèöèàëèçàöèÿ âåêòîðà ðåçóëüòàòîâ
	//data1>>b[i];	//ñ÷èòûâàíèå âåêòîðà èç ôàéëà
    }
    
    for(i=0;i<(n+1);i++)
        ANL[i]=0;
    
cout<<"Point2"<<endl;
    l=-1;
    //ñ÷èòûâàíèå ìàòðèöû èç ôàéëà â ôîðìàò CSR
    for(i=0;i<m;i++)
    {
	data>>ik>>jk>>k;
	AV[i]=k;
	ANC[i]=jk-1;
                ANR[i]=ik-1;
//		if (l!=(ik-1))
//		{
//			ANL[++l]=i;
//		}
    }
    cout<<"Check1"<<endl;
    for(i=0;i<m;i++)
        for(int j=0;j<m-1;j++)
            if(ANR[j]>ANR[j+1])
            {
                ik=ANR[j];
                ANR[j]=ANR[j+1];
                ANR[j+1]=ik;
                
                jk=ANC[j];
                ANC[j]=ANC[j+1];
                ANC[j+1]=jk;
                
                k=AV[j];
                AV[j]=AV[j+1];
                AV[j+1]=k;
            }
    
/*ofstream jjj;
jjj.open("matrix1.mtx");
jjj<<"%aaa"<<endl;
jjj<<n<<" "<<n<<" "<<m<<endl<<endl;
for(i=0;i<m;i++)
    jjj<<ANR[i]<<" "<<ANC[i]<<" "<<AV[i]<<endl;
jjj.close();
*/
    cout<<"Check2"<<endl;
    for(i=0;i<m;i++)
    {
        ANL[ANR[i]+1]++;
    }

/*        ik=0; //t
    jk=0; //s
        unsigned int t1 = 0, s1 = 0;
        
    for(i=0 ; i <= n ; i++)
        {
            ik = ANL[i];
        ANL[i]=jk;
        jk+=ik;
    }
*/

    for(i=1; i<=n ; i++)
        ANL[i]+=ANL[i-1];

        
        for(i=0;i<n;i++)
        {
        b[i]=0;
          for (int j=ANL[i]; j<ANL[i+1]; j++)
            b[i]+=AV[j];
        }
cout<<"Point3"<<endl;
    
    
    for(i=ANL[0]; i<ANL[1]; i++)
        cout<<AV[i]<<"   "<<ANC[i]<<endl;
    cout<<"&&&&"<<endl;
                                         
    for(i=ANL[1]; i<ANL[2]; i++)
        cout<<AV[i]<<"   "<<ANC[i]<<endl;
        
    ANL[n]=m;
    *dim=n;
    *NumEl=m;
    *AV1=AV;
    *ANC1=ANC;
    *ANL1=ANL;
    *res1=res;
    data.close();
//	data1.close();
}



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
void JG(double *AV, int *IA, int *JA, int *dim, double ***AV_11, int *NumEl) //, int **IA_11, int **JA_11,int str, double **mass, double **&LU,double **&A_1 )
		{
	int n, m;
	n = *dim;
	double ttime;
	double **JG,**E, buf, **check;
	m = *NumEl;

	JG = new double*[n];
	for (int i = 0; i < n; i++)
		{JG[i] = new double[n];}

	E = new double*[n];
	for (int i = 0; i < n; i++)
		E[i] = new double[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}
	cout << "check in JG " << n << " " << m << endl;


		for(int i = 0; i<n; i++ )
	        for (int j = 0;j<n;j++)
		    JG[i][j]=0;

	for(int i = 0; i<n; i++ )
	    for (int j = JA[i];j<JA[i+1];j++)
		    JG[i][IA[j]]=AV[j];

	timeval tim;
	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
//JORDANO_GAUSSA_METHOD
//PRIVODIM K TREUGOLNOMY VIDY
	for (int k = 0; k < n; k++)
	    {
	        buf = JG[k][k];
	        for (int j = 0; j < n; j++)
	            {
	                JG[k][j] /= buf;
	                E[k][j] /= buf;
	            }
	        for (int i = k + 1; i < n; i++)
	        {
	            buf = JG[i][k];
	            for (int j = 0; j < n; j++)
	            {
	                JG[i][j] -= JG[k][j]*buf;
	                E[i][j] -= E[k][j]*buf;
	            }
	        }
	    }
	gettimeofday(&tim, NULL);
	double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
//OBRATII HOD
// obrazuem nuli nad glavnoi diagonal'u
	for (int k = n - 1; k > 0; k--)
	{
	        for (int i = k - 1; i >= 0; i--)
	        {
	            buf = JG[i][k];
	                for (int j = 0; j < n; j++)
	                {
	                    JG[i][j] -= JG[k][j]*buf;
	                    E[i][j] -= E[k][j]*buf;
	                }
	        }
        }
//
//END JORDANO_GAUSSA_METHOD
//
	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	ttime = t2 - t1;
	cout << "time = " << ttime << endl;
	double time_tr,time_gl_diog;
	time_tr = t3-t1;
	time_gl_diog = t2-t3;
	cout<<"privodim k treugolnomy vidy = "<<time_tr<<endl;
	cout<<"obrazuem nuli nad gl diog = "<<time_gl_diog<<endl;
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
				check[i][k] += E[i][IA[j]] * AV[j];
		}
	}


double  sum_buff,buff=0;
//vichitaem 1 i vichislaem normy
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
//cout<<"norma= "<<sum_buff<<endl;
}


cout<<"norma!!!= "<<buff<<endl;

	double Norma_invA=0;double Sum_norm = 0;
	for(int i=0;i<n;i++)
	{    for(int j=0; j<n; j++)
	   {
		Sum_norm+=fabs(check[i][j])*fabs(check[i][j]);
	    }
	}
	    Norma_invA = sqrt(Sum_norm);

cout<<"norma sum= "<<Norma_invA<<endl;


	ofstream write;
	write.open("check.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			write << check[i][j] << " ";
		}
		write << endl;
	}
	write.close();
	delete[] JG;
	*AV_11 = E;
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
	read(&AV, &IA, &JA, &buf, &dim, &NumEl);
//	read1(&AV, &IA, &JA, &buf, &dim, &NumEl);
//	read3(&AV, &IA, &JA, &buf, &dim, &NumEl);
//	read2(&AV, &IA, &JA, &buf, &dim, &NumEl);
	cout << "read2" << endl;
	JG(AV, IA, JA, &dim, &AV_1, &NumEl);
	cout << "read3" << endl;
	write(AV_1, dim);
	cout << "read4" << endl;
	delete[] IA;
	delete[] AV;
	delete[] JA;
	delete[] buf;
	return 0;
}
