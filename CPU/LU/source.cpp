#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sys/time.h>

using namespace std;

void read1(double **AV1, int **IA1, int **JA1,double **buf1,int *dim,int *NumEl)
{
	int col=0,num=0,N=0,n,m,nc,i,ik,jk,k,l;//,
	double *AV,*buf;
	int *IA,*JA;
	ifstream cout1;
	cout1.open("new.mtx");
	char s[255];
	cout1.getline(s,255);
	cout1>>n>>nc>>m;
	AV = new double [m];
	IA = new int[m];
	JA = new int[n+1];
	buf = new double [n];
	l=-1;
	for(i=0;i<m;i++)
	{
		cout1>>ik>>jk>>k;
		AV[i]=k;
		IA[i]=jk-1;
		if (l!=(ik-1))
		{
			JA[++l]=i;
		}
	}
        for(i=0;i<n;i++)
        {
        buf[i]=0;
          for (int j=JA[i]; j<JA[i+1]; j++)
            buf[i]+=AV[j];
        }
cout<<"check2"<<endl;
	JA[n]=m;
	*dim=n;
	*NumEl=m;
	*AV1=AV;
	*IA1=IA;
	*JA1=JA;
	*buf1=buf;
	cout1.close();

}
void read3(double **AV1,int **IA1,int **JA1,double **buf1,int *dim,int *NumEl)
{
    double *AV, *V, *buf;
    int *IA, *NC, *JA, *NR;

    ifstream data, data1;
    int i,j,ik,jk,l,n,nc,m,m1;
    double k, *e;
    data.open("new.mtx");
    char s[255], c[255];
    data.getline(s,255);
        short flag=1;
    while ((flag)&&(!data.eof()))
    {
        data.getline(c,255);
        if (c[0]!='%')
            flag=0;
    }

        cout<<"check2"<<endl;
        int len = strlen(c);
        cout<<len<<"  "<<c<<endl;
        int size[3];
        j=0;
        i=0;
        while (i<len)
        {
            n=0;
            while( (c[i]!=' ') && (i<len) )
            {
                n=n*10 + c[i] - 0x30;
                i++;
            }
            size[j++]=n;
            i++;
        }

        n=size[0];
        nc=size[1];
        m1=size[2];
        m = 2*m1-n;

        cout<<"!!"<<n<<"   "<<nc<<"   "<<m1<<"   "<<m<<endl;
        
    V=new double [m1];
    NC=new int [m1];
    NR=new int [m1];
    AV=new double[m];
    IA=new int [m];
    JA=new int [n+1];
    buf=new double [n];
    for(i=0;i<n;i++)
    {
	JA[i]=0;
    }
    JA[n]=0;
        i=0;
    while(!data.eof())
    {
	data>>ik>>jk>>k;
	V[i]=k;
	NC[i]=jk-1;
	NR[i]=ik-1;
            JA[ik]++;
	if (ik!=jk)
	    JA[jk]++;
	i++;
    }

    cout<<"check2"<<endl;

    ik=0;jk=0;
    for(i=1;i<=n;i++)
    {
        ik=JA[i];
        JA[i]=jk;
        jk+=ik;
    }

        cout<<"check3"<<endl;

    for(i=0;i<m1;i++)
    {
        int col = NC[i];
        int row = JA[NR[i]+1];
        double val = V[i];
        AV[row]=val;
        IA[row]=col;
        JA[NR[i]+1]++;
        if(NR[i]!=col)
        {
            AV[JA[col+1]]=val;
            IA[JA[col+1]]=NR[i];
            JA[col+1]++;
        }
    }

    cout<<"check4"<<endl;

        for(i=0;i<n;i++)
        {
        buf[i]=0;
          for (j=JA[i]; j<JA[i+1]; j++)
            buf[i]+=AV[j];
        }
//        delete [] V;
//        delete [] NC;
//        delete [] NR;

    *dim=n;
    *NumEl=m;
    *AV1=AV;
    *IA1=IA;
    *JA1=JA;
    *buf1=buf;
    data.close();
}
void LUfact(double *AV, int *IA, int *JA,int *dim,double ***AV_11,int *NumEl)//, int **IA_11, int **JA_11,int str, double **mass, double **&LU,double **&A_1 )
{
	int n,m;
	n = *dim;
	double ttime;
	double **AV_1, **LU,**sum2, sum,**check, *AV_1_csr;
	m = *NumEl;
	sum2 = new double*[n];
	for(int i=0;i<n;i++)
	{sum2[i]=new double[n];}


	LU = new double*[n];
	for(int i=0;i<n;i++)
		LU[i]=new double[n];
	for(int i=0;i<n;i++)
	{
	    for(int j=0; j<n;j++)
	{
	    if (i==j)LU[i][j]=1;
	else
	LU[i][j]=0;
	}
	}
	for(int i=0;i<n;i++)
	{
	    for(int j=0; j<n;j++)
	{
	    sum2[i][j]=0;
	}
	}
	cout<<"check in LUfact 1 "<<n<<" "<<m<<endl;

timeval tim;
gettimeofday(&tim, NULL);
double t1=tim.tv_sec+(tim.tv_usec/1000000.0);

for (int i = 0; i < n; i++)
{
    for (int j = JA[i];j < JA[i+1]; j++) 
    {
	
	sum=0;
	if(i<=JA[i])// (i<=IA[a])
	{
	for(int k = 0 ; k < i; k++)
	{
//	sum2[i][IA[j]]+= LU[i][k]*LU[k][IA[j]];
	sum+= LU[i][k]*LU[k][IA[j]];
	if(i==0)LU[0][IA[j]] = AV[IA[j]];
	else	LU[i][IA[j]] = AV[IA[j]]-sum;//2[i][IA[j]];
	}
	}
	}
	for (int j = JA[i]+1;j < JA[i+1]; j++)
        {	if(i>JA[i])
	{
	for(int k = 0 ; k < i; k++)
	{
	sum2[i][k]+= LU[IA[j]][k]*LU[k][i];
	if(j==0) LU[IA[j]][i] = LU[i][IA[j]]/LU[j][j];
	else     LU[IA[j]][i] = (AV[j]-sum2[i][k])/LU[IA[j]][IA[j]];
	}
	}
    }
}

AV_1 = new double*[n];
for (int i = 0; i < n; i++)
AV_1[i]=new double [n];

	for (int i = n-1; i >= 0; i--)
	{
		for (int j = n-1; j >= 0; j--)
		{
			sum=0;
			if (i == j)
		{
				for(int k=j+1;k<n;k++)
				{
					sum+= LU[j][k]*AV_1[k][j];
					AV_1[j][j] = (1-sum)/LU[j][j];
				}
				
			}
			else
				if (i < j)
				{
					for (int k = i+1; k < n; k++)
					{
						sum+= LU[i][k]*AV_1[k][j];
						AV_1[i][j] =-sum/LU[i][i];
					}
				}
				else
				{
					for (int k = j+1; k <n; k++)
					{
						sum+= LU[k][j]*AV_1[i][k];
						AV_1[i][j] = -sum;
					}
				}
		}
	}
 gettimeofday(&tim, NULL);
 double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
 ttime=t2-t1;

cout<<"time = "<<ttime<<endl;

check = new double*[n];
for (int i = 0; i < n; i++)
check[i]=new double [n];



//first_martix = new double*[n];//pervonachalnaia matrica zapisannaia v obichnom formate
//for (int i = 0; i < n; i++)
//AV_1_csr[i]=new double [n];

	for(int i = 0; i < n; i ++)
	{
	    for (int j=JA[i];j <JA[i+1]; j++) 
	{
	    check[i][j]=0;
	    for(int k=0;k<n;k++)
	    check[i][j]+=AV_1[k][IA[j]]*AV[j];
}
	}

ofstream wr;
	wr.open("LU.txt");
	for(int i = 0; i < n; i ++)
	{
	    for (int j = 0; j < n; j++)
	    {
		wr<<check[i][j]<<" ";
	    }
	    wr<<endl;
	}
wr.close();

delete[]LU;
delete[] sum2;
	*AV_11 = AV_1;
}

void write(double **AV_1, int dim)
{
	ofstream wr;
	cout<<"check9"<<endl;
	wr.open("result.txt");
	for(int i = 0; i < dim; i ++)
	{
	    for (int j = 0; j < dim; j++)
	    {
		wr<<AV_1[i][j]<<" ";	//запись результата в файл
	    }
	    wr<<endl;
	}
	//	*AV_11 = *AV_1;
cout<<"check8"<<endl;
	wr.close();
cout<<"check20"<<endl;
}
int main()
{
	double *AV,*buf,**AV_1;
		int *IA,*JA, dim, NumEl;
		//read1(&AV, &IA, &JA, &buf, &dim, &NumEl);
		read3(&AV, &IA, &JA, &buf, &dim, &NumEl);
cout<<"read2"<<endl;
		LUfact(AV, IA, JA, &dim, &AV_1, &NumEl);
cout<<"read3"<<endl;
		write(AV_1, dim);
cout<<"read4"<<endl;
	delete [] IA;
	delete [] AV;
	delete [] JA;
	delete [] buf;
	return 0;
}