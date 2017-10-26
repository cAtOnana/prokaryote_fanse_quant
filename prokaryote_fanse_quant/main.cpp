#include<iostream>
#include<fstream>
#include<string>
#include<climits>
#include<algorithm>
#include<vector>
using namespace std;

const double million = 1000000;

struct reflat
{
	string gname;
	string gid;
	string gtype;
	int Start;
	int End;
	double length = 0;
	int readcount = 0;
	double rpkm = 0.0;
	char strand;
};
bool cmp_Start(const reflat &a, const reflat&b)
{
	return a.Start < b.Start;
}

void read_reflat(string fname, vector<reflat>& line);
void fill_rc(string fanse3, vector<reflat>& line);


int main(int agvc, char *agvr[])//��mian��������ɽ��������в�������ʽ��������������һ��Ϊreflat�ļ���ַ���ڶ���Ϊfanse3�ļ���ַ
{
	string arga =agvr[1];
	string argb = agvr[2];
	if (agvc != 3)
	{
		cout << "�����˶���Ĳ���";
		exit(EXIT_FAILURE);
	}

	vector<reflat> line;
	read_reflat(arga, line);
	fill_rc(argb, line);//��readcount
	int total_rc = 0;
	for (int i = 0; i < line.size(); i++)//����readcount
		total_rc += line[i].readcount;
	for (int i = 0; i < line.size(); i++)//��rpkm
		line[i].rpkm = line[i].readcount / (total_rc / million*line[i].length / 1000);
	//��ӡ�������
	sort(line.begin(), line.end(), cmp_Start);
	string outname1 = argb;
	outname1 = outname1.erase(outname1.size() - 7) + "-rpkm.txt";
	ofstream fout1(outname1);
	cout << "outputing the file..." << endl;
	fout1.precision(7);
	fout1 << "GeneName" <<"	" <<"Strand" <<"	"<< "Start" <<"	"<< "End"<<"	"<< "Length"
		<< "	" << "Read count" << "	" << "rpkM" << endl;
	for (int i = 0; i < line.size(); i++)
	{
		fout1 << line[i].gname << "	" << line[i].strand<< "	"<<line[i].Start<< "	"<<line[i].End<< "	"<<line[i].length << "	" << line[i].readcount << "	" << line[i].rpkm << endl;
	}
	fout1.close();
	cout << "Done!";

	return 0;
}

void read_reflat(string fname, vector<reflat>& line)
{
	ifstream fin1(fname);
	if (!fin1.is_open())
	{
		cout << "reFLAT��ȡʧ�ܣ������ж�";
		exit(EXIT_FAILURE);
	}
	cout << "reFLAT is in!" << endl;
	char s;
	for (int i = 0; i < 7; i++)
	{
		while (cin.get() == '\n')
			break;
	}
	string waste;
	reflat temp;
	while (fin1 >> waste) {//ѭ�����룬ֻ��ȡ������Ϊ gene �ֶε��е���Ϣ
		fin1 >> waste;
		fin1 >> waste;
		if (waste != "gene") {
			getline(fin1, waste);
			continue;
		}
		fin1 >> temp.Start;
		fin1 >> temp.End;
		fin1 >> waste;
		fin1 >> temp.strand;
		fin1 >> waste;
		fin1 >> waste;//�������һ�в����в�֣��õ�gname,gid,gtype����Ϣ
		int pos = waste.find("GeneID:");
		for (int i = pos + string("GeneID:").length(); waste[i] != ';'; i++)
			temp.gid += waste[i];
		pos = waste.find("Name=");
		for (int i = pos + string("Name=").length(); waste[i] != ';'; i++)
			temp.gname += waste[i];
		pos = waste.find("gene_biotype=");
		for (int i = pos + string("gene_biotype=").length(); waste[i] != ';'; i++)
			temp.gtype += waste[i];
		line.push_back(temp);
	}
}

void fill_rc(string fanse3, vector<reflat>& line)
{
	ifstream fin2(fanse3);
	if (!fin2.is_open())
	{
		cout << "fanse3��ȡʧ�ܣ������ж�";
		exit(EXIT_FAILURE);
	}
	cout << "fanse3 is in!" << endl;
	string waste;
	string id;

	char test;
	while (fin2.get(test))
	{

		getline(fin2, waste);
		if (waste[0] =='\0')
			break;
		fin2 >> waste;
		fin2 >> id;
		getline(fin2, waste);

		if (id == "")
			break;
		int end = line.size() - 1, front = 0, mid = (front + end) / 2;
		while (front < end&&line[mid].gname != id)//���ַ���nm��
		{
			if (line[mid].gname< id)
				front = mid + 1;
			else if (line[mid].gname> id)
				end = mid - 1;
			mid = (front + end) / 2;
		}
		if (line[mid].gname== id)
			line[mid].readcount += 1;
		else if
	}1
}


