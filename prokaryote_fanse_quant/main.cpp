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
	string gname = "";
	string gid = "";
	string gtype = "";
	int Start = 0;
	int End = 0;
	double length = 0;
	int readcount = 0;
	double rpkm = 0.0;
	char strand = '.';
};
bool cmp_Start(const reflat &a, const reflat&b)
{
	return a.Start < b.Start;
}

void read_reflat(string fname, vector<reflat>& line);
int fill_rc(string fanse3, vector<reflat>& line,ofstream& log);


int main()//int agvc, char *agvr[])//��mian��������ɽ��������в�������ʽ��������������һ��Ϊreflat�ļ���ַ���ڶ���Ϊfanse3�ļ���ַ
{
	string gff = "GCF_000005845.2_ASM584v2_genomic.gff";//agvr[1];
	string fanse3 = "1.fanse3";//agvr[2];
	/*if (agvc != 3)
	{
		cout << "�����˶���Ĳ���";
		exit(EXIT_FAILURE);
	}*/
	string logname = fanse3;
	logname = logname.erase(logname.find('.')) + ".log";
	ofstream log(logname);


	vector<reflat> line;
	read_reflat(gff, line);
	sort(line.begin(), line.end(), cmp_Start);
	int unquant=fill_rc(fanse3, line,log);//��readcount

	int total_rc = 0;
	for (int i = 0; i < line.size(); i++)//����readcount
		total_rc += line[i].readcount;

	log << "�� " << total_rc + unquant << " ��reads���� " << unquant << " ��reads�޷���������ռ��Ϊ��" << double(unquant) / (total_rc + unquant) * 100 << "%\n";

	for (int i = 0; i < line.size(); i++)//��rpkm
		line[i].rpkm = line[i].readcount / (total_rc / million*line[i].length / 1000);

	//��ӡ�������
	string outname1 = fanse3;
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

	log.close();
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

	for (int i = 0; i < 7; i++)
	{
		while (fin1.get() != '\n')
			continue;
	}
	string waste;
	reflat temp = {"","","",0,0,0,0,0,'.'};
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
		temp.length = temp.End - temp.Start;
		line.push_back(temp);
		temp = { "","","",0,0,0,0,0,'.' };
	}
}

int fill_rc(string fanse3, vector<reflat>& line,ofstream& log)
{
	ifstream fin2(fanse3);
	if (!fin2.is_open())
	{
		cout << "fanse3��ȡʧ�ܣ������ж�";
		exit(EXIT_FAILURE);
	}
	cout << "fanse3 is in!" << endl;

	string waste;
	int id;
	int unquant = 0;//�����޷�������reads��

	char test;
	while (fin2.get(test))
	{

		getline(fin2, waste);
		if (waste[0] =='\0')
			break;
		fin2 >> waste>>waste>>waste;
		fin2 >> id;
		getline(fin2, waste);

		int end = line.size() - 1, front = 0, mid = (front + end) / 2;
		while (end-front>1&&line[mid].Start != id)//���ַ���nm��
		{
			if (line[mid].Start< id)
				front = mid;
			else if (line[mid].Start> id)
				end = mid;
			mid = (front + end) / 2;
		}
		if (line[mid].Start <= id&&id<=line[mid].End)
			line[mid].readcount += 1;
		else
		{
			unquant++;
			//log << "��read�޷�������ע���ļ��У���λ��Ϊ��" << id <<"��������Ļ�����:"<<line[mid].gname << "����ʼλ��Ϊ��" << line[mid].Start << "��ֹλ��Ϊ��"<<line[mid].End << endl;
		}
	}
	return unquant;
}


