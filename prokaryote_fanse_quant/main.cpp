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


int main()//int agvc, char *agvr[])//将mian（）改造成接受命令行参数的形式，两个参数，第一个为reflat文件地址，第二个为fanse3文件地址
{
	string gff = "GCF_000005845.2_ASM584v2_genomic.gff";//agvr[1];
	string fanse3 = "1.fanse3";//agvr[2];
	/*if (agvc != 3)
	{
		cout << "输入了多余的参数";
		exit(EXIT_FAILURE);
	}*/
	string logname = fanse3;
	logname = logname.erase(logname.find('.')) + ".log";
	ofstream log(logname);


	vector<reflat> line;
	read_reflat(gff, line);
	sort(line.begin(), line.end(), cmp_Start);
	int unquant=fill_rc(fanse3, line,log);//算readcount

	int total_rc = 0;
	for (int i = 0; i < line.size(); i++)//算总readcount
		total_rc += line[i].readcount;

	log << "在 " << total_rc + unquant << " 条reads中有 " << unquant << " 条reads无法被定量，占比为：" << double(unquant) / (total_rc + unquant) * 100 << "%\n";

	for (int i = 0; i < line.size(); i++)//算rpkm
		line[i].rpkm = line[i].readcount / (total_rc / million*line[i].length / 1000);

	//打印定量结果
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
		cout << "reFLAT读取失败，程序中断";
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
	while (fin1 >> waste) {//循环读入，只留取第三行为 gene 字段的行的信息
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
		fin1 >> waste;//读入最后一列并进行拆分，得到gname,gid,gtype等信息
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
		cout << "fanse3读取失败，程序中断";
		exit(EXIT_FAILURE);
	}
	cout << "fanse3 is in!" << endl;

	string waste;
	int id;
	int unquant = 0;//计算无法定量的reads数

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
		while (end-front>1&&line[mid].Start != id)//二分法找nm号
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
			//log << "有read无法定量到注释文件中，其位点为：" << id <<"离它最近的基因是:"<<line[mid].gname << "其起始位点为：" << line[mid].Start << "终止位点为："<<line[mid].End << endl;
		}
	}
	return unquant;
}


