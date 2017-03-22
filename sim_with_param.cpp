//
//  main.cpp
//  Methylation_Server
//
//  Created by 任红雷 on 2016/12/29.
//  Copyright (c) 2016 ___renhonglei___. All rights reserved.
//

#include<iostream>
#include<cstdio>
#include<cmath>
#include<stdlib.h>
#include<sstream>
#include<string>
#include <cstdlib>
#include<string.h>
#include<vector>
#include<algorithm>
#include<functional>
#include<cstdlib>
#include<ctime>
#include<time.h>
#include<regex>
#include<numeric>
#include<functional>
#include<fstream>
#include<iomanip>
#include<map>
#include<random>

#ifdef WIN32
#include <io.h>
#include <direct.h>
#else

#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>

#endif


#define pb push_back

using namespace std;
const float EPSINON = 1e-6;

vector<int> nucleosome_start_pos_list;
vector<int> nucleosome_end_pos_list;
string nucleosome_status;

double beta = 0.195;
double gamma_param=0.04;
double alpha = 12.0;
double power=1.0;
double c_off=0.13;
double pi = 3.1415926535898;

int exist_or_make_dir(string dir_path)
{
    if (access(dir_path.c_str(), 0) == -1)
    {
        cout<<dir_path<<" is not existing"<<endl;
        cout<<"now make it"<<endl;
        int flag = mkdir(dir_path.c_str(),0777);
        if (flag == 0)
        {
            cout<<"make "<< dir_path.c_str() <<" successfully!"<<endl;
            return 0;
        } else {
            cout<<"make "<< dir_path.c_str() <<" errorly!"<<endl;
            return 1;
        }
    }
    return 0;
}

//CpG对应反应编号所需的原CpG状态
char right_status_of_reaction[] = {'U','H','M','H'};

//核小体对应反应编号所需的原核小体状态
char right_status_of_nuc_reaction[] = {'U','H','A','H'};

//对应反应编号，反应后的CpG状态
char right_status_hash[] = {'H','M','H','U'};

//对应反应编号，反应后的核小体状态
char right_nuc_status_hash[] = {'H','A','H','U'};

//string to_string(int snum)
//{
//    ostringstream stream;
//    stream<<snum;
//    return stream.str();
//}

//对0-range采sample_num个样本编号
vector<int> index_random(int sample_num, int range)
{
    vector<int> index,return_index;
    srand(time(0));
    for(int i=0; i<range; i++)
    {
        index.pb(0);
    }
    for(int i=1; i<=sample_num; i++)
    {
        int idx = rand()%range;
        if(index[idx]==0)
        {
            index[idx] = 1;
        }
        else
        {
            while(1)
            {
                idx = rand()%range;
                if(index[idx]==0)
                {
                    index[idx] = 1;
                    break;
                }
            }
        }
    }
    for(int i=0; i<range; i++)
    {
        if(index[i]==1)
        {
            return_index.pb(i);
        }
    }
    return return_index;
}

//构建一个最大CpG数量为max_cpg_sites，位置符合超几何分布G(p)的CpG位置向量
vector<int> construct_n_cpg_sites_for_exp_distribution(int max_cpg_sites, float p)
{
    vector<int> index_pos_list;
    int first_cpg_pos = rand() % 9 + 1;
    index_pos_list.pb(first_cpg_pos);

    default_random_engine generator(time(NULL));
    geometric_distribution<int> distribution(p);
    int number;
    for (int i=1; i< max_cpg_sites; ++i)
    {
        number= distribution(generator)+2;
        first_cpg_pos=first_cpg_pos+number;
        //        cout << i << " : " << first_cpg_pos <<endl;
        index_pos_list.pb(first_cpg_pos);
    }

    return index_pos_list;
}

//从输入input_bed_file_path的bed文件中，获取max_cpg_sites个CpG的位置列表
vector<int> get_pos_list_from_bed_file(string input_bed_file_path,int max_cpg_sites=1000)
{
    ifstream bed_file(input_bed_file_path);
    int  buf_size=100;
    int cnt_umh=0;
    int u_count=0;
    int m_count=0;
    int h_count=0;

    char buffer[buf_size];
    vector<int> pos_list;

    if (!bed_file)
    {
        cout << "Unable to open " << input_bed_file_path << endl;
    }
    else
    {
        int pos=0;
        float methy_level=0.0;
        int cnt=0;

        while (! bed_file.eof() )
        {
            if (max_cpg_sites < 0 or (cnt < max_cpg_sites and max_cpg_sites > 0))
            {
                bed_file.getline(buffer,buf_size);
                sscanf(buffer,"%d %f",&pos,&methy_level);
                if (methy_level<0.2)
                {
                    u_count=u_count+1;
                }
                else if(methy_level>0.8)
                {
                    m_count=m_count+1;
                }
                else
                {
                    h_count=h_count+1;
                }
                pos_list.pb(pos);
                cnt++;
            }
            else
            {
                break;
            }
        }
        cnt_umh=m_count+h_count+u_count;
        //        float u_ratio=float(u_count)/float(cnt_umh);
        //        float m_ratio=float(m_count)/float(cnt_umh);

        bed_file.close();
        printf("get pos list from %s successfully!\n",input_bed_file_path.c_str() );
    }
    return pos_list;
}

//输入CpG位点的数量,M,U状态所占的比例,输出一条符合UMH比例分布的string链，每个状态用'U','H','M'代表
string generate_CpG_in_methylation_percent_UHM(int CpG_sites_counts,double m_ratio, double u_ratio)
{
    string CpG_str;

    for (int i = 0; i < CpG_sites_counts; i++)
    {
        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;

        if (random_num < m_ratio)
        {
            CpG_str.append(1, 'M');
        }
        else if (random_num > m_ratio and random_num < m_ratio + u_ratio)
        {
            CpG_str.append(1, 'U');
        }
        else
        {
            CpG_str.append(1, 'H');
        }
    }

    return CpG_str;// 返回生成的状态字符串, 长度=CpG_sites_counts
}

//生成核小体的位置文件，输入: 起始核小体在染色体上的位置, 结束位置, 输出文件路径, 每个核小体长度, 每个gap的长度, 每种状态核小体的比例(分别为AA,AU,UU三种状态的浮点列表，它们的大小比值即为比例)
void generate_nucleosome_pos_list(int start_pos, int end_pos, string output_file_path, int constant_len, int constant_gap, double ratio_list[])
{
    double sum_ratio = 0.0;
    for(int i=0; i < 3; i++)
    {
        sum_ratio += ratio_list[i];
    }

    char nucleosome_status_list[] = {'A','H','U'};

    srand(time(0));

    ofstream outfile(output_file_path);
    int buffer_size = 50;
    if(!outfile)
    {
        cout << "Unable to open file: " << output_file_path << endl;
        exit(1);
    }

    int now_start = start_pos;
    int now_end = now_start;

    int nucleosome_count = 0;
    while (now_end + constant_len < end_pos)
    {
        nucleosome_count += 1;
        if (nucleosome_count % 10000 == 0)
        {
            printf("Now write %d nucleosomes\n", nucleosome_count);
        }

        now_end = now_start + constant_len;

        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
        double tmp_sum_ratio = 0.0;
        random_num = random_num * sum_ratio;
        char n_status='U';

        for(int i = 0; i< 3; i++)
        {
            tmp_sum_ratio = tmp_sum_ratio + ratio_list[i];
            if(random_num < tmp_sum_ratio)
            {
                n_status = nucleosome_status_list[i];
                break;
            }
        }

        char buffer[buffer_size];
        sprintf(buffer, "%d\t%d\t%c\n",now_start, now_end, n_status);
        outfile << buffer;

        now_start = now_end + constant_gap;
    }
    outfile.close();
    printf("Write successful\n");
}

//产生核小体位置文件
void prepare_nucleosome_pos_file()
{
    string base_dir = "//Mac/Home/Documents/Codeblocks/Nucleosome/output/";

    int constant_len = 160 ;

    int start_pos = 3001600;

    int end_pos = 126590304 + constant_len;

    string output_file_path = base_dir + "nucleosome_positions.np";

    int constant_gap = 5;

    double ratio_list[] = {0.5, 0.3, 0.2};

    generate_nucleosome_pos_list(start_pos, end_pos, output_file_path, constant_len, constant_gap, ratio_list);
}

vector<int> get_neucleosome_id_list_from_pos_list(vector<int> index_pos_list, string nucleosome_pos_file_path)
{
    vector<int> nucleo_id_list;
    ifstream nucleosome_pos_file(nucleosome_pos_file_path);

    if(!nucleosome_pos_file)
    {
        cout << "Unable to open read " << nucleosome_pos_file_path << endl;
        exit(1);
    }
    int buffer_size = 50;
    char buff[buffer_size];
    int now_nuceo_start = -1;
    int now_nuceo_end = -1;
    char status = 'N';

    int cpg_pos_index = 0;
    int nuceo_pos_index = 0;

    nucleosome_pos_file.getline(buff,buffer_size);
    sscanf(buff,"%d\t%d\t%c\n",&now_nuceo_start,&now_nuceo_end,&status);
    nucleosome_start_pos_list.push_back(now_nuceo_start);
    nucleosome_end_pos_list.push_back(now_nuceo_end);
    nucleosome_status.append(1,status);

    for (int i=0; i< index_pos_list.size(); i++)
    {
        int now_cpg_pos = index_pos_list[i];
        if ( now_cpg_pos <= now_nuceo_end)
        {
            if(now_cpg_pos >= now_nuceo_start)
            {
                //在组蛋白的范围内,则把该组蛋白的id加入到nucleo_id_list中
                nucleo_id_list.push_back(nuceo_pos_index);
            }
            else
            {
                //在gap中的CpG对应的核小体编号设置为-1
                nucleo_id_list.push_back(-1);
            }
        }
        else
        {
            while(now_cpg_pos > now_nuceo_end && ! nucleosome_pos_file.eof())
            {
                //当前CpG的位置>核小体的范围,继续读核小体的文件
                nucleosome_pos_file.getline(buff,buffer_size);
                sscanf(buff,"%d\t%d\t%c\n",&now_nuceo_start,&now_nuceo_end,&status);

                nucleosome_start_pos_list.push_back(now_nuceo_start);
                nucleosome_end_pos_list.push_back(now_nuceo_end);
                nucleosome_status.append(1,status);

                nuceo_pos_index += 1;
            }
            i--;
        }
    }
    nucleosome_pos_file.close();
    return nucleo_id_list;
}

void simulate(int round_no,int generation, int time_step, string init_cell,vector<int> neucleosome_id_list,vector<vector<int>> cpg_id_list, string detail_file_dir, string ratio_file_dir,string nuc_detail_dir, string nuc_ratio_dir,string nucleosome_pos_file_path, vector<double> propensity_list,vector<double> nucleosome_propensity_list, vector<double> nuceo_to_cpg_efficiency,vector<double> nearby_promote_efficiency,vector<double> k_range, int update_nuc_status_frequency, vector<int> index_pos_list, int max_cells,int out_start_gen,int out_end_gen,int round_start)
{
    double U_plus=propensity_list[0];
    double H_plus=propensity_list[1];
    double M_minus=propensity_list[2];
    double H_minus=propensity_list[3];
    double H_p_H=propensity_list[4];
    double H_p_M=propensity_list[5];
    double U_p_M=propensity_list[6];
    double H_m_U=propensity_list[7];
    double M_m_U=propensity_list[8];
    double pij = 1.0;

    random_device rd;
    mt19937 eng(rd());

    srand(time(0));
    //CpG链的细胞列表
    vector<string> cell_collection;
    //核小体的细胞列表，和CpG链的细胞列表保持同步
    vector<string> nucleo_collection;
    cell_collection.pb(init_cell);
    nucleo_collection.push_back(nucleosome_status);
    uniform_int_distribution<int> nuc_rand(0,nucleosome_status.length()-1);
    uniform_int_distribution<int> cpg_rand(1,init_cell.length()-2);


    //拷贝全局初始化的核小体状态以便于模拟的时候做更改

    for(int i = 1; i<= generation; i++)
    {
        //计时
        clock_t gen_start_time=clock();
        time_t raw_time_start;
        struct tm * start_time_info;
        time ( &raw_time_start );
        start_time_info = localtime ( &raw_time_start );
        printf ( "round %d, generataion %d start time is: %s", round_no, i,asctime (start_time_info) );

        //对本轮模拟的细胞采样
        vector<string> cells_wait_to_add;
        vector<string> cells_of_nuc_wait_to_add;
        unsigned long cc_size = cell_collection.size();
        //        unsigned long nc_size = nucleo_collection.size();

        if(cc_size > max_cells/2)
        {
            vector<string>::iterator it = cell_collection.begin();
            vector<string>::iterator it2 = nucleo_collection.begin();
            vector<int> index_vec = index_random(max_cells/2, max_cells);
            int pos=0;

            for(int j = 0 ; it != cell_collection.end() && it2 != nucleo_collection.end() ; j++)
            {
                if(( pos >= index_vec.size()) or (j != index_vec[pos]))
                {
                    cell_collection.erase(it);
                    nucleo_collection.erase(it2);
                }
                else
                {
                    pos++;
                    it++;
                    it2++;
                }
            }
        }

        vector<vector<int> > M_count_statistics(time_step);
        vector<vector<int> > H_count_statistics(time_step);
        vector<vector<int> > U_count_statistics(time_step);
        vector<vector<string> > out_detail_seq_arr(cell_collection.size());

        vector<vector<int> > AA_count_statistics(time_step);
        vector<vector<int> > AU_count_statistics(time_step);
        vector<vector<int> > UU_count_statistics(time_step);
        vector<vector<string> > out_nuc_detail_seq_arr(nucleo_collection.size());

        //迭代所有细胞
        for(int idx=0; idx < cell_collection.size(); idx++)
        {

            vector<string> vect_tmp;

            //迭代所有时间步
            for(int j=0 ; j< time_step ; j++)
            {

                int cell_len = (int)cell_collection[idx].length();
                int nuc_len = (int)nucleo_collection[idx].length();

                if ((j % update_nuc_status_frequency ==0) && (j != 0))
                {
                    for(int tk = 0; tk < nuc_len; tk++)
                    {
                        int k = nuc_rand(eng);  //0~nuc_len-1
                        double uu_promote_ratio = 0.0;
                        double au_promote_ratio = 0.0;
                        double aa_inhibit_ratio = 0.0;
                        double au_inhibit_ratio = 0.0;

                        //计算核小体内所有CpG的影响系数
                        if (! cpg_id_list[k].empty())
                        {
                            //含有CpG的核小体
                            //核小体内总共的CpG的个数
                            int tot = 0;
                            int m_or_h_cnt = 0;

                            int cpg_id_list_size = (int)cpg_id_list[k].size();

                            for(int nt = 0; nt < cpg_id_list_size; nt++)
                            {
                                tot++;
                                int cpg_id = cpg_id_list[k][nt];
                                char CpG_status_tmp = cell_collection[idx][cpg_id];
                                if ( CpG_status_tmp=='H' || CpG_status_tmp =='M')
                                {
                                    m_or_h_cnt++;
                                }
                            }
                            uu_promote_ratio = m_or_h_cnt/float(tot);
                            au_promote_ratio = uu_promote_ratio;
                            aa_inhibit_ratio = uu_promote_ratio;
                            au_inhibit_ratio = uu_promote_ratio;
                        }

                        double near_by_nuc_uu_promote_ratio = 0.0;
                        double near_by_nuc_au_promote_ratio = 0.0;
                        double near_by_nuc_aa_inhibit_ratio = 0.0;
                        double near_by_nuc_au_inhibit_ratio = 0.0;

                        //计算邻居核小体的影响系数
                        if (k == 0 || k == cpg_id_list.size()-1)
                        {
                            int target_k = 1;
                            if (k == cpg_id_list.size()-1)
                            {
                                target_k = (int)cpg_id_list.size()-2;
                            }
                            if (nucleo_collection[idx][target_k] == 'A' || nucleo_collection[idx][target_k] =='H')
                            {
                                near_by_nuc_uu_promote_ratio = 1.0;
                                near_by_nuc_au_promote_ratio = 1.0;
                                near_by_nuc_aa_inhibit_ratio = 1.0;
                                near_by_nuc_au_inhibit_ratio = 1.0;
                            }
                        }
                        else
                        {
                            char status_left = nucleo_collection[idx][k-1];
                            char status_right = nucleo_collection[idx][k+1];
                            if ( status_left == 'A' || status_left =='H' || status_right == 'A' || status_right =='H')
                            {
                                near_by_nuc_uu_promote_ratio = 1.0;
                                near_by_nuc_au_promote_ratio = 1.0;
                                near_by_nuc_aa_inhibit_ratio = 1.0;
                                near_by_nuc_au_inhibit_ratio = 1.0;
                            }
                        }

                        //核小体的四种反应速率
                        double uu_k_plus = nucleosome_propensity_list[0];
                        uu_k_plus = uu_k_plus + uu_promote_ratio * k_range[0] + near_by_nuc_uu_promote_ratio * nearby_promote_efficiency[0];

                        double au_k_plus = nucleosome_propensity_list[1];
                        au_k_plus = au_k_plus + au_promote_ratio * k_range[1] + near_by_nuc_au_promote_ratio * nearby_promote_efficiency[1];

                        double aa_k_minus=nucleosome_propensity_list[2];
                        aa_k_minus = aa_k_minus + aa_inhibit_ratio* k_range[2] + near_by_nuc_aa_inhibit_ratio * nearby_promote_efficiency[2];

                        double au_k_minus=nucleosome_propensity_list[3];
                        au_k_minus = au_k_minus + au_inhibit_ratio* k_range[3] + near_by_nuc_au_inhibit_ratio * nearby_promote_efficiency[3];

                        vector<double> propensity_tmp = {uu_k_plus,au_k_plus,aa_k_minus,au_k_minus};
                        int num_of_reactions = (int)propensity_tmp.size();
                        double sum_propensity = 0.0;

                        for(int l=0; l<num_of_reactions; l++)
                        {
                            sum_propensity += propensity_tmp[l];
                        }

                        double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
                        double tmp_sum_propencity = 0.0;
                        int reaction_id = -1;
                        random_num = random_num * sum_propensity;

                        //选反应
                        for(int l = 0; l < num_of_reactions; l++)
                        {
                            tmp_sum_propencity = tmp_sum_propencity + propensity_tmp[l];
                            if(random_num < tmp_sum_propencity)
                            {
                                reaction_id = l;
                                break;
                            }
                        }

                        char status_of_target_nuc_site = nucleo_collection[idx][k];

                        if(right_status_of_nuc_reaction[reaction_id] != status_of_target_nuc_site)
                        {
                            continue;
                        }

                        nucleo_collection[idx][k] = right_nuc_status_hash[reaction_id];
                    }
                }
                //迭代所有CpG位点
                //                float cnt_tmp=0;
                //                vector<int> reaction_list_cnt = {0,0,0,0};
                for(int k = 0; k < cell_len; k++)
                {
                    int target_reaction_CpG_site = cpg_rand(eng);  //1~cell_len-2
                    int col_CpG_site_index = target_reaction_CpG_site + ((rand()%2)?1:-1);
                    char status_col_site = cell_collection[idx][col_CpG_site_index];

                    char xj_status=status_col_site;

                    int pos_target = index_pos_list[target_reaction_CpG_site];
                    int col_site_pos = index_pos_list[col_CpG_site_index];
                    int distance = int(fabs(pos_target - col_site_pos));

                    if(distance > 1000){
                        continue;
                    }

                    double sum_propensity = 0.0;
                    int reaction_id=-1;
                    //double histone_effect = (alpha / (30.0 * sqrt(2 * pi))) * exp(-(distance - 160.0) * (distance - 160.0) / (1800.0));
                    double phi_d = (1.0 / pow((beta + gamma_param * distance),power)) + c_off;



                    //核小体促进u+的比率
                    double u_add = 0.0;
                    //核小体促进h+的比率
                    double h_add = 0.0;
                    //核小体减小m-的比率
                    double m_minus = 0.0;
                    //核小体减小h-的比率
                    double h_minus = 0.0;

                    int nucleo_id = neucleosome_id_list[target_reaction_CpG_site];

                    if (nucleo_id != -1)
                    {
                        char status_of_its_nucleosome = nucleo_collection[idx][nucleo_id];
                        if(status_of_its_nucleosome == 'H' or status_of_its_nucleosome == 'A')
                        {
                            u_add = nuceo_to_cpg_efficiency[0];
                            h_add = nuceo_to_cpg_efficiency[1];
                            m_minus = nuceo_to_cpg_efficiency[2];
                            h_minus = nuceo_to_cpg_efficiency[3];
                        }
                    }

                    double u_i_plus = U_plus + pij * int(xj_status == 'M') * phi_d * (U_p_M - U_plus);
                    double h_i_plus = H_plus + pij * int(xj_status == 'M') * phi_d * (H_p_M - H_plus) + pij * int(xj_status == 'H') * phi_d * (H_p_H - H_plus);
                    double m_i_minus = M_minus + pij * int(xj_status == 'U') * phi_d * (M_m_U - M_minus);
                    double h_i_minus = H_minus + pij * int(xj_status == 'U') * phi_d * (H_m_U - H_minus);

                    u_i_plus = u_i_plus + u_add;
                    h_i_plus = h_i_plus + h_add;
                    m_i_minus = max(m_i_minus + m_minus,0.0);
                    h_i_minus = max(h_i_minus + h_minus,0.0);

                    vector<double> propensity_tmp = {u_i_plus,h_i_plus,m_i_minus,h_i_minus};
                    int num_of_reactions = (int)propensity_tmp.size();

                    for(int l=0; l<num_of_reactions; l++)
                    {
                        sum_propensity += propensity_tmp[l];
                    }

                    double random_num = (double)(rand()%RAND_MAX)/RAND_MAX;
                    double tmp_sum_propencity = 0.0;
                    random_num = random_num * sum_propensity;
                    for(int l = 0; l < num_of_reactions; l++)
                    {
                        tmp_sum_propencity = tmp_sum_propencity + propensity_tmp[l];
                        if(random_num < tmp_sum_propencity)
                        {
                            reaction_id = l;
                            break;
                        }
                    }

                    char status_of_target_site = cell_collection[idx][target_reaction_CpG_site];

                    if(right_status_of_reaction[reaction_id] != status_of_target_site)
                    {
                        continue;
                    }

                    cell_collection[idx][target_reaction_CpG_site] = right_status_hash[reaction_id];
                    //                    cnt_tmp+=1;
                    //                    reaction_list_cnt[reaction_id]+=1;
                }
                //                printf("prabability %f , up: %.2f, hp: %.2f, mm: %.2f, hm: %.2f \n",cnt_tmp/cell_len,reaction_list_cnt[0]/cnt_tmp,reaction_list_cnt[1]/cnt_tmp,reaction_list_cnt[2]/cnt_tmp,reaction_list_cnt[3]/cnt_tmp);

                int m_count=0,h_count=0,u_count=0;
                for (int it=0; it < cell_len; it++)
                {
                    switch (cell_collection[idx][it])
                    {
                        case 'M':
                            m_count++;
                            break;
                        case 'H':
                            h_count++;
                            break;
                        case 'U':
                            u_count++;
                            break;
                        default:
                            break;
                    }
                }
                M_count_statistics[j].pb(m_count);
                H_count_statistics[j].pb(h_count);
                U_count_statistics[j].pb(u_count);
                //                float tot_t = m_count+h_count+u_count;
                //                printf("m: %.2f h:%.2f u:%.2f\n",m_count/tot_t,h_count/tot_t,u_count/tot_t);
                int aa_count=0, au_count=0, uu_count=0;
                for (int it=0; it < nuc_len; it++)
                {
                    switch (nucleo_collection[idx][it])
                    {
                        case 'A':
                            aa_count++;
                            break;
                        case 'H':
                            au_count++;
                            break;
                        case 'U':
                            uu_count++;
                            break;
                        default:
                            break;
                    }
                }
                AA_count_statistics[j].pb(aa_count);
                AU_count_statistics[j].pb(au_count);
                UU_count_statistics[j].pb(uu_count);

                out_detail_seq_arr[idx].pb(cell_collection[idx]);
                out_nuc_detail_seq_arr[idx].pb(nucleo_collection[idx]);
            }

            //处理细胞分裂
            string cell_1, cell_2;
            for(int j=0; j< cell_collection[idx].length(); j++)
            {
                char ch = cell_collection[idx][j];
                if(ch == 'M')
                {
                    cell_1.append(1, 'H');
                    cell_2.append(1, 'H');
                }
                else if(ch == 'U')
                {
                    cell_1.append(1, 'U');
                    cell_2.append(1, 'U');
                }
                else if(ch == 'H')
                {
                    cell_1.append(1, 'H');
                    cell_2.append(1, 'U');
                }
            }
            cells_wait_to_add.pb(cell_1);
            cells_wait_to_add.pb(cell_2);

            //处理核小体的细胞分裂
            string nuc_cell_1, nuc_cell_2;
            for(int k = 0; k < nucleo_collection[idx].length(); k++)
            {
                char ch = nucleo_collection[idx][k];
                if(ch == 'A')
                {
                    nuc_cell_1.append(1, 'H');
                    nuc_cell_2.append(1, 'H');
                }
                else if(ch == 'U')
                {
                    nuc_cell_1.append(1, 'U');
                    nuc_cell_2.append(1, 'U');
                }
                else if(ch == 'H')
                {
                    nuc_cell_1.append(1, 'H');
                    nuc_cell_2.append(1, 'U');
                }
            }
            cells_of_nuc_wait_to_add.pb(nuc_cell_1);
            cells_of_nuc_wait_to_add.pb(nuc_cell_2);
        }
        vector<double> m_means_ratio;
        vector<double> h_means_ratio;
        vector<double> u_means_ratio;

        vector<double> aa_means_ratio;
        vector<double> au_means_ratio;
        vector<double> uu_means_ratio;

        for (int j = 0; j < time_step; j++)
        {

            double m_sum = accumulate(M_count_statistics[j].begin(), M_count_statistics[j].end(), 0.0);
            double m_mean =  m_sum / init_cell.length(); //按照M_count_statistics的行求均值，即每个cell_collection内所有细胞在特定time_step的均值
            m_means_ratio.pb(m_mean);

            double h_sum = accumulate(H_count_statistics[j].begin(), H_count_statistics[j].end(), 0.0);
            double h_mean =  h_sum / init_cell.length();
            h_means_ratio.pb(h_mean);

            double u_sum = accumulate(U_count_statistics[j].begin(), U_count_statistics[j].end(), 0.0);
            double u_mean =  u_sum / init_cell.length();
            u_means_ratio.pb(u_mean);

            double aa_sum = accumulate(AA_count_statistics[j].begin(), AA_count_statistics[j].end(), 0.0);
            double aa_mean = aa_sum / nucleosome_status.length();
            aa_means_ratio.pb(aa_mean);

            double au_sum = accumulate(AU_count_statistics[j].begin(), AU_count_statistics[j].end(), 0.0);
            double au_mean =  au_sum / nucleosome_status.length();
            au_means_ratio.pb(au_mean);

            double uu_sum = accumulate(UU_count_statistics[j].begin(), UU_count_statistics[j].end(), 0.0);
            double uu_mean =  uu_sum / nucleosome_status.length();
            uu_means_ratio.pb(uu_mean);
        }


        int buf_ratio_size=50;

        if ((i >= out_start_gen) and (i<= out_end_gen))
        {
            for(int t = 0; t < time_step; t++)
            {
                string ratio_file_path = ratio_file_dir+to_string(i)+"_"+to_string(t)+".csv";
                ofstream out_ratio(ratio_file_path,(round_no == round_start)? ios::trunc : ios::app);
                if(!out_ratio)
                {
                    cout << "Unable To Open Write " << ratio_file_path << endl;
                    exit(1);
                }
                else
                {
                    char buffer[buf_ratio_size];
                    sprintf(buffer, "%d,%.4f,%.4f,%.4f\n",round_no, m_means_ratio[t],h_means_ratio[t],u_means_ratio[t]);
                    out_ratio<<buffer;
                    out_ratio.close();
                }

                string nuc_ratio_file_path = nuc_ratio_dir+to_string(i)+"_"+to_string(t)+".csv";
                ofstream out_nuc_ratio(nuc_ratio_file_path,(round_no == round_start)? ios::trunc : ios::app);
                if(!out_nuc_ratio)
                {
                    cout << "Unable To Open Write " << nuc_ratio_file_path << endl;
                    exit(1);
                }
                else
                {
                    char buffer[buf_ratio_size];
                    sprintf(buffer, "%d,%.4f,%.4f,%.4f\n",round_no, aa_means_ratio[t],au_means_ratio[t],uu_means_ratio[t]);
                    out_nuc_ratio<<buffer;
                    out_nuc_ratio.close();
                }
            }
            int len_init_cell = (int)init_cell.length();
            int buf_detail_size = len_init_cell + 50;
            int buf_nuc_detail_size = ((int) nucleosome_status.length()) + 50;
            for(int idx = 0; idx < cell_collection.size(); idx++)
            {
                for(int j= 0; j < time_step; j++)
                {
                    string detail_file_path=detail_file_dir+to_string(i)+"_"+to_string(j)+".csv";
                    ofstream out_detail(detail_file_path,(round_no == round_start)?ios::trunc:ios::app);

                    if(!out_detail)
                    {
                        cout << "Unable To Open Write " << detail_file_path << endl;
                        exit(1);
                    }
                    else
                    {
                        char buffer[buf_detail_size];//在index的size不等于1时,得到的结果是多细胞的
                        sprintf(buffer, "%d,%s\n", round_no, out_detail_seq_arr[idx][j].c_str());
                        out_detail<<buffer;
                        out_detail.close();
                    }

                    string nuc_detail_file_path=nuc_detail_dir + to_string(i) + "_" + to_string(j) + ".csv";
                    ofstream out_nuc_detail(nuc_detail_file_path,(round_no == round_start)?ios::trunc:ios::app);

                    if(!out_nuc_detail)
                    {
                        cout << "Unable To Open Write " << nuc_detail_file_path << endl;
                        exit(1);
                    }
                    else
                    {
                        char buffer[buf_nuc_detail_size];//在index的size不等于1时,得到的结果是多细胞的
                        sprintf(buffer, "%d,%s\n", round_no, out_nuc_detail_seq_arr[idx][j].c_str());
                        out_nuc_detail << buffer;
                        out_nuc_detail.close();
                    }

                }
            }
        }

        cell_collection = cells_wait_to_add;
        nucleo_collection = cells_of_nuc_wait_to_add;

        clock_t gen_end_time=clock();
        cout<< "generation time:"<<static_cast<double>(gen_end_time-gen_start_time)/CLOCKS_PER_SEC<<"s"<<endl;//输出1代运行时间

    }
}

vector<vector<int>> get_corresponding_cpg_id_list(int cpg_id_list_size,vector<int> nucleosome_id_list)
{
    vector<vector<int>> cpg_id_list(cpg_id_list_size);
    for (int i = 0; i < nucleosome_id_list.size(); i++)
    {
        if(nucleosome_id_list[i] != -1)
        {
            cpg_id_list[nucleosome_id_list[i]].push_back(i);
        }
    }
    return cpg_id_list;
}

void sort_to_bed(int round_size,vector<int> pos_list,string sort_detail_dir,string bed_files_dir,int generation_start,int generation_end,int time_steps,int max_cpg_sites)
{
    int  buf_size=max_cpg_sites+50;
    char buffer[buf_size];

    for(int gen=generation_start; gen<=generation_end; gen++)
    {
        for (int t=0; t<time_steps; t++)
        {
            int round=0;
            unsigned long len_methy_seq=0;
            vector<int> methy_status;

            string detail_file_path=sort_detail_dir+to_string(gen)+"_"+to_string(t)+".csv";
            ifstream detail_file(detail_file_path);

            if(!detail_file)
            {
                cout << "Unable to open read " << detail_file_path << endl;
                exit(1);
            }
            else
            {
                string bed_file_path=bed_files_dir+to_string(gen)+"_"+to_string(t)+".bed";
                ofstream bed_file(bed_file_path);

                char methy_seq[buf_size];
                int line_cnt=0;
                while (! detail_file.eof() && line_cnt< round_size)
                {
                    line_cnt=line_cnt+1;
                    detail_file.getline(buffer,buf_size);
                    sscanf(buffer,"%d,%s\n",&round,methy_seq);
                    len_methy_seq=strlen(methy_seq);
                    for(unsigned long i=0; i<len_methy_seq; i++)
                    {
                        char ch=methy_seq[i];
                        if (ch=='M')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(2);
                            }
                            else
                            {
                                methy_status[i]=methy_status[i]+2;
                            }
                        }
                        else if (ch=='H')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(1);
                            }
                            else
                            {
                                methy_status[i]=methy_status[i]+1;
                            }
                        }
                        else if (ch=='U')
                        {
                            if(line_cnt==1)
                            {
                                methy_status.push_back(0);
                            }
                        }
                    }
                }
                //                vector<float> methy_mean;
                char wrt_buffer[100];
                float round_cnt=2*(float)round;
                for(unsigned long i=0; i<len_methy_seq; i++)
                {
                    //                    methy_mean.pb(methy_status[i]/round_cnt);
                    sprintf(wrt_buffer,"%d %.6f\n",pos_list[i],methy_status[i]/round_cnt);
                    bed_file<<wrt_buffer;
                }
                detail_file.close();
                bed_file.close();
            }
            printf( "finished convert %d.%d to bed\n",gen,t);
        }

    }
}
map<int,float> read_bed_file_and_store_pos_to_a_struct(string bed_file_path,bool ignore_d)
{
    map<int,float> struct_to_store;
    char buffer[100];
    int pos;
    float methy_level;
    int index=0;
    ifstream bed_file(bed_file_path);
    if(!bed_file)
    {
        cout << "Unable to open read " << bed_file_path << endl;
        exit(1);
    }
    else
    {
        while (!bed_file.eof() )
        {
            index=index+1;
            bed_file.getline(buffer,100);
            sscanf(buffer,"%d %f",&pos,&methy_level);
            if (ignore_d)
            {
                pos=index;
            }
            struct_to_store[pos]=methy_level;
        }
        bed_file.close();
    }
    return struct_to_store;
}

vector<vector<float>> filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(int d,vector<int> keys,vector<float> vals)
{
    int key_size=(int)keys.size();
    int pre_key,post_key;
    float pre_val,post_val;
    vector<vector<float>> array_to_store_pairs(2);

    for (int i=0; i< key_size; i++)
    {
        pre_key = keys[i];
        post_key = keys[i+1];

        if (pre_key + d == post_key)
        {
            pre_val=vals[i];
            post_val=vals[i+1];
            //            printf("<%d,%d> : %.2f, %.2f\n",pre_key,post_key,pre_val,post_val);
            array_to_store_pairs[0].push_back(pre_val);
            array_to_store_pairs[1].push_back(post_val);
        }
    }

    return array_to_store_pairs;
}
vector<vector<float>> filter_d_length_to_generate_CpG_pairs(int d,map<int,float> CpG_pos_and_methy_struct)
{
    int map_size=(int)CpG_pos_and_methy_struct.size();
    int pre_key,post_key;
    float pre_val,post_val;
    vector<vector<float>> array_to_store_pairs(2);
    for(map<int,float>::iterator it = CpG_pos_and_methy_struct.begin(); it != CpG_pos_and_methy_struct.end(); ++it)
    {
        pre_key=it->first;
        pre_val=it->second;
        post_key=pre_key+d;
        if (CpG_pos_and_methy_struct.find(post_key) != CpG_pos_and_methy_struct.end())
        {
            post_val=CpG_pos_and_methy_struct[post_key];
            array_to_store_pairs[0].push_back(pre_val);
            array_to_store_pairs[1].push_back(post_val);
        }
    }
    return array_to_store_pairs;
}

float calc_C_d_by_pearson_correlation(vector<vector<float>> CpG_pairs)
{
    float sum_pre=0.0;
    float sum_post=0.0;

    int length = (int)CpG_pairs[0].size();
    for (int i=0; i<length; i++)
    {
        sum_pre=sum_pre+CpG_pairs[0][i];
        sum_post=sum_post+CpG_pairs[1][i];
    }

    float mean1=sum_pre/float(length);
    float mean2=sum_post/float(length);

    float sum_up=0.0;
    float sum_down_left=0.0;
    float sum_down_right=0.0;

    float xi,yi;
    for (int i=0; i<length; i++)
    {
        xi=CpG_pairs[0][i];
        yi=CpG_pairs[1][i];

        sum_up=sum_up+(xi-mean1)*(yi-mean2);
        sum_down_left=sum_down_left+(xi-mean1)*(xi-mean1);
        sum_down_right=sum_down_right+(yi-mean2)*(yi-mean2);
    }

    float sum_down=sqrt(sum_down_left*sum_down_right);

    if((sum_down >= - EPSINON) && (sum_down <= EPSINON))
    {
        return -5;
    }
    float rd=sum_up/sum_down;
    return rd;
}

//读bed文件，计算距离相关性
void calc_correlation(string bed_file_path,string rd_file_path,int d_max,bool is_inter_with_other_cpg,bool ignore_d=false)
{
    ofstream rd_file(rd_file_path);

    if (!rd_file)
    {
        cout << "Unable to open read " << rd_file_path << endl;
        exit(1);
    }

    map<int,float> CpG_pos_and_methy_struct = read_bed_file_and_store_pos_to_a_struct(bed_file_path, ignore_d);

    vector<int> keys;
    vector<float> vals;
    for(map<int,float>::iterator it = CpG_pos_and_methy_struct.begin(); it != CpG_pos_and_methy_struct.end(); ++it)
    {
        keys.push_back(it->first);
        vals.push_back(it->second);
        //        cout << it->first<<" "<< it->second << endl;
    }

    vector<vector<float>> CpG_pairs;
    int d_count=0;
    float rd=0.0;
    char ltw[100];
    for (int d=2; d<d_max; d++)
    {
        if (is_inter_with_other_cpg)
        {
            CpG_pairs=filter_d_length_to_generate_CpG_pairs(d,CpG_pos_and_methy_struct);
        }
        else
        {
            CpG_pairs=filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(d,keys,vals);
        }
        d_count = (int)CpG_pairs.size();

        if (d_count)
        {
            rd=calc_C_d_by_pearson_correlation(CpG_pairs);
            if (rd > -2)
            {
                //printf("chr%d d=%d, rd=%f\n",1,d,rd);
                sprintf(ltw,"%d,%f\n",d,rd);
                rd_file<<ltw;
            }
        }
        else
        {
            printf("chr%d passed d=%d",1,d);
        }
    }
    rd_file.close();

    printf("\n\n");
}

void calc_mean_rd_from_rd_dir(string rd_dir_name,string out_file_path,int d_max,int generation_start,int generation_end,int time_steps)
{

    string rd_file_path;
    ofstream out_file(out_file_path);
    if (!out_file)
    {
        cout << "Unable to open out " << out_file_path << endl;
        exit(1);
    }
    vector<vector<float>> rd_vect;
    for (int i=0; i<d_max; i++)
    {
        vector<float> tmp_vect;
        rd_vect.push_back(tmp_vect);
    }

    int d=-1;
    float methy_level=0.0;
    char buffer[50];
    float m_sum,m_mean;
    for(int gen=generation_start; gen <= generation_end; gen++)
    {
        for(int step = 0; step < time_steps; step++)
        {
            rd_file_path=rd_dir_name+to_string(gen)+"_"+to_string(step)+".csv";
            ifstream rd_file(rd_file_path);
            if (!rd_file)
            {
                cout << "Unable to open read " << rd_file_path << endl;
                exit(1);
            }
            while (!rd_file.eof() )
            {
                rd_file.getline(buffer,50);
                sscanf(buffer,"%d,%f",&d,&methy_level);
                rd_vect[d].pb(methy_level);
            }
            rd_file.close();
        }
    }
    vector<float> rd_list;
    char ltw[50];
    for (d=2; d<d_max; d++)
    {
        if (!rd_vect[d].size())
            continue;
        m_sum = accumulate(rd_vect[d].begin(), rd_vect[d].end(), 0.0);
        m_mean =  m_sum / (float) rd_vect[d].size(); //按照M_count_statistics的行求均值，即每个cell_collection内所有细胞在特定time_step的均值
        sprintf(ltw,"%d,%f\n",d,m_mean);
        out_file<<ltw;
    }
    out_file.close();
    printf("writing mean rd finished!\n");
}

void calc_correlation_for_generations(int generation_start,int generation_end,int time_steps,string bed_dir,string rd_without_dir,int d_max,bool calc_interval=false,bool ignore_d=false)
{
    for (int gen=generation_start; gen <= generation_end; gen++)
    {
        for(int t=0; t< time_steps; t++)
        {
            string input_bed_file_path=bed_dir+to_string(gen)+"_"+to_string(t)+".bed";
            string out_rd_file_path=rd_without_dir+to_string(gen)+"_"+to_string(t)+".csv";
            calc_correlation(input_bed_file_path, out_rd_file_path, d_max, calc_interval,ignore_d);
        }
    }

}

void start_simulation()
{
    //prepare_nucleosome_pos_file();
    int round_start = 1;
    int round_end = 4;

    int generations = 10;
    int max_cpg_sites = 100000;
    // max_cpg_sites = 224481;
    string init_cell;
    double m_ratio = 0.181214;
    double u_ratio = 0.391004;

    int out_start_gen = generations;
    int out_end_gen = generations;

    int update_nuc_status_frequency = 1;

    //CpG的趋向性函数,即9种反应各自的比例,顺序:u+,h+,m-,h-,h+h,h+m,u+m,h-u,m-u
    vector<double> propensity_list = {0.008,0.008, 0.04, 0.04,0.24,0.24,0.24,0.05,0.05};

    //核小体的趋向性函数,4种反应各自的比例
    vector<double> nucleosome_propensity_list = {0.04,0.04,0.2,0.2};

    //相邻核小体对当前核小体的促进效率,按照uu+,au+,aa-,au-的顺序
    vector<double> nearby_promote_efficiency = {0.2,0.2,0.0,0.0};

    //UU+的Kmax-Kmin,和AU+的Kmax-Kmin
    vector<double> k_range = {1.0,1.0,0.0,0.0};

    //核小体对u+以及m-两个反应速率的促进程度，按照其u+,m-原始速率的倍数计算
    //vector<double> nuceo_to_cpg_efficiency = {0.032,0.032, -0.0016,-0.0016};
    vector<double> nuceo_to_cpg_efficiency = {0.0,0.0, -0.00,-0.00};

    vector<int> index_pos_list;

    //当前文件夹的全路径
    string now_dir="/Users/Ren/XCodeProjects/Cpp_Test/Cpp_Test/";

    string output_dir_name = "output/";
    string input_dir_name = "input/";
    string output_dir = now_dir + output_dir_name;
    string input_bed_file_path = now_dir + input_dir_name + "chr1.bed";
    string nucleosome_pos_file_path = output_dir + "nucleosome_positions.np";
    string ratio_file_dir = output_dir+"ratio/";//CpG ratio文件夹
    string nuc_ratio_dir = output_dir + "nuc_ratio/"; //核小体 ratio 文件夹
    string detail_file_dir = output_dir+"detail/";//CpG detail文件夹
    string nuc_detail_dir = output_dir + "nuc_detail/"; //核小体 detail 文件夹
    string bed_file_dir;
    string rd_without_dir;

    int time_step=100;
    int corr_d_max=1000;
    bool calc_interval=false;
    bool ignore_d=false;

    int max_cells = 2;

    bool simulation=true;
    bool calc_corr=true;

    bool real_chr_pos=true;

    if(!real_chr_pos)
    {
        //超几何分布参数
        float geometric_p = 0.3;
        index_pos_list = construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p);
    }
    else
    {
        index_pos_list = get_pos_list_from_bed_file(input_bed_file_path,max_cpg_sites);
    }


    //根据CpG在数组index_pos_list中的下标,查找对应的核小体下标
    vector<int> neucleosome_id_list = get_neucleosome_id_list_from_pos_list(index_pos_list,nucleosome_pos_file_path);

    int nucleosome_list_size = (int) nucleosome_end_pos_list.size();
    vector<vector<int>> cpg_id_list = get_corresponding_cpg_id_list(nucleosome_list_size,neucleosome_id_list);

    init_cell = generate_CpG_in_methylation_percent_UHM(max_cpg_sites,m_ratio,u_ratio);

    if (simulation)
    {
        for(int round_i=round_start; round_i<=round_end; round_i++)
        {
            simulate(round_i, generations,time_step,init_cell,neucleosome_id_list,cpg_id_list,detail_file_dir,ratio_file_dir,nuc_detail_dir,nuc_ratio_dir,nucleosome_pos_file_path,propensity_list,nucleosome_propensity_list,nuceo_to_cpg_efficiency,nearby_promote_efficiency,k_range,update_nuc_status_frequency,index_pos_list,max_cells,out_start_gen,out_end_gen,round_start);
        }
    }
    if(calc_corr)
    {
        bed_file_dir = output_dir + "bed/";//detail的生成序列

        rd_without_dir=output_dir+"rd_without/";

        int round_size=round_end-round_start+1;
        sort_to_bed(round_size,index_pos_list,detail_file_dir,bed_file_dir,out_start_gen,out_end_gen,time_step,max_cpg_sites);
        calc_correlation_for_generations(out_start_gen,out_end_gen,time_step,bed_file_dir,rd_without_dir,corr_d_max,calc_interval,ignore_d=false);
        string out_mean_rd_file = output_dir+"rd_mean.csv";
        calc_mean_rd_from_rd_dir(rd_without_dir,out_mean_rd_file,corr_d_max,out_start_gen,out_end_gen,time_step);
    }
}

void preprocess_file(string input_file_path,string output_file_path)
{
    ifstream csv_file(input_file_path);
    ofstream outfile(output_file_path);
    if(!outfile)
    {
        cout << "Unable to open file: " << output_file_path << endl;
        exit(1);
    }

    if (!csv_file)
    {
        cout << "Unable to open " << input_file_path << endl;
        exit(1);
    }
    int  buf_size=100100;
    char buffer[buf_size];
    char methy_seq[buf_size];

    csv_file.getline(buffer,buf_size);
    int round=-1;
    sscanf(buffer,"%d,%s\n",&round,methy_seq);
    csv_file.close();
    int len_methy_seq=(int)strlen(methy_seq);

    string ret_str = "";

    for(unsigned long i=0; i< len_methy_seq; i++)
    {
        char ch=methy_seq[i];
        string now_ch_to_ap = "";
        if (ch=='M')
        {
             now_ch_to_ap= "2";
        }
        else if (ch=='H')
        {
            now_ch_to_ap  = "1";
        }
        else if (ch=='U')
        {
            now_ch_to_ap  = "0";
        }
        if (i != len_methy_seq - 1) {
            ret_str += now_ch_to_ap;
            ret_str += ' ';
        }
        else
        {
            ret_str += now_ch_to_ap;
        }
    }
    int buf_size_wt = 2 * buf_size;
    char buffer_wt[buf_size_wt];//在index的size不等于1时,得到的结果是多细胞的
    sprintf(buffer_wt, "%s\n", ret_str.c_str());
    outfile<<buffer_wt;
    outfile.close();
    printf("Write Success!");

}

//将参数文件里面的所有参数名和参数值读取并存入hashmap，返回得到的hashmap

map<string,string> read_param_file_into_map(string param_file_path)
{
    map<string,string> k_v_hash;
    ifstream param_file(param_file_path);
    int  buf_size=100;

    char buffer[buf_size];

    if (!param_file)
    {
        cout << "param file path: " << param_file_path << "does not exist!"<< endl;
        exit(0);
    }

    char param_name[buf_size];
    char param_val[buf_size];
    while (! param_file.eof() )
    {
        param_file.getline(buffer,buf_size);
        sscanf(buffer,"%s %s",param_name,param_val);
//        printf("%s : %s \n",param_name,param_val);
        k_v_hash[param_name] = param_val;
    }
    return k_v_hash;
}

void start_simualation_from(vector<string> param_file_path_list)
{
    //遍历参数文件路径列表
    for(int i = 0; i < param_file_path_list.size(); i++)
    {
        //当前参数文件路径
        string param_file_path = param_file_path_list[i];
        //获取参数和参数值
        map<string,string> param_map = read_param_file_into_map(param_file_path);

        int round_start = std::atoi(param_map["round_start"].c_str());
        printf("round_start : %d \n",round_start);

        int round_end = std::atoi(param_map["round_end"].c_str());
        printf("round_end : %d \n",round_end);

        int generations = std::atoi(param_map["generations"].c_str());
        printf("generations : %d \n",generations);

        int max_cpg_sites = std::atoi(param_map["max_cpg_sites"].c_str());
        printf("max_cpg_sites : %d \n",max_cpg_sites);

        string init_cell;
        double m_ratio = std::atof(param_map["m_ratio"].c_str());
        printf("m_ratio : %.6f \n",m_ratio);

        double u_ratio = std::atof(param_map["u_ratio"].c_str());
        printf("u_ratio : %.6f \n",u_ratio);

        int out_start_gen = std::atoi(param_map["out_start_gen"].c_str());
        printf("out_start_gen : %d \n",out_start_gen);

        int out_end_gen = std::atoi(param_map["out_end_gen"].c_str());
        printf("out_end_gen : %d \n",out_end_gen);

        int update_nuc_status_frequency = std::atoi(param_map["update_nuc_status_frequency"].c_str());
        printf("update_nuc_status_frequency : %d \n",update_nuc_status_frequency);

        double U_plus = std::atof(param_map["U_plus"].c_str());
        printf("U_plus : %.6f \n",U_plus);

        double H_plus = std::atof(param_map["H_plus"].c_str());
        printf("H_plus : %.6f \n",H_plus);

        double M_minus = std::atof(param_map["M_minus"].c_str());
        printf("M_minus : %.6f \n",M_minus);

        double H_minus = std::atof(param_map["H_minus"].c_str());
        printf("H_minus : %.6f \n",H_minus);

        double H_p_H = std::atof(param_map["H_p_H"].c_str());
        printf("H_p_H : %.6f \n",H_p_H);

        double H_p_M = std::atof(param_map["H_p_M"].c_str());
        printf("H_p_M : %.6f \n",H_p_M);

        double U_p_M = std::atof(param_map["U_p_M"].c_str());
        printf("U_p_M : %.6f \n",U_p_M);

        double H_m_U = std::atof(param_map["H_m_U"].c_str());
        printf("H_m_U : %.6f \n",H_m_U);

        double M_m_U = std::atof(param_map["M_m_U"].c_str());
        printf("M_m_U : %.6f \n",M_m_U);
        
        
        beta = std::atof(param_map["beta"].c_str());
        printf("beta : %.6f \n",beta);
        
        gamma_param = std::atof(param_map["gamma_param"].c_str());
        printf("gamma_param : %.6f \n",gamma_param);
        
        alpha = std::atof(param_map["alpha"].c_str());
        printf("alpha : %.6f \n",alpha);
        
        power = std::atof(param_map["power"].c_str());
        printf("power : %.6f \n",power);
        
        c_off = std::atof(param_map["c_off"].c_str());
        printf("c_off : %.6f \n",c_off);
        
        

        //CpG的趋向性函数,即9种反应各自的比例,顺序:u+,h+,m-,h-,h+h,h+m,u+m,h-u,m-u
        vector<double> propensity_list = {U_plus, H_plus, M_minus, H_minus,H_p_H, H_p_M,U_p_M,H_m_U,M_m_U};

        double uu_k_plus = std::atof(param_map["uu_k_plus"].c_str());
        printf("uu_k_plus : %.6f \n",uu_k_plus);

        double au_k_plus = std::atof(param_map["au_k_plus"].c_str());
        printf("au_k_plus : %.6f \n",au_k_plus);

        double aa_k_minus = std::atof(param_map["aa_k_minus"].c_str());
        printf("aa_k_minus : %.6f \n",aa_k_minus);

        double au_k_minus = std::atof(param_map["au_k_minus"].c_str());
        printf("au_k_minus : %.6f \n",au_k_minus);

        //核小体的趋向性函数,4种反应各自的比例
        vector<double> nucleosome_propensity_list = {uu_k_plus, au_k_plus, aa_k_minus,au_k_minus};

        double nearby_nucleosome_to_uu_plus = std::atof(param_map["nearby_nucleosome_to_uu_plus"].c_str());
        printf("nearby_nucleosome_to_uu_plus : %.6f \n",nearby_nucleosome_to_uu_plus);

        double nearby_nucleosome_to_au_plus = std::atof(param_map["nearby_nucleosome_to_au_plus"].c_str());
        printf("nearby_nucleosome_to_au_plus : %.6f \n",nearby_nucleosome_to_au_plus);

        double nearby_nucleosome_to_aa_minus = std::atof(param_map["nearby_nucleosome_to_aa_minus"].c_str());
        printf("nearby_nucleosome_to_aa_minus : %.6f \n",nearby_nucleosome_to_aa_minus);

        double nearby_nucleosome_to_au_minus = std::atof(param_map["nearby_nucleosome_to_au_minus"].c_str());
        printf("nearby_nucleosome_to_au_minus : %.6f \n",nearby_nucleosome_to_au_minus);


        //相邻核小体对当前核小体的促进效率,按照uu+,au+,aa-,au-的顺序
        vector<double> nearby_promote_efficiency = {nearby_nucleosome_to_uu_plus,nearby_nucleosome_to_au_plus,nearby_nucleosome_to_aa_minus,nearby_nucleosome_to_au_minus};

        double uu_k_plus_k_range = std::atof(param_map["uu_k_plus_k_range"].c_str());
        printf("uu_k_plus_k_range : %.6f \n",uu_k_plus_k_range);

        double au_k_plus_k_range = std::atof(param_map["au_k_plus_k_range"].c_str());
        printf("au_k_plus_k_range : %.6f \n",au_k_plus_k_range);

        double aa_k_minus_k_range = std::atof(param_map["aa_k_minus_k_range"].c_str());
        printf("aa_k_minus_k_range : %.6f \n",aa_k_minus_k_range);

        double au_k_minus_k_range = std::atof(param_map["au_k_minus_k_range"].c_str());
        printf("au_k_minus_k_range : %.6f \n",au_k_minus_k_range);

        //UU+的Kmax-Kmin,和AU+的Kmax-Kmin
        vector<double> k_range = {uu_k_plus_k_range,au_k_plus_k_range,aa_k_minus_k_range,au_k_minus_k_range};

        double nuceo_to_cpg_u_add = std::atof(param_map["nuceo_to_cpg_u_add"].c_str());
        printf("nuceo_to_cpg_u_add : %.6f \n",nuceo_to_cpg_u_add);

        double nuceo_to_cpg_h_add = std::atof(param_map["nuceo_to_cpg_h_add"].c_str());
        printf("nuceo_to_cpg_h_add : %.6f \n",nuceo_to_cpg_h_add);

        double nuceo_to_cpg_m_minus = std::atof(param_map["nuceo_to_cpg_m_minus"].c_str());
        printf("nuceo_to_cpg_m_minus : %.6f \n",nuceo_to_cpg_m_minus);

        double nuceo_to_cpg_h_minus = std::atof(param_map["nuceo_to_cpg_h_minus"].c_str());
        printf("nuceo_to_cpg_h_minus : %.6f \n",nuceo_to_cpg_h_minus);


        //核小体对u+以及m-两个反应速率的促进程度，按照其u+,m-原始速率的倍数计算
        //vector<double> nuceo_to_cpg_efficiency = {0.032,0.032, -0.0016,-0.0016};
        vector<double> nuceo_to_cpg_efficiency = {nuceo_to_cpg_u_add,nuceo_to_cpg_h_add,nuceo_to_cpg_m_minus,nuceo_to_cpg_h_minus};

        vector<int> index_pos_list;

        //当前文件夹的全路径
        string now_dir = param_map["now_dir"];
        printf("now_dir : %s \n",now_dir.c_str());

        //输出文件夹的名字
        string output_dir_name = param_map["output_dir_name"]+"/";
        printf("output_dir_name : %s \n",param_map["output_dir_name"].c_str());

        //输入文件夹的名字
        string input_dir_name = param_map["input_dir_name"]+"/";
        printf("input_dir_name : %s \n",param_map["input_dir_name"].c_str());

        string output_dir = now_dir + output_dir_name;

        int rtn_code=exist_or_make_dir(output_dir);

        string input_bed_file_path = now_dir + input_dir_name + "chr1.bed";
        string nucleosome_pos_file_path = now_dir + input_dir_name + "nucleosome_positions.np";
        string ratio_file_dir = output_dir+"ratio/";//CpG ratio文件夹
        rtn_code=exist_or_make_dir(ratio_file_dir);
        string nuc_ratio_dir = output_dir + "nuc_ratio/"; //核小体 ratio 文件夹
        rtn_code=exist_or_make_dir(nuc_ratio_dir);
        string detail_file_dir = output_dir+"detail/";//CpG detail文件夹
        rtn_code=exist_or_make_dir(detail_file_dir);
        string nuc_detail_dir = output_dir + "nuc_detail/"; //核小体 detail 文件夹
        rtn_code=exist_or_make_dir(nuc_detail_dir);
        string bed_file_dir;
        string rd_without_dir;
        //时间步数
        int time_step = std::atoi(param_map["time_step"].c_str());
        printf("time_step : %d \n",time_step);

        int corr_d_max = std::atoi(param_map["corr_d_max"].c_str());
        printf("corr_d_max : %d \n",corr_d_max);

        bool calc_interval = std::atoi(param_map["calc_interval"].c_str());
        printf("calc_interval : %d \n",calc_interval);

        bool ignore_d = std::atoi(param_map["ignore_d"].c_str());
        printf("ignore_d : %d \n",ignore_d);

        int max_cells = std::atoi(param_map["max_cells"].c_str());
        printf("max_cells : %d \n",max_cells);

        bool simulation = std::atoi(param_map["simulation"].c_str());
        printf("simulation : %d \n",simulation);

        bool calc_corr = std::atoi(param_map["calc_corr"].c_str());
        printf("calc_corr : %d \n",calc_corr);

        bool real_chr_pos = std::atoi(param_map["real_chr_pos"].c_str());
        printf("real_chr_pos : %d \n",real_chr_pos);

        if(!real_chr_pos)
        {
            //超几何分布参数
            float geometric_p = 0.3;
            index_pos_list = construct_n_cpg_sites_for_exp_distribution(max_cpg_sites, geometric_p);
        }
        else
        {
            index_pos_list = get_pos_list_from_bed_file(input_bed_file_path,max_cpg_sites);
        }


        //根据CpG在数组index_pos_list中的下标,查找对应的核小体下标
        vector<int> neucleosome_id_list = get_neucleosome_id_list_from_pos_list(index_pos_list,nucleosome_pos_file_path);

        int nucleosome_list_size = (int) nucleosome_end_pos_list.size();
        vector<vector<int>> cpg_id_list = get_corresponding_cpg_id_list(nucleosome_list_size,neucleosome_id_list);

        init_cell = generate_CpG_in_methylation_percent_UHM(max_cpg_sites,m_ratio,u_ratio);

        if (simulation)
        {
            for(int round_i=round_start; round_i<=round_end; round_i++)
            {
                simulate(round_i, generations,time_step,init_cell,neucleosome_id_list,cpg_id_list,detail_file_dir,ratio_file_dir,nuc_detail_dir,nuc_ratio_dir,nucleosome_pos_file_path,propensity_list,nucleosome_propensity_list,nuceo_to_cpg_efficiency,nearby_promote_efficiency,k_range,update_nuc_status_frequency,index_pos_list,max_cells,out_start_gen,out_end_gen,round_start);
            }
        }
        if(calc_corr)
        {
            bed_file_dir = output_dir + "bed/";//detail的生成序列
            rtn_code=exist_or_make_dir(bed_file_dir);
            rd_without_dir=output_dir+"rd_without/";
            rtn_code=exist_or_make_dir(rd_without_dir);
            int round_size=round_end-round_start+1;
            sort_to_bed(round_size,index_pos_list,detail_file_dir,bed_file_dir,out_start_gen,out_end_gen,time_step,max_cpg_sites);
            calc_correlation_for_generations(out_start_gen,out_end_gen,time_step,bed_file_dir,rd_without_dir,corr_d_max,calc_interval,ignore_d=false);
            string out_mean_rd_file = output_dir+"rd_mean.csv";
            calc_mean_rd_from_rd_dir(rd_without_dir,out_mean_rd_file,corr_d_max,out_start_gen,out_end_gen,time_step);
        }
    }
}



int main(int argc, const char * argv[])
{
    //    start_simulation();
    //    system("say Mission completed!");
    vector<string> param_file_path_list = {"/Users/Ren/XCodeProjects/Cpp_Test/Cpp_Test/config.txt"};
    start_simualation_from(param_file_path_list);
    return 0;
}
