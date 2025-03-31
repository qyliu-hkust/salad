#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
ofstream logStream;

// 加载数据返回data vector
template <typename K> vector<K> load_data(const string &filename) {
//    cout << "Test1 Dataset:" << filename << endl;
    /* Open file. */
    ifstream in(filename, ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char*>(&n_keys), sizeof(K));
//    cout << "Number of keys: " << n_keys << " " << sizeof (K) << endl;

    /* Initialize vector. */
    vector<K> data;
    data.resize(n_keys);
    // cout << data.size() << endl;

    /* Read keys. */
    in.read(reinterpret_cast<char*>(data.data()), n_keys * sizeof(K));  
    in.close();
    // set<K> s(data.begin(), data.end());
    // data = vector<K>(s.begin(), s.end());
    // return data;
//    for (int i = 0; i < 50; i++)
//    {
//        cout << data[i] << " " << endl;
//    }
    /* Sort the data in increasing order. */

    // sort(data.begin(), data.end());
    // for (int i = 0; i < 50; i++)
    // {
    //     cout << data[i] << " " << endl;
    // }
    // auto newEnd = unique(data.begin(), data.end());
    // data.erase(newEnd, data.end());
    sort(data.begin(), data.end());
    return data;
}

string to_string_int128(__int128_t n) {
    if (n == 0) return "0";  
    string res;
    bool negative = n < 0;  
    if (negative) n = -n;  

    while (n > 0) {  
        res += (n % 10) + '0'; // 把数字转换为字符  
        n /= 10;  
    }  

    if (negative) res += '-';  
    reverse(res.begin(), res.end());
    return res;  
}

string formatted_time(int hours_offset = 8){
    // set current time as log file name
    chrono::seconds offset_seconds(hours_offset * 3600);
    auto now = chrono::system_clock::now();
    now = now + offset_seconds;
    auto in_time_t = chrono::system_clock::to_time_t(now);
    // decode time
    tm buf;
    localtime_r(&in_time_t, &buf); // Linux
//     localtime_s(&buf, &in_time_t); // Windows

    // get time
    ostringstream ss;
    ss << put_time(&buf, "%Y%m%d_%H%M%S"); // 格式化为YYYYMMDD_HHMMSS
    return ss.str();
}

void output_message(string message){
    cout << message << endl;
    logStream << message << endl;
}

// const string create_time = formatted_time();
inline void delete_files_in_directory(string code_out_path) {
    filesystem::path directoryPath = code_out_path; // 替换为实际的目录路径
    if (!filesystem::exists(directoryPath)) {
        cerr << "Directory does not exist: " << directoryPath << endl;
        return;
    }
    try {
        for (const auto& entry : filesystem::directory_iterator(directoryPath)) {
            if (filesystem::is_regular_file(entry.path())) {
                filesystem::remove(entry.path());
            }
        }
        output_message("Files in " + code_out_path + " have been cleared.");
    } catch (const filesystem::filesystem_error& ex) {
        cerr << "Error: " << ex.what() << endl;
    }
}

void create_log_path(const string home_path, string dataset_name, string create_time=formatted_time(), string exp_name=""){
    // delete_files_in_directory(home_path + "code/out/");
    // cout << "create_time " << create_time << endl;
    error_code ec;
    filesystem::path dir_path(home_path + "result/" + exp_name + "/" + create_time);
    filesystem::create_directories(dir_path, ec);
    if(ec){
        cerr << "create log directory failed: " << ec.message() << endl;
        return;
    }
    string log_path = home_path + "result/" +  exp_name + "/" + create_time + "/" + exp_name + "_" + dataset_name + "_" + create_time  + ".txt";
    output_message("——————————————Dataset:" + dataset_name + "——————————————");
    logStream.open(log_path);
    if(!logStream.is_open()){
        cerr << "open log file failed" << endl;
    }
}

std::vector<std::string> split_str(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        // 去除前后空格
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        tokens.push_back(token);
    }
    return tokens;
}

// 转换为 double 的函数
std::vector<double> parseDoubles(const std::string &line, char delimiter) {
    std::vector<double> values;
    std::vector<std::string> tokens = split_str(line, delimiter);

    for (const auto &token : tokens) {
        try {
            // 尝试将每个 token 转换为 double
            values.push_back(std::stod(token));
        } catch (const std::exception &e) {
            // 如果无法转换，跳过或处理异常
            std::cerr << "Error parsing token: " << token << " (" << e.what() << ")\n";
        }
    }
    return values;
}