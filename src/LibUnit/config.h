#pragma once
//
// Created by chenc on 2020/12/16.
//

#ifndef IRTKLIB_CONFIG_H
#define IRTKLIB_CONFIG_H

#include <iomanip>
#include <iostream>
#include <regex>
#include <fstream>
#include "rtklib.h"

#include <map>


using namespace std;

extern string StrTrim(string s) {
    if (s.empty()) return s;
    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
}

extern void SplitString(const string& s, vector<string>& v, string c) {
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (string::npos != pos2) {
        v.push_back(s.substr(pos1, pos2 - pos1));
        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

extern vector<string> MultiSplitStr(const string& s, const string& seperator) {
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x) {
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
            }
        }
        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x) {
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            }
            if (flag == 0) ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }
    return result;
}

extern vector<string> TextSplit(const string& in, const string& delim) {
    std::vector<std::string> ret;
    try {
        std::regex re{ delim };
        return std::vector<std::string>{std::sregex_token_iterator(in.begin(), in.end(), re, -1),
            std::sregex_token_iterator()};
    }
    catch (const std::exception& e) {
        std::cout << "error:" << e.what() << std::endl;
    }
    return ret;
}

class Config {
public:
    using Ptr_ = std::shared_ptr<Config>;

private:
    Config() = default;

    Config(Config&&) = delete;

    Config(const Config&) = delete;

    Config& operator=(Config&&) = delete;

    Config& operator=(const Config&) = delete;

    static Ptr_ config_info_;
    std::map<std::string, std::string> storage_;

public:
    ~Config() = default;

public:
    static Ptr_ GetInstance();

    bool Open(std::string config_file);

    template<typename T>
    T Get(std::string key) {
        transform(key.begin(), key.end(), key.begin(), ::tolower);
        if (storage_.count(key) > 0) {
            try {
                double value = stod(storage_[key]);
                return static_cast<T>(value);
            }
            catch (const std::exception& e) {
                std::cerr << e.what() << '\n';
            }
        }
        else {
            //            LOG(ERROR)<<"The key of "<<key<<" does not exist";
            //            getchar();
            return T(0x0);
        }
    }

    template<typename T>
    std::vector<T> GetArray(std::string key) {
        std::vector<T> data;
        transform(key.begin(), key.end(), key.begin(), ::tolower);
        if (storage_.count(key) > 0) {
            try {
                auto text = TextSplit(storage_[key], ",");
                for (auto index : text) {
                    double value = stod(index);
                    data.emplace_back(static_cast<T>(value));
                }
            }
            catch (const std::exception& e) {
                std::cerr << e.what() << '\n';
            }
        }
        else {
            //            LOG(ERROR)<<"The key of "<<key<<" does not exist";
            //            getchar();
        }
        return data;
    }
};

template<>
inline std::string Config::Get<std::string>(std::string key) {
    transform(key.begin(), key.end(), key.begin(), ::tolower);
    if (storage_.count(key) > 0) {
        try {
            return std::string(storage_[key]);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }
    }
    else {
        //        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
        //            getchar();
        return "";
    }
}

template<>
inline std::vector<std::string> Config::GetArray<std::string>(std::string key) {
    std::vector<std::string> data;
    transform(key.begin(), key.end(), key.begin(), ::tolower);
    if (storage_.count(key) > 0) {
        try {
            data = TextSplit(storage_[key], ",");
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }
    }
    else {
        //        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
        getchar();
    }
    return data;
}

//extern int ReadMatConf(prcopt_t *popt,solopt_t *sopt,filopt_t *fopt,MATFile *pmat);

#endif //IRTKLIB_CONFIG_H
