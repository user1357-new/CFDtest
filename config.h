#pragma once
#include <string>

class Config {
public:
    // 物理参数
    double gamma;          // 绝热指数
    double xMin;           // 计算域左边界
    double xMax;           // 计算域右边界
    double tEnd;           // 终止时间
    double cfl;            // CFL数
    double mu;             // 人工粘性系数
    double mu4;            // 四阶人工粘性系数
    
    // 数值参数
    int defaultNx;         // 默认网格数
    
    // 输出参数
    std::string outputDir; // 输出目录
    
    Config();
    
    // 从命令行参数解析
    void parseCommandLine(int argc, char* argv[]);
    
    // 打印配置信息
    void print() const;
};