#ifndef _SERVER_H
#define _SERVER_H

#include "header.h"

class Server {
 public:
  Server() = default;
  explicit Server(int id);
  Server(const Server &server);
  Server(double o, double s, double g, double w) : o(o), s(s), g(g), w(w) {}

  int getId();
  double getO();
  double getS();
  double getG();
  double getW();

  void setID(int value);
  void setO(double value);
  void setS(double value);
  void setG(double value);
  void setW(double value);

 private:
  int idx = 0;   // 处理机id
  double o = 0;  // 传输启动开销
  double s = 0;  // 计算启动开销
  double g = 0;  // 传输单位任务量所需要时间
  double w = 0;  // 计算单位任务量所需要时间
};

#endif  //_SERVER_H
