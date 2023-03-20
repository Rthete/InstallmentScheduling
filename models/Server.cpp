/*
 * @FilePath: \InstallmentScheduling\models\Server.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-13 22:38:44
 * @LastEditTime: 2023-03-15 16:02:31
 */

#include "Server.h"

Server::Server(int id) {
    this->idx = id;
}

Server::Server(const Server &server) {
    this->idx = server.idx;
    this->o = server.o;
    this->s = server.s;
    this->g = server.g;
    this->w = server.w;
}

/**
 * @brief Set server's ID.
 * 
 * @param value 
 */
void Server::setID(int value) {
    this->idx = value;
}

/**
 * @brief Get server's ID.
 * 
 * @return int 
 */
int Server::getId() {
    return this->idx;
}

/**
 * @brief Set server's communication time g.
 * 
 * @param value 
 */
void Server::setG(double value) {
    this->g = value;
}

/**
 * @brief Get server's communication time g.
 * 
 * @return double 
 */
double Server::getG() {
    return this->g;
}

/**
 * @brief Set server's communication startup time o.
 * 
 * @param value 
 */
void Server::setO(double value) {
    this->o = value;
}

/**
 * @brief Get server's communication startup time o.
 * 
 * @return double 
 */
double Server::getO() {
    return this->o;
}

/**
 * @brief Set server's computation startup overhead s.
 * 
 * @param value 
 */
void Server::setS(double value) {
    this->s = value;
}

/**
 * @brief Get server's computation startup overhead s.
 * 
 * @return double 
 */
double Server::getS() {
    return this->s;
}

/**
 * @brief Get server's computation time w.
 * 
 * @param value 
 */
void Server::setW(double value) {
    this->w = value;
}

/**
 * @brief Get server's computation time w.
 * 
 * @return double 
 */
double Server::getW() {
    return this->w;
}
