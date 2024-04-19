#include <iostream>

#include "RK.h"

const double GRAV_PARAMETER = 3.986e14;

int main() {

    const double rho = 7e6;
    const double velocity = std::sqrt(GRAV_PARAMETER / rho);

    const Eigen::Vector<double, 6> RV = {rho, 0, 0, 0, velocity, 0};
    const double initialTime = 0;
    const double step = 1;

    const State initialState = {RV, initialTime};

    /*Проверяем работу РК*/

    const double T = 2 * M_PI * rho * std::sqrt(rho / GRAV_PARAMETER);
    const double endTime = 5 * T;

    const std::vector<State> solution = rungeKutta(initialState, step, endTime, GRAV_PARAMETER);

    std::cout << (solution[solution.size() - 1].time - initialState.time) - T << std::endl;
    std::cout << solution[solution.size() - 1].RV - initialState.RV;
}
