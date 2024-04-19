//
// Created by ivankhripunov on 19.04.24.
//

#ifndef RK_RK_H
#define RK_RK_H

#include "eigen/Eigen/Dense"
#include <vector>

struct State {
    Eigen::Vector<double, 6> RV;
    double time;
};

[[nodiscard]] Eigen::Vector<double, 3>
centralFieldAcceleration(const Eigen::Vector3d &relativePosition, const double gravParameter) {

    const double positionNorm = relativePosition.norm();
    const double positionNormCube = positionNorm * positionNorm * positionNorm;

    return -gravParameter * relativePosition / positionNormCube;
}

[[nodiscard]] Eigen::Vector<double, 6>
centralFieldRHS(const Eigen::Vector<double, 6> &RV, const double time, const double gravParameter) {

    const Eigen::Vector3d position = RV.segment<3>(0);
    const Eigen::Vector3d velocity = RV.segment<3>(3);
    const Eigen::Vector3d acceleration = centralFieldAcceleration(position, gravParameter);

    Eigen::Vector<double, 6> result;

    result.segment<3>(0) = velocity;
    result.segment<3>(3) = acceleration;

    return result;
}

[[nodiscard]] State rungeKuttaStep(const State &state, const double step, const double gravParameter) {

    const Eigen::Vector<double, 6> k1 = centralFieldRHS(state.RV, state.time, gravParameter);
    const Eigen::Vector<double, 6> k2 = centralFieldRHS(state.RV + step * k1 / 2, state.time + step / 2, gravParameter);
    const Eigen::Vector<double, 6> k3 = centralFieldRHS(state.RV + step * k2 / 2, state.time + step / 2, gravParameter);
    const Eigen::Vector<double, 6> k4 = centralFieldRHS(state.RV + step * k3, state.time + step, gravParameter);

    const Eigen::Vector<double, 6> updatedRV = state.RV + step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    const double updatedTime = state.time + step;

    return {updatedRV, updatedTime};
}

[[nodiscard]] std::vector<State>
rungeKutta(const State &initialState, const double step, const double endTime, const double gravParameter) {

    const auto solutionSize = static_cast<std::size_t>(((endTime - initialState.time) / step)) + 2;

    std::vector<State> solution;
    solution.reserve(solutionSize);
    solution.push_back(initialState);

    for (std::size_t i = 1; i < solutionSize - 1; ++i) {
        solution.push_back(rungeKuttaStep(solution[i - 1], step, gravParameter));
    }

    const double finalStep = endTime - solution[solutionSize - 2].time;

    solution.push_back(rungeKuttaStep(solution[solutionSize - 2], finalStep, gravParameter));

    return solution;
}

#endif //RK_RK_H
