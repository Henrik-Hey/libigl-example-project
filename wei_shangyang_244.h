#include <Eigen/Core>

void WT_Lifting_1(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
);

void WT_Lifting_2(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
);

void WT_Lifting_3(
    const std::vector<Eigen::Vector3d> vertices,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_prime
);

void WT_Lifting_5(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v1_prime,
    Eigen::Vector3d& v2_prime,
    Eigen::Vector3d& v3_prime,
    Eigen::Vector3d& v4_prime
);

void WT_Solve_Weights(
    const double n_0,
    const double n_1,
    const double n_2,
    const double n_3,
    Eigen::Vector4d& W
);

void WT_Lifting_6(
    const Eigen::Vector3d v_i,
    const double omega_i,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_i_prime
);