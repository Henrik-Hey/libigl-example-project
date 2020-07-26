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