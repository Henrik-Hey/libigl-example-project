#include <Eigen/Core>

void WT_Lifting_1(
	const Eigen::Vector3d v,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
		  	Eigen::Vector3d &v_prime
);

void WT_Lifting_1_inverse(
	const Eigen::Vector3d v_prime,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
				Eigen::Vector3d &v
);

void WT_Lifting_2(
	const Eigen::Vector3d v,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
		  	Eigen::Vector3d &v_prime
);

void WT_Lifting_2_inverse(
	const Eigen::Vector3d v_prime,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
				Eigen::Vector3d &v
);

void WT_Lifting_3(
	const std::vector<Eigen::Vector3d> vertices,
	const Eigen::Vector3d v,
				Eigen::Vector3d &v_prime
);

void WT_Lifting_3_inverse(
	const std::vector<Eigen::Vector3d> vertices,
	const Eigen::Vector3d v_prime,
				Eigen::Vector3d &v
);

void WT_Scaling(
	const Eigen::Vector3d v,
	const double n,
				Eigen::Vector3d &v_prime
);

void WT_Scaling_inverse(
	const Eigen::Vector3d v_prime,
	const double n,
				Eigen::Vector3d &v
);

void WT_Lifting_4(
	const Eigen::Vector3d v,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
	const Eigen::Vector3d v3,
	const Eigen::Vector3d v4,
				Eigen::Vector3d &v_prime
);

void WT_Lifting_4_inverse(
	const	Eigen::Vector3d v_prime,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
	const Eigen::Vector3d v3,
	const Eigen::Vector3d v4,
				Eigen::Vector3d &v
);

void WT_Lifting_5(
	const Eigen::Vector3d v,
	const Eigen::Vector3d v1,
	const Eigen::Vector3d v2,
	const Eigen::Vector3d v3,
	const Eigen::Vector3d v4,
				Eigen::Vector3d &v1_prime,
				Eigen::Vector3d &v2_prime,
				Eigen::Vector3d &v3_prime,
				Eigen::Vector3d &v4_prime
);

void WT_Lifting_5_inverse(
	const Eigen::Vector3d v,
	const Eigen::Vector3d v1_prime,
	const Eigen::Vector3d v2_prime,
	const Eigen::Vector3d v3_prime,
	const Eigen::Vector3d v4_prime,
				Eigen::Vector3d &v1,
				Eigen::Vector3d &v2,
				Eigen::Vector3d &v3,
				Eigen::Vector3d &v4
);

void WT_Solve_Weights(
	const double n_0,
	const double n_1,
	const double n_2,
	const double n_3,
				Eigen::Vector4d &W
);

void WT_Lifting_6(
	const Eigen::Vector3d v_i,
	const double omega_i,
	const Eigen::Vector3d v,
				Eigen::Vector3d &v_i_prime
);

void WT_Lifting_6_inverse(
	const Eigen::Vector3d v_i_prime,
	const double omega_i,
	const Eigen::Vector3d v,
				Eigen::Vector3d &v_i
);
